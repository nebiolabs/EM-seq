#!/usr/bin/env python
"""Merge per-library C tables and compute cross-library depth-dependent statistics.

The merge walks the genome one contig at a time. For each contig it builds a
site x library matrix over every C covered in at least one library, updates
running sums, and discards the matrix. Peak memory therefore follows the largest
contig, not the genome. Every statistic is exact; nothing is subsampled.

The two Cs of a CpG  are counted
separately. A C counts as "covered" in a
library when its depth is at least --min-depth. Methylation level is nmeth / cov at
covered Cs.

Outputs (prefix depth_dependent_). Tables marked * hold one block of rows per minimum depth
in --report-depths (column min_depth), so the report can switch between depths.
  pairwise.tsv *   one row per library pair. Gives shared and private counts, Jaccard,
                   methylation Pearson r on shared Cs, the mean methylation of shared
                   and single-library Cs, and Pearson r of log1p(depth) over all
                   reference Cs (the same at every depth).
  pair_hist2d.tsv * 2-D methylation histogram of each pair. Bin -1 marks Cs the
                   library does not cover, so the bin_a = -1 row holds Cs only b
                   covers and the bin_b = -1 column Cs only a covers.
  meth_hist.tsv *  per-library methylation histogram of all covered Cs and of private
                   Cs, which are covered in that library alone.
  summary.tsv *    per-library totals.
  upset.tsv *      number of reference Cs per coverage pattern. The empty pattern
                   counts Cs covered in no library.
  depth_hist2d.tsv 2-D depth histogram of each pair over every reference C, uncovered
                   ones included. depth_a and depth_b are the lower edge of 2x bins; the
                   last bin (100) holds every C at 100x or deeper.
  core.tsv         the core set: libraries chosen at --core-depth (whole name groups when
                   the names group) as the largest set whose set Jaccard of covered Cs
                   reaches --core-jaccard, with that Jaccard and the common and union counts.
  core_pairwise.tsv * methylation r between core libraries on the Cs that every core
                   library covers at that depth (n_common, the same for every pair).
  core_pair_hist2d.tsv * 2-D methylation histograms of core pairs on those common Cs.
  core_meth_hist.tsv * each core library's methylation histogram on those common Cs.
  matrix.parquet   (with --matrix) contig, pos, strand, then cov_<lib> and nmeth_<lib>
                   for every C covered in at least one library.
"""
import argparse
import itertools
import os
import re
import sys
from collections import Counter

import numpy as np
import polars as pl

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from library_order import natural_key  # noqa: E402
from feature_depth_bokeh import analyze_name_components, get_component  # noqa: E402

MAX_UPSET_LIBRARIES = 64  # patterns are bitmasks in a uint64
DEPTH_BIN_WIDTH = 2
DEPTH_BINS = 51           # 0-1, 2-3, ..., 98-99, then >= 100


def library_name(path):
    return re.sub(r'\.cs\.parquet$', '', os.path.basename(path))


def pattern_codes(mask):
    """One uint64 per row with bit j set when column j is True (up to 64 columns)."""
    packed = np.packbits(mask, axis=1, bitorder='little')
    padded = np.zeros((mask.shape[0], 8), dtype=np.uint8)
    padded[:, :packed.shape[1]] = packed
    return padded.view('<u8').ravel()


def choose_core(patterns, libs, threshold):
    """Largest set of libraries whose set Jaccard of covered Cs reaches threshold.

    Set Jaccard = Cs covered by every library in the set / Cs covered by any of them,
    both read from the coverage-pattern counts (pattern = bitmask of covering
    libraries). Libraries are dropped in whole name groups (0.1ng_*, 1ng_*, ...) when
    the names group, so replicates stay together, and one at a time otherwise. Ties in
    size go to the set with more common Cs. The search is exhaustive up to 16 units and
    greedy (drop the unit whose removal raises the Jaccard most) beyond that.

    When no set of whole groups reaches the threshold, the search runs again on single
    libraries. Mixed naming schemes can group by the wrong name part, such as protocol
    instead of input amount, and then every group holds a shallow library.
    """
    group_idx, _ = analyze_name_components(libs)
    groups = {}
    for j, lib in enumerate(libs):
        groups.setdefault(get_component(lib, group_idx) or lib, []).append(j)
    singles = [[j] for j in range(len(libs))]
    grouped = list(groups.values())
    nonzero = [(code, n) for code, n in patterns.items() if code]
    pats = (np.array([c for c, _ in nonzero], dtype=np.uint64),
            np.array([n for _, n in nonzero], dtype=np.int64))
    if len(grouped) >= 2 and not all(len(u) == 1 for u in grouped):
        found = _search_core(grouped, pats, libs, threshold)
        if found is not None:
            return found
    return _search_core(singles, pats, libs, threshold)


def _search_core(units, pats, libs, threshold):
    """choose_core's search over one way of splitting the libraries into units."""
    unit_masks = [sum(1 << j for j in u) for u in units]

    codes, counts = pats

    def score(mask):
        hit = codes & np.uint64(mask)
        return int(counts[hit == np.uint64(mask)].sum()), int(counts[hit != 0].sum())

    def result(mask):
        inter, union = score(mask)
        return {'members': [j for j in range(len(libs)) if mask >> j & 1],
                'jaccard': inter / union if union else 0.0, 'inter': inter, 'union': union}

    best = None
    if len(units) <= 16:
        # Units can hold different numbers of libraries, so fewer units can still mean more
        # libraries. Search every combination and rank by library count, then common Cs.
        for k in range(len(units), 0, -1):
            for combo in itertools.combinations(range(len(units)), k):
                mask = sum(unit_masks[u] for u in combo)
                n_libs = bin(mask).count('1')
                if n_libs < 2:
                    continue
                inter, union = score(mask)
                if union and inter / union >= threshold and (best is None or (n_libs, inter) > best[0]):
                    best = ((n_libs, inter), mask)
        return result(best[1]) if best else None
    mask = sum(unit_masks)
    while bin(mask).count('1') >= 2:
        inter, union = score(mask)
        if union and inter / union >= threshold:
            return result(mask)
        drops = [mask & ~m for m in unit_masks if mask & m and bin(mask & ~m).count('1') >= 2]
        if not drops:
            break
        mask = max(drops, key=lambda m: (lambda iu: iu[0] / iu[1] if iu[1] else 0)(score(m)))
    return None


def add_pair_hist(out, bins, nb):
    """Add every library pair's joint histogram of bins to out.

    bins is (n, L) with values in [0, nb). Block (a, b) of out, rows a*nb.. and columns
    b*nb.., gains the count of rows where library a is in bin i and b in bin j. Both
    triangles are filled, so block (b, a) is the transpose of block (a, b).

    One bincount per pair replaces the product of an (n, L*nb) one-hot matrix with itself.
    That product is mostly zeros and costs (L*nb)^2 per row; this costs L(L+1)/2 per row.
    """
    L = bins.shape[1]
    # intp, because bincount converts any other integer type to it with a copy.
    cols = np.ascontiguousarray(bins.T, dtype=np.intp)
    scaled = cols * nb
    for a in range(L):
        for b in range(a, L):
            h = np.bincount(scaled[a] + cols[b], minlength=nb * nb).reshape(nb, nb)
            out[a * nb:(a + 1) * nb, b * nb:(b + 1) * nb] += h
            if b != a:
                out[b * nb:(b + 1) * nb, a * nb:(a + 1) * nb] += h.T


def pearson(n, sx, sy, sxx, syy, sxy):
    if n < 3:
        return float('nan')
    cov = sxy - sx * sy / n
    vx = sxx - sx * sx / n
    vy = syy - sy * sy / n
    if vx <= 0 or vy <= 0:
        return float('nan')
    return float(cov / np.sqrt(vx * vy))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('reference_cpgs', help='output of count_reference_cpgs.py')
    ap.add_argument('parquets', nargs='+', help='<library>.cs.parquet from depth_dependent_prep.py')
    ap.add_argument('--min-depth', type=int, default=10)
    ap.add_argument('--bins', type=int, default=20, help='methylation-level histogram bins')
    ap.add_argument('--report-depths', default='1,3,5,10,20,30',
                    help='minimum depths written to every depth-dependent table, for the '
                         'report\'s depth toggle (default 1,3,5,10,20,30; --min-depth is always added)')
    ap.add_argument('--matrix', action='store_true', help='write depth_dependent_matrix.parquet')
    ap.add_argument('--core-jaccard', type=float, default=0.85,
                    help='core set: the largest set of libraries (whole name groups when the names '
                         'group) whose set Jaccard of covered Cs reaches this (default 0.85; 0 disables)')
    ap.add_argument('--core-depth', type=int, default=1,
                    help='depth at which a C counts as covered when choosing the core set (default 1)')
    args = ap.parse_args()

    d = args.min_depth
    meth_depths = sorted({d} | {int(x) for x in args.report_depths.split(',') if x.strip()})
    # A depth of 0 would count every uncovered reference C as covered.
    if meth_depths[0] < 1 or args.core_depth < 1:
        sys.exit(f'--min-depth, --report-depths and --core-depth must be at least 1 '
                 f'(got depths {meth_depths}, core depth {args.core_depth})')
    upset_depths = meth_depths
    B = args.bins

    paths = sorted(args.parquets, key=lambda p: natural_key(library_name(p)))
    libs = [library_name(p) for p in paths]
    if len(set(libs)) != len(libs):
        sys.exit(f'duplicate library names: {libs}')
    N = len(libs)
    if N < 2:
        sys.exit('need at least two libraries')
    do_upset = N <= MAX_UPSET_LIBRARIES
    if not do_upset:
        print(f'[depth_dependent_merge] {N} libraries exceeds {MAX_UPSET_LIBRARIES}; skipping upset',
              file=sys.stderr)

    ref = pl.read_csv(args.reference_cpgs, separator='\t', schema_overrides={'contig': pl.Utf8})
    contigs = ref['contig'].to_list()
    ref_len = ref['length'].to_list()
    ref_n = ref['n_c'].to_list()
    n_ref_total = int(sum(ref_n))
    pairs = list(itertools.combinations(range(N), 2))

    # Every pairwise sum is a column inner product, so each accumulator below is an
    # N x N matrix filled by one matrix product per contig. With M the covered mask and
    # X the methylation level at covered sites (0 elsewhere):
    #   MM[a,b] = Cs covered in both a and b    (diagonal: covered in a)
    #   XM[a,b] = sum of a's level over those     (diagonal: over all of a's covered Cs)
    #   QM[a,b] = sum of a's squared level over those
    #   XX[a,b] = sum of a's level times b's level
    # One set per minimum depth; the --min-depth set feeds the main outputs.
    W = B + 1
    acc = {t: {'MM': np.zeros((N, N)), 'XM': np.zeros((N, N)), 'QM': np.zeros((N, N)),
               'XX': np.zeros((N, N)), 'hist2d': np.zeros((N * W, N * W), dtype=np.int64),
               'private': np.zeros((N, B), dtype=np.int64)}
           for t in meth_depths}
    # Coverage correlation over every reference C; uncovered ones contribute zeros.
    LL = np.zeros((N, N))
    L_sum = np.zeros(N)
    n_cov_sites = 0
    # hist2d[a*W + i, b*W + j] counts Cs where library a's level falls in bin i-1 and b's
    # in bin j-1, with W = B + 1. Bin -1 (index 0) means not covered, so the block for
    # (a, b) holds the joint histogram over shared Cs plus each library's Cs that the
    # other misses. The diagonal block (a, a) holds a's own histogram.
    upset = {t: Counter() for t in upset_depths}
    # 2-D depth histograms over every reference C, built the same way. Bin k holds depths
    # k*DEPTH_BIN_WIDTH up to the next bin; the last bin holds everything deeper.
    D = DEPTH_BINS
    depth2d = np.zeros((N * D, N * D), dtype=np.int64)

    writer = None
    if args.matrix:
        import pyarrow as pa
        import pyarrow.parquet as pq
        fields = [pa.field('contig', pa.string()), pa.field('pos', pa.uint32()),
                  pa.field('strand', pa.string())]
        for lib in libs:
            fields += [pa.field(f'cov_{lib}', pa.uint32()), pa.field(f'nmeth_{lib}', pa.uint32())]
        writer = pq.ParquetWriter('depth_dependent_matrix.parquet', pa.schema(fields), compression='zstd')

    scans = [pl.scan_parquet(p) for p in paths]

    def load_contig(ci, with_nmeth=True):
        """(keys, cov, nmeth) for every C covered in any library on contig ci, or None.

        nmeth is None when with_nmeth is False. Each library is held as plain numpy arrays
        and freed once it is copied into the n x N matrices, so the peak stays near one
        copy of the contig's data.
        """
        # One key per C: the two Cs of a CpG sit on opposite strands. A uint32 key holds
        # contigs up to 2 Gb, which halves the memory of the per-library keys.
        key_type = pl.UInt32 if 2 * ref_len[ci] + 2 < 2**32 else pl.Int64
        cols = ['cov', 'nmeth'] if with_nmeth else ['cov']
        per_lib = []
        for scan in scans:
            df = (scan.filter(pl.col('contig_idx') == ci)
                      .select((pl.col('pos').cast(key_type) * 2 + pl.col('strand')).alias('key'), *cols)
                      .collect())
            per_lib.append([df[c].to_numpy() for c in ['key', *cols]])
            del df
        if all(arrs[0].size == 0 for arrs in per_lib):
            return None
        # Union of keys through a bitmap over the contig: one pass, where np.unique on the
        # concatenation sorts every library's keys together.
        seen = np.zeros(2 * ref_len[ci] + 2, dtype=bool)
        for lib, (k, *_rest) in zip(libs, per_lib):
            # A position past the contig end means the methylKit came from another build.
            if k.size and k.max() >= seen.size:
                sys.exit(f'{lib}: position {k.max() // 2} on {contigs[ci]} is past the contig '
                         f'length {ref_len[ci]}; was it aligned to a different reference?')
            seen[k] = True
        keys = np.flatnonzero(seen)
        del seen
        cov = np.zeros((keys.size, N), dtype=np.uint32)
        nmeth = np.zeros((keys.size, N), dtype=np.uint32) if with_nmeth else None
        for j in range(N):
            k, *vals = per_lib[j]
            per_lib[j] = None
            if k.size:
                idx = np.searchsorted(keys, k)
                cov[idx, j] = vals[0]
                if with_nmeth:
                    nmeth[idx, j] = vals[1]
        return keys, cov, nmeth

    # Pass 2 works through each contig in chunks of rows, so its float64 temporaries (about
    # ten n x N arrays per depth) scale with the chunk rather than the contig. Every statistic
    # is a sum over rows, so the chunk size does not change the results.
    chunk_rows = max(10_000, 4_000_000 // N)

    # Pass 1: coverage patterns at --core-depth, to choose the core set before pass 2
    # restricts its statistics to the Cs every core library covers.
    core = None
    if args.core_jaccard > 0 and N <= MAX_UPSET_LIBRARIES:
        core_patterns = Counter()
        for ci in range(len(contigs)):
            loaded = load_contig(ci, with_nmeth=False)
            if loaded is None:
                continue
            cov = loaded[1]
            for start in range(0, cov.shape[0], chunk_rows):
                vals, counts = np.unique(pattern_codes(cov[start:start + chunk_rows] >= args.core_depth),
                                         return_counts=True)
                for v, c in zip(vals.tolist(), counts.tolist()):
                    core_patterns[v] += c
            del loaded, cov
        core = choose_core(core_patterns, libs, args.core_jaccard)
        if core is None:
            print(f'[depth_dependent_merge] no core set: no two libraries reach Jaccard '
                  f'{args.core_jaccard} at >= {args.core_depth}x', file=sys.stderr)
        else:
            print(f'[depth_dependent_merge] core set at >= {args.core_depth}x: '
                  f'{[libs[j] for j in core["members"]]} (Jaccard {core["jaccard"]:.3f})', file=sys.stderr)
    core_idx = np.array(core['members']) if core else None
    K = len(core_idx) if core else 0
    core_acc = {t: {'n': 0, 's': np.zeros(K), 'q': np.zeros(K), 'xx': np.zeros((K, K)),
                    'hist2d': np.zeros((K * B, K * B), dtype=np.int64)} for t in meth_depths}

    def add_chunk(cov, nmeth):
        """Add one block of rows (Cs) of the current contig to every accumulator."""
        level = nmeth / np.maximum(cov, 1)
        for t in meth_depths:
            a_t = acc[t]
            covered_t = cov >= t
            M = covered_t.astype(np.float64)
            X_t = np.where(covered_t, level, 0.0)
            a_t['MM'] += M.T @ M
            a_t['XM'] += X_t.T @ M
            a_t['QM'] += (X_t * X_t).T @ M
            a_t['XX'] += X_t.T @ X_t
            # Bin index per site and library; -1 where the library does not cover the site.
            bins_t = np.where(covered_t, np.minimum((X_t * B).astype(np.intp), B - 1), -1)
            # Bin -1 (not covered) shifts to index 0, so W = B + 1 bins per library.
            add_pair_hist(a_t['hist2d'], bins_t + 1, W)
            if core:
                # Core-set view: only Cs that every core library covers at this depth.
                rows = covered_t[:, core_idx].all(axis=1)
                Xc = X_t[rows][:, core_idx]
                c_t = core_acc[t]
                c_t['n'] += Xc.shape[0]
                c_t['s'] += Xc.sum(axis=0)
                c_t['q'] += (Xc * Xc).sum(axis=0)
                c_t['xx'] += Xc.T @ Xc
                bins_c = bins_t[rows][:, core_idx]
                add_pair_hist(c_t['hist2d'], bins_c, B)
            # Private Cs: covered at this depth in one library and no other.
            private_t = covered_t & (covered_t.sum(axis=1, keepdims=True) == 1)
            for a in range(N):
                a_t['private'][a] += np.bincount(bins_t[private_t[:, a], a], minlength=B)

        L = np.log1p(cov, dtype=np.float64)
        # In-place adds, because LL and L_sum belong to main().
        np.add(LL, L.T @ L, out=LL)
        np.add(L_sum, L.sum(axis=0), out=L_sum)

        depth_bin = np.minimum(cov // DEPTH_BIN_WIDTH, D - 1).astype(np.intp)
        add_pair_hist(depth2d, depth_bin, D)

        if do_upset:
            for t in upset_depths:
                codes = pattern_codes(cov >= t)
                vals, counts = np.unique(codes, return_counts=True)
                for v, c in zip(vals.tolist(), counts.tolist()):
                    upset[t][v] += c

    for ci, contig in enumerate(contigs):
        loaded = load_contig(ci)
        if loaded is None:
            # No reads anywhere: every reference C here is uncovered in every library.
            for t in upset_depths:
                upset[t][0] += ref_n[ci]
            n_cov_sites += ref_n[ci]
            depth2d[0::D, 0::D] += ref_n[ci]
            continue

        keys, cov_all, nmeth_all = loaded
        n = keys.size
        n_uncovered_ref = max(ref_n[ci] - n, 0)
        for start in range(0, n, chunk_rows):
            rows = slice(start, start + chunk_rows)
            add_chunk(cov_all[rows], nmeth_all[rows])
            if writer is not None:
                k = keys[rows]
                cols = {'contig': pa.array([contig] * k.size, pa.string()),
                        'pos': pa.array(k // 2, pa.uint32()),
                        'strand': pa.array(np.where(k % 2, 'R', 'F'), pa.string())}
                for j, lib in enumerate(libs):
                    cols[f'cov_{lib}'] = pa.array(cov_all[rows, j])
                    cols[f'nmeth_{lib}'] = pa.array(nmeth_all[rows, j])
                writer.write_table(pa.table(cols, schema=writer.schema))
        n_cov_sites += n + n_uncovered_ref
        # Reference Cs no library covers sit at depth 0 in every library.
        depth2d[0::D, 0::D] += n_uncovered_ref
        if do_upset:
            for t in upset_depths:
                upset[t][0] += n_uncovered_ref

        print(f'[depth_dependent_merge] {contig}: {n} Cs covered in >=1 library', file=sys.stderr)

    if writer is not None:
        writer.close()

    # ---- write outputs ----
    # Every depth-dependent table has one block of rows per minimum depth (column min_depth),
    # so the report can switch between them.
    def ratio(num, den):
        return float(num / den) if den else float('nan')

    cov_r = {(a, b): pearson(n_cov_sites, L_sum[a], L_sum[b], LL[a, a], LL[b, b], LL[a, b])
             for a, b in pairs}
    pw_rows, h2_rows, mh_rows, sum_rows = [], [], [], []
    for t in meth_depths:
        MM, XM, QM, XX, hist2d = (acc[t][k] for k in ('MM', 'XM', 'QM', 'XX', 'hist2d'))
        MMi = np.rint(MM).astype(np.int64)

        def block(a, b):
            return hist2d[a * W:(a + 1) * W, b * W:(b + 1) * W]

        for a, b in pairs:
            shared = int(MMi[a, b])
            only_a = int(MMi[a, a] - shared)
            only_b = int(MMi[b, b] - shared)
            pw_rows.append({
                'min_depth': t, 'lib_a': libs[a], 'lib_b': libs[b],
                'n_shared': shared, 'n_only_a': only_a, 'n_only_b': only_b,
                'jaccard': ratio(shared, shared + only_a + only_b),
                'meth_r': pearson(MM[a, b], XM[a, b], XM[b, a], QM[a, b], QM[b, a], XX[a, b]),
                'meth_mean_a_shared': ratio(XM[a, b], shared),
                'meth_mean_b_shared': ratio(XM[b, a], shared),
                'meth_mean_a_only': ratio(XM[a, a] - XM[a, b], only_a),
                'meth_mean_b_only': ratio(XM[b, b] - XM[b, a], only_b),
                'cov_log1p_r': cov_r[(a, b)],
            })
            # Bin -1 means not observed (below this depth) in that library. The (-1, -1)
            # cell counts Cs that only other libraries cover, so it is left out.
            h = block(a, b)
            for i, j in zip(*np.nonzero(h)):
                if i or j:
                    h2_rows.append((t, libs[a], libs[b], int(i) - 1, int(j) - 1, int(h[i, j])))
        private = acc[t]['private']
        for j, lib in enumerate(libs):
            for cls, h in (('covered', np.diag(block(j, j))[1:]), ('private', private[j])):
                for i in range(B):
                    mh_rows.append((t, lib, cls, i / B, (i + 1) / B, int(h[i])))
        n_covered = np.diag(MMi)
        for j, lib in enumerate(libs):
            sum_rows.append({
                'min_depth': t, 'library': lib, 'n_reference_c': n_ref_total,
                'n_covered': int(n_covered[j]), 'frac_covered': n_covered[j] / n_ref_total,
                'n_private': int(private[j].sum()),
                'mean_meth_covered': ratio(XM[j, j], n_covered[j]),
            })
    pl.DataFrame(pw_rows).write_csv('depth_dependent_pairwise.tsv', separator='\t')
    pl.DataFrame(h2_rows, orient='row', schema=['min_depth', 'lib_a', 'lib_b', 'bin_a', 'bin_b', 'n_c']
                 ).write_csv('depth_dependent_pair_hist2d.tsv', separator='\t')
    pl.DataFrame(mh_rows, orient='row',
                 schema=['min_depth', 'library', 'site_class', 'bin_lo', 'bin_hi', 'n_c']
                 ).write_csv('depth_dependent_meth_hist.tsv', separator='\t')
    pl.DataFrame(sum_rows).write_csv('depth_dependent_summary.tsv', separator='\t')

    # Core set: membership, then every core statistic on the Cs all core libraries cover.
    pl.DataFrame({
        'library': libs,
        'in_core': [bool(core) and j in core['members'] for j in range(N)],
        'core_depth': [args.core_depth] * N,
        'core_jaccard_threshold': [args.core_jaccard] * N,
        'core_jaccard': [core['jaccard'] if core else float('nan')] * N,
        'core_common_c': [core['inter'] if core else 0] * N,
        'core_union_c': [core['union'] if core else 0] * N,
    }).write_csv('depth_dependent_core.tsv', separator='\t')
    cpw, ch2, cmh = [], [], []
    if core:
        names = [libs[j] for j in core_idx]
        for t in meth_depths:
            c_t = core_acc[t]
            n_t = c_t['n']
            for a, b in itertools.combinations(range(K), 2):
                cpw.append((t, names[a], names[b], n_t,
                            pearson(n_t, c_t['s'][a], c_t['s'][b], c_t['q'][a], c_t['q'][b], c_t['xx'][a, b])))
                blk = c_t['hist2d'][a * B:(a + 1) * B, b * B:(b + 1) * B]
                for i, j in zip(*np.nonzero(blk)):
                    ch2.append((t, names[a], names[b], int(i), int(j), int(blk[i, j])))
            for a in range(K):
                own = np.diag(c_t['hist2d'][a * B:(a + 1) * B, a * B:(a + 1) * B])
                for i in range(B):
                    cmh.append((t, names[a], i / B, (i + 1) / B, int(own[i])))
    pl.DataFrame(cpw, orient='row', schema=['min_depth', 'lib_a', 'lib_b', 'n_common', 'meth_r']
                 ).write_csv('depth_dependent_core_pairwise.tsv', separator='\t')
    pl.DataFrame(ch2, orient='row', schema=['min_depth', 'lib_a', 'lib_b', 'bin_a', 'bin_b', 'n_c']
                 ).write_csv('depth_dependent_core_pair_hist2d.tsv', separator='\t')
    pl.DataFrame(cmh, orient='row', schema=['min_depth', 'library', 'bin_lo', 'bin_hi', 'n_c']
                 ).write_csv('depth_dependent_core_meth_hist.tsv', separator='\t')

    rows = []
    for a, b in pairs:
        h = depth2d[a * D:(a + 1) * D, b * D:(b + 1) * D]
        for i, j in zip(*np.nonzero(h)):
            rows.append((libs[a], libs[b], int(i) * DEPTH_BIN_WIDTH, int(j) * DEPTH_BIN_WIDTH,
                         int(h[i, j])))
    pl.DataFrame(rows, orient='row', schema=['lib_a', 'lib_b', 'depth_a', 'depth_b', 'n_c']
                 ).write_csv('depth_dependent_depth_hist2d.tsv', separator='\t')

    rows = []
    if do_upset:
        for t in upset_depths:
            for code, cnt in sorted(upset[t].items(), key=lambda kv: (-kv[1], kv[0])):
                members = [libs[j] for j in range(N) if code >> j & 1]
                rows.append((t, len(members), ','.join(members), cnt))
    pl.DataFrame(rows, orient='row', schema=['min_depth', 'n_libraries', 'libraries', 'n_c'],
                 ).write_csv('depth_dependent_upset.tsv', separator='\t')


if __name__ == '__main__':
    main()
