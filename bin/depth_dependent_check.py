#!/usr/bin/env python
"""Decide which libraries are deep enough to trigger the depth-dependent analyses.

Read one mosdepth summary per library and find the shallowest checked contig in each. A library
passes when the mean depth of every checked contig is greater than --min-contig-depth. The
workflow runs the depth-dependent analyses when at least --min-libraries libraries pass.

Checked contigs are chosen in one of two ways. With --contig-groups, keep contigs whose
contig_type is in --contig-types (default autosomal). Sex chromosomes are left out by default
because chrY has no coverage in female samples and chrX has half coverage in male samples.
Without --contig-groups, keep every contig whose name does not match --exclude-regex, which
drops controls, organelles, sex chromosomes and unplaced or alternate contigs.

Library names come from the summary file names: <library>.mosdepth.summary.txt.

Outputs:
  depth_dependent_check.tsv  one row per library: the shallowest checked contig, its mean
                               depth, the number of checked contigs and whether the library passes.
  depth_dependent_contig_depth.tsv   mean depth of every contig in every library, with a column
                               marking the checked contigs.
"""
import argparse
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from library_order import natural_sorted  # noqa: E402

DEFAULT_CONTIG_TYPES = 'autosomal'
DEFAULT_EXCLUDE_REGEX = (
    r'(?i)^((chr)?(M|MT|X|Y|W|Z|EBV|Un.*)|.*_(random|alt|fix|decoy)|phage_.*|plasmid_.*|'
    r'lambda.*|puc19.*|T4|Xp12|chrC|Pt|chloroplast|mitochondri.*)$'
)
SUFFIX = '.mosdepth.summary.txt'


def read_contig_groups(path, contig_types):
    keep = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        try:
            name_i = header.index('contig_name')
            type_i = header.index('contig_type')
        except ValueError:
            sys.exit(f'{path}: needs contig_name and contig_type columns, got {header}')
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) > type_i and fields[type_i] in contig_types:
                keep.append(fields[name_i])
    if not keep:
        sys.exit(f'no contigs in {path} have contig_type in {sorted(contig_types)}')
    return keep


def read_summary(path):
    """Return {contig: (length, mean depth)} in file order, without the total and region rows."""
    depths = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        if header[:4] != ['chrom', 'length', 'bases', 'mean']:
            sys.exit(f'{path}: not a mosdepth summary (header {header})')
        for line in fh:
            chrom, length, _bases, mean = line.rstrip('\n').split('\t')[:4]
            if chrom == 'total' or chrom.endswith('_region'):
                continue
            depths[chrom] = (int(length), float(mean))
    return depths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('summaries', nargs='+', help=f'<library>{SUFFIX} files')
    ap.add_argument('--min-contig-depth', type=float, required=True,
                    help='a library passes when every checked contig has a mean depth above this')
    ap.add_argument('--contig-groups', help='contig_groups.tsv with a contig_type column')
    ap.add_argument('--contig-types', default=DEFAULT_CONTIG_TYPES,
                    help=f'comma-separated contig_type values to check (default {DEFAULT_CONTIG_TYPES})')
    ap.add_argument('--exclude-regex', default=DEFAULT_EXCLUDE_REGEX,
                    help='contig names to leave out when no --contig-groups is given')
    ap.add_argument('--check-out', default='depth_dependent_check.tsv')
    ap.add_argument('--depth-out', default='depth_dependent_contig_depth.tsv')
    args = ap.parse_args()

    by_library = {}
    for path in args.summaries:
        name = os.path.basename(path)
        if not name.endswith(SUFFIX):
            sys.exit(f'{path}: expected a name ending in {SUFFIX}')
        by_library[name[:-len(SUFFIX)]] = read_summary(path)
    libraries = natural_sorted(by_library)

    if args.contig_groups:
        checked_contigs = read_contig_groups(args.contig_groups, set(args.contig_types.split(',')))
    else:
        exclude = re.compile(args.exclude_regex)
        first = by_library[libraries[0]]
        checked_contigs = [c for c in first if not exclude.match(c)]
        if not checked_contigs:
            sys.exit('every contig matches --exclude-regex; give --contig-groups instead')
    checked_set = set(checked_contigs)

    # mosdepth lists every contig in the BAM header, even one with no reads. So a checked
    # contig missing from a summary means the names differ, such as chr1 against 1.
    seen = set().union(*(depths.keys() for depths in by_library.values()))
    missing = [c for c in checked_contigs if c not in seen]
    if len(missing) == len(checked_contigs):
        sys.exit(f'none of the {len(checked_contigs)} checked contigs (e.g. {checked_contigs[:3]}) '
                 f'appears in any mosdepth summary (e.g. {sorted(seen)[:3]}). Check that '
                 f'contig_groups.tsv uses the same contig names as the BAMs.')
    if missing:
        print(f'[depth_dependent_check] WARNING: {len(missing)} checked contigs are in no mosdepth '
              f'summary and count as depth 0: {", ".join(missing[:10])}'
              f'{" ..." if len(missing) > 10 else ""}', file=sys.stderr)

    with open(args.check_out, 'w') as check_fh, open(args.depth_out, 'w') as depth_fh:
        check_fh.write('library\tn_checked_contigs\tshallowest_contig\tshallowest_mean_depth\t'
                      'min_contig_depth\tpasses\n')
        depth_fh.write('library\tcontig\tlength\tmean_depth\tchecked_contig\n')
        for lib in libraries:
            depths = by_library[lib]
            for contig, (length, mean) in depths.items():
                depth_fh.write(f'{lib}\t{contig}\t{length}\t{mean:.4f}\t'
                               f'{"yes" if contig in checked_set else "no"}\n')
            # A checked contig missing from the summary has no reads, so it counts as depth 0.
            shallowest = min(checked_contigs, key=lambda c: depths.get(c, (0, 0.0))[1])
            low = depths.get(shallowest, (0, 0.0))[1]
            passes = low > args.min_contig_depth
            check_fh.write(f'{lib}\t{len(checked_contigs)}\t{shallowest}\t{low:.4f}\t'
                          f'{args.min_contig_depth:g}\t{"yes" if passes else "no"}\n')
    print(f'[depth_dependent_check] {len(libraries)} libraries, {len(checked_contigs)} checked contigs',
          file=sys.stderr)


if __name__ == '__main__':
    main()
