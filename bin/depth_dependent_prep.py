#!/usr/bin/env python
"""Convert one library's CpG methylKit file into per-C counts and a depth histogram.

MethylDackel writes one row per cytosine: the C of each CpG on the forward strand and
the C on the reverse strand are separate rows. This script keeps them separate. Each C
is keyed by its contig, its 1-based position and its strand, and its depth is the
reads on that strand alone.

Rows on contigs missing from the reference table are dropped. That table already
excludes spike-in controls and organelles (see count_reference_cpgs.py).

Outputs:
  <library>.cs.parquet     contig_idx, pos, strand (0 = F, 1 = R), cov, nmeth,
                           sorted by contig_idx, pos, strand
  <library>.depth_hist.tsv library, depth, n_c. Depth 0 counts reference Cs with no
                           reads. The last row holds every C at or above --max-depth.
"""
import argparse
import sys

import polars as pl
import pyarrow as pa
import pyarrow.csv as pacsv


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('methylkit', help='MethylDackel CpG methylKit file (.gz ok)')
    ap.add_argument('reference_cpgs', help='output of count_reference_cpgs.py')
    ap.add_argument('--library', required=True)
    ap.add_argument('--max-depth', type=int, default=1000,
                    help='histogram top bin; Cs at or above it share one row (default 1000)')
    args = ap.parse_args()

    ref = pl.read_csv(args.reference_cpgs, separator='\t',
                      schema_overrides={'contig': pl.Utf8})
    contig_idx = {c: i for i, c in enumerate(ref['contig'].to_list())}
    n_ref = int(ref['n_c'].sum())

    # Stream the file through pyarrow, which decompresses gzip as it reads, so the
    # whole file is never held as text in memory.
    reader = pacsv.open_csv(
        args.methylkit,
        read_options=pacsv.ReadOptions(block_size=16 << 20),
        parse_options=pacsv.ParseOptions(delimiter='\t'),
        convert_options=pacsv.ConvertOptions(
            include_columns=['chr', 'base', 'strand', 'coverage', 'freqC'],
            # freqC is right-aligned with leading spaces, which the float parser accepts.
            column_types={'chr': pa.string(), 'base': pa.uint32(), 'strand': pa.string(),
                          'coverage': pa.uint32(), 'freqC': pa.float64()}),
    )
    contigs = list(contig_idx)
    n_rows_in = 0
    parts = []
    for batch in reader:
        n_rows_in += batch.num_rows
        parts.append(
            pl.from_arrow(batch)
            .filter(pl.col('chr').is_in(contigs))
            .select(
                pl.col('chr').replace_strict(contig_idx, return_dtype=pl.UInt32).alias('contig_idx'),
                pl.col('base').alias('pos'),
                (pl.col('strand') == 'R').cast(pl.UInt8).alias('strand'),
                pl.col('coverage').alias('cov'),
                # freqC has two decimals, so this rounding is exact below 10,000x depth.
                (pl.col('coverage') * pl.col('freqC') / 100.0).round(0).cast(pl.UInt32)
                  .alias('nmeth'),
            )
        )
    # A library with no usable reads writes a methylKit file with only a header, so the
    # reader yields no batches. It still gets an empty table and an all-zero histogram.
    schema = {'contig_idx': pl.UInt32, 'pos': pl.UInt32, 'strand': pl.UInt8,
              'cov': pl.UInt32, 'nmeth': pl.UInt32}
    cs = (pl.concat(parts) if parts else pl.DataFrame(schema=schema)).sort(['contig_idx', 'pos', 'strand'])
    del parts
    n_dup = cs.select(pl.struct(['contig_idx', 'pos', 'strand']).is_duplicated().sum()).item()
    if n_dup:
        sys.exit(f'{args.methylkit}: {n_dup} rows repeat a contig/position/strand')

    cs.write_parquet(f'{args.library}.cs.parquet', row_group_size=1_000_000)

    hist = (
        cs.select(pl.col('cov').clip(upper_bound=args.max_depth).alias('depth'))
        .group_by('depth').len('n_c')
    )
    n_covered = cs.height
    n_zero = n_ref - n_covered
    if n_zero < 0:
        print(f'[depth_dependent_prep] WARNING: {n_covered} covered Cs exceed {n_ref} reference '
              f'Cs; depth-0 count clamped to 0', file=sys.stderr)
        n_zero = 0
    hist = (
        pl.concat([
            pl.DataFrame({'depth': [0], 'n_c': [n_zero]},
                         schema={'depth': pl.UInt32, 'n_c': pl.UInt32}),
            hist.with_columns(pl.col('n_c').cast(pl.UInt32)),
        ])
        .sort('depth')
        .with_columns(pl.lit(args.library).alias('library'))
        .select(['library', 'depth', 'n_c'])
    )
    hist.write_csv(f'{args.library}.depth_hist.tsv', separator='\t')

    mean_depth = float(cs['cov'].sum()) / n_ref if n_ref else 0.0
    # Rows on contigs outside the reference table are dropped on purpose: spike-in controls
    # and organelles are excluded there. When every row is dropped, the contig names most
    # likely do not match the reference, so say so loudly.
    n_dropped = n_rows_in - n_covered
    print(f'[depth_dependent_prep] {args.library}: {n_rows_in} rows -> {n_covered} Cs '
          f'({n_covered / n_ref:.1%} of {n_ref} reference Cs; {n_dropped} rows on other contigs '
          f'dropped), mean depth {mean_depth:.1f}', file=sys.stderr)
    if n_rows_in and not n_covered:
        print(f'[depth_dependent_prep] WARNING: {args.library}: none of the {n_rows_in} rows is on a '
              f'contig in {args.reference_cpgs}. Check that the methylKit file and the reference '
              f'use the same contig names.', file=sys.stderr)


if __name__ == '__main__':
    main()
