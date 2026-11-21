#!/usr/bin/env python
"""Plot coverage and methylation per feature category, one dot per category per library.

This replaces the "Feature Coverage vs Methylation" Tableau sheet. Each feature
category (exon, CDS, CpG island, promoter, ...) gets one value per library:

  relative coverage   mean per-base read depth over the category's features divided by
                      the genome-wide mean depth. 1.0 means the category is covered as
                      deeply as the genome on average.
                        category depth = sum(region mean depth x region length) /
                                         sum(region length), from mosdepth --by
                        genome depth   = bases / length over the contigs that carry
                                         features, from the mosdepth summary; spike-in
                                         controls carry no features and are excluded
  mean methylation    mean of the per-feature methylation levels from
                      feature_methylation (CpGs with depth >= 5 only)

Per-base depth, unlike a count of overlapping fragments, does not inflate short features:
a 300 bp fragment overlaps a 60 bp promoter from 360 start positions but adds at most
one read to the depth of any base in it.

Outputs:
  feature_cov_meth.html         three views: coverage vs methylation, coverage by
                                library, methylation by library
  feature_cov_meth_summary.tsv  the plotted values, one row per library, context and category
"""
import argparse
import os
import re
import sys

import polars as pl
from bokeh.layouts import column
from bokeh.models import (ColumnDataSource, Div, FactorRange, HoverTool, Legend, LegendItem,
                          Span, TabPanel, Tabs)
from bokeh.plotting import figure

from feature_depth_bokeh import analyze_name_components, get_component, library_sort_key
from report_overlay import save_report
from report_style import OTHER, SERIES, TEXT, TEXT_2, note, style

# feature_methylation names two categories in the plural; the regions BED (and
# combine_counts) use the singular. Methylation names are mapped onto the BED's.
FEATURE_ALIASES = {'cpg_islands': 'cpg_island', 'epd_promoters': 'epd_promoter'}
CONTEXT_SUFFIX = re.compile(r'(_(CpG|CHG|CHH))+$')
MARKERS = ['circle', 'square', 'triangle', 'diamond', 'inverted_triangle', 'hex', 'star', 'plus']


def methylation_library_and_context(path):
    """'<sample_id>_<context>_combined_methylation.tsv' -> (library, context).

    sample_id is the methylKit basename, which already ends in _<context>, so the
    context suffix appears twice. Trailing lowercase dot-tags (.ds, .md) are removed
    so the name matches the BAM-derived one.
    """
    stem = re.sub(r'_combined_methylation\.tsv$', '', os.path.basename(path))
    m = re.search(r'_(CpG|CHG|CHH)$', stem)
    context = m.group(1) if m else 'CpG'
    library = re.sub(r'(\.[a-z]+)+$', '', CONTEXT_SUFFIX.sub('', stem))
    return library, context


def mosdepth_library(path, suffix):
    """'<id>.regions.bed.gz' -> library name, stripping trailing dot-tags as for BAMs."""
    stem = re.sub(re.escape(suffix) + '$', '', os.path.basename(path))
    return re.sub(r'(\.[a-z]+)+$', '', stem)


def coverage_table(regions_files, summary_files):
    """Mean depth per feature category relative to the genome mean, from mosdepth output."""
    summaries = {mosdepth_library(f, '.mosdepth.summary.txt'): f for f in summary_files}
    frames = []
    for f in regions_files:
        library = mosdepth_library(f, '.regions.bed.gz')
        if library not in summaries:
            sys.exit(f'{f}: no mosdepth summary for library {library}')
        regions = pl.read_csv(f, separator='\t', has_header=False,
                              new_columns=['chrom', 'start', 'end', 'feature', 'mean'],
                              schema_overrides={'chrom': pl.Utf8, 'feature': pl.Utf8})
        # Genome mean over the contigs that carry features, so spike-in controls and
        # unplaced scaffolds with no annotation do not shift it.
        contigs = set(regions['chrom'].unique().to_list())
        summary = pl.read_csv(summaries[library], separator='\t')
        genome = summary.filter(pl.col('chrom').is_in(list(contigs)))
        genome_mean = float(genome['bases'].sum()) / float(genome['length'].sum())
        frames.append(
            regions.with_columns((pl.col('end') - pl.col('start')).alias('len'))
            .group_by('feature')
            .agg(pl.len().alias('n_features'), pl.col('len').sum().alias('total_length'),
                 (pl.col('mean') * pl.col('len')).sum().alias('depth_bases'))
            .with_columns((pl.col('depth_bases') / pl.col('total_length')).alias('mean_depth'),
                          pl.lit(library).alias('library'), pl.lit(genome_mean).alias('genome_mean_depth'))
            .with_columns((pl.col('mean_depth') / pl.col('genome_mean_depth')).alias('relative_coverage'),
                          pl.col('feature').replace(FEATURE_ALIASES))
            .select(['library', 'feature', 'n_features', 'total_length', 'mean_depth',
                     'genome_mean_depth', 'relative_coverage'])
        )
    return pl.concat(frames)


def methylation_table(paths):
    frames = []
    for p in paths:
        library, context = methylation_library_and_context(p)
        df = pl.read_csv(p, separator='\t', schema_overrides={'Feature': pl.Utf8, 'Locus': pl.Utf8,
                                                               'Meth': pl.Float64})
        if df.is_empty():
            continue
        frames.append(
            df.group_by('Feature')
            .agg(pl.col('Meth').mean().alias('mean_methylation'),
                 pl.col('Meth').len().alias('n_methylation_features'))
            .with_columns(pl.lit(library).alias('library'), pl.lit(context).alias('context'),
                          pl.col('Feature').replace(FEATURE_ALIASES).alias('feature'))
            .drop('Feature')
        )
    if not frames:
        sys.exit('no methylation rows in any input file')
    return pl.concat(frames)


def ordered_libraries(libs):
    """(libraries in feature_depth_bokeh's order, index of the name part that groups them)."""
    group_idx, subgroup_idx = analyze_name_components(libs)
    return sorted(libs, key=lambda lib: library_sort_key(lib, group_idx, subgroup_idx)), group_idx


def build_context(df, libs, group_idx, features, ctx):
    color = {f: (SERIES[i] if i < len(SERIES) else OTHER) for i, f in enumerate(features)}
    marker = {f: MARKERS[i % len(MARKERS)] for i, f in enumerate(features)}
    df = df.with_columns(
        pl.col('feature').replace_strict(color, return_dtype=pl.Utf8).alias('color'),
        pl.col('feature').replace_strict(marker, return_dtype=pl.Utf8).alias('marker'),
    )
    tips = [('library', '@library'), ('feature', '@feature'),
            ('relative coverage', '@relative_coverage{0.000}'),
            ('mean methylation', '@mean_methylation{0.000}'),
            ('features counted', '@n_features{0,0}'),
            ('features with methylation', '@n_methylation_features{0,0}')]

    def scatter_by_feature(p, data, x, y):
        """One marker series per feature category, with a legend that hides categories."""
        items = []
        for f in features:
            # A category with no methylation rows has a null mean; Bokeh would draw it at 0.
            sub = data.filter((pl.col('feature') == f) & pl.col(x).is_not_null() & pl.col(y).is_not_null())
            if sub.is_empty():
                continue
            cols = sub.to_dict(as_series=False)
            # A grouped x axis takes (group, library) tuples; polars hands back lists.
            if 'x' in cols:
                cols['x'] = [tuple(v) if isinstance(v, list) else v for v in cols['x']]
            r = p.scatter(x, y, source=ColumnDataSource(cols), size=10, marker=marker[f],
                          fill_color=color[f], line_color='#ffffff', line_width=1.5, fill_alpha=0.9)
            items.append(LegendItem(label=f, renderers=[r]))
        p.add_tools(HoverTool(tooltips=tips))
        p.add_layout(Legend(items=items, click_policy='hide', label_text_font_size='9pt'), 'right')

    # 1. coverage vs methylation
    scatter = style(figure(width=1000, height=520, x_axis_type='log',
                           title=f'{ctx}: mean methylation vs relative coverage, '
                                 'one dot per feature category per library',
                           x_axis_label='relative coverage (1 = genome-wide average, log scale)',
                           y_axis_label='mean methylation', tools='pan,box_zoom,wheel_zoom,reset,save'))
    scatter.add_layout(Span(location=1.0, dimension='height', line_color=TEXT_2,
                            line_dash='dashed', line_width=1))
    scatter_by_feature(scatter, df, 'relative_coverage', 'mean_methylation')

    # 2/3. by library: group libraries on the x axis
    if group_idx is not None:
        group_of = {lib: get_component(lib, group_idx) for lib in libs}
        x_range = FactorRange(*[(group_of[lib], lib) for lib in libs])
        df = df.with_columns(pl.concat_list(
            pl.col('library').replace_strict(group_of, return_dtype=pl.Utf8), 'library').alias('x'))
    else:
        x_range = FactorRange(*libs)
        df = df.with_columns(pl.col('library').alias('x'))

    def by_library(y, label, title, log=False):
        p = style(figure(width=1000, height=440, x_range=x_range, y_axis_type='log' if log else 'linear',
                         title=title, y_axis_label=label, tools='pan,box_zoom,wheel_zoom,reset,save'))
        p.xaxis.major_label_orientation = 0.9
        p.xgrid.grid_line_color = None
        scatter_by_feature(p, df, 'x', y)
        return p

    cov_lib = by_library('relative_coverage', 'relative coverage (log scale)',
                         f'{ctx}: relative coverage per feature category', log=True)
    cov_lib.add_layout(Span(location=1.0, dimension='width', line_color=TEXT_2,
                            line_dash='dashed', line_width=1))
    meth_lib = by_library('mean_methylation', 'mean methylation',
                          f'{ctx}: mean methylation per feature category')
    return Tabs(tabs=[TabPanel(title='Coverage vs methylation', child=scatter),
                      TabPanel(title='Coverage by library', child=cov_lib),
                      TabPanel(title='Methylation by library', child=meth_lib)])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--regions', nargs='+', required=True,
                    help='<library>.regions.bed.gz from mosdepth --by feature_regions.bed')
    ap.add_argument('--summaries', nargs='+', required=True,
                    help='<library>.mosdepth.summary.txt from the same mosdepth runs')
    ap.add_argument('--methylation', nargs='+', required=True,
                    help='<sample>_<context>_combined_methylation.tsv files')
    ap.add_argument('-o', '--out', default='feature_cov_meth.html')
    ap.add_argument('--table', default='feature_cov_meth_summary.tsv')
    args = ap.parse_args()

    cov = coverage_table(args.regions, args.summaries)
    meth = methylation_table(args.methylation)

    unmatched = sorted(set(meth['library'].unique()) - set(cov['library'].unique()))
    if unmatched:
        print(f'[feature_cov_meth_bokeh] WARNING: methylation libraries with no mosdepth output: '
              f'{unmatched}; mosdepth libraries are {sorted(cov["library"].unique())}', file=sys.stderr)

    contexts = sorted(meth['context'].unique().to_list())
    table = (cov.join(meth.select('context').unique(), how='cross')
             .join(meth, on=['library', 'context', 'feature'], how='left')
             .sort(['context', 'library', 'feature']))
    table.write_csv(args.table, separator='\t')

    layout = build_from_table(table)
    header = Div(text=f'<h2 style="color:{TEXT};margin-bottom:0">Feature coverage and methylation</h2>')
    save_report(column(header, layout), args.out, 'Feature coverage and methylation', 'feature coverage plots')
    print(f'[feature_cov_meth_bokeh] wrote {args.out}: {table["library"].n_unique()} libraries, '
          f'{table["feature"].n_unique()} categories, contexts {contexts}', file=sys.stderr)


def build_from_table(table):
    """Build the views from feature_cov_meth_summary.tsv; depth_dependent_report.py embeds them too."""
    features = sorted(table['feature'].unique().to_list())
    if len(features) > len(SERIES):
        print(f'[feature_cov_meth_bokeh] {len(features)} feature categories; those past '
              f'{len(SERIES)} share the grey "other" colour and are told apart by marker',
              file=sys.stderr)
    libs, group_idx = ordered_libraries(table['library'].unique().to_list())

    panels = []
    for ctx in sorted(table['context'].drop_nulls().unique().to_list()):
        sub = table.filter(pl.col('context') == ctx)
        panels.append(TabPanel(title=ctx, child=build_context(sub, libs, group_idx, features, ctx)))

    intro = note(
        'Relative coverage is the mean per-base read depth over a category\u2019s features (mosdepth) '
        'divided by the genome-wide mean depth; 1 (dashed) is the genome average. Mean methylation averages the '
        'per-feature methylation of CpGs with depth \u2265 5. Click a legend entry to hide a category.')
    body = panels[0].child if len(panels) == 1 else Tabs(tabs=panels)
    return column(intro, body)


if __name__ == '__main__':
    main()
