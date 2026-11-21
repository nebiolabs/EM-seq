#!/usr/bin/env python
"""Build an interactive Bokeh report of per-C coverage and cross-library agreement.

Counts are CpG-context Cs, one per strand: each CpG holds two.

The report has these tabs:
  Depth        fraction of reference Cs at depth >= X for each library, and each
               library's depth check result when --depth-check is given.
               Always present.
  Upset        reference Cs per coverage pattern (which libraries cover them).
  Correlation  pairwise heatmaps of methylation r, coverage r and Jaccard, plus a
               pair explorer. The explorer shows the 2-D methylation histogram of
               shared Cs and the methylation of Cs covered by one library only.
  Clustering   average-linkage dendrograms over 1 - r and 1 - Jaccard.
  Library-specific Cs
               Cs covered in exactly one library, and how their methylation
               differs from all covered Cs in that library.

Only the Depth tab is built when no --merge-dir is given. That happens when there is
only one library.
"""
import argparse
import html
import math
import re
import os
import sys
from collections import Counter

import numpy as np
import polars as pl
from bokeh.layouts import column, row
from bokeh.models import (BasicTicker, BoxAnnotation, ColorBar, ColumnDataSource, CustomJS,
                          CustomJSHover, CustomJSTickFormatter, CheckboxGroup, HTMLTemplateFormatter,
                          DataTable, Div, FixedTicker, HoverTool, LabelSet, Legend, LegendItem,
                          LinearAxis, LinearColorMapper, LogColorMapper, NumberFormatter,
                          NumeralTickFormatter, RadioButtonGroup, Range1d, Select, Spacer, Span,
                          TableColumn, TabPanel, Tabs)
from bokeh.palettes import Greys256, Viridis256
from bokeh.plotting import figure

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from feature_depth_bokeh import RE_REPLICATE_LIKE, analyze_name_components, get_component  # noqa: E402
from report_overlay import save_report  # noqa: E402
from library_order import masses, natural_sorted  # noqa: E402
from depth_dependent_merge import DEPTH_BIN_WIDTH, DEPTH_BINS  # noqa: E402

# The All pairs depth axis stops at this percentile of covered depth, pooled over every
# library. Deeper Cs share the top bin, so a few very deep Cs cannot squeeze the rest into
# a corner of each cell.
DEPTH_AXIS_PERCENTILE = 99

# Categorical slots in fixed order (dataviz reference palette, light mode).
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#6250d6', '#e34948']
OTHER = '#9b9a95'
TEXT = '#0b0b0b'
TEXT_2 = '#52514e'
GRID = '#e6e5e0'
# Sequential blue ramp, light to dark.
BLUES = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7', '#3987e5',
         '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
DASHES = ['solid', 'dashed', 'dotted', 'dotdash', 'dashdot']
FIG_W = 1100


def label_px(labels):
    """Left border wide enough for 9pt tick labels.

    Plots stacked in a column align only when their left borders match, so each
    stack sets min_border_left to this one value. Bokeh's gridplot would do the
    alignment, but it lays out wrongly inside a tab that is hidden at load.
    """
    return max(90, 7 * max((len(x) for x in labels), default=0) + 24)


def style(p):
    p.outline_line_color = None
    p.xgrid.grid_line_color = GRID
    p.ygrid.grid_line_color = GRID
    p.axis.axis_line_color = TEXT_2
    p.axis.major_tick_line_color = TEXT_2
    p.axis.minor_tick_line_color = None
    p.axis.major_label_text_color = TEXT_2
    p.axis.axis_label_text_color = TEXT
    p.axis.axis_label_text_font_style = 'normal'
    p.title.text_color = TEXT
    p.toolbar.logo = None
    return p


def note(html):
    return Div(text=f'<div style="color:{TEXT_2};max-width:{FIG_W}px">{html}</div>')


def library_styles(libs):
    """Colour by library group (fixed slot order) and dash by position within the group."""
    group_idx, _ = analyze_name_components(libs)
    groups = natural_sorted({get_component(lib, group_idx) or lib for lib in libs})
    color = {g: (SERIES[i] if i < len(SERIES) else OTHER) for i, g in enumerate(groups)}
    seen = {}
    out = {}
    for lib in natural_sorted(libs):
        g = get_component(lib, group_idx) or lib
        k = seen.get(g, 0)
        seen[g] = k + 1
        out[lib] = {'group': g, 'color': color[g], 'dash': DASHES[k % len(DASHES)]}
    return out


# ---------------------------------------------------------------- overview tab

MASS_PART = re.compile(r'^\d+(?:\.\d+)?\s*(?:pg|ng|ug|µg)$', re.IGNORECASE)


def name_columns(libs):
    """Split library names on '_' into labelled columns.

    A position whose every value is a mass (0.1ng, 100pg) is the input amount; a last
    position whose every value looks like a replicate (1, rep2, B) is the replicate;
    the rest are numbered name parts.
    """
    parts = [lib.split('_') for lib in libs]
    width = max(len(p) for p in parts)
    padded = [p + [''] * (width - len(p)) for p in parts]
    cols = []
    for i in range(width):
        values = [p[i] for p in padded]
        present = [v for v in values if v]
        if present and all(MASS_PART.match(v) for v in present):
            name = 'input'
        elif i == width - 1 and width > 1 and present and all(RE_REPLICATE_LIKE.match(v) for v in present):
            name = 'replicate'
        else:
            name = f'name part {i + 1}'
        cols.append((name, values))
    return cols


def overview_tab(hist, read_stats, min_depth):
    libs = natural_sorted(hist['library'].unique().to_list())
    data = {'library': libs}
    columns = [TableColumn(field='library', title='library', width=200)]
    for name, values in name_columns(libs):
        field = name.replace(' ', '_')
        data[field] = values
        columns.append(TableColumn(field=field, title=name, width=90))

    if read_stats is not None:
        rs = {r['library']: r for r in read_stats.iter_rows(named=True)}
        get = lambda lib, k: (rs.get(lib) or {}).get(k)  # noqa: E731

        def ratio(a, b):
            return a / b if a is not None and b else float('nan')
        data['primary_reads'] = [get(lib, 'primary_reads') for lib in libs]
        data['pct_mapped'] = [ratio(get(lib, 'primary_mapped'), get(lib, 'primary_reads')) for lib in libs]
        data['pct_dups'] = [ratio(get(lib, 'primary_duplicates'), get(lib, 'primary_mapped')) for lib in libs]
        data['pct_usable'] = [ratio(get(lib, 'usable_reads'), get(lib, 'primary_reads')) for lib in libs]
        pct = NumberFormatter(format='0.0%')
        if 'flowcell' in read_stats.columns:
            data['flowcell'] = [get(lib, 'flowcell') or '' for lib in libs]
            columns.append(TableColumn(field='flowcell', title='flowcell', width=110))
        columns += [TableColumn(field='primary_reads', title='primary reads', width=110,
                                formatter=NumberFormatter(format='0,0')),
                    # Two decimals: aligners leave few reads unmapped, so 99.96% would show as 100.0%.
                    TableColumn(field='pct_mapped', title='% mapped', width=80,
                                formatter=NumberFormatter(format='0.00%')),
                    TableColumn(field='pct_dups', title='% duplicates', width=90, formatter=pct),
                    TableColumn(field='pct_usable', title='% usable', width=80, formatter=pct)]

    # Effective coverage: mean depth over every reference C (the per-C depth the other
    # tabs use), drawn as a labelled bar scaled to the deepest library.
    cov = []
    for lib in libs:
        h = hist.filter(pl.col('library') == lib)
        total = float(h['n_c'].sum())
        cov.append(float((h['depth'] * h['n_c']).sum()) / total if total else 0.0)
    top = max(cov) or 1.0
    data['cov'] = cov
    data['cov_pct'] = [100 * c / top for c in cov]
    data['cov_label'] = [f'{c:.1f}×' for c in cov]
    bar = HTMLTemplateFormatter(template=f"""
        <div style="position:relative;height:18px;margin-top:2px;background:{GRID};border-radius:3px">
          <div style="width:<%= cov_pct %>%;height:100%;background:{SERIES[0]};border-radius:3px"></div>
          <span style="position:absolute;left:6px;top:1px;font-size:11px;font-weight:600;color:{TEXT};
                       text-shadow:0 0 3px #fff,0 0 3px #fff"><%= cov_label %></span>
        </div>""")
    columns.append(TableColumn(field='cov', title='effective coverage (mean depth per C)',
                               width=260, formatter=bar))

    width = sum(c.width for c in columns) + 20
    table = DataTable(source=ColumnDataSource(data), columns=columns, width=max(width, FIG_W),
                      height=min(40 + 28 * len(libs), 800), index_position=None, row_height=26)
    notes = ['Library names are split on “_”: a mass (such as 0.1ng) is the input amount and '
             'a trailing number the replicate.']
    if read_stats is not None:
        q = int(read_stats['min_mapq'][0]) if 'min_mapq' in read_stats.columns else None
        notes.append('Primary reads, % mapped and % duplicates come from samtools flagstat '
                     '(duplicates as a share of mapped primary reads). % usable is the share of '
                     'primary reads MethylDackel uses with this pipeline’s settings: mapped, not '
                     'duplicates, QC failures or supplementary, '
                     + (f'MAPQ ≥ {q}, ' if q is not None else '')
                     + 'and properly paired when paired.')
    else:
        notes.append('Read counts need the BAMs: run with --bam_files (depth_dependent.nf) or from main.nf.')
    notes.append('Effective coverage is the strand-specific mean depth over every reference C')
    return TabPanel(title='Overview', child=column(note(' '.join(notes)), table))


# ---------------------------------------------------------------- depth tab

def depth_tab(hist, min_depth, depth_check, styles):
    libs = natural_sorted(hist['library'].unique().to_list())
    n_ref = hist.group_by('library').agg(pl.col('n_c').sum()).row(0)[1]
    max_bin = int(hist['depth'].max())

    # x range: the depth below which 99.5% of covered Cs fall in the deepest library
    x_max = 10
    stats = []
    curves = {}
    for lib in libs:
        h = hist.filter(pl.col('library') == lib).sort('depth')
        depth = h['depth'].to_numpy()
        n = h['n_c'].to_numpy().astype(np.float64)
        total = n.sum()
        dense = np.zeros(max_bin + 1)
        dense[depth] = n
        ge = np.cumsum(dense[::-1])[::-1] / total  # ge[x] = fraction with depth >= x
        curves[lib] = (dense / total, ge)
        covered = dense[1:].sum()
        if covered:
            cum = np.cumsum(dense[1:]) / covered
            x_max = max(x_max, int(np.searchsorted(cum, 0.995)) + 1)
        median = int(np.searchsorted(np.cumsum(dense) / total, 0.5))
        frac_min_depth = float(ge[min_depth]) if min_depth <= max_bin else 0.0
        stats.append({
            'library': lib,
            'mean_depth': float((dense * np.arange(max_bin + 1)).sum() / total),
            'median_depth': median,
            'frac_ge1': float(ge[1]) if max_bin >= 1 else 0.0,
            'frac_min_depth': frac_min_depth,
        })
        if depth_check is not None:
            g = depth_check.get(lib)
            stats[-1].update({
                'shallowest_contig': g['shallowest_contig'] if g else '',
                'shallowest_contig_depth': float(g['shallowest_mean_depth']) if g else float('nan'),
                'passes_depth_check': g['passes'] if g else 'not measured',
            })
    x_max = min(x_max, max_bin)
    xs = np.arange(1, x_max + 1)

    p_cum = style(figure(width=FIG_W, height=420, x_range=Range1d(0, x_max), y_range=Range1d(0, 1.02),
                         title=f'Fraction of {n_ref:,} reference Cs covered at depth ≥ X '
                               f'(controls and organelles excluded)',
                         tools='pan,box_zoom,wheel_zoom,reset,save',
                         x_axis_label='minimum depth X', y_axis_label='fraction of reference Cs'))
    p_dist = style(figure(width=FIG_W, height=300, x_range=p_cum.x_range,
                          title='Depth distribution (fraction of reference Cs at exactly depth X)',
                          tools='pan,box_zoom,wheel_zoom,reset,save',
                          x_axis_label='depth', y_axis_label='fraction of reference Cs'))
    items_cum = []
    for lib in libs:
        st = styles[lib]
        dist, ge = curves[lib]
        # Fractions and the matching C counts; the y-axis switch picks which one is drawn.
        src = ColumnDataSource(dict(x=xs, ge=ge[1:x_max + 1], at=dist[1:x_max + 1],
                                    ge_n=np.rint(ge[1:x_max + 1] * n_ref), at_n=np.rint(dist[1:x_max + 1] * n_ref),
                                    lib=[lib] * len(xs)))
        r1 = p_cum.line('x', 'ge', source=src, line_color=st['color'], line_dash=st['dash'], line_width=2)
        p_dist.line('x', 'at', source=src, line_color=st['color'], line_dash=st['dash'], line_width=2)
        items_cum.append(LegendItem(label=lib, renderers=[r1]))
    p_cum.add_tools(HoverTool(tooltips=[('library', '@lib'), ('depth ≥', '@x'),
                                        ('fraction', '@ge{0.000}'), ('Cs', '@ge_n{0,0}')], mode='mouse'))
    p_dist.add_tools(HoverTool(tooltips=[('library', '@lib'), ('depth', '@x'),
                                         ('fraction', '@at{0.0000}'), ('Cs', '@at_n{0,0}')]))
    p_cum.add_layout(Span(location=min_depth, dimension='height', line_color=TEXT_2,
                          line_dash='dashed', line_width=1))

    # y-axis switch: fraction of reference Cs, or the number of Cs.
    y_units = RadioButtonGroup(labels=['fraction of reference Cs', 'number of Cs'], active=0)
    # Millions as "M": NumeralTickFormatter's "a" format writes a lowercase "m".
    count_fmt = CustomJSTickFormatter(code="""
        const a = Math.abs(tick);
        if (a >= 1e6) return (tick / 1e6).toFixed(a % 1e6 ? 1 : 0) + 'M';
        if (a >= 1e3) return (tick / 1e3).toFixed(a % 1e3 ? 1 : 0) + 'k';
        return String(tick);
    """)
    y_units.js_on_change('active', CustomJS(
        args=dict(cum=[r for r in p_cum.renderers], dist=[r for r in p_dist.renderers],
                  cum_range=p_cum.y_range, cum_axis=p_cum.yaxis[0], dist_axis=p_dist.yaxis[0],
                  frac_fmt_cum=p_cum.yaxis[0].formatter, frac_fmt_dist=p_dist.yaxis[0].formatter,
                  count_fmt=count_fmt, n_ref=n_ref,
                  cum_title=p_cum.title, dist_title=p_dist.title),
        code="""
            const counts = cb_obj.active === 1;
            for (const r of cum) r.glyph.y = {field: counts ? 'ge_n' : 'ge'};
            for (const r of dist) r.glyph.y = {field: counts ? 'at_n' : 'at'};
            cum_range.end = counts ? n_ref * 1.02 : 1.02;
            cum_axis.formatter = counts ? count_fmt : frac_fmt_cum;
            dist_axis.formatter = counts ? count_fmt : frac_fmt_dist;
            cum_axis.axis_label = counts ? 'number of Cs' : 'fraction of reference Cs';
            dist_axis.axis_label = counts ? 'number of Cs' : 'fraction of reference Cs';
            const unit = counts ? 'Number of Cs' : 'Fraction of ' + n_ref.toLocaleString('en-US') + ' reference Cs';
            cum_title.text = unit + ' covered at depth \u2265 X (controls and organelles excluded)';
            dist_title.text = 'Depth distribution (' + (counts ? 'number of Cs' : 'fraction of reference Cs') +
                              ' at exactly depth X)';
        """))
    legend = Legend(items=items_cum, click_policy='hide', label_text_font_size='9pt')
    p_cum.add_layout(legend, 'right')
    p_dist.add_layout(Legend(items=[LegendItem(label=it.label, renderers=[r])
                                    for it, r in zip(items_cum, p_dist.renderers)],
                             click_policy='hide', label_text_font_size='9pt'), 'right')

    tsrc = ColumnDataSource(pl.DataFrame(stats).to_dict(as_series=False))
    pct = NumberFormatter(format='0.0%')
    table = DataTable(source=tsrc, width=FIG_W, height=min(40 + 28 * len(libs), 600), index_position=None,
                      columns=[TableColumn(field='library', title='library', width=320),
                               TableColumn(field='mean_depth', title='mean depth',
                                           formatter=NumberFormatter(format='0.0')),
                               TableColumn(field='median_depth', title='median depth'),
                               TableColumn(field='frac_ge1', title='covered ≥ 1', formatter=pct),
                               TableColumn(field='frac_min_depth', title=f'covered ≥ {min_depth}',
                                           formatter=pct)]
                      + ([TableColumn(field='shallowest_contig', title='shallowest checked contig'),
                          TableColumn(field='shallowest_contig_depth', title='its mean depth',
                                      formatter=NumberFormatter(format='0.0')),
                          TableColumn(field='passes_depth_check', title='passes depth check')]
                         if depth_check is not None else []))
    check_text = ''
    if depth_check is not None:
        n_pass = sum(s['passes_depth_check'] == 'yes' for s in stats)
        any_row = next(iter(depth_check.values()))
        check_text = (
            f'These analyses ran because enough libraries passed the depth check. A library passes '
            f'when each of its {any_row["n_checked_contigs"]} checked contigs has a mean depth above '
            f'{float(any_row["min_contig_depth"]):g}x. {n_pass} of {len(libs)} libraries pass. ')
    intro = note(
        'Every count is of CpG-context Cs. A CpG has one C on each strand, and each C counts '
        'separately. A C\'s depth is the number of reads on its own strand. '
        f'The {n_ref:,} reference Cs make up {n_ref // 2:,} CpGs. '
        + check_text +
        'The cross-library tabs include every library, even the shallow ones. The toggle at the '
        'top sets the minimum depth for a C to count as covered in those tabs. The dashed line '
        f'marks the depth the report opens at ({min_depth}x).')
    switch = row(Div(text=f'<b style="color:{TEXT}">y-axis:</b>'), y_units)
    return TabPanel(title='Depth', child=column(intro, switch, p_cum, p_dist, table)), stats


# ---------------------------------------------------------------- upset tab

def upset_panel(up, libs, n_ref, top_k=30):
    empty = up.filter(pl.col('n_libraries') == 0)['n_c'].sum()
    pats = up.filter(pl.col('n_libraries') > 0).sort('n_c', descending=True)
    shown = pats.head(top_k)
    hidden = pats.height - shown.height

    # Set sizes: Cs each library covers at this depth.
    size = {lib: 0 for lib in libs}
    for members, n in pats.select(['libraries', 'n_c']).iter_rows():
        for m in members.split(','):
            size[m] += n
    order = natural_sorted(libs)[::-1]  # y grows upward, so reverse to put the first library on top
    y_of = {lib: i for i, lib in enumerate(order)}
    K = shown.height
    xs = list(range(K))
    left = label_px(libs)
    frame_w = max(420, min(26 * K, FIG_W - 260 - left))
    sz_w = 240

    bar_src = ColumnDataSource(dict(
        x=xs, n=shown['n_c'].to_list(),
        frac=[n / n_ref for n in shown['n_c'].to_list()],
        libs=[s.replace(',', ', ') for s in shown['libraries'].to_list()],
        k=shown['n_libraries'].to_list()))
    # Log scale: the all-library pattern usually dwarfs the rest.
    bars = style(figure(frame_width=frame_w, height=300, x_range=Range1d(-0.6, K - 0.4),
                        y_axis_type='log', min_border_left=left, toolbar_location=None,
                        title='Cs per coverage pattern (largest first)',
                        tools='', y_axis_label='Cs (log scale)'))
    r = bars.vbar(x='x', top='n', bottom=1, width=0.7, source=bar_src, fill_color=BLUES[7],
                  line_color=None)
    bars.add_tools(HoverTool(renderers=[r], tooltips=[('libraries', '@libs'), ('# libraries', '@k'),
                                                      ('Cs', '@n{0,0}'), ('of reference', '@frac{0.0%}')]))
    bars.xaxis.visible = False
    bars.xgrid.grid_line_color = None
    bars.y_range.start = 1

    dot_h = max(160, 28 * len(libs) + 40)
    dots = style(figure(frame_width=frame_w, height=dot_h, x_range=bars.x_range,
                        y_range=Range1d(-0.6, len(libs) - 0.4), tools='',
                        min_border_left=left, toolbar_location=None))
    gx, gy = [], []
    on_x, on_y, seg_x0, seg_x1, seg_y0, seg_y1 = [], [], [], [], [], []
    for i, members in enumerate(shown['libraries'].to_list()):
        ys = [y_of[m] for m in members.split(',')]
        for lib in libs:
            gx.append(i)
            gy.append(y_of[lib])
        on_x += [i] * len(ys)
        on_y += ys
        seg_x0.append(i)
        seg_x1.append(i)
        seg_y0.append(min(ys))
        seg_y1.append(max(ys))
    dots.scatter(gx, gy, size=10, fill_color=GRID, line_color=None)
    dots.segment(seg_x0, seg_y0, seg_x1, seg_y1, line_color=TEXT, line_width=2)
    dots.scatter(on_x, on_y, size=10, fill_color=TEXT, line_color=None)
    dots.xaxis.visible = False
    dots.xgrid.grid_line_color = None
    dots.yaxis.ticker = FixedTicker(ticks=list(range(len(order))))
    dots.yaxis.major_label_overrides = {i: lib for i, lib in enumerate(order)}
    dots.yaxis.major_label_text_font_size = '9pt'

    max_size = max(size.values()) if size else 1
    sz = style(figure(width=sz_w, height=dot_h, y_range=dots.y_range, tools='',
                      x_range=Range1d(max_size * 1.05, 0), toolbar_location=None,
                      x_axis_label='Cs covered'))
    ssrc = ColumnDataSource(dict(y=list(range(len(order))), n=[size[lib] for lib in order],
                                 lib=order, frac=[size[lib] / n_ref for lib in order]))
    rs = sz.hbar(y='y', right='n', height=0.6, source=ssrc, fill_color=OTHER, line_color=None)
    sz.add_tools(HoverTool(renderers=[rs], tooltips=[('library', '@lib'), ('Cs', '@n{0,0}'),
                                                     ('of reference', '@frac{0.0%}')]))
    sz.xaxis.formatter = NumeralTickFormatter(format='0a')
    sz.yaxis.visible = False
    sz.ygrid.grid_line_color = None

    msg = (f'{empty:,} reference Cs ({empty / n_ref:.1%}) are covered in no library at this depth.')
    if hidden:
        msg += f' {hidden} smaller patterns are not drawn; they are in depth_dependent_upset.tsv.'
    return column(note(msg), row(Spacer(width=sz_w), bars), row(sz, dots))


def upset_view(up, libs, n_ref):
    intro = note('Each column is a coverage pattern: the libraries marked with dots cover those Cs '
                 'and the others do not. The bar above gives the number of reference Cs with that '
                 'exact pattern. The bars at the left give each library’s total.')
    return column(intro, upset_panel(up, libs, n_ref))


# ---------------------------------------------------------------- correlation tab

METRICS = [('meth_r', 'Methylation Pearson r (shared Cs)'),
           ('cov_log1p_r', 'Coverage Pearson r, log(1+depth), over all reference Cs'),
           ('jaccard', 'Jaccard index of covered Cs')]
# Coverage r is computed over every reference C, so it is the same at every minimum depth.
# Its heatmaps and clustering are built once and shown whatever depth is chosen.
DEPTH_INDEPENDENT = {'cov_log1p_r'}


def metric_name(title):
    return title.split(' (')[0].split(',')[0]


def square(pw, libs, metric):
    n = len(libs)
    m = np.eye(n)
    idx = {lib: i for i, lib in enumerate(libs)}
    for r in pw.iter_rows(named=True):
        a, b = idx[r['lib_a']], idx[r['lib_b']]
        m[a, b] = m[b, a] = r[metric]
    return m


def count_label(v):
    """Compact C count: 1.2M, 340k or 512."""
    return f'{v / 1e6:.1f}M' if v >= 1e6 else f'{v / 1e3:.0f}k' if v >= 1e3 else f'{v:,}'


def cell_px(n):
    return max(28, min(70, 640 // max(n, 1)))


def heatmap(pw, libs, metric, title, order=None):
    order = order or libs
    n = len(order)
    m = square(pw, libs, metric)
    pos = {lib: libs.index(lib) for lib in order}
    info = {}
    for r in pw.iter_rows(named=True):
        info[(r['lib_a'], r['lib_b'])] = info[(r['lib_b'], r['lib_a'])] = r
    xs, ys, vals, la, lb, ns, oa, ob = [], [], [], [], [], [], [], []
    for i, a in enumerate(order):
        for j, b in enumerate(order):
            xs.append(i)
            ys.append(n - 1 - j)
            vals.append(m[pos[a], pos[b]])
            la.append(a)
            lb.append(b)
            r = info.get((a, b))
            if r is None:
                ns.append(None), oa.append(None), ob.append(None)
            elif r['lib_a'] == a:
                ns.append(r['n_shared']), oa.append(r['n_only_a']), ob.append(r['n_only_b'])
            else:
                ns.append(r['n_shared']), oa.append(r['n_only_b']), ob.append(r['n_only_a'])
    off = [v for v, a, b in zip(vals, la, lb) if a != b and v == v]
    low = min(off) if off else 0.0
    mapper = LinearColorMapper(palette=Viridis256, low=low, high=1.0)
    cell = cell_px(n)
    p = style(figure(frame_width=cell * n, frame_height=cell * n, title=title, tools='save',
                     min_border_left=label_px(order), min_border_bottom=label_px(order),
                     x_range=Range1d(-0.5, n - 0.5), y_range=Range1d(-0.5, n - 0.5)))
    # The Cs behind each value, shown under it: shared Cs for methylation r, Cs covered in
    # either library for Jaccard. Coverage r uses every reference C, which the title states.
    def behind(n_sh, o_a, o_b):
        if n_sh is None:
            return ''
        if metric == 'meth_r':
            return '\n' + count_label(n_sh) + ' Cs'
        if metric == 'jaccard':
            return '\n' + count_label(n_sh + o_a + o_b) + ' Cs'
        return ''
    src = ColumnDataSource(dict(x=xs, y=ys, v=vals, a=la, b=lb, n=ns, oa=oa, ob=ob,
                                label=[f'{v:.3f}' + behind(n_sh, o_a, o_b) if a != b else ''
                                       for v, a, b, n_sh, o_a, o_b in zip(vals, la, lb, ns, oa, ob)]))
    r = p.rect('x', 'y', 0.94, 0.94, source=src, fill_color={'field': 'v', 'transform': mapper},
               line_color=None)
    if n <= 12:
        mid = (low + 1.0) / 2
        # Viridis runs dark to light, so low values get white text.
        src.data['tc'] = [('#ffffff' if v < mid else TEXT) for v in vals]
        p.text('x', 'y', text='label', source=src, text_align='center', text_baseline='middle',
               text_font_size='8pt', text_color='tc')
    tips = [('A', '@a'), ('B', '@b'), ('value', '@v{0.0000}')]
    # Shared and single-library counts depend on the minimum depth. A heatmap built once
    # for every depth leaves them out; the other heatmaps show them for the chosen depth.
    if metric not in DEPTH_INDEPENDENT:
        tips += [('shared Cs', '@n{0,0}'), ('only A', '@oa{0,0}'), ('only B', '@ob{0,0}')]
    p.add_tools(HoverTool(renderers=[r], tooltips=tips))
    p.xaxis.ticker = FixedTicker(ticks=list(range(n)))
    p.xaxis.major_label_overrides = {i: lib for i, lib in enumerate(order)}
    p.xaxis.major_label_orientation = math.pi / 3
    p.axis.major_label_text_font_size = '9pt'
    p.yaxis.ticker = FixedTicker(ticks=list(range(n)))
    p.yaxis.major_label_overrides = {n - 1 - j: lib for j, lib in enumerate(order)}
    p.grid.grid_line_color = None
    p.add_layout(ColorBar(color_mapper=mapper, width=12, ticker=BasicTicker(desired_num_ticks=6),
                          major_label_text_color=TEXT_2), 'right')
    return p


NOT_OBSERVED = -0.09   # axis position of the "not observed" strip in the pair explorer
NOT_OBSERVED_W = 0.06


def pair_views(pw, ph, ph2, libs, B):
    """Precompute every pair explorer view, for both orders of each library pair."""
    centers = [(i + 0.5) / B for i in range(B)]

    def pos(i):
        return NOT_OBSERVED if i < 0 else (i + 0.5) / B

    def size(i):
        return NOT_OBSERVED_W if i < 0 else 1.0 / B

    def label(i):
        return 'not observed' if i < 0 else f'{i / B:.2f}–{(i + 1) / B:.2f}'

    def fmt(v):
        return 'n/a' if v is None or (isinstance(v, float) and math.isnan(v)) else f'{v:.3f}'

    views = {}
    for r in pw.iter_rows(named=True):
        lo, hi = r['lib_a'], r['lib_b']
        h2 = ph2.filter((pl.col('lib_a') == lo) & (pl.col('lib_b') == hi))
        cells = list(zip(h2['bin_a'].to_list(), h2['bin_b'].to_list(), h2['n_c'].to_list()))
        for a, b, swap in ((lo, hi, False), (hi, lo, True)):
            ab = [(j, i, n) if swap else (i, j, n) for i, j, n in cells]
            rect = dict(na=[a] * len(ab), nb=[b] * len(ab),
                        x=[pos(i) for i, _, _ in ab], y=[pos(j) for _, j, _ in ab],
                        w=[size(i) for i, _, _ in ab], h=[size(j) for _, j, _ in ab],
                        n=[n for _, _, n in ab],
                        la=[label(i) for i, _, _ in ab], lb=[label(j) for _, j, _ in ab])

            def margin(axis):
                # Histogram of one library's methylation, split by whether the other library
                # covers the C. The not-observed bin holds Cs that only the other covers.
                both, only, other = [0] * B, [0] * B, 0
                for cell in ab:
                    mine, theirs, n = (cell[0], cell[1], cell[2]) if axis == 0 else (cell[1], cell[0], cell[2])
                    if mine < 0:
                        other += n
                    elif theirs < 0:
                        only[mine] += n
                    else:
                        both[mine] += n
                return dict(pos=centers + [NOT_OBSERVED], w=[1.0 / B] * B + [NOT_OBSERVED_W],
                            both=both + [0], only=only + [0], other=[0] * B + [other],
                            top=[b + o for b, o in zip(both, only)] + [other],
                            # The not-observed bar is often a sliver, so its count is printed.
                            note=[''] * B + [f'{other:,}'],
                            label=[label(i) for i in range(B)] + ['not observed'])

            k = (lambda x, y: r[y] if swap else r[x])
            stats = (
                '<table style="border-spacing:14px 2px">'
                f'<tr><td><b>shared Cs</b></td><td>{r["n_shared"]:,}</td>'
                f'<td><b>methylation r</b></td><td>{fmt(r["meth_r"])}</td></tr>'
                f'<tr><td><b>{a} only</b></td><td>{k("n_only_a", "n_only_b"):,}</td>'
                f'<td><b>mean methylation in {a}, shared / {a} only</b></td>'
                f'<td>{fmt(k("meth_mean_a_shared", "meth_mean_b_shared"))} / '
                f'{fmt(k("meth_mean_a_only", "meth_mean_b_only"))}</td></tr>'
                f'<tr><td><b>{b} only</b></td><td>{k("n_only_b", "n_only_a"):,}</td>'
                f'<td><b>mean methylation in {b}, shared / {b} only</b></td>'
                f'<td>{fmt(k("meth_mean_b_shared", "meth_mean_a_shared"))} / '
                f'{fmt(k("meth_mean_b_only", "meth_mean_a_only"))}</td></tr>'
                f'<tr><td><b>Jaccard</b></td><td>{fmt(r["jaccard"])}</td>'
                f'<td><b>coverage r, log(1+depth)</b></td><td>{fmt(r["cov_log1p_r"])}</td></tr></table>')
            views[f'{a}|{b}'] = {
                'rect': rect, 'below': margin(0), 'left': margin(1), 'stats': stats,
            }
    return views


def pair_explorer(pw_by_depth, ph2_by_depth, libs, pair_bins, depths, t0):
    """Methylation of library A against B, with each library's histogram along its axis.

    Both axes carry a "not observed" strip: Cs one library covers and the other does
    not. Their methylation in the covering library spreads along the strip, and the
    marginal histograms count them in a not-observed bar.

    One explorer serves every minimum depth. It returns the layout and a hidden depth
    selector, which the report's depth toggle sets.
    """
    views = pair_views(pw_by_depth[t0], None, ph2_by_depth[t0], libs, pair_bins)
    # Each marginal histogram scales to its own pair, so a shallow pair is not squashed by
    # the deepest one; the callback below resets the ranges when the selection changes.
    for v in views.values():
        for side in ('below', 'left'):
            m = v[side]
            v[side + '_max'] = max(max(b + o + t for b, o, t in zip(m['both'], m['only'], m['other'])), 1)
    a0, b0 = libs[0], libs[1]
    v0 = views[f'{a0}|{b0}']
    # The page embeds only each pair's (B+1) x (B+1) count matrix, index 0 meaning not
    # observed, and rebuilds the cells and histograms from it in the browser. Embedding the
    # built views for every ordered pair grew as libraries squared times depths (21 MB for
    # 8 libraries on one chromosome).
    W = pair_bins + 1
    keys = [str(t) for t in depths]
    mats, stats, max_n = {}, {}, {}
    for t, key in zip(depths, keys):
        mats[key] = {}
        for la, lb, ia, ib, cnt in ph2_by_depth[t].select(['lib_a', 'lib_b', 'bin_a', 'bin_b',
                                                           'n_c']).iter_rows():
            m = mats[key].setdefault(f'{la}|{lb}', [0] * (W * W))
            m[(ia + 1) * W + (ib + 1)] = cnt
        depth_views = views if t == t0 else pair_views(pw_by_depth[t], None, ph2_by_depth[t], libs,
                                                         pair_bins)
        stats[key] = {k: v['stats'] for k, v in depth_views.items()}
        # Each depth keeps its own colour scale, as when every depth had its own explorer.
        max_n[key] = max(max((max(m) for m in mats[key].values()), default=10), 10)
    depth_sel = RadioButtonGroup(labels=keys, active=depths.index(t0), visible=False)
    sel_a = Select(title='x-axis library', value=a0, options=libs, width=360)
    sel_b = Select(title='y-axis library', value=b0, options=libs, width=360)
    stats_div = Div(width=FIG_W, text=v0['stats'])

    frame = 440
    left_border = 80
    side_w = 210
    axis_range = (NOT_OBSERVED - NOT_OBSERVED_W / 2 - 0.01, 1.01)
    ticks = [NOT_OBSERVED, 0, 0.2, 0.4, 0.6, 0.8, 1.0]
    tick_labels = {NOT_OBSERVED: 'not obs.'}

    mapper = LogColorMapper(palette=Viridis256, low=1, high=max_n[str(t0)])
    src_rect = ColumnDataSource(v0['rect'])
    center = style(figure(frame_width=frame, frame_height=frame, min_border_left=left_border,
                          min_border_top=10, x_range=Range1d(*axis_range), y_range=Range1d(*axis_range),
                          tools='save,box_zoom,reset', x_axis_label=f'methylation, {a0}',
                          y_axis_label=f'methylation, {b0}'))
    strip_edge = NOT_OBSERVED + NOT_OBSERVED_W / 2 + 0.01
    center.add_layout(BoxAnnotation(right=strip_edge, fill_color=GRID, fill_alpha=0.6, level='underlay'))
    center.add_layout(BoxAnnotation(top=strip_edge, fill_color=GRID, fill_alpha=0.6, level='underlay'))
    r = center.rect('x', 'y', 'w', 'h', source=src_rect, line_color=None,
                    fill_color={'field': 'n', 'transform': mapper})
    center.add_tools(HoverTool(renderers=[r], tooltips=[('x', '@na: @la'), ('y', '@nb: @lb'),
                                                        ('Cs', '@n{0,0}')]))
    center.add_layout(ColorBar(color_mapper=mapper, width=12, title='Cs per cell',
                               major_label_text_color=TEXT_2), 'right')
    for ax in (center.xaxis[0], center.yaxis[0]):
        ax.ticker = FixedTicker(ticks=ticks)
        ax.major_label_overrides = tick_labels
    center.grid.grid_line_color = None

    stack = ['both', 'only', 'other']
    colors = [SERIES[0], SERIES[1], OTHER]
    tips = [('methylation', '@label'), ('covered in both', '@both{0,0}'),
            ('this library only', '@only{0,0}'), ('not observed here', '@other{0,0}')]

    src_below = ColumnDataSource(v0['below'])
    below = style(figure(frame_width=frame, frame_height=150, min_border_left=left_border,
                         x_range=center.x_range, toolbar_location=None, tools='',
                         y_axis_label=f'Cs, {a0}'))
    rb = below.vbar_stack(stack, x='pos', width='w', source=src_below, color=colors, line_color=None)
    # Left-aligned with the bar's left edge: the not-observed bar sits at the frame edge, so a
    # centred count would be clipped.
    below.text('pos', 'top', text='note', source=src_below, y_offset=4, x_offset=-11,
               text_align='left', text_baseline='top', text_font_size='8pt', text_color=TEXT_2)
    below.add_tools(HoverTool(renderers=rb, tooltips=tips))
    below.xaxis.visible = False
    below.xgrid.grid_line_color = None
    below.yaxis.formatter = NumeralTickFormatter(format='0.[0]a')

    src_left = ColumnDataSource(v0['left'])
    left = style(figure(width=side_w, frame_height=frame, min_border_top=10, y_range=center.y_range,
                        toolbar_location=None, tools='', x_axis_label=f'Cs, {b0}'))
    rl = left.hbar_stack(stack, y='pos', height='w', source=src_left, color=colors, line_color=None)
    left.text('top', 'pos', text='note', source=src_left, text_align='right', x_offset=-4,
              text_baseline='middle', text_font_size='8pt', text_color=TEXT_2)
    left.add_tools(HoverTool(renderers=rl, tooltips=tips))
    left.add_layout(Legend(items=[LegendItem(label='covered in both', renderers=[rl[0]]),
                                  LegendItem(label='this library only', renderers=[rl[1]]),
                                  LegendItem(label='not observed here', renderers=[rl[2]])],
                           location='center_left', label_text_font_size='8pt'))
    left.x_range = Range1d(v0['left_max'] * 1.45, 0)
    left.yaxis.visible = False
    left.ygrid.grid_line_color = None
    left.xaxis.formatter = NumeralTickFormatter(format='0.[0]a')
    left.xaxis.major_label_orientation = 0.8
    # Flipped, so the bars hang down from the 2-D plot as the left ones extend out from it.
    below.y_range = Range1d(v0['below_max'] * 1.3, 0)

    title = Div(text='')

    def title_html(a, b):
        return (f'<b style="color:{TEXT}">Methylation in {a} (x) vs {b} (y)</b>'
                f'<span style="color:{TEXT_2}"> — the grey strips hold Cs one library '
                'covers and the other does not</span>')
    title.text = title_html(a0, b0)
    cb = CustomJS(args=dict(sa=sel_a, sb=sel_b, depth_sel=depth_sel, keys=keys, mats_by_depth=mats,
                            stats_by_depth=stats, max_n=max_n, mapper=mapper, B=pair_bins,
                            NO=NOT_OBSERVED, NOW=NOT_OBSERVED_W, rect=src_rect, below=src_below,
                            left=src_left, div=stats_div, title=title,
                            cx=center.xaxis[0], cy=center.yaxis[0], by=below.yaxis[0], lx=left.xaxis[0],
                            lr=left.x_range, br=below.y_range,
                            tmpl=title_html('__A__', '__B__')), code="""
        const key = keys[depth_sel.active];
        const mats = mats_by_depth[key], stats = stats_by_depth[key];
        mapper.high = max_n[key];
        const a = sa.value, b = sb.value;
        if (a === b) { div.text = '<i>Pick two different libraries.</i>'; return; }
        let M = mats[a + '|' + b], swap = false;
        if (!M) { M = mats[b + '|' + a]; swap = true; }
        const W = B + 1;
        // val(i, j): Cs with x-library bin i-1 and y-library bin j-1 (bin -1 = not observed)
        const val = (i, j) => M ? (swap ? M[j * W + i] : M[i * W + j]) : 0;
        const pos = (k) => k < 0 ? NO : (k + 0.5) / B;
        const size = (k) => k < 0 ? NOW : 1.0 / B;
        const lab = (k) => k < 0 ? 'not observed' : (k / B).toFixed(2) + '\u2013' + ((k + 1) / B).toFixed(2);
        const r = {na: [], nb: [], x: [], y: [], w: [], h: [], n: [], la: [], lb: []};
        for (let i = 0; i < W; i++) for (let j = 0; j < W; j++) {
            const n = val(i, j);
            if ((i === 0 && j === 0) || !n) continue;
            r.na.push(a); r.nb.push(b); r.n.push(n);
            r.x.push(pos(i - 1)); r.y.push(pos(j - 1)); r.w.push(size(i - 1)); r.h.push(size(j - 1));
            r.la.push(lab(i - 1)); r.lb.push(lab(j - 1));
        }
        // Each library's histogram, split by whether the other library covers the C; the
        // not-observed bar counts Cs only the other library covers.
        const margin = (mine) => {
            const both = [], only = [], other = [], top = [], note = [], label = [], p = [], w = [];
            let missed = 0;
            for (let k = 0; k < B; k++) {
                let s = 0;
                for (let o = 1; o < W; o++) s += mine(k + 1, o);
                const x = mine(k + 1, 0);
                both.push(s); only.push(x); other.push(0); top.push(s + x); note.push('');
                label.push(lab(k)); p.push(pos(k)); w.push(size(k));
            }
            for (let o = 1; o < W; o++) missed += mine(0, o);
            both.push(0); only.push(0); other.push(missed); top.push(missed);
            note.push(missed.toLocaleString('en-US')); label.push('not observed'); p.push(NO); w.push(NOW);
            return {pos: p, w: w, both: both, only: only, other: other, top: top, note: note, label: label};
        };
        const mb = margin((i, j) => val(i, j)), ml = margin((j, i) => val(i, j));
        const peak = (m) => Math.max(1, ...m.top);
        rect.data = r; below.data = mb; left.data = ml;
        lr.start = peak(ml) * 1.45; br.start = peak(mb) * 1.3;
        div.text = stats[a + '|' + b];
        title.text = tmpl.replace('__A__', a).replace('__B__', b);
        cx.axis_label = 'methylation, ' + a; cy.axis_label = 'methylation, ' + b;
        by.axis_label = 'Cs, ' + a; lx.axis_label = 'Cs, ' + b;
    """)
    sel_a.js_on_change('value', cb)
    sel_b.js_on_change('value', cb)
    depth_sel.js_on_change('active', cb)

    return column(row(sel_a, sel_b), depth_sel, stats_div, title,
                  row(left, center), row(Spacer(width=side_w), below)), depth_sel


def correlation_tab(pw_by_depth, h2_by_depth, libs, depths, t0, pair_bins):
    """Correlation tab, the layouts to show at each depth, and the explorer's depth selector.

    Each depth gets its own intro and depth-dependent heatmaps. The depth-independent
    heatmaps and the pair explorer are built once.
    """
    per_depth = []
    intros = []
    heat = {m: [] for m, _ in METRICS if m not in DEPTH_INDEPENDENT}
    for t in depths:
        pw = pw_by_depth[t]
        shown = [correlation_intro(t)]
        intros.append(shown[0])
        for m, title in METRICS:
            if m in heat:
                heat[m].append(heatmap(pw, libs, m, title))
                shown.append(heat[m][-1])
        per_depth.append(shown)
    explorer, explorer_depth = pair_explorer(pw_by_depth, h2_by_depth, libs, pair_bins, depths, t0)
    heat_tabs = Tabs(tabs=[
        TabPanel(title=metric_name(title),
                 child=column(*heat[m]) if m in heat else heatmap(pw_by_depth[depths[0]], libs, m, title))
        for m, title in METRICS])
    panel = TabPanel(title='Correlation', child=column(column(*intros), heat_tabs,
                                                      Div(text='<h3>Pair explorer</h3>'),
                                                      explorer))
    return panel, per_depth, explorer_depth


def correlation_intro(min_depth):
    return note(
        f'A C is considerd covered in a library if the number of strand-specific overlapping reads reaches {min_depth}. '
        'Methylation r uses Cs covered in both libraries. Coverage r uses every reference C, '
        'uncovered ones included, so it measures how alike the coverage profiles are. '
        'Jaccard is shared Cs divided by Cs covered in either library. '
        'In the pair explorer, each library\u2019s histogram runs along its own axis. Orange '
        'counts Cs that library covers and the other does not; the grey bar counts Cs it misses '
        'that the other covers. If the orange part has a different methylation profile from the '
        'blue (shared) part, coverage dropout is not random with respect to methylation.')


# ---------------------------------------------------------------- all-pairs tab

def all_pairs_tab(pw, dh2, by_depth, order, min_depth, upset=None, n_ref=0, core=None, depth_top=None):
    """One figure with a cell per library pair, with a minimum-depth toggle.

    Below the diagonal: 2-D depth histogram over every reference C, in greyscale, with
    Pearson r of log(1 + depth). On the diagonal: the library's methylation histogram.
    Above the diagonal: 2-D methylation histogram of shared Cs, with Pearson r. The
    diagonal bars and the methylation cells share one viridis C-count scale. The toggle
    swaps the methylation half between minimum depths; the depth half covers every
    reference C and does not change.

    by_depth maps each minimum depth to (pair 2-D histogram, per-library histogram,
    pairwise r) tables. In every cell x is the column library and y the row library. Each
    cell maps its values into [pad, 1 - pad] of a unit square, so all cells share one
    coordinate system; the outer axes carry the within-cell scales (depth bottom and left,
    methylation top and right).

    depth_top is the depth where the depth axis stops; deeper Cs join the top bin. Without
    it the axis runs over the merge's full grid.
    """
    n = len(order)
    pad = 0.06
    span = 1 - 2 * pad
    pos = {lib: i for i, lib in enumerate(order)}
    depths = sorted(by_depth)
    B = int(by_depth[depths[0]][1].group_by('library').len()['len'].max())
    # The merge writes depth_hist2d.tsv on a fixed grid of DEPTH_BINS bins, each
    # DEPTH_BIN_WIDTH wide. The axis keeps the bins up to depth_top and folds every deeper
    # bin into the last one, which is open-ended. At least 3 bins keep the 0, middle and
    # top ticks apart.
    depth_w = DEPTH_BIN_WIDTH
    D = DEPTH_BINS if depth_top is None else max(3, min(DEPTH_BINS, depth_top // depth_w + 1))
    top_lo = (D - 1) * depth_w
    dh2 = (dh2.with_columns(pl.col('depth_a').clip(upper_bound=top_lo),
                            pl.col('depth_b').clip(upper_bound=top_lo))
              .group_by(['lib_a', 'lib_b', 'depth_a', 'depth_b']).agg(pl.col('n_c').sum())
              .sort(['lib_a', 'lib_b', 'depth_a', 'depth_b']))
    r_cov = {}
    for r in pw.iter_rows(named=True):
        r_cov[(r['lib_a'], r['lib_b'])] = r_cov[(r['lib_b'], r['lib_a'])] = r['cov_log1p_r']

    def origin(row_lib, col_lib):
        return pos[col_lib], n - 1 - pos[row_lib]

    def cells(df, a_col, b_col, lower, nbins, scale, coords=True):
        # Library names and bins go in as small integer codes that the hover formatters turn
        # back into text, and every column as a numpy array, which Bokeh embeds as binary.
        # As JSON text and strings these cells made up most of the page. The render callback
        # places cells from the codes, so x and y are only needed for the initial view
        # (coords=True). Every cell in a source has the same width, which the glyph holds.
        out = dict(x=[], y=[], n=[], lx=[], ly=[], vx=[], vy=[])
        for la, lb, ia, ib, cnt in df.select(['lib_a', 'lib_b', a_col, b_col, 'n_c']).iter_rows():
            ia, ib = min(ia // scale, nbins - 1), min(ib // scale, nbins - 1)
            if ia < 0 or ib < 0 or la not in pos or lb not in pos:
                continue
            # The column library goes on x. In the lower triangle the row sits below the
            # diagonal, so the column is the library earlier in the order.
            col_is_a = (pos[la] < pos[lb]) == lower
            col_lib, row_lib = (la, lb) if col_is_a else (lb, la)
            bx, by = (ia, ib) if col_is_a else (ib, ia)
            x0, y0 = origin(row_lib, col_lib)
            out['x'].append(x0 + pad + (bx + 0.5) / nbins * span)
            out['y'].append(y0 + pad + (by + 0.5) / nbins * span)
            out['n'].append(cnt)
            out['lx'].append(pos[col_lib])
            out['ly'].append(pos[row_lib])
            out['vx'].append(bx)
            out['vy'].append(by)
        types = dict(x=np.float32, y=np.float32, n=np.int32, lx=np.int16,
                     ly=np.int16, vx=np.int16, vy=np.int16)
        return {k: np.asarray(v, dtype=types[k]) for k, v in out.items()
                if coords or k not in ('x', 'y')}

    def diag_data(h1):
        out = dict(left=[], right=[], bottom=[], top=[], lib=[], bin=[], n=[], code=[])
        for lib in order:
            h = h1.filter(pl.col('library') == lib).sort('bin_lo')
            counts = h['n_c'].to_list()
            peak = max(counts) if counts and max(counts) else 1
            x0, y0 = origin(lib, lib)
            for lo, hi, c in zip(h['bin_lo'].to_list(), h['bin_hi'].to_list(), counts):
                out['left'].append(x0 + pad + lo * span)
                out['right'].append(x0 + pad + hi * span)
                out['bottom'].append(y0 + pad)
                out['top'].append(y0 + pad + c / peak * span)
                out['lib'].append(lib)
                out['bin'].append(f'{lo:.2f}–{hi:.2f}')
                out['n'].append(c)
                out['code'].append(pos[lib])
        return out

    cell = max(90, min(190, 1000 // max(n, 1)))
    box_h = 14 / cell

    def count(v):
        return f'{v / 1e6:.1f}M' if v >= 1e6 else f'{v / 1e3:.0f}k' if v >= 1e3 else f'{v:,}'

    def label_data(texts, centred=False):
        """Label sources for one cell per text; centred labels sit above the diagonal bars."""
        # rc/cc are the row and column codes, al 1 for centred labels; the library
        # checkboxes use them to move labels when the matrix is rebuilt.
        lab = dict(x=[], y=[], text=[], bx=[], by=[], bw=[], bh=[], rc=[], cc=[], al=[])
        for (row_lib, col_lib), text in texts.items():
            x0, y0 = origin(row_lib, col_lib)
            w = (len(text) * 5.4 + 8) / cell
            x = x0 + 0.5 if centred else x0 + pad
            lab['x'].append(x)
            lab['y'].append(y0 + 1 - pad)
            lab['text'].append(text)
            lab['bx'].append(x if centred else x + w / 2)
            lab['by'].append(y0 + 1 - pad - box_h / 2)
            lab['bw'].append(w)
            lab['bh'].append(box_h)
            lab['rc'].append(pos[row_lib])
            lab['cc'].append(pos[col_lib])
            lab['al'].append(1 if centred else 0)
        return lab

    t0 = str(min_depth if min_depth in by_depth else depths[-1])
    views = {}
    for t in depths:
        h2, h1, rt = by_depth[t]
        upper_text = {}
        for la, lb, n_sh, v in rt.select(['lib_a', 'lib_b', 'n_shared', 'meth_r']).iter_rows():
            if la not in pos or lb not in pos or v is None or v != v:
                continue
            row_lib, col_lib = (la, lb) if pos[la] < pos[lb] else (lb, la)
            upper_text[(row_lib, col_lib)] = f'r = {v:.3f} \u00b7 {count(n_sh)} Cs'
        diag_text = {(lib, lib): f'{count(c)} Cs at \u2265 {t}x'
                     for lib, c in h1.group_by('library').agg(pl.col('n_c').sum()).iter_rows()
                     if lib in pos}
        views[str(t)] = {'upper': cells(h2, 'bin_a', 'bin_b', False, B, 1, coords=str(t) == t0),
                         'diag': diag_data(h1), 'labels': label_data(upper_text),
                         'diag_labels': label_data(diag_text, centred=True)}
    # Core-set views: every statistic on the Cs that all core libraries cover at that depth.
    core_views = {}
    core_codes = []
    if core is not None:
        core_codes = [pos[x] for x in core['members'] if x in pos]
        for t in depths:
            ch2, cmh, cpw = core['by_depth'].get(t, (None, None, None))
            if ch2 is None:
                continue
            n_common = int(cpw['n_common'][0]) if cpw.height else 0
            upper_text = {}
            for la, lb, v in cpw.select(['lib_a', 'lib_b', 'meth_r']).iter_rows():
                if la in pos and lb in pos and v is not None and v == v:
                    row_lib, col_lib = (la, lb) if pos[la] < pos[lb] else (lb, la)
                    upper_text[(row_lib, col_lib)] = f'r = {v:.3f} \u00b7 {count(n_common)} Cs'
            diag_text = {(lib, lib): f'{count(n_common)} common Cs at \u2265 {t}x'
                         for lib in core['members'] if lib in pos}
            core_views[str(t)] = {'upper': cells(ch2, 'bin_a', 'bin_b', False, B, 1, coords=False),
                                  'diag': diag_data(cmh), 'labels': label_data(upper_text),
                                  'diag_labels': label_data(diag_text, centred=True),
                                  'n_common': n_common}
    lower = cells(dh2, 'depth_a', 'depth_b', True, D, depth_w)

    border = label_px(order)
    p = style(figure(frame_width=cell * n, frame_height=cell * n, x_range=Range1d(0, n),
                     y_range=Range1d(0, n), min_border_left=border, min_border_bottom=border,
                     tools='save,box_zoom,wheel_zoom,pan,reset'))
    p.grid.grid_line_color = None
    src_grid = ColumnDataSource(dict(
        left=[j for j in range(n) for _ in range(n)], right=[j + 1 for j in range(n) for _ in range(n)],
        bottom=[i for _ in range(n) for i in range(n)], top=[i + 1 for _ in range(n) for i in range(n)]))
    p.quad(left='left', right='right', bottom='bottom', top='top', source=src_grid,
           fill_color=None, line_color=GRID, line_width=1)

    # Depth counts in greys (light to dark, starting above white so single Cs show).
    greys = list(reversed(Greys256))[48:]
    depth_mapper = LogColorMapper(palette=greys, low=1, high=max(int(lower['n'].max(initial=10)), 10))
    src_low = ColumnDataSource(lower)
    r_low = p.rect('x', 'y', span / D, span / D, source=src_low, line_color=None,
                   fill_color={'field': 'n', 'transform': depth_mapper})
    lib_name = CustomJSHover(args=dict(names=order), code='return names[value];')
    depth_bin = CustomJSHover(args=dict(w=depth_w, last=D - 1), code="""
        const lo = value * w;
        return value === last ? '\u2265 ' + lo + 'x' : lo + '\u2013' + (lo + w - 1) + 'x';
    """)
    meth_bin = CustomJSHover(args=dict(B=B), code="""
        return (value / B).toFixed(2) + '\u2013' + ((value + 1) / B).toFixed(2);
    """)
    p.add_tools(HoverTool(renderers=[r_low], tooltips=[
        ('x library', '@lx{lib}'), ('y library', '@ly{lib}'), ('x depth', '@vx{bin}'),
        ('y depth', '@vy{bin}'), ('reference Cs', '@n{0,0}')],
        formatters={'@lx': lib_name, '@ly': lib_name, '@vx': depth_bin, '@vy': depth_bin}))

    # One viridis C-count scale for the methylation cells and the diagonal bars, fixed
    # across depths so colours stay comparable when the toggle changes.
    meth_max = max(max(int(v['upper']['n'].max(initial=0)), max(v['diag']['n'], default=0))
                   for v in list(views.values()) + list(core_views.values()))
    meth_mapper = LogColorMapper(palette=Viridis256, low=1, high=max(meth_max, 10))
    src_up = ColumnDataSource(views[t0]['upper'])
    r_up = p.rect('x', 'y', span / B, span / B, source=src_up, line_color=None,
                  fill_color={'field': 'n', 'transform': meth_mapper})
    p.add_tools(HoverTool(renderers=[r_up], tooltips=[
        ('x library', '@lx{lib}'), ('y library', '@ly{lib}'), ('x methylation', '@vx{bin}'),
        ('y methylation', '@vy{bin}'), ('shared Cs', '@n{0,0}')],
        formatters={'@lx': lib_name, '@ly': lib_name, '@vx': meth_bin, '@vy': meth_bin}))
    src_diag = ColumnDataSource(views[t0]['diag'])
    r_diag = p.quad(left='left', right='right', bottom='bottom', top='top', source=src_diag,
                    fill_color={'field': 'n', 'transform': meth_mapper}, line_color='#ffffff',
                    line_width=0.5)
    p.add_tools(HoverTool(renderers=[r_diag], tooltips=[('library', '@lib'), ('methylation', '@bin'),
                                                        ('Cs', '@n{0,0}')]))

    # r labels on white boxes: LabelSet's own background does not cover the full text,
    # and labels otherwise sit on dark cells.
    def add_labels(src, align='left'):
        p.rect('bx', 'by', 'bw', 'bh', source=src, fill_color='#ffffff', fill_alpha=0.85,
               line_color=None)
        p.add_layout(LabelSet(x='x', y='y', text='text', source=src, text_align=align,
                              text_baseline='top', text_font_size='8pt', text_color=TEXT,
                              x_offset=3 if align == 'left' else 0, y_offset=-1))
    lower_text = {}
    for (la, lb), v in r_cov.items():
        if la in pos and lb in pos and pos[la] > pos[lb] and v == v:
            lower_text[(la, lb)] = f'r = {v:.3f}'
    src_llab = ColumnDataSource(label_data(lower_text))
    add_labels(src_llab)
    src_lab = ColumnDataSource(views[t0]['labels'])
    add_labels(src_lab)
    src_dlab = ColumnDataSource(views[t0]['diag_labels'])
    add_labels(src_dlab, align='center')

    # Within-cell scales on the outer edges. Depth cells fill the lower left, so their
    # scale sits along the bottom (columns 0..n-2) and the left (rows 1..n-1). Methylation
    # cells and the diagonal fill the upper right, so methylation sits along the top (every
    # column, the diagonal included) and the right (rows 0..n-2).
    max_depth = (D - 1) * depth_w
    depth_ticks = {0: '0', max_depth // 2: f'{max_depth // 2}', max_depth: f'≥{max_depth}'}
    meth_ticks = {0: '0', 0.5: '0.5', 1: '1'}

    def depth_pos(c, v):
        # depth bins are centred at (k + 0.5) / D; place the tick at its bin centre
        k = min(int(v // depth_w), D - 1)
        return c + pad + (k + 0.5) / D * span

    def meth_pos(c, v):
        return c + pad + v * span

    def scale_ticks(ticks, cells_at, to_cell):
        return {round(to_cell(c, v), 6): label for c in cells_at for v, label in ticks.items()}

    def styled_axis(axis, ticks, label):
        axis.ticker = FixedTicker(ticks=sorted(ticks))
        axis.major_label_overrides = ticks
        axis.major_label_text_font_size = '7pt'
        axis.major_label_text_color = TEXT_2
        axis.axis_label = label
        axis.axis_label_text_font_size = '9pt'
        axis.axis_label_text_font_style = 'normal'
        axis.minor_tick_line_color = None
        return axis

    # Keep the depth axes now. Once the top and right axes are added, p.xaxis[0] and
    # p.yaxis[0] can name those instead.
    ax_bottom = styled_axis(p.xaxis[0], scale_ticks(depth_ticks, range(n - 1), depth_pos),
                            'depth (lower-left cells)')
    ax_left = styled_axis(p.yaxis[0], scale_ticks(depth_ticks, [n - 1 - i for i in range(1, n)], depth_pos),
                      'depth (lower-left cells)')
    ax_top = styled_axis(LinearAxis(), scale_ticks(meth_ticks, range(n), meth_pos),
                         'methylation (upper-right cells and diagonal)')
    ax_right = styled_axis(LinearAxis(), scale_ticks(meth_ticks, [n - 1 - i for i in range(n - 1)], meth_pos),
                           'methylation (upper-right cells)')
    p.add_layout(ax_top, 'above')
    p.add_layout(ax_right, 'right')
    names_x = LinearAxis(ticker=FixedTicker(ticks=[j + 0.5 for j in range(n)]),
                            major_label_overrides={j + 0.5: lib for j, lib in enumerate(order)},
                            major_label_orientation=math.pi / 3, major_label_text_font_size='9pt',
                            major_tick_line_color=None, axis_line_color=None,
                            axis_label='column library (x within each cell)',
                            axis_label_text_font_style='normal')
    p.add_layout(names_x, 'below')
    names_y = LinearAxis(ticker=FixedTicker(ticks=[n - 1 - i + 0.5 for i in range(n)]),
                            major_label_overrides={n - 1 - i + 0.5: lib for i, lib in enumerate(order)},
                            major_label_text_font_size='9pt', major_tick_line_color=None,
                            axis_line_color=None, axis_label='row library (y within each cell)',
                            axis_label_text_font_style='normal')
    p.add_layout(names_y, 'left')
    p.add_layout(ColorBar(color_mapper=depth_mapper, width=12, title='Cs per cell, depth (lower left)',
                          title_text_font_size='9pt', major_label_text_color=TEXT_2), 'right')
    p.add_layout(ColorBar(color_mapper=meth_mapper, width=12,
                          title='Cs per cell or bar, methylation (upper right and diagonal)',
                          title_text_font_size='9pt', major_label_text_color=TEXT_2), 'right')

    labels = [f'≥ {t}x' for t in depths]
    toggle = RadioButtonGroup(labels=labels, active=depths.index(int(t0)))
    # Coverage patterns per depth (library codes and C counts), for the count of Cs that
    # every shown library covers.
    patterns = None
    if upset is not None and upset.height:
        patterns = {}
        for t in depths:
            members, counts = [], []
            for libs_str, cnt in upset.filter(pl.col('min_depth') == t).select(['libraries', 'n_c']).iter_rows():
                codes = [pos[x] for x in (libs_str or '').split(',') if x in pos]
                if codes:
                    members.append(codes)
                    counts.append(int(cnt))
            patterns[str(t)] = {'members': members, 'n': counts}
    common_div = Div(text='', width=FIG_W)
    mode = None
    core_note = ''
    if core_views:
        dropped = [x for x in order if x not in core['members']]
        core_note = (f'Core set, chosen at \u2265 {core["depth"]}x as the largest set of libraries with '
                     f'set Jaccard \u2265 {core["threshold"]:g}: {len(core["members"])} libraries '
                     f'(Jaccard {core["jaccard"]:.3f}); dropped: {", ".join(dropped) or "none"}.')
        mode = RadioButtonGroup(labels=['all libraries (pairwise Cs)',
                                        'core set (Cs common to every core library)'], active=0)
    boxes = CheckboxGroup(labels=list(order), active=list(range(n)), inline=True)
    # Rebuild the matrix from the full data for the chosen depth and libraries. Every
    # element was placed at (column code + offset, row offset + within-cell offset), so a
    # subset only renumbers rows and columns; label boxes are sized in screen pixels and
    # shrink in data units as cells grow.
    # The sources start with the full depth matrix and the t0 view, and the callback keeps
    # those first contents as its master copy. They are therefore left out of the args, so
    # the page embeds them once.
    render = CustomJS(
        args=dict(toggle=toggle, boxes=boxes, views={k: v for k, v in views.items() if k != t0},
                  t0=t0, keys=[str(t) for t in depths], B=B, names=list(order), N=n, pad=pad,
                  span=span, D=D, dw=depth_w, maxd=max_depth, low=src_low, up=src_up, diag=src_diag,
                  lab=src_lab, dlab=src_dlab, llab=src_llab, grid=src_grid, xr=p.x_range,
                  yr=p.y_range, axb=ax_bottom, axl=ax_left, axt=ax_top, axr=ax_right,
                  nx=names_x, ny=names_y, common_div=common_div, patterns=patterns, n_ref=n_ref,
                  mode=mode, core_views=core_views or None, core_codes=core_codes, core_note=core_note),
        code="""
          const busy = window.reportBusy || (() => {});
          const done = window.reportBusyDone || (() => {});
          busy(true);
          setTimeout(() => { try {
            // On the first call the sources still hold the full matrix as built in Python.
            // Keep those arrays; every later rebuild starts from them.
            const M = low._master || (low._master = {
                lower: low.data, lower_labels: llab.data,
                view: {upper: up.data, diag: diag.data, labels: lab.data, diag_labels: dlab.data}});
            const coreMode = !!(mode && mode.active === 1 && core_views);
            const sel = boxes.active.filter(c => !coreMode || core_codes.includes(c)).sort((a, b) => a - b);
            const m = sel.length;
            const at = new Map(sel.map((c, i) => [c, i]));
            const key = keys[toggle.active];
            const v = coreMode ? core_views[key] : (key === t0 ? M.view : views[key]);
            const shrink = m / N;
            const rowY = (code) => m - 1 - at.get(code);
            const keep = (master, ok, place) => {
                const out = {};
                const cols = Object.keys(master);
                for (const k of cols) out[k] = [];
                const len = master[cols[0]].length;
                for (let i = 0; i < len; i++) {
                    if (!ok(master, i)) continue;
                    for (const k of cols) out[k].push(master[k][i]);
                    place(out, out[cols[0]].length - 1, master, i);
                }
                return out;
            };
            // Place each kept cell from its codes: column and row library, then x and y bin
            // out of nb bins.
            const cell = (d, nb) => {
                const o = {x: [], y: [], n: [], lx: [], ly: [], vx: [], vy: []};
                for (let i = 0; i < d.n.length; i++) {
                    if (!at.has(d.lx[i]) || !at.has(d.ly[i])) continue;
                    o.x.push(at.get(d.lx[i]) + pad + (d.vx[i] + 0.5) / nb * span);
                    o.y.push(rowY(d.ly[i]) + pad + (d.vy[i] + 0.5) / nb * span);
                    o.n.push(d.n[i]); o.lx.push(d.lx[i]); o.ly.push(d.ly[i]);
                    o.vx.push(d.vx[i]); o.vy.push(d.vy[i]);
                }
                return o;
            };
            const bars = (master) => keep(master, (d, i) => at.has(d.code[i]),
                (o, j, d, i) => {
                    const dx = at.get(d.code[i]) - d.code[i];
                    const dy = rowY(d.code[i]) - (N - 1 - d.code[i]);
                    o.left[j] += dx; o.right[j] += dx; o.bottom[j] += dy; o.top[j] += dy;
                });
            const labels = (master) => keep(master,
                (d, i) => at.has(d.rc[i]) && at.has(d.cc[i]),
                (o, j, d, i) => {
                    const x = d.x[i] - d.cc[i] + at.get(d.cc[i]);
                    const y = d.y[i] - (N - 1 - d.rc[i]) + rowY(d.rc[i]);
                    o.x[j] = x; o.y[j] = y;
                    o.bw[j] = d.bw[i] * shrink; o.bh[j] = d.bh[i] * shrink;
                    o.bx[j] = d.al[i] ? x : x + o.bw[j] / 2;
                    o.by[j] = y - o.bh[j] / 2;
                });
            low.data = cell(M.lower, D); up.data = cell(v.upper, B); diag.data = bars(v.diag);
            llab.data = labels(M.lower_labels); lab.data = labels(v.labels); dlab.data = labels(v.diag_labels);

            const g = {left: [], right: [], bottom: [], top: []};
            for (let j = 0; j < m; j++) for (let i = 0; i < m; i++) {
                g.left.push(j); g.right.push(j + 1); g.bottom.push(i); g.top.push(i + 1);
            }
            grid.data = g;
            xr.end = Math.max(m, 1); yr.end = Math.max(m, 1);

            const r6 = (x) => Math.round(x * 1e6) / 1e6;
            const depthPos = (c, d) => c + pad + (Math.min(Math.floor(d / dw), D - 1) + 0.5) / D * span;
            const methPos = (c, f) => c + pad + f * span;
            const depthTicks = [[0, '0'], [Math.floor(maxd / 2), String(Math.floor(maxd / 2))],
                                [maxd, '\u2265' + maxd]];
            const methTicks = [[0, '0'], [0.5, '0.5'], [1, '1']];
            const setTicks = (axis, cells, vals, posf) => {
                const t = [], o = new Map();
                for (const c of cells) for (const [val, label] of vals) {
                    const q = r6(posf(c, val)); t.push(q); o.set(q, label);
                }
                axis.ticker.ticks = t; axis.major_label_overrides = o;
            };
            const range = (a, b) => Array.from({length: Math.max(b - a, 0)}, (_, i) => a + i);
            setTicks(axb, range(0, m - 1), depthTicks, depthPos);
            setTicks(axl, range(1, m).map(k => m - 1 - k), depthTicks, depthPos);
            setTicks(axt, range(0, m), methTicks, methPos);
            setTicks(axr, range(0, m - 1).map(k => m - 1 - k), methTicks, methPos);
            nx.ticker.ticks = range(0, m).map(j => j + 0.5);
            nx.major_label_overrides = new Map(sel.map((c, j) => [j + 0.5, names[c]]));
            ny.ticker.ticks = range(0, m).map(i => m - 1 - i + 0.5);
            ny.major_label_overrides = new Map(sel.map((c, i) => [m - 1 - i + 0.5, names[c]]));

            // Cs every shown library covers at this depth, summed over the upset patterns
            // that include all of them. Each cell's r and count stay pairwise.
            if (coreMode) {
                common_div.text = core_note + ' At \u2265 ' + keys[toggle.active] + 'x, <b>' +
                    v.n_common.toLocaleString('en-US') + '</b> Cs are covered by every core library; ' +
                    'every cell and histogram in this view uses exactly those Cs.';
            } else if (common_div && patterns) {
                const pat = patterns[keys[toggle.active]];
                let common = 0;
                if (pat && m > 0) {
                    for (let i = 0; i < pat.n.length; i++) {
                        const mem = new Set(pat.members[i]);
                        if (sel.every(c => mem.has(c))) common += pat.n[i];
                    }
                }
                const pct = n_ref ? (100 * common / n_ref).toFixed(1) : '0';
                common_div.text = '<b>' + common.toLocaleString('en-US') + '</b> Cs (' + pct +
                    '% of reference Cs) are covered at \u2265 ' + keys[toggle.active] + 'x in all ' + m +
                    ' shown librar' + (m === 1 ? 'y' : 'ies') + '. Each cell\u2019s r and C count ' +
                    'use only the Cs covered in both of its two libraries, so they do not change ' +
                    'when other libraries are hidden.';
            }
          } finally {
            done();
          } }, 30);
        """)
    toggle.js_on_change('active', render)
    boxes.js_on_change('active', render)
    if mode is not None:
        mode.js_on_change('active', render)
    # The report's global depth toggle sets this selector; it stays hidden.
    toggle.visible = False
    intro = note(
        'Each cell compares the column library (x) with the row library (y). Lower left, in grey: '
        f'depth of every reference C, from 0 to ≥ {max_depth}x along each axis, with Pearson r '
        'of log(1 + depth). The top depth bin holds every C at or above about the '
        f'{DEPTH_AXIS_PERCENTILE}th percentile of covered depth. '
        'Diagonal: the library’s methylation histogram. Upper right: '
        'methylation of Cs both libraries cover at the minimum depth chosen at the top of the '
        'report, from 0 to 1 along '
        'each axis, with Pearson r and the number of Cs the cell shows. The diagonal bars and the methylation cells share one viridis '
        'scale of C counts. The scales for each cell are on the outer axes: depth along the bottom '
        'and left, methylation along the top and right. Hover over a cell for its values.')
    pick = row(Div(text=f'<b style="color:{TEXT}">libraries shown:</b>'), boxes)
    if mode is not None:
        pick = column(row(Div(text=f'<b style="color:{TEXT}">view:</b>'), mode), pick)
    elif core is not None:
        pick = column(note(f'No core set: no two libraries reach set Jaccard \u2265 {core["threshold"]:g} '
                           f'at \u2265 {core["depth"]}x.'), pick)
    if patterns:
        pat0 = patterns[t0]
        common0 = sum(c for mem, c in zip(pat0['members'], pat0['n']) if len(set(mem)) == n)
        common_div.text = (f'<b>{common0:,}</b> Cs ({100 * common0 / n_ref:.1f}% of reference Cs) '
                           f'are covered at \u2265 {t0}x in all {n} shown libraries. Each cell\u2019s r '
                           'and C count use only the Cs covered in both of its two libraries, so they '
                           'do not change when other libraries are hidden.')
    return TabPanel(title='All pairs', child=column(intro, pick, common_div, toggle, p)), toggle


def read_core(rd, merge_dir):
    """Core-set membership and per-depth core tables from the merge, or None if absent."""
    if not os.path.exists(os.path.join(merge_dir, 'depth_dependent_core.tsv')):
        return None
    info = rd('depth_dependent_core.tsv')
    members = info.filter(pl.col('in_core'))['library'].to_list()
    core = {'members': members, 'depth': int(info['core_depth'][0]),
            'threshold': float(info['core_jaccard_threshold'][0]),
            'jaccard': float(info['core_jaccard'][0]), 'by_depth': {}}
    if not members:
        return core
    cpw = rd('depth_dependent_core_pairwise.tsv')
    ch2 = rd('depth_dependent_core_pair_hist2d.tsv')
    cmh = rd('depth_dependent_core_meth_hist.tsv')
    for t in sorted(cpw['min_depth'].unique().to_list()):
        core['by_depth'][int(t)] = (ch2.filter(pl.col('min_depth') == t),
                                    cmh.filter(pl.col('min_depth') == t),
                                    cpw.filter(pl.col('min_depth') == t))
    return core


def matrix_by_depth(pw, h2, mh):
    """{min depth: (shared 2-D histogram, per-library histogram, pairwise r)} for All pairs."""
    out = {}
    for t in sorted(pw['min_depth'].unique().to_list()):
        out[int(t)] = (
            h2.filter((pl.col('min_depth') == t) & (pl.col('bin_a') >= 0) & (pl.col('bin_b') >= 0)),
            mh.filter((pl.col('min_depth') == t) & (pl.col('site_class') == 'covered')),
            pw.filter(pl.col('min_depth') == t).select(['lib_a', 'lib_b', 'n_shared', 'meth_r']),
        )
    return out


# ---------------------------------------------------------------- clustering tab

def newick(node, labels):
    if node.is_leaf():
        return labels[node.id]
    left = newick(node.left, labels)
    right = newick(node.right, labels)
    return f'({left}:{node.dist - node.left.dist:.6f},{right}:{node.dist - node.right.dist:.6f})'


def cluster(pw, libs, metric):
    """Average-linkage clustering on 1 - metric: (dendrogram, clusters.tsv row)."""
    from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
    from scipy.spatial.distance import squareform

    m = square(pw, libs, metric)
    dist = np.clip(1.0 - np.nan_to_num(m, nan=0.0), 0, None)
    np.fill_diagonal(dist, 0.0)
    Z = linkage(squareform(dist, checks=False), method='average')
    dn = dendrogram(Z, no_plot=True, labels=libs)
    return dn, {'metric': metric, 'leaf_order': ','.join(dn['ivl']),
                'newick': newick(to_tree(Z), libs) + ';'}


def cluster_view(pw, libs, metric, title, dn):
    order = dn['ivl']
    n = len(libs)
    dg = style(figure(frame_width=cell_px(n) * n, height=200, x_range=Range1d(-0.5, n - 0.5),
                      tools='save', min_border_left=label_px(order),
                      title=f'Average linkage on 1 − {metric_name(title)}',
                      y_axis_label='distance'))
    xs = [[(x - 5) / 10 for x in ic] for ic in dn['icoord']]
    dg.multi_line(xs, dn['dcoord'], line_color=TEXT, line_width=1.5)
    dg.xaxis.visible = False
    dg.xgrid.grid_line_color = None
    dg.y_range.start = 0
    return column(dg, heatmap(pw, libs, metric, '', order=order))


def clustering_tab(pw_by_depth, libs, depths):
    """Clustering tab, the layouts to show at each depth, and the clusters.tsv rows.

    Depth-independent metrics get one view. clusters.tsv still has a row for every depth
    and metric.
    """
    views = {m: [] for m, _ in METRICS}
    per_depth, rows = [], []
    for t in depths:
        shown = []
        for metric, title in METRICS:
            dn, r = cluster(pw_by_depth[t], libs, metric)
            rows.append({'min_depth': t, **r})
            if metric not in DEPTH_INDEPENDENT or not views[metric]:
                views[metric].append(cluster_view(pw_by_depth[t], libs, metric, title, dn))
                if metric not in DEPTH_INDEPENDENT:
                    shown.append(views[metric][-1])
        per_depth.append(shown)
    panels = [TabPanel(title=metric_name(title), child=column(*views[m])) for m, title in METRICS]
    intro = note('Libraries that join low in the tree have the most similar profiles. Methylation '
                 'clusters on agreement in methylation at shared Cs. Coverage and Jaccard cluster '
                 'on where reads land, which tracks library preparation and input amount.')
    return TabPanel(title='Clustering', child=column(intro, Tabs(tabs=panels))), per_depth, rows


# ---------------------------------------------------------------- library-specific tab

def private_tab(summary, mh, libs, styles, min_depth):
    # Categorical y runs bottom to top, so reverse to put the first library on top.
    rank = {lib: i for i, lib in enumerate(natural_sorted(summary['library'].to_list())[::-1])}
    s = summary.sort(pl.col('library').replace_strict(rank, return_dtype=pl.Int64))
    src = ColumnDataSource(dict(
        lib=s['library'].to_list(), n=s['n_private'].to_list(),
        frac=[p / c if c else 0 for p, c in zip(s['n_private'].to_list(), s['n_covered'].to_list())],
        covered=s['n_covered'].to_list(),
        color=[styles[lib]['color'] for lib in s['library'].to_list()]))
    p = style(figure(width=FIG_W, height=max(200, 30 * len(libs) + 80), y_range=s['library'].to_list(),
                     title=f'Cs covered (depth ≥ {min_depth}) in this library and no other',
                     tools='save', x_axis_label='Cs'))
    r = p.hbar(y='lib', right='n', height=0.6, source=src, fill_color='color', line_color=None)
    p.add_tools(HoverTool(renderers=[r], tooltips=[('library', '@lib'), ('private Cs', '@n{0,0}'),
                                                   ('of its covered Cs', '@frac{0.00%}'),
                                                   ('covered Cs', '@covered{0,0}')]))
    p.x_range.start = 0
    p.ygrid.grid_line_color = None

    figs = []
    for lib in libs:
        h = mh.filter(pl.col('library') == lib)
        cov = h.filter(pl.col('site_class') == 'covered').sort('bin_lo')
        prv = h.filter(pl.col('site_class') == 'private').sort('bin_lo')
        x = ((cov['bin_lo'] + cov['bin_hi']) / 2).to_list()

        def norm(v):
            t = sum(v)
            return [a / t if t else 0 for a in v]
        f = style(figure(width=270, height=220, title=lib, x_range=Range1d(0, 1), tools='',
                         toolbar_location=None,
                         x_axis_label='methylation', y_axis_label='fraction'))
        f.title.text_font_size = '8pt'
        r1 = f.line(x, norm(cov['n_c'].to_list()), line_color=OTHER, line_width=2)
        r2 = f.line(x, norm(prv['n_c'].to_list()), line_color=styles[lib]['color'], line_width=2)
        if not figs:
            f.add_layout(Legend(items=[LegendItem(label='all covered', renderers=[r1]),
                                       LegendItem(label='private', renderers=[r2])],
                                location='top_center', label_text_font_size='8pt'))
        f.add_tools(HoverTool(renderers=[r1, r2], tooltips=[('methylation', '$x{0.00}'),
                                                            ('fraction', '$y{0.000}')]))
        figs.append(f)
    grid = column(*[row(*figs[i:i + 4]) for i in range(0, len(figs), 4)])
    intro = note('Private Cs reach the depth threshold in one library only. Comparing their '
                 'methylation distribution with all of the library’s covered Cs shows whether '
                 'its extra coverage falls in a distinct methylation class, such as CpG islands.')
    return TabPanel(title='Library-specific Cs', child=column(intro, p, grid))


# ---------------------------------------------------------------- feature coverage tab

def feature_coverage_tab(cov_summary, counts):
    """Embed the feature_cov_meth.nf views, when their outputs are given."""
    panels = []
    if cov_summary:
        from feature_cov_meth_bokeh import build_from_table
        table = pl.read_csv(cov_summary, separator='\t',
                            schema_overrides={'library': pl.Utf8, 'feature': pl.Utf8, 'context': pl.Utf8})
        panels.append(TabPanel(title='By category', child=build_from_table(table)))
    if counts:
        from feature_depth_bokeh import build_layout, make_parser
        layout = build_layout(make_parser().parse_args([counts]))
        panels.append(TabPanel(title='Depth per feature', child=layout))
    if not panels:
        return None
    return TabPanel(title='Feature coverage', child=Tabs(tabs=panels))


# ---------------------------------------------------------------- main

def report_label(libs, n_atoms=3):
    """(file name tag, heading summary) that describe the libraries from their names.

    Each name splits on _ and - into atoms. An atom counts once per library that uses it. The
    tag is up to n_atoms of the most frequent atoms that at least half the libraries share
    (and at least two, unless there is only one library), most frequent first, with ties in order of first appearance. Masses (10ng), replicate
    numbers (1, rep2) and single letters (A, B) mark the design rather than the experiment,
    so they never count. Without such atoms the tag is the range of input amounts in the
    names (0.1ng-100ng), and without either it is ''. The summary adds the number of
    libraries and the input range.
    """
    first_seen, counts = {}, Counter()
    for lib in libs:
        atoms = [a for a in re.split(r'[_\-]', lib) if a]
        for a in dict.fromkeys(atoms):
            if masses(a) or RE_REPLICATE_LIKE.match(a):
                continue
            counts[a] += 1
            first_seen.setdefault(a, len(first_seen))
    # An atom must be in at least two libraries, so one library's own label cannot name the set.
    need = 1 if len(libs) == 1 else max(2, math.ceil(len(libs) / 2))
    common = [a for a, c in counts.items() if c >= need]
    common.sort(key=lambda a: (-counts[a], first_seen[a]))
    atoms = common[:n_atoms]

    found = [m for lib in libs for m in masses(lib)[:1]]
    amounts = ''
    if found:
        lo, hi = min(found)[1], max(found)[1]
        amounts = lo if lo == hi else f'{lo}-{hi}'

    tag = re.sub(r'[^A-Za-z0-9._+-]+', '-', '_'.join(atoms) or amounts).strip('-')
    summary = [', '.join(atoms)] if atoms else []
    summary.append(f'{len(libs)} {"library" if len(libs) == 1 else "libraries"}')
    if found:
        summary.append(f'{amounts.replace("-", " to ")} input')
    return tag, ', '.join(summary)


def depth_percentile(hist, q):
    """Depth at percentile q of covered Cs (depth > 0), pooled over every library, or None."""
    covered = (hist.filter(pl.col('depth') > 0).group_by('depth').agg(pl.col('n_c').sum())
                   .sort('depth'))
    if covered.is_empty() or covered['n_c'].sum() == 0:
        return None
    cum = covered['n_c'].cum_sum()
    return int(covered['depth'][int((cum >= q / 100 * cum[-1]).arg_max())])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('depth_hists', nargs='+', help='<library>.depth_hist.tsv files')
    ap.add_argument('--merge-dir', help='directory holding depth_dependent_merge.py outputs')
    ap.add_argument('--min-depth', type=int, required=True)
    ap.add_argument('--depth-check',
                    help='depth_dependent_check.tsv; adds each library\'s shallowest checked contig '
                         'and depth check result to the depth tab')
    ap.add_argument('-o', '--out', help='output HTML (default <tag>_depth_dependent_report.html, where the '
                                        'tag comes from the library names; see report_label)')
    ap.add_argument('--depth-summary', default='depth_dependent_depth_summary.tsv')
    ap.add_argument('--clusters', default='depth_dependent_clusters.tsv')
    ap.add_argument('--read-stats', nargs='+',
                    help='<library>.read_stats.tsv files (depth_dependent_read_stats); adds read counts to '
                         'the overview')
    ap.add_argument('--feature-counts',
                    help='combined_feature_counts.tsv from feature_cov_meth.nf; adds per-feature depth')
    ap.add_argument('--feature-cov-summary',
                    help='feature_cov_meth_summary.tsv from feature_cov_meth.nf; adds coverage and '
                         'methylation per feature category')
    args = ap.parse_args()

    hist = pl.concat([pl.read_csv(f, separator='\t', schema_overrides={'library': pl.Utf8})
                      for f in args.depth_hists])
    all_libs = natural_sorted(hist['library'].unique().to_list())
    styles = library_styles(all_libs)
    depth_top = depth_percentile(hist, DEPTH_AXIS_PERCENTILE)

    depth_check = None
    if args.depth_check:
        depth_check = {r['library']: r for r in pl.read_csv(
            args.depth_check, separator='\t',
            schema_overrides={'library': pl.Utf8, 'shallowest_contig': pl.Utf8}).iter_rows(named=True)}
    depth_panel, stats = depth_tab(hist, args.min_depth, depth_check, styles)
    summary_df = pl.DataFrame(stats)
    # At a 1x requirement the frac_min_depth column is the same as frac_ge1.
    summary_df = (summary_df.drop('frac_min_depth') if args.min_depth == 1
                  else summary_df.rename({'frac_min_depth': f'frac_ge{args.min_depth}'}))
    summary_df.write_csv(args.depth_summary, separator='\t')
    # Flowcell IDs can be all digits (Element), so read them as text.
    read_stats = (pl.concat([pl.read_csv(f, separator='\t',
                                         schema_overrides={'library': pl.Utf8, 'flowcell': pl.Utf8})
                             for f in args.read_stats], how='diagonal') if args.read_stats else None)
    tag, summary = report_label(all_libs)
    out = args.out or (f'{tag}_' if tag else '') + 'depth_dependent_report.html'
    tabs = [overview_tab(hist, read_stats, args.min_depth), depth_panel]

    depth_toggle = None
    if args.merge_dir:
        def rd(name):
            path = os.path.join(args.merge_dir, name)
            with open(path) as fh:
                header = fh.readline().rstrip('\n').split('\t')
            text_cols = {'lib_a', 'lib_b', 'library', 'libraries', 'site_class'}
            return pl.read_csv(path, separator='\t',
                               schema_overrides={c: pl.Utf8 for c in header if c in text_cols})
        pw_all = rd('depth_dependent_pairwise.tsv')
        sum_all = rd('depth_dependent_summary.tsv')
        h2_all = rd('depth_dependent_pair_hist2d.tsv')
        mh_all = rd('depth_dependent_meth_hist.tsv')
        up_all = rd('depth_dependent_upset.tsv')
        libs = natural_sorted(sum_all['library'].unique().to_list())
        # Normally the same libraries as the depth histograms; style any extras too.
        styles = library_styles(natural_sorted(set(all_libs) | set(libs)))
        n_ref = int(sum_all['n_reference_c'][0])
        depths = sorted(int(t) for t in pw_all['min_depth'].unique().to_list())
        t0 = args.min_depth if args.min_depth in depths else depths[-1]
        pair_bins = mh_all.group_by(['min_depth', 'library', 'site_class']).len()['len'].max()

        def at(df, t):
            return df.filter(pl.col('min_depth') == t).drop('min_depth')

        # Each depth gets its own copy of every depth-dependent view; the toggle shows one.
        # Views that do not depend on depth are built once and always shown.
        pw_by_depth = {t: at(pw_all, t) for t in depths}
        upset_views, private_views = [], []
        for t in depths:
            if up_all.height:
                upset_views.append(upset_view(at(up_all, t), libs, n_ref))
            private_views.append(private_tab(at(sum_all, t), at(mh_all, t), libs, styles, t).child)
        corr_panel, corr_by_depth, explorer_depth = correlation_tab(
            pw_by_depth, {t: at(h2_all, t) for t in depths}, libs, depths, t0, pair_bins)
        clus_panel, clus_by_depth, cluster_rows = clustering_tab(pw_by_depth, libs, depths)
        groups = [([upset_views[i]] if upset_views else []) + corr_by_depth[i] + clus_by_depth[i]
                  + [private_views[i]] for i in range(len(depths))]
        for t, layouts in zip(depths, groups):
            for layout in layouts:
                layout.visible = (t == t0)
        pl.DataFrame(cluster_rows).write_csv(args.clusters, separator='\t')

        if upset_views:
            tabs.append(TabPanel(title='Upset', child=column(*upset_views)))
        tabs.append(corr_panel)
        pairs_panel, pairs_toggle = all_pairs_tab(pw_all.filter(pl.col('min_depth') == t0),
                                                  rd('depth_dependent_depth_hist2d.tsv'),
                                                  matrix_by_depth(pw_all, h2_all, mh_all), libs, t0,
                                                  upset=up_all, n_ref=n_ref, core=read_core(rd, args.merge_dir),
                                                  depth_top=depth_top)
        tabs.append(pairs_panel)
        tabs.append(clus_panel)
        tabs.append(TabPanel(title='Library-specific Cs', child=column(*private_views)))

        depth_toggle = RadioButtonGroup(labels=[f'\u2265 {t}x' for t in depths], active=depths.index(t0))
        depth_toggle.js_on_change('active', CustomJS(
            args=dict(groups=groups, pairs=pairs_toggle, explorer=explorer_depth),
            code="""
                const k = cb_obj.active;
                const busy = window.reportBusy || (() => {});
                const done = window.reportBusyDone || (() => {});
                busy(true);
                // Yield a frame so the busy panel paints before the switch blocks the page.
                setTimeout(() => {
                    try {
                        groups.forEach((layouts, i) => layouts.forEach(l => { l.visible = (i === k); }));
                        pairs.active = k;
                        explorer.active = k;
                    } finally {
                        done();
                    }
                }, 30);
            """))

    feature_tab = feature_coverage_tab(args.feature_cov_summary, args.feature_counts)
    if feature_tab is not None:
        tabs.append(feature_tab)

    title = f'Depth-dependent analyses \u2013 {summary}'
    header = Div(text=f'<h2 style="color:{TEXT};margin-bottom:0">{html.escape(title)}</h2>'
                      f'<div style="color:{TEXT_2}">C coverage and methylation, CpG context</div>')
    body = [header]
    if depth_toggle is not None:
        body.append(row(Div(text=f'<b style="color:{TEXT}">minimum depth per C for the cross-library '
                                 f'tabs:</b>'), depth_toggle))
    body.append(Tabs(tabs=tabs))
    save_report(column(*body), out, title, 'depth-dependent analyses')
    print(f'[depth_dependent_report] wrote {out} with {len(tabs)} tabs', file=sys.stderr)


if __name__ == '__main__':
    main()
