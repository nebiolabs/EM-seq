#!/usr/bin/env python
"""Hex-density plot of per-feature library depth.

Dense regions render as hex tiles colored by point count; sparse hexes
(<= --outlier-threshold) drop their points back in as individual markers so
outliers stay visible. Within each library column points are spread along x
by chromosome so they don't collapse on a single line.
"""
import argparse
import math
import re
import sys

import numpy as np
import polars as pl
from bokeh.io import output_file, save
from bokeh.layouts import column as bk_column, row as bk_row
from bokeh.models import (
    BasicTicker,
    ColorBar,
    ColumnDataSource,
    CustomJS,
    CustomJSTickFormatter,
    Div,
    FixedTicker,
    HoverTool,
    Label,
    LabelSet,
    LinearColorMapper,
    RadioButtonGroup,
    Range1d,
    Select,
    Span,
)
from bokeh.palettes import Category10, Viridis256
from bokeh.plotting import figure
from bokeh.util.hex import cartesian_to_axial, hexbin

# Depth threshold options available in the HTML UI dropdown.
# Each value counts features with depth_orig >= that threshold.
# The largest value (9) is the highest selectable option; nothing is excluded above it.
DEPTH_THRESHOLDS = [1, 2, 5, 7, 9]
# Strips any trailing combination of lowercase dot-tags before .bam
# e.g. sample.ds.md.bam → sample, sample.sorted.dedup.bam → sample, sample.bam → sample
DEFAULT_STRIP_SUFFIX = r'(\.[a-z]+)*\.bam$'
ID_COLS = ['File', 'Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length']


def analyze_name_components(libraries):
    """Choose which name-component positions to use for group and subgroup labels.

    All library names are tokenised on '_' and '-'.  At each component position
    the number of distinct values across all libraries is counted.  Positions
    where every library shares the same value carry no information and are
    skipped.  The remaining positions are ranked by ascending distinct-value
    count:

    * fewest distinct values  → ``group``    (top bold separator + label)
    * second fewest distinct  → ``subgroup`` (angled annotation above replicates)

    The last component  serves as ``replicate`` (x-axis tick label)
    if it looks like a replicate
        r'^(?:rep(?:licate)?|r)?\d+$'
        r'|^[A-Za-z]$'

    Returns ``(group_idx, subgroup_idx)``.  Either may be ``None`` when there
    are not enough informative positions.
    """
    if not libraries:
        return None, None

    all_parts = [re.split(r'[_\-]', lib) for lib in libraries]
    max_len = max(len(p) for p in all_parts)
    # Pad shorter splits so positions align across all libraries
    padded = [p + [''] * (max_len - len(p)) for p in all_parts]

    # Only exclude the last position from frequency ranking if every library's
    # last component looks like a replicate identifier (digits, rep1, r2, A, …).
    # If the last component is something meaningful like 'control' or 'treated'
    # it should compete normally for the group/subgroup roles.
    last_vals = [re.split(r'[_\-]', lib)[-1] for lib in libraries]
    exclude_last = all(RE_REPLICATE_LIKE.match(v) for v in last_vals)
    search_end = max_len - 1 if exclude_last else max_len

    informative = []
    for i in range(search_end):
        n_distinct = len({row[i] for row in padded})
        if n_distinct > 1:          # skip positions identical across all libraries
            informative.append((n_distinct, i))

    informative.sort()              # ascending: fewest distinct values first
    group_idx    = informative[0][1] if len(informative) >= 1 else None
    subgroup_idx = informative[1][1] if len(informative) >= 2 else None
    return group_idx, subgroup_idx


def get_component(name, idx):
    """Return the *idx*-th component of *name* (split on '_' or '-'), or ''."""
    if idx is None:
        return ''
    parts = re.split(r'[_\-]', name)
    return parts[idx] if idx < len(parts) else ''


RE_METACHAR = re.compile(r'[\\^$.|?*+(){}[\]]')

# Matches component values that look like replicate identifiers:
#   pure digits (1, 2, 12), optional rep/r prefix + digits (rep1, r2, REP3,
#   replicate1), or a single letter (A, B, a, b).
RE_REPLICATE_LIKE = re.compile(
    r'^(?:rep(?:licate)?|r)?\d+$'
    r'|^[A-Za-z]$',
    re.IGNORECASE,
)


def build_strip_re(raw):
    """Return a regex pattern that strips *raw* from library names.

    If *raw* contains regex metacharacters it is used verbatim as a pattern.
    Otherwise it is treated as a literal suffix and anchored to end-of-string,
    so plain strings like '.ds.md.bam' work without requiring the user to
    escape dots or add a '$'.
    """
    if RE_METACHAR.search(raw):
        return raw
    return re.escape(raw) + '$'


def chr_sort_key(c):
    s = str(c)
    m = re.match(r'(?:chr)?(\d+|X|Y|M|MT|Mt|Mito)$', s, re.IGNORECASE)
    if m:
        v = m.group(1)
        if v.isdigit():
            return (0, int(v), '')
        return (1, 0, v)
    return (2, 0, s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('counts_tsv')
    ap.add_argument('-o', '--out', default='feature_depth_bokeh.html')
    ap.add_argument('--hex-size-x', type=float, default=0.06,
                    help='hex tile size in x data units (default 0.06)')
    ap.add_argument('--outlier-threshold', type=int, default=2,
                    help='hexes with this count or fewer render as points (default 2)')
    ap.add_argument('--y-scale', choices=['log', 'linear'], default='log',
                    help="initial vertical axis scale: 'log' (depth log10-transformed, "
                         "default) or 'linear' (raw depth, includes zeros). Switchable in HTML UI.")
    ap.add_argument('--plot-type', choices=['hex', 'violin'], default='hex',
                    help="initial plot type: 'hex' (hexbin density + outliers, default) "
                         "or 'violin' (per-library kernel density). Switchable in HTML UI.")
    ap.add_argument('--y-range', choices=['per-feature', 'common'], default='per-feature',
                    help="initial y-range mode: 'per-feature' (each panel autoscales, default) "
                         "or 'common' (all panels share one scale). Switchable in HTML UI.")
    ap.add_argument('--violin-bins', type=int, default=40,
                    help='histogram bins per library for violin density estimate (default 40)')
    ap.add_argument('--depth-cap', type=float, default=None,
                    help='cap raw depth at this value (rows above stack at the cap, '
                         'preventing single outliers from compressing the visible scale).')
    ap.add_argument('--depth-cap-percentile', type=float, default=99.0,
                    help='cap at this percentile of non-zero depths (default 99). '
                         'Ignored if --depth-cap is given. Pass 100 to disable.')
    ap.add_argument('--strip-suffix', default=DEFAULT_STRIP_SUFFIX,
                    help='Suffix to strip from BAM column names when building library labels. '
                         'Plain text (e.g. ".ds.md.bam") is matched literally at end-of-string. '
                         'A value containing regex metacharacters is used as a raw pattern '
                         '(e.g. r"(\\.[a-z]+)*\\.bam$"). '
                         'Default: strip any trailing lowercase dot-tags before .bam.')
    args = ap.parse_args()
    strip_re = build_strip_re(args.strip_suffix)

    with open(args.counts_tsv) as fh:
        header_cols = fh.readline().rstrip('\n').split('\t')
    sample_cols = [c for c in header_cols if c not in ID_COLS]
    if not sample_cols:
        sys.exit('no sample columns found in counts table')

    overrides = {'Geneid': pl.Utf8, 'Chr': pl.Utf8}
    overrides.update({c: pl.Float64 for c in sample_cols})

    lf = pl.scan_csv(
        args.counts_tsv,
        separator='\t',
        schema_overrides=overrides,
    )

    long_lf = (
        lf.select(['File', 'Geneid', 'Chr', *sample_cols])
          .unpivot(
              index=['File', 'Geneid', 'Chr'],
              on=sample_cols,
              variable_name='library',
              value_name='depth',
          )
          .with_columns([
              pl.col('library').str.replace(strip_re, ''),
              pl.col('depth').cast(pl.Float64, strict=False),
          ])
          .filter(pl.col('depth').is_not_null())
    )
    long = long_lf.collect(engine='streaming')
    if long.is_empty():
        sys.exit('no rows to plot')

    # Snapshot original (pre-cap) depths for the per-library threshold counter.
    long = long.with_columns(pl.col('depth').alias('depth_orig'))

    # Apply optional depth cap. Rows above the cap collapse onto the cap value
    # so a few extreme outliers can't compress the visible distribution.
    depth_cap = None
    if args.depth_cap is not None:
        depth_cap = float(args.depth_cap)
    elif args.depth_cap_percentile is not None:
        nonzero = long.filter(pl.col('depth') > 0)['depth']
        if nonzero.len() > 0:
            depth_cap = float(nonzero.quantile(args.depth_cap_percentile / 100.0))
    if depth_cap is not None and depth_cap > 0:
        n_over = long.filter(pl.col('depth') > depth_cap).height
        print(
            f'[feature_depth_bokeh] depth cap = {depth_cap:g} '
            f'(clipping {n_over} of {long.height} rows; '
            f'depth max before clip = {float(long["depth"].max()):g})',
            file=sys.stderr,
        )
        long = long.with_columns(
            pl.when(pl.col('depth') > depth_cap).then(depth_cap).otherwise(pl.col('depth')).alias('depth')
        )
    else:
        print('[feature_depth_bokeh] no depth cap applied', file=sys.stderr)

    libraries = long['library'].unique().to_list()
    group_idx, subgroup_idx = analyze_name_components(libraries)

    meta_rows = []
    for lib in libraries:
        parts = re.split(r'[_\-]', lib)
        replicate = parts[-1] if len(parts) > 1 else lib
        group    = get_component(lib, group_idx) or lib   # fall back to full name if no grouping
        subgroup = get_component(lib, subgroup_idx)
        meta_rows.append({'library': lib, 'group': group, 'subgroup': subgroup, 'replicate': replicate})
    meta = pl.DataFrame(meta_rows)

    # Sort by group → subgroup → replicate for a consistent, predictable x-axis ordering
    meta = (
        meta.sort(['group', 'subgroup', 'replicate'])
            .with_columns(pl.int_range(0, pl.len()).alias('idx'))
    )

    groups = meta['group'].unique(maintain_order=False).sort().to_list()
    n_groups = len(groups)
    palette_size = max(3, min(n_groups, 10))
    palette = Category10[palette_size]
    color_map = {g: palette[i % palette_size] for i, g in enumerate(groups)}
    meta = meta.with_columns(
        pl.col('group').replace_strict(color_map, default='#888888', return_dtype=pl.Utf8).alias('color')
    )

    long = long.join(
        meta.select(['library', 'group', 'subgroup', 'replicate', 'idx', 'color']),
        on='library',
        how='inner',
    )

    chroms = sorted(long['Chr'].drop_nulls().unique().to_list(), key=chr_sort_key)
    n_chr = len(chroms)
    if n_chr > 1:
        chr_offset = {c: -0.4 + 0.8 * i / (n_chr - 1) for i, c in enumerate(chroms)}
    elif chroms:
        chr_offset = {chroms[0]: 0.0}
    else:
        chr_offset = {}

    long = long.with_columns(
        pl.col('Chr').replace_strict(chr_offset, default=0.0, return_dtype=pl.Float64).alias('chr_off')
    ).with_columns(
        (pl.col('idx') + pl.col('chr_off')).alias('x'),
    )

    N = meta.height
    replicate_by_idx = {int(i): str(r) for i, r in zip(meta['idx'].to_list(), meta['replicate'].to_list())}
    ab_x = meta['idx'].to_list()
    ab_text = meta['subgroup'].to_list()

    # Precompute per-(feature, library, threshold) counts of rows with
    # depth_orig >= threshold. Used for the in-HTML "show count above" widget.
    threshold_counts = {}  # {(feature, lib_idx, ti): count}
    for ti, t in enumerate(DEPTH_THRESHOLDS):
        g = long.filter(pl.col('depth_orig') >= t).group_by(['File', 'idx']).len()
        for row in g.iter_rows(named=True):
            threshold_counts[(row['File'], int(row['idx']), ti)] = int(row['len'])

    thr = args.outlier_threshold
    size = args.hex_size_x

    FIG_W, FIG_H = 1500, 600

    # cap value in y_raw space for both scales (None if no cap)
    cap_y_raw = {
        'log': np.log10(depth_cap) if depth_cap is not None and depth_cap > 0 else None,
        'linear': depth_cap if depth_cap is not None else None,
    }

    lib_color = dict(zip(meta['idx'].to_list(), meta['color'].to_list()))
    n_bins = args.violin_bins
    # per-column metadata, populated inside build_column and consumed by the JS callback.
    # Keyed by (scale, plot_type, y_mode). Storing only the visible column's items
    # gets touched per JS callback so dropdown changes stay fast.
    col_labels = {}            # {key: [LabelSet, ...]}
    col_spans = {}             # {key: [Span, ...]}
    col_span_y_min = {}        # {key: [float, ...]}
    col_span_y_to_norm = {}    # {key: [float, ...]}
    col_is_log = {}            # {key: bool}

    def render_violin(p, sub):
        """Per-library histogram-density violin polygons."""
        smooth_kernel = np.array([1, 2, 4, 2, 1], dtype=float)
        smooth_kernel /= smooth_kernel.sum()
        # one partition pass instead of N filter scans
        parts = sub.partition_by('idx', as_dict=True)
        for key, lib_sub in parts.items():
            lib_idx = key[0] if isinstance(key, tuple) else int(key)
            y_vals = lib_sub['y'].to_numpy()
            if len(y_vals) < 3:
                continue
            ymin_l, ymax_l = float(y_vals.min()), float(y_vals.max())
            if ymax_l <= ymin_l:
                continue
            edges = np.linspace(ymin_l, ymax_l, n_bins + 1)
            counts, _ = np.histogram(y_vals, bins=edges)
            smoothed = np.convolve(counts.astype(float), smooth_kernel, mode='same')
            d_max = smoothed.max()
            if d_max <= 0:
                continue
            half_width = smoothed / d_max * 0.4
            centers = (edges[:-1] + edges[1:]) / 2
            xs = np.concatenate([lib_idx - half_width, (lib_idx + half_width)[::-1]])
            ys = np.concatenate([centers, centers[::-1]])
            color = lib_color.get(lib_idx, '#888888')
            p.patch(
                xs.tolist(), ys.tolist(),
                fill_color=color, fill_alpha=0.45,
                line_color=color, line_alpha=0.9, line_width=1,
            )
            # median tick
            median = float(np.median(y_vals))
            p.line(
                [lib_idx - 0.4, lib_idx + 0.4], [median, median],
                line_color=color, line_width=2, line_alpha=0.9,
            )

    def build_column(scale, plot_type, y_mode):
        col_key = (scale, plot_type, y_mode)
        col_labels[col_key] = []
        col_spans[col_key] = []
        col_span_y_min[col_key] = []
        col_span_y_to_norm[col_key] = []
        col_is_log[col_key] = (scale == 'log')

        use_log = scale == 'log'
        if use_log:
            y_axis_label = 'depth (log scale)'
            raw_expr = pl.when(pl.col('depth') > 0).then(pl.col('depth').log10())
        else:
            y_axis_label = 'depth (linear)'
            raw_expr = pl.col('depth')

        base = long.with_columns(raw_expr.alias('y_raw'))
        plottable_all = base.filter(pl.col('y_raw').is_not_null())
        if plottable_all.is_empty():
            return None

        # Pixel-matched y normalization so hex tiles render regular regardless of
        # data range. y_raw is rescaled into [0, y_norm_scale] using either the
        # whole dataset's min/max (y_mode == 'common') or each feature's own
        # min/max (y_mode == 'per_feature'). Original-depth ticks are restored
        # via the JS tick formatter below.
        y_norm_scale = N * FIG_H / FIG_W
        if y_mode == 'common':
            common_y_min = float(plottable_all['y_raw'].min())
            common_y_max = float(plottable_all['y_raw'].max())
            common_y_span = max(common_y_max - common_y_min, 1e-9)

        plots = []
        for feature in sorted(base['File'].unique().to_list()):
            sub_all = base.filter(pl.col('File') == feature)
            if sub_all.is_empty():
                continue
            sub_nonnull = sub_all.filter(pl.col('y_raw').is_not_null())

            if sub_nonnull.is_empty():
                # placeholder panel using a nominal y range
                y_top = y_norm_scale * 1.55
                y_bot = -y_norm_scale * 0.05
                p = figure(
                    x_range=Range1d(start=-0.6, end=N - 0.4),
                    y_range=Range1d(start=y_bot, end=y_top),
                    title=f'{feature}: depth per library  (no non-zero counts; '
                          f'{sub_all.height} rows total, all zero)',
                    width=FIG_W, height=300,
                    tools='pan,box_zoom,wheel_zoom,reset,save',
                )
                p.add_layout(Label(
                    x=(N - 1) / 2, y=(y_bot + y_top) / 2,
                    text='no non-zero read counts in this feature category',
                    text_align='center', text_font_style='italic',
                    text_color='#888888',
                ))
                p.xaxis[0].ticker = FixedTicker(ticks=list(range(N)))
                p.xaxis[0].major_label_overrides = replicate_by_idx
                p.xgrid.grid_line_color = None
                plots.append(p)
                continue

            if y_mode == 'common':
                y_raw_min = common_y_min
                y_raw_max = common_y_max
                y_raw_span = common_y_span
            else:
                y_raw_min = float(sub_nonnull['y_raw'].min())
                y_raw_max = float(sub_nonnull['y_raw'].max())
                y_raw_span = max(y_raw_max - y_raw_min, 1e-9)
            y_to_norm = y_norm_scale / y_raw_span

            sub_all = sub_all.with_columns(
                ((pl.col('y_raw') - y_raw_min) * y_to_norm).alias('y')
            )
            sub = sub_all.filter(pl.col('y').is_not_null())

            y_top = y_norm_scale * 1.55
            y_bot = -y_norm_scale * 0.05
            ab_label_y = y_norm_scale * 1.05
            count_label_y = y_norm_scale * 1.30
            enz_label_y = y_norm_scale * 1.45

            if plot_type == 'hex':
                title_suffix = '(hex density + outliers; x-jitter within each library = chromosome)'
            else:
                title_suffix = '(per-library histogram density; median line at horizontal tick)'

            p = figure(
                x_range=Range1d(start=-0.6, end=N - 0.4),
                y_range=Range1d(start=y_bot, end=y_top),
                title=f'{feature}: depth per library  {title_suffix}',
                width=FIG_W,
                height=FIG_H,
                tools='pan,box_zoom,wheel_zoom,reset,save',
            )

            if plot_type == 'hex':
                x = sub['x'].to_numpy()
                y = sub['y'].to_numpy()
                hb = hexbin(x, y, size)
                bins = pl.DataFrame({
                    'q': np.asarray(hb.q).astype(np.int64),
                    'r': np.asarray(hb.r).astype(np.int64),
                    'counts': np.asarray(hb.counts).astype(np.int64),
                })
                qs, rs = cartesian_to_axial(x, y, size, 'pointytop')
                sub_qr = sub.with_columns([
                    pl.Series('q', qs.astype(np.int64)),
                    pl.Series('r', rs.astype(np.int64)),
                ])
                # hexbin() already returns per-bin counts; reuse them instead of a fresh groupby.
                outlier_qr = bins.filter(pl.col('counts') <= thr).select(['q', 'r'])
                outliers = sub_qr.join(outlier_qr, on=['q', 'r'], how='inner')
                keep_bins = bins.filter(pl.col('counts') > thr)

                if keep_bins.height > 0:
                    mapper = LinearColorMapper(
                        palette=Viridis256,
                        low=float(keep_bins['counts'].min()),
                        high=float(keep_bins['counts'].max()),
                    )
                    p.hex_tile(
                        q='q', r='r', size=size,
                        source=ColumnDataSource(dict(
                            q=keep_bins['q'].to_list(),
                            r=keep_bins['r'].to_list(),
                            counts=keep_bins['counts'].to_list(),
                        )),
                        fill_color={'field': 'counts', 'transform': mapper},
                        line_color=None, alpha=0.85,
                    )
                    p.add_layout(ColorBar(
                        color_mapper=mapper,
                        ticker=BasicTicker(desired_num_ticks=6),
                        title='points per hex',
                        title_text_font_size='9pt',
                        major_label_text_font_size='8pt',
                        width=12, padding=4,
                    ), 'right')

                if outliers.height > 0:
                    osrc = ColumnDataSource(dict(
                        x=outliers['x'].to_list(),
                        y=outliers['y'].to_list(),
                        depth=outliers['depth'].to_list(),
                        library=outliers['library'].to_list(),
                        group=outliers['group'].to_list(),
                        subgroup=outliers['subgroup'].to_list(),
                        chrom=outliers['Chr'].to_list(),
                        geneid=outliers['Geneid'].cast(pl.Utf8).to_list(),
                    ))
                    r = p.scatter(
                        x='x', y='y', source=osrc,
                        size=4, alpha=0.55, color='#444444', line_color=None,
                    )
                    p.add_tools(HoverTool(renderers=[r], tooltips=[
                        ('library', '@library'),
                        ('group', '@group'),
                        ('subgroup', '@subgroup'),
                        ('chr', '@chrom'),
                        ('geneid', '@geneid'),
                        ('depth', '@depth'),
                    ]))
            else:
                render_violin(p, sub)

            cap_raw = cap_y_raw['log' if use_log else 'linear']
            if cap_raw is not None and y_raw_min <= cap_raw <= y_raw_max:
                cap_y = (cap_raw - y_raw_min) * y_to_norm
                p.add_layout(Span(
                    location=cap_y, dimension='width',
                    line_color='#cc3333', line_width=1,
                    line_dash='dashed', line_alpha=0.6,
                ))
                p.add_layout(Label(
                    x=N - 0.4, y=cap_y,
                    text=f'cap @ depth={depth_cap:g}',
                    text_align='right', text_baseline='bottom',
                    text_font_size='8pt', text_color='#cc3333',
                ))

            # threshold marker line, position updated via JS on dropdown change.
            # level='overlay' forces it above hex tiles / violins.
            initial_t = DEPTH_THRESHOLDS[2]  # 5: a middle value visually obvious by default
            initial_y_raw = math.log10(initial_t) if use_log else float(initial_t)
            initial_loc = (initial_y_raw - y_raw_min) * y_to_norm
            tspan = Span(
                location=initial_loc, dimension='width',
                line_color='#ff6600', line_width=2.0,
                line_dash='dashed', line_alpha=0.95,
                level='overlay',
            )
            p.add_layout(tspan)
            col_spans[col_key].append(tspan)
            col_span_y_min[col_key].append(float(y_raw_min))
            col_span_y_to_norm[col_key].append(float(y_to_norm))

            p.yaxis.axis_label = y_axis_label
            p.yaxis.formatter = CustomJSTickFormatter(
                args=dict(y_min=y_raw_min, slope=1.0 / y_to_norm, is_log=use_log),
                code="""
                    const v = tick * slope + y_min;
                    if (is_log) {
                        const raw = Math.pow(10, v);
                        if (raw >= 1000) return raw.toExponential(0);
                        return raw.toFixed(raw >= 10 ? 0 : 1);
                    }
                    if (Math.abs(v) >= 1000) return v.toExponential(1);
                    return v.toFixed(0);
                """,
            )
            p.xaxis[0].ticker = FixedTicker(ticks=list(range(N)))
            p.xaxis[0].major_label_overrides = replicate_by_idx
            p.xaxis[0].axis_label = 'replicate'
            p.xgrid.grid_line_color = None

            p.add_layout(LabelSet(
                x='x', y='y', text='text',
                source=ColumnDataSource(dict(
                    x=ab_x, text=ab_text, y=[ab_label_y] * N,
                )),
                text_font_size='8pt', angle=0.6, text_align='left',
            ))

            count_cds_data = {
                'x': list(range(N)),
                'y': [count_label_y] * N,
            }
            for ti in range(len(DEPTH_THRESHOLDS)):
                count_cds_data[f'count_t{ti}'] = [
                    (str(c) if (c := threshold_counts.get((feature, lib_idx, ti), 0)) > 0 else '')
                    for lib_idx in range(N)
                ]
            count_label = LabelSet(
                x='x', y='y', text='count_t2',  # initial = '5' (DEPTH_THRESHOLDS[2])
                source=ColumnDataSource(count_cds_data),
                text_font_size='8pt', text_align='center', text_color='#555555',
            )
            p.add_layout(count_label)
            col_labels[col_key].append(count_label)

            for grp_name in groups:
                grp = meta.filter(pl.col('group') == grp_name)
                if grp.is_empty():
                    continue
                x_start = float(grp['idx'].min())
                x_end = float(grp['idx'].max())
                p.add_layout(Span(
                    location=x_end + 0.5, dimension='height',
                    line_color='black', line_width=1,
                    line_dash='dotted', line_alpha=0.4,
                ))
                p.add_layout(Label(
                    x=(x_start + x_end) / 2, y=enz_label_y,
                    text=grp_name,
                    text_font_size='11pt', text_font_style='bold',
                    text_align='center', text_color='#333333',
                ))

            plots.append(p)

        return bk_column(*plots) if plots else None

    cols = {}
    for scale in ('log', 'linear'):
        for ptype in ('hex', 'violin'):
            for ymode in ('per_feature', 'common'):
                cols[(scale, ptype, ymode)] = build_column(scale, ptype, ymode)

    if all(c is None for c in cols.values()):
        sys.exit('no plottable depth rows')

    initial_scale = args.y_scale
    initial_type = args.plot_type
    initial_ymode = 'common' if args.y_range == 'common' else 'per_feature'
    for (scale, ptype, ymode), col in cols.items():
        if col is None:
            continue
        col.visible = (scale == initial_scale and ptype == initial_type and ymode == initial_ymode)

    js_args = {}
    for (scale, ptype, ymode), col in cols.items():
        if col is not None:
            js_args[f'col_{scale}_{ptype}_{ymode}'] = col

    scale_toggle = RadioButtonGroup(
        labels=['log10', 'linear'],
        active=0 if initial_scale == 'log' else 1,
    )
    type_toggle = RadioButtonGroup(
        labels=['hex', 'violin'],
        active=0 if initial_type == 'hex' else 1,
    )
    range_toggle = RadioButtonGroup(
        labels=['per-feature', 'common'],
        active=0 if initial_ymode == 'per_feature' else 1,
    )
    js_args['scale_toggle'] = scale_toggle
    js_args['type_toggle'] = type_toggle
    js_args['range_toggle'] = range_toggle

    vis_lines = []
    for (scale, ptype, ymode), col in cols.items():
        if col is None:
            continue
        name = f'col_{scale}_{ptype}_{ymode}'
        vis_lines.append(
            f"{name}.visible = (scale === '{scale}' && ptype === '{ptype}' && ymode === '{ymode}');"
        )

    threshold_select = Select(
        value=str(DEPTH_THRESHOLDS[2]),  # default = 5, matches initial label/span positions
        options=[str(t) for t in DEPTH_THRESHOLDS],
        width=80,
    )

    # serialize per-column meta into flat dicts keyed by 'scale_ptype_ymode' for JS lookup
    def _ckey(s, p, m):
        return f'{s}_{p}_{m}'
    labels_by_col_js = {_ckey(*k): v for k, v in col_labels.items()}
    spans_by_col_js = {_ckey(*k): v for k, v in col_spans.items()}
    ymin_by_col_js = {_ckey(*k): v for k, v in col_span_y_min.items()}
    ynorm_by_col_js = {_ckey(*k): v for k, v in col_span_y_to_norm.items()}
    islog_by_col_js = {_ckey(*k): v for k, v in col_is_log.items()}

    js_args['threshold_select'] = threshold_select
    js_args['threshold_options'] = [str(t) for t in DEPTH_THRESHOLDS]
    js_args['labels_by_col'] = labels_by_col_js
    js_args['spans_by_col'] = spans_by_col_js
    js_args['ymin_by_col'] = ymin_by_col_js
    js_args['ynorm_by_col'] = ynorm_by_col_js
    js_args['islog_by_col'] = islog_by_col_js

    js_code = """
        const want_log = scale_toggle.labels[scale_toggle.active] === 'log10';
        const want_hex = type_toggle.labels[type_toggle.active] === 'hex';
        const want_perf = range_toggle.labels[range_toggle.active] === 'per-feature';
        const scale = want_log ? 'log' : 'linear';
        const ptype = want_hex ? 'hex' : 'violin';
        const ymode = want_perf ? 'per_feature' : 'common';
    """ + '\n        '.join(vis_lines) + """
        // only touch the now-visible column's labels and spans
        const ckey = scale + '_' + ptype + '_' + ymode;
        const ti = threshold_options.indexOf(threshold_select.value);
        const t_val = parseFloat(threshold_select.value);
        const labels = labels_by_col[ckey];
        if (labels && ti >= 0) {
            const field_name = 'count_t' + ti;
            for (const ls of labels) {
                ls.text = {field: field_name};
            }
        }
        const spans = spans_by_col[ckey];
        const ymins = ymin_by_col[ckey];
        const ynorms = ynorm_by_col[ckey];
        const is_log = islog_by_col[ckey];
        if (spans) {
            const y_raw = is_log ? Math.log10(t_val) : t_val;
            for (let i = 0; i < spans.length; i++) {
                spans[i].location = (y_raw - ymins[i]) * ynorms[i];
            }
        }
    """
    callback = CustomJS(args=js_args, code=js_code)
    scale_toggle.js_on_change('active', callback)
    type_toggle.js_on_change('active', callback)
    range_toggle.js_on_change('active', callback)
    threshold_select.js_on_change('value', callback)

    scale_div = Div(text='<b>y-axis scale:</b>')
    type_div = Div(text='<b>plot type:</b>')
    range_div = Div(text='<b>y-axis range:</b>')
    threshold_div = Div(text='<b>show count of points with depth &ge;</b>')
    controls = bk_row(
        bk_column(scale_div, scale_toggle),
        bk_column(type_div, type_toggle),
        bk_column(range_div, range_toggle),
        bk_column(threshold_div, threshold_select),
    )

    children = [controls]
    for col in cols.values():
        if col is not None:
            children.append(col)
    layout = bk_column(*children)

    output_file(args.out, title='Feature depth per library')
    save(layout)


if __name__ == '__main__':
    main()
