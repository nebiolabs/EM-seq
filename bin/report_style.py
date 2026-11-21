"""Colours and figure styling shared by the Bokeh reports.

depth_dependent_report.py and feature_cov_meth_bokeh.py both draw with these, so their pages
look alike. Keeping them here also lets feature_cov_meth_bokeh.py avoid importing the whole
depth-dependent report.
"""
from bokeh.models import Div

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


def style(p):
    """Apply the shared look to a figure and return it."""
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
    """A block of secondary text, as wide as a standard figure."""
    return Div(text=f'<div style="color:{TEXT_2};max-width:{FIG_W}px">{html}</div>')
