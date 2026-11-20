"""Save Bokeh HTML reports with a loading overlay and a busy indicator.

A large report (tens of MB of embedded data) leaves the browser blank for several
seconds while BokehJS parses the data and draws the plots. save_report renders the
page with a "loading" panel ahead of the plots and the data, so the panel shows as
soon as the HTML arrives. A script removes it once Bokeh has rendered.

The page also gets window.reportBusy(on), for slow in-page updates such as a depth
switch. CustomJS callbacks call reportBusy(true), yield a frame so the panel paints,
do their work, and call window.reportBusyDone(). Calls nest: the panel hides when
every caller has finished.

BokehJS reports no progress while it works, so the bar shows activity, not a
percentage. It is a CSS animation, which the browser keeps running while the page's
JavaScript is busy building the plots.
"""
from bokeh.embed import file_html
from bokeh.resources import CDN

# A Jinja template body. Bokeh prepends "{% extends base %}" to string templates, where
# base is its standard page. Overriding inner_body puts the overlay before the plot
# divs and the data script; super() then renders Bokeh's own content.
TEMPLATE = """
{% block inner_body %}
<div id="report-loading" role="status" aria-live="polite"
     style="position:fixed;inset:0;z-index:2147483647;display:flex;flex-direction:column;
            align-items:center;justify-content:center;gap:12px;background:rgba(252,252,251,0.97);
            font-family:system-ui,-apple-system,'Segoe UI',sans-serif;color:#0b0b0b">
  <div style="font-size:16px;font-weight:600">Loading {{ loading_label | e }}</div>
  <div style="width:320px;height:6px;border-radius:3px;background:#e6e5e0;overflow:hidden">
    <div class="report-loading-bar"
         style="width:40%;height:100%;border-radius:3px;background:#2a78d6"></div>
  </div>
  <div id="report-loading-stage" style="font-size:13px;color:#52514e">reading the data</div>
</div>
<style>
  .report-loading-bar { animation: report-loading-slide 1.2s ease-in-out infinite; }
  @keyframes report-loading-slide {
    0%   { transform: translateX(-100%); }
    100% { transform: translateX(250%); }
  }
  @media (prefers-reduced-motion: reduce) { .report-loading-bar { animation: none; } }
</style>
<div id="report-busy" role="status" aria-live="polite"
     style="position:fixed;inset:0;z-index:2147483646;display:none;flex-direction:column;
            align-items:center;justify-content:center;gap:12px;background:rgba(252,252,251,0.75);
            font-family:system-ui,-apple-system,'Segoe UI',sans-serif;color:#0b0b0b">
  <div style="font-size:15px;font-weight:600">Updating the plots</div>
  <div style="width:320px;height:6px;border-radius:3px;background:#e6e5e0;overflow:hidden">
    <div class="report-loading-bar"
         style="width:40%;height:100%;border-radius:3px;background:#2a78d6"></div>
  </div>
</div>
<script>
  (function () {
    var pending = 0;
    window.reportBusy = function (on) {
      var el = document.getElementById('report-busy');
      pending = Math.max(0, pending + (on ? 1 : -1));
      if (el) el.style.display = pending ? 'flex' : 'none';
    };
    // Hide once Bokeh has drawn the change: two animation frames after the work. The
    // timeout covers a browser that pauses animation frames (a background tab, headless
    // rendering); whichever fires first hides the panel, exactly once.
    window.reportBusyDone = function () {
      var fired = false;
      var finish = function () { if (!fired) { fired = true; window.reportBusy(false); } };
      requestAnimationFrame(function () { requestAnimationFrame(finish); });
      setTimeout(finish, 1500);
    };
  })();
  (function () {
    var stage = document.getElementById('report-loading-stage');
    var timer = setInterval(function () {
      var B = window.Bokeh;
      if (B && B.documents && B.documents.length && stage) {
        stage.textContent = 'drawing the plots';
      }
      // Bokeh.index holds one view per rendered root; it fills once drawing finishes.
      if (B && B.index && Object.keys(B.index).length) {
        clearInterval(timer);
        var el = document.getElementById('report-loading');
        if (el) el.remove();
      }
    }, 200);
  })();
</script>
{{ super() }}
{% endblock %}
"""


def save_report(layout, path, title, loading_label='report'):
    """Write layout to path as a standalone HTML page with the loading overlay.

    The page loads BokehJS from the CDN, as bokeh.io.output_file does by default.
    """
    html = file_html(layout, CDN, title, template=TEMPLATE,
                     template_variables={'loading_label': loading_label})
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(html)
