"""
A single plotly theme for every figure HRIBO produces.

Having one place to change colour, typography and grid styling is what keeps the
figures looking like one set rather than a pile of defaults, and it is why the
plotting code below carries no literal colours.

The palette is colour-vision-deficiency safe. The categorical order is fixed:
slot 1 is always the first series, never cycled by rank, so a colour keeps
meaning the same thing when a series is filtered out. The sequential ramp is a
single hue, light to dark, as magnitude encodings require.

The figures commit to a light surface. They are written to SVG and PDF for
publication as well as to HTML, and a single deliberate look travels better
across those than an automatic dark flip would.

Author: Rick Gelhausen
"""

from __future__ import annotations

import plotly.graph_objects as go

# --------------------------------------------------------------------------
# Palette
# --------------------------------------------------------------------------

SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8984"
GRID = "#e8e7e3"
AXIS = "#c9c8c3"

# Fixed categorical order. Validated for colour-vision deficiency; the first
# three slots are safe for all-pairs comparison, which is why reading frames
# (the only place three colours sit side by side) use exactly those.
CATEGORICAL = [
    "#2a78d6",  # blue
    "#eb6834",  # orange
    "#1baf7a",  # aqua
    "#eda100",  # yellow
    "#e87ba4",  # magenta
    "#008300",  # green
    "#4a3aa7",  # violet
    "#e34948",  # red
]

FRAME_COLORS = CATEGORICAL[:3]

# Single-hue sequential ramp for magnitude, light to dark.
SEQUENTIAL_BLUE = [
    [0.00, "#fcfcfb"],
    [0.08, "#cde2fb"],
    [0.20, "#9ec5f4"],
    [0.35, "#6da7ec"],
    [0.50, "#3987e5"],
    [0.65, "#2a78d6"],
    [0.80, "#1c5cab"],
    [0.90, "#184f95"],
    [1.00, "#0d366b"],
]

# Reserved for state, never for a series.
STATUS = {
    "good": "#1baf7a",
    "warning": "#eda100",
    "serious": "#eb6834",
    "critical": "#e34948",
}

CONFIDENCE_COLORS = {
    "high": STATUS["good"],
    "medium": STATUS["warning"],
    "low": STATUS["serious"],
    "none": STATUS["critical"],
}

FONT_FAMILY = "Helvetica Neue, Helvetica, Arial, sans-serif"


def template() -> go.layout.Template:
    """The HRIBO plotly template."""
    return go.layout.Template(
        layout=go.Layout(
            font=dict(family=FONT_FAMILY, size=12, color=INK_PRIMARY),
            title=dict(font=dict(size=16, color=INK_PRIMARY), x=0.01, xanchor="left"),
            paper_bgcolor=SURFACE,
            plot_bgcolor=SURFACE,
            colorway=CATEGORICAL,
            margin=dict(l=70, r=30, t=70, b=60),
            xaxis=dict(
                showgrid=True,
                gridcolor=GRID,
                gridwidth=1,
                zeroline=False,
                linecolor=AXIS,
                ticks="outside",
                tickcolor=AXIS,
                ticklen=4,
                tickfont=dict(size=10, color=INK_SECONDARY),
                title=dict(font=dict(size=12, color=INK_SECONDARY)),
            ),
            yaxis=dict(
                showgrid=True,
                gridcolor=GRID,
                gridwidth=1,
                zeroline=False,
                linecolor=AXIS,
                ticks="outside",
                tickcolor=AXIS,
                ticklen=4,
                tickfont=dict(size=10, color=INK_SECONDARY),
                title=dict(font=dict(size=12, color=INK_SECONDARY)),
            ),
            legend=dict(
                font=dict(size=10, color=INK_SECONDARY),
                bgcolor="rgba(0,0,0,0)",
                borderwidth=0,
            ),
            hoverlabel=dict(font=dict(family=FONT_FAMILY, size=11), namelength=-1),
        )
    )


def apply(fig: go.Figure, title: str = "", subtitle: str = "") -> go.Figure:
    """Apply the template, and set a title with an optional muted subtitle."""
    fig.update_layout(template=template())
    if title:
        text = title
        if subtitle:
            text = f"{title}<br><span style='font-size:11px;color:{INK_MUTED}'>{subtitle}</span>"
        fig.update_layout(title=dict(text=text))
    return fig


def series_color(index: int) -> str:
    """Categorical colour for a series, by fixed position rather than by rank."""
    return CATEGORICAL[index % len(CATEGORICAL)]


# --------------------------------------------------------------------------
# Page shell for the interactive reports
# --------------------------------------------------------------------------

PAGE_CSS = f"""
:root {{
  --surface: {SURFACE};
  --ink: {INK_PRIMARY};
  --ink-secondary: {INK_SECONDARY};
  --ink-muted: {INK_MUTED};
  --grid: {GRID};
}}
* {{ box-sizing: border-box; }}
body {{
  margin: 0;
  padding: 0 0 4rem;
  background: var(--surface);
  color: var(--ink);
  font-family: {FONT_FAMILY};
  font-size: 15px;
  line-height: 1.55;
}}
.page {{ max-width: 1180px; margin: 0 auto; padding: 0 1.5rem; }}
header {{ padding: 2.5rem 0 1.25rem; border-bottom: 1px solid var(--grid); margin-bottom: 2rem; }}
h1 {{ font-size: 1.6rem; font-weight: 600; margin: 0 0 .35rem; letter-spacing: -0.01em; }}
h2 {{ font-size: 1.15rem; font-weight: 600; margin: 2.5rem 0 .5rem; }}
h3 {{ font-size: .95rem; font-weight: 600; margin: 1.75rem 0 .5rem; color: var(--ink-secondary); }}
p  {{ margin: .5rem 0; max-width: 74ch; }}
.subtitle {{ color: var(--ink-muted); font-size: .95rem; margin: 0; }}
.plot {{ margin: 1rem 0 2rem; }}
table {{ border-collapse: collapse; font-size: .85rem; margin: 1rem 0; width: 100%; }}
th, td {{ text-align: left; padding: .4rem .7rem; border-bottom: 1px solid var(--grid); }}
th {{ font-weight: 600; color: var(--ink-secondary); white-space: nowrap; }}
td.num {{ text-align: right; font-variant-numeric: tabular-nums; }}
.table-wrap {{ overflow-x: auto; }}
.verdict {{ border-left: 3px solid var(--grid); padding: .6rem 0 .6rem 1rem; margin: 1rem 0; }}
.verdict strong {{ display: block; font-size: 1.05rem; }}
.pill {{
  display: inline-block; padding: .1rem .55rem; border-radius: 999px;
  font-size: .75rem; font-weight: 600; color: #fff; vertical-align: middle;
}}
.warn {{ color: var(--ink-secondary); font-size: .9rem; margin: .35rem 0; }}
.warn::before {{ content: "!"; display: inline-block; width: 1.1em; font-weight: 700; color: {STATUS['warning']}; }}
code {{ background: #f2f1ee; padding: .1rem .35rem; border-radius: 3px; font-size: .85em; }}
pre {{ background: #f2f1ee; padding: .9rem 1.1rem; border-radius: 5px; overflow-x: auto; font-size: .82rem; }}
pre code {{ background: none; padding: 0; }}
"""


def page(title: str, subtitle: str, body: str) -> str:
    """Wrap rendered figures and tables in a standalone HTML page."""
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>{PAGE_CSS}</style>
</head>
<body>
<div class="page">
<header>
<h1>{title}</h1>
<p class="subtitle">{subtitle}</p>
</header>
{body}
</div>
</body>
</html>
"""
