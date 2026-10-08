"""
Metagene and P-site figures.

The previous version drew every read length as a line on one pair of shared
axes, cycling ten colours and six dash patterns. With the ten read lengths the
default configuration asks for, that is twenty overlapping traces and nothing
can be read off it.

The forms here are chosen by what each figure has to answer:

  which read lengths pile up where     -> heatmap, read length against position
  what one read length looks like      -> small multiples, one panel each
  is the signal in frame                -> grouped bars, three reading frames
  how deep is each read length          -> bar chart

Colour comes from lib.theme; nothing here contains a literal colour.

Author: Rick Gelhausen
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

from lib import theme

ENRICHMENT_LABEL = "Enrichment<br>(x own background)"

# Colour saturates at a high quantile of the data rather than a fixed value, so
# that both a lone tall initiation peak and a weaker periodic body stay
# distinguishable. Bounded so a flat profile cannot produce a misleadingly
# sensitive scale, nor one peak an insensitive one.
ENRICHMENT_CEILING_BOUNDS = (4.0, 25.0)
ENRICHMENT_CEILING_QUANTILE = 0.995


def _enrichment_ceiling(*matrices) -> float:
    values = [m for m in matrices if m.size]
    if not values:
        return ENRICHMENT_CEILING_BOUNDS[0]
    pooled = np.concatenate([m.ravel() for m in values])
    pooled = pooled[np.isfinite(pooled)]
    if not pooled.size:
        return ENRICHMENT_CEILING_BOUNDS[0]
    ceiling = float(np.quantile(pooled, ENRICHMENT_CEILING_QUANTILE))
    return float(np.clip(ceiling, *ENRICHMENT_CEILING_BOUNDS))


def _row_enrichment(matrix: np.ndarray) -> np.ndarray:
    """Each row as a multiple of its own median.

    Normalising each row by its maximum, the obvious choice, inverts the figure:
    a read length that is pure noise has its noise stretched to full scale, while
    a read length with a genuine initiation peak has everything else compressed
    to white by that peak. Dividing by the row median instead measures the thing
    that actually distinguishes them -- how far anything rises above that read
    length's own background.
    """
    matrix = np.asarray(matrix, dtype=float)
    backgrounds = np.median(matrix, axis=1, keepdims=True)
    # A row that is mostly zeros has a zero median; fall back to its mean so the
    # division stays finite, and to 1 if the row is empty outright.
    fallback = matrix.mean(axis=1, keepdims=True)
    backgrounds = np.where(backgrounds > 0, backgrounds, fallback)
    backgrounds = np.where(backgrounds > 0, backgrounds, 1.0)
    return matrix / backgrounds


def _profile_matrix(frame, read_lengths):
    """Rows of a metagene dataframe as a matrix, in read length order."""
    columns = [str(length) for length in read_lengths if str(length) in frame.columns]
    lengths = [int(column) for column in columns]
    if not columns:
        return np.zeros((0, len(frame))), []
    return np.vstack([frame[column].to_numpy(dtype=float) for column in columns]), lengths


def plot_metagene_heatmap(
    df_start,
    df_stop,
    read_lengths,
    title,
    subtitle="",
    offsets=None,
    read_end="fiveprime",
    anchor="start",
    site="P",
    background_references=None,
):
    """Read length against position, for the start and stop codon windows.

    This is the figure that actually answers "which read lengths should I use":
    a read length with a clean initiation signal shows as a bright vertical band
    a fixed distance from the start codon, and the distance is the offset.
    """
    start_matrix, lengths = _profile_matrix(df_start, read_lengths)
    stop_matrix, _ = _profile_matrix(df_stop, read_lengths)

    start_enrichment = _row_enrichment(start_matrix)
    stop_enrichment = _row_enrichment(stop_matrix)
    if background_references is not None:
        calibrated_matrix = stop_matrix if anchor == "stop" else start_matrix
        references = np.array([background_references.get(length, 0.0) for length in lengths])[:, None]
        calibrated_enrichment = np.divide(
            calibrated_matrix, references, out=np.full_like(calibrated_matrix, np.nan), where=references > 0
        )
        if anchor == "stop":
            stop_enrichment = calibrated_enrichment
        else:
            start_enrichment = calibrated_enrichment
    ceiling = _enrichment_ceiling(start_enrichment, stop_enrichment)

    fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=tuple(
            f"Relative to {boundary} codon" + (
                " (scored)" if boundary == anchor else " (diagnostic)"
            ) if background_references is not None else f"Relative to {boundary} codon"
            for boundary in ("start", "stop")
        ),
        horizontal_spacing=0.08,
        shared_yaxes=True,
    )

    for column, (matrix, frame, panel_anchor) in enumerate(
        ((start_enrichment, df_start, "start"), (stop_enrichment, df_stop, "stop")), start=1
    ):
        if matrix.size == 0:
            continue
        raw_matrix = start_matrix if panel_anchor == "start" else stop_matrix
        background_label = "upstream site background" if (
            background_references is not None and panel_anchor == anchor
        ) else "whole-window median (diagnostic)" if background_references is not None else "this length's background"
        fig.add_trace(
            go.Heatmap(
                z=matrix,
                x=frame["coordinates"],
                y=lengths,
                customdata=raw_matrix,
                colorscale=theme.SEQUENTIAL_BLUE,
                zmin=1,
                zmax=ceiling,
                showscale=column == 2,
                colorbar=dict(
                    title=dict(text=ENRICHMENT_LABEL if background_references is None else
                               "Enrichment<br>(x panel background)", font=dict(size=10)),
                    thickness=12,
                    len=0.85,
                    tickfont=dict(size=9),
                ),
                hovertemplate=(
                    f"read length %{{y}} nt<br>%{{x}} nt from {panel_anchor} codon"
                    f"<br>%{{z:.1f}}x {background_label}"
                    "<br>%{customdata:.0f} aligned read ends<extra></extra>"
                ),
            ),
            row=1,
            col=column,
        )
        # The codon itself, as a reference line rather than a gridline.
        fig.add_vline(
            x=0,
            line=dict(color=theme.INK_SECONDARY, width=1, dash="dot"),
            row=1,
            col=column,
        )

    # Mark the calibrated site offset on the corresponding boundary panel.
    if offsets and read_end in {"fiveprime", "threeprime"}:
        marked = [(length, offset) for length, offset in sorted(offsets.items()) if offset is not None]
        if marked:
            # The mapped end lies upstream for 5' offsets and downstream for 3'.
            direction = -1 if read_end == "fiveprime" else 1
            fig.add_trace(
                go.Scatter(
                    x=[direction * offset for _, offset in marked],
                    y=[length for length, _ in marked],
                    mode="markers",
                    marker=dict(
                        symbol="line-ns",
                        size=11,
                        line=dict(color=theme.STATUS["serious"], width=2),
                    ),
                    name=f"estimated {site}-site offset",
                    hovertemplate="read length %{y} nt<br>offset %{customdata} nt<extra></extra>",
                    customdata=[offset for _, offset in marked],
                ),
                row=1,
                col=1 if anchor == "start" else 2,
            )

    fig.update_xaxes(title_text="Distance from start codon (nt)", row=1, col=1)
    fig.update_xaxes(title_text="Distance from stop codon (nt)", row=1, col=2)
    fig.update_yaxes(title_text="Read length (nt)", dtick=1, row=1, col=1)
    fig.update_annotations(font_size=12)

    theme.apply(fig, title, subtitle)
    fig.update_layout(
        height=max(320, 26 * len(lengths) + 190),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def plot_read_length_profiles(
    df_start,
    read_lengths,
    title,
    subtitle="",
    offsets=None,
    read_end="fiveprime",
    color_list=None,
    max_panels=12,
    candidate_support=None,
    value_label="Reads",
):
    """One panel per read length, sharing an x axis.

    Small multiples rather than overlaid lines: comparing shapes across panels is
    what the eye is good at, and it does not need a colour per read length.

    User-supplied series colours are assigned in read-length order and cycle
    deterministically when fewer colours than panels are supplied.  An empty
    list falls back to the fixed colour-vision-deficiency-safe theme palette.
    """
    matrix, lengths = _profile_matrix(df_start, read_lengths)
    if matrix.size == 0:
        return None

    lengths = lengths[:max_panels]
    matrix = matrix[: len(lengths)]
    palette = list(color_list) if color_list else theme.CATEGORICAL
    colors_by_length = {
        int(read_length): palette[index % len(palette)]
        for index, read_length in enumerate(read_lengths)
    }

    fig = make_subplots(
        rows=len(lengths),
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.012,
        subplot_titles=[
            f"{length} nt" + (
                f" · {candidate_support[length]['contributing_cds']} contributing CDSs"
                f" · {candidate_support[length]['raw_count_contributions']} raw count contributions"
                if candidate_support is not None and length in candidate_support else ""
            )
            for length in lengths
        ],
    )

    for index, (length, row_values) in enumerate(zip(lengths, matrix), start=1):
        color = colors_by_length[length]
        fig.add_trace(
            go.Scatter(
                x=df_start["coordinates"],
                y=row_values,
                mode="lines",
                line=dict(color=color, width=1.5),
                fill="tozeroy",
                name=f"{length} nt",
                showlegend=False,
                hovertemplate=f"%{{x}} nt from start<br>%{{y:.3g}} {value_label}<extra></extra>",
            ),
            row=index,
            col=1,
        )
        fig.add_vline(
            x=0, line=dict(color=theme.INK_SECONDARY, width=1, dash="dot"), row=index, col=1
        )
        if (
            offsets
            and read_end in {"fiveprime", "threeprime"}
            and offsets.get(length) is not None
        ):
            direction = -1 if read_end == "fiveprime" else 1
            fig.add_vline(
                x=direction * offsets[length],
                line=dict(color=theme.STATUS["serious"], width=1.5),
                row=index,
                col=1,
            )

    fig.update_xaxes(title_text="Distance from start codon (nt)", row=len(lengths), col=1)
    fig.update_yaxes(title_text=value_label, row=max(1, len(lengths) // 2), col=1)
    fig.update_annotations(font_size=10)

    theme.apply(fig, title, subtitle)
    fig.update_layout(height=max(300, 95 * len(lengths)))
    return fig


def plot_frame_composition(scores, title, subtitle="", site="P"):
    """Share of calibrated sites in each reading frame, per read length.

    Grouped rather than stacked: the question is whether one frame stands above
    the 1/3 line, and grouped bars put every frame on the same baseline. The
    reference line at 1/3 is what "no periodicity" looks like.
    """
    usable = [s for s in scores if s.offset is not None and any(s.frame_fractions)]
    if not usable:
        return None

    # Categorical x: only read lengths with an estimated offset appear, and they
    # sit side by side instead of leaving gaps that read as missing data.
    lengths = [str(s.read_length) for s in usable]
    fig = go.Figure()

    for frame in range(3):
        values = [s.frame_fractions[frame] for s in usable]
        fig.add_trace(
            go.Bar(
                x=lengths,
                y=values,
                name=f"frame {frame}",
                marker=dict(
                    color=theme.FRAME_COLORS[frame],
                    line=dict(color=theme.SURFACE, width=2),
                ),
                # Direct labels: the aqua slot sits below 3:1 against the
                # surface, so the values must be readable without the colour.
                text=[f"{value:.0%}" for value in values],
                textposition="outside",
                textfont=dict(size=9, color=theme.INK_SECONDARY),
                hovertemplate=f"read length %{{x}} nt<br>frame {frame}: %{{y:.1%}}<extra></extra>",
            )
        )

    fig.add_hline(
        y=1 / 3,
        line=dict(color=theme.INK_MUTED, width=1, dash="dash"),
        annotation=dict(
            text="no frame preference",
            font=dict(size=9, color=theme.INK_MUTED),
            xanchor="left",
            yanchor="bottom",
        ),
        annotation_position="top left",
    )

    fig.update_xaxes(title_text="Read length (nt)", type="category")
    fig.update_yaxes(title_text=f"Share of {site}-sites", tickformat=".0%", range=[0, 1.08])

    theme.apply(fig, title, subtitle)
    fig.update_layout(
        barmode="group",
        bargap=0.25,
        bargroupgap=0.05,
        height=380,
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def plot_read_length_distribution(scores, title, subtitle="", anchor="start"):
    """Library composition by read length, with the usable ones picked out."""
    if not scores:
        return None

    # Two traces rather than one with per-bar colours, so that usability is
    # carried by a legend entry and not by colour alone.
    groups = (
        ("usable for TIS calling" if anchor == "start" else "usable for TTS peaks",
         True, theme.CATEGORICAL[0]),
        ("not usable", False, theme.INK_MUTED),
    )

    fig = go.Figure()
    for name, wanted, color in groups:
        selected = [s for s in scores if s.usable is wanted]
        if not selected:
            continue
        fig.add_trace(
            go.Bar(
                x=[s.read_length for s in selected],
                y=[s.abundance for s in selected],
                name=name,
                marker=dict(color=color, line=dict(color=theme.SURFACE, width=2)),
                hovertemplate=(
                    f"read length %{{x}} nt<br>%{{y:.1%}} of reads<br>{name}<extra></extra>"
                ),
            )
        )

    fig.update_xaxes(title_text="Read length (nt)", dtick=1)
    fig.update_yaxes(title_text="Share of reads", tickformat=".0%")

    theme.apply(fig, title, subtitle)
    fig.update_layout(
        height=320,
        barmode="overlay",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
    )
    return fig


def plot_pooled_profile(pooled, coordinates, title, subtitle="", anchor="start", site="P"):
    """The offset-corrected, pooled profile of the selected boundary."""
    fig = go.Figure(
        go.Scatter(
            x=coordinates,
            y=pooled,
            mode="lines",
            line=dict(color=theme.CATEGORICAL[0], width=2),
            fill="tozeroy",
            fillcolor="rgba(42,120,214,0.12)",
            name=f"pooled {site}-sites",
            showlegend=False,
            hovertemplate=f"%{{x}} nt from {anchor}<br>%{{y:.0f}} {site}-sites<extra></extra>",
        )
    )
    fig.add_vline(
        x=0,
        line=dict(color=theme.STATUS["serious"], width=1.5),
        annotation=dict(text=f"{anchor} codon", font=dict(size=9, color=theme.STATUS["serious"])),
        annotation_position="top right",
    )
    fig.update_xaxes(title_text=f"Distance from {anchor} codon (nt), after offset correction")
    fig.update_yaxes(title_text=f"Pooled {site}-sites")

    theme.apply(fig, title, subtitle)
    fig.update_layout(height=340)
    return fig
