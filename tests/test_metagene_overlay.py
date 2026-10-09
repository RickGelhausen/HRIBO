"""Overlaid metagenes preserve the reported count scale and every read length."""

import numpy as np
import pandas as pd
import pytest

import metagene_profiling
from lib import misc, plotting


@pytest.mark.parametrize("normalization", ["raw", "cpm", "window"])
@pytest.mark.parametrize("start_only", [False, True], ids=["general", "sorfs"])
def test_restored_overlay_uses_exact_normalized_counts_for_each_displayed_anchor(
    tmp_path, normalization, start_only,
):
    start_values = {25: np.arange(1, 7, dtype=float), 30: np.arange(12, 6, -1, dtype=float)}
    stop_values = {25: np.arange(7, 13, dtype=float), 30: np.arange(6, 0, -1, dtype=float)}
    start = {"chr": {**start_values, 99: np.full(6, 99.0)}}
    stop = {"chr": {**stop_values, 99: np.full(6, 99.0)}}
    if normalization == "cpm":
        start = misc.normalize_coverage(start, {"chr": 2, "other": 8})
        stop = misc.normalize_coverage(stop, {"chr": 2, "other": 8})
    figures = metagene_profiling.create_metagene_figures(
        start, stop, [25, 27, 30], tmp_path, "fiveprime", normalization,
        2, 4, ["#112233", "#abcdef"], start_only=start_only,
    )
    overlay = next(
        figure for name, _, figure in figures if name.endswith("(overlaid read lengths)")
    )
    lines = [trace for trace in overlay.data if trace.mode == "lines"]
    assert len(lines) == (3 if start_only else 6)
    assert {trace.name for trace in lines} == {"25 nt", "27 nt", "30 nt"}
    expected_label = {"raw": "Reads", "cpm": "CPM", "window": "Window-normalized counts"}[normalization]
    assert overlay.layout.yaxis.title.text == expected_label

    for anchor, xaxis, raw_values in (
        ("start", "x", start_values), ("stop", "x2", stop_values),
    ):
        if start_only and anchor == "stop":
            assert not any(trace.xaxis == xaxis for trace in lines)
            assert not (tmp_path / "fiveprime_readcounts_stop.xlsx").exists()
            continue
        workbook = pd.read_excel(
            tmp_path / f"fiveprime_readcounts_{anchor}.xlsx", sheet_name="chr"
        )
        axis_layout = overlay.layout.xaxis if anchor == "start" else overlay.layout.xaxis2
        assert anchor in axis_layout.title.text.lower()
        panel = {trace.name: trace for trace in lines if trace.xaxis == xaxis}
        assert list(panel) == ["25 nt", "27 nt", "30 nt"]
        assert [panel[f"{length} nt"].line.color for length in (25, 27, 30)] == [
            "#112233", "#abcdef", "#112233",
        ]
        for length in (25, 27, 30):
            trace = panel[f"{length} nt"]
            assert list(trace.x) == workbook["coordinates"].tolist()
            assert np.asarray(trace.y) == pytest.approx(workbook[str(length)].to_numpy())
            expected = raw_values.get(length, np.zeros(6)).copy()
            if normalization == "cpm":
                expected *= 100_000
            elif normalization == "window" and expected.any():
                expected /= expected.sum() / 6
            assert np.asarray(trace.y) == pytest.approx(expected)
            assert expected_label in trace.hovertemplate
            assert anchor in trace.hovertemplate.lower()

    if not start_only:
        # Both codons use a common count scale so their amplitudes can be compared.
        assert (overlay.layout.yaxis2.matches == "y"
                or overlay.layout.yaxis.range == overlay.layout.yaxis2.range)
        assert overlay.layout.yaxis.range[0] <= 0
        assert overlay.layout.yaxis.range[1] >= max(np.max(trace.y) for trace in lines)


def test_general_overlay_and_small_multiples_include_all_nineteen_selected_lengths(tmp_path):
    lengths = list(range(22, 41))
    coverage = {"chr": {length: np.ones(6) for length in lengths}}
    figures = metagene_profiling.create_metagene_figures(
        coverage, coverage, lengths, tmp_path, "fiveprime", "raw", 2, 4, [],
    )
    overlay = next(
        figure for name, _, figure in figures if name.endswith("(overlaid read lengths)")
    )
    multiples = next(
        figure for name, _, figure in figures if name.endswith("(per read length)")
    )
    assert len(overlay.data) == 2 * len(lengths)
    for axis in ("x", "x2"):
        assert [trace.name for trace in overlay.data if trace.xaxis == axis] == [
            f"{length} nt" for length in lengths
        ]
    assert [trace.name for trace in multiples.data] == [f"{length} nt" for length in lengths]


@pytest.mark.parametrize("start_only", [False, True], ids=["general", "sorfs"])
def test_zero_evidence_overlay_preserves_zero_lines_on_a_usable_count_axis(start_only):
    frame = pd.DataFrame({"coordinates": [-2, -1, 0, 1], "30": np.zeros(4)})
    figure = plotting.plot_metagene_profiles(
        frame, frame, [30], "No evidence", value_label="CPM", start_only=start_only,
    )
    assert figure is not None
    assert len(figure.data) == (1 if start_only else 2)
    assert all(not np.asarray(trace.y).any() for trace in figure.data)
    assert figure.layout.yaxis.range[0] <= 0 < figure.layout.yaxis.range[1]


def assert_comparable_panel_scales(figure, shared_upper, local_uppers):
    shared, independent = figure.layout.updatemenus[0].buttons
    assert shared.method == independent.method == "relayout"
    shared_args, independent_args = shared.args[0], independent.args[0]
    for index, local_upper in enumerate(local_uppers, start=1):
        axis_name = "yaxis" if index == 1 else f"yaxis{index}"
        axis = figure.layout[axis_name]
        matching_axis = None if index == 1 else "y"
        assert axis.type in (None, "linear")
        assert axis.range == pytest.approx([0, shared_upper])
        assert axis.matches == matching_axis
        assert shared_args[f"{axis_name}.range"] == pytest.approx([0, shared_upper])
        assert shared_args[f"{axis_name}.matches"] == matching_axis
        assert independent_args[f"{axis_name}.range"] == pytest.approx([0, local_upper])
        assert independent_args[f"{axis_name}.matches"] is None


def test_read_length_panels_preserve_varied_counts_with_shared_and_independent_scales():
    frame = pd.DataFrame({
        "coordinates": [-2, -1, 0, 1],
        "30": [0, 100, 50, 0], "31": [0, 2, 1, 0], "32": [0, 0, 0, 0],
    })
    figure = plotting.plot_read_length_profiles(
        frame, [30, 31, 32], "Comparable read lengths", value_label="CPM",
    )
    assert [trace.name for trace in figure.data] == ["30 nt", "31 nt", "32 nt"]
    for trace, length in zip(figure.data, [30, 31, 32]):
        assert list(trace.y) == frame[str(length)].tolist()
        assert "CPM" in trace.hovertemplate
    assert_comparable_panel_scales(figure, 110, [110, 2.2, 1])


def test_capped_read_length_panels_ignore_hidden_extreme_counts_in_all_scale_controls():
    frame = pd.DataFrame({
        "coordinates": [-2, -1, 0, 1],
        "30": [0, 100, 50, 0], "31": [0, 2, 1, 0], "32": [1_000_000] * 4,
    })
    figure = plotting.plot_read_length_profiles(
        frame, [30, 31, 32], "Two displayed lengths", max_panels=2,
    )
    assert [trace.name for trace in figure.data] == ["30 nt", "31 nt"]
    assert_comparable_panel_scales(figure, 110, [110, 2.2])
    for button in figure.layout.updatemenus[0].buttons:
        assert not any(key.startswith("yaxis3.") for key in button.args[0])
