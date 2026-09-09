"""Regression tests for metagene parsing, orientation, and profile handling."""

import numpy as np
import pandas as pd
import pytest

import lib.io as io
import lib.annotation as annotation
import lib.metagene as metagene
import lib.misc as misc
import metagene_profiling
import tis_advisor


class IntervalCollection:
    """Small inclusive-interval index matching the ``InterLap.find`` interface."""

    def __init__(self, intervals):
        self.intervals = intervals

    def find(self, bounds):
        start, stop = bounds
        return [
            interval
            for interval in self.intervals
            if interval[0] <= stop and interval[1] >= start
        ]


def test_metagene_cpm_uses_the_complete_two_contig_library_total():
    coverage = {
        "chrA": {30: np.array([1.0, 2.0])},
        "chrB": {30: np.array([3.0, 4.0])},
    }

    normalized = misc.normalize_coverage(
        coverage, {"chrA": 2, "chrB": 8}
    )

    assert normalized["chrA"][30] == pytest.approx([100_000, 200_000])
    assert normalized["chrB"][30] == pytest.approx([300_000, 400_000])


def test_shared_metagene_tis_rpkm_filter_uses_the_complete_library(tmp_path):
    annotation_path = tmp_path / "annotation.gff"
    annotation_path.write_text(
        "chrA\ttest\tCDS\t101\t200\t.\t+\t0\tID=gene\n"
    )
    reads = {
        ("chrA", "+"): IntervalCollection([(100, 120, 21)]),
    }

    starts, stops = annotation.retrieve_annotation_positions(
        annotation_path,
        reads,
        {"chrA": 1, "chrB": 99},
        {"chrA": 1000, "chrB": 1000},
        ["rpkm"],
        "fiveprime",
        200_000,
        0,
        10,
        20,
    )

    # 1 read / (100 nt * 100 library reads) = 100,000 RPKM.  The old
    # contig-local denominator was 1 and incorrectly retained this feature.
    assert starts == {"-": {}, "+": {}}
    assert stops == {"-": {}, "+": {}}


@pytest.mark.parametrize(
    "mapping_method, expected",
    [
        ("fiveprime", [0, 1, 0, 0, 0, 0]),
        ("threeprime", [0, 0, 0, 1, 0, 0]),
        ("centered", [0, 0, 1, 0, 0, 0]),
        ("global", [0, 1, 1, 1, 0, 0]),
    ],
)
@pytest.mark.parametrize("strand", ["+", "-"])
@pytest.mark.parametrize(
    "mapper, codon, reads",
    [
        (
            metagene.metagene_mapping_start,
            (100, 102),
            {"+": (99, 101, 3), "-": (101, 103, 3)},
        ),
        (
            metagene.metagene_mapping_stop,
            (300, 302),
            {"+": (300, 302, 3), "-": (300, 302, 3)},
        ),
    ],
    ids=["start", "stop"],
)
def test_metagene_mapping_is_transcript_oriented(
    mapper, codon, reads, strand, mapping_method, expected
):
    """Every anchor/strand combination must produce the same oriented profile."""
    positions_out_ORF = 2
    positions_in_ORF = 4
    chromosome = "chr"

    coverage = mapper(
        {strand: {chromosome: [codon]}},
        {(chromosome, strand): IntervalCollection([reads[strand]])},
        positions_out_ORF,
        positions_in_ORF,
        mapping_method,
    )

    assert coverage[chromosome][3].tolist() == expected


def test_tis_advisor_uses_the_stop_profile_coordinate_axis():
    start, stop = tis_advisor.metagene_coordinates(2, 4)

    assert start.tolist() == [-2, -1, 0, 1, 2, 3]
    assert stop.tolist() == [-4, -3, -2, -1, 0, 1]


@pytest.mark.parametrize(
    "mapper",
    [metagene.metagene_mapping_start, metagene.metagene_mapping_stop],
    ids=["start", "stop"],
)
@pytest.mark.parametrize("strand", ["+", "-"])
def test_metagene_mapping_tolerates_a_missing_read_index(mapper, strand):
    coverage = mapper(
        {strand: {"chr": [(100, 102)]}},
        {},
        2,
        4,
        "global",
    )

    assert coverage == {"chr": {}}


def test_equalize_dictionary_keys_accepts_fully_empty_evidence():
    start, stop = misc.equalize_dictionary_keys({}, {}, 2, 4)

    assert start == {}
    assert stop == {}


@pytest.mark.parametrize("evidence_side", ["start", "stop"])
def test_equalize_dictionary_keys_fills_the_empty_evidence_side(evidence_side):
    evidence = {"chr": {30: np.ones(6, dtype=np.intp)}}
    start = evidence if evidence_side == "start" else {}
    stop = evidence if evidence_side == "stop" else {}

    start, stop = misc.equalize_dictionary_keys(start, stop, 2, 4)

    assert set(start) == set(stop) == {"chr"}
    assert set(start["chr"]) == set(stop["chr"]) == {30}
    empty_side = stop if evidence_side == "start" else start
    assert empty_side["chr"][30].tolist() == [0, 0, 0, 0, 0, 0]


def test_empty_metagene_evidence_has_an_explicit_dataframe():
    frames = misc.create_data_frame({}, 2, 4, "start")

    assert set(frames) == {"no_evidence"}
    assert frames["no_evidence"]["coordinates"].tolist() == [-2, -1, 0, 1, 2, 3]


def test_empty_metagene_evidence_writes_valid_workbooks(tmp_path):
    figures = metagene_profiling.create_metagene_figures(
        {},
        {},
        [30],
        tmp_path,
        "fiveprime",
        "raw",
        2,
        4,
        [],
    )

    start = pd.read_excel(
        tmp_path / "fiveprime_readcounts_start.xlsx",
        sheet_name="no_evidence",
    )
    stop = pd.read_excel(
        tmp_path / "fiveprime_readcounts_stop.xlsx",
        sheet_name="no_evidence",
    )

    assert len(figures) == 2
    assert start["coordinates"].tolist() == [-2, -1, 0, 1, 2, 3]
    assert stop["coordinates"].tolist() == [-4, -3, -2, -1, 0, 1]
    assert list(start.columns) == ["coordinates", "30", "sum"]
    assert list(stop.columns) == ["coordinates", "30", "sum"]
    assert start["30"].tolist() == [0, 0, 0, 0, 0, 0]
    assert stop["30"].tolist() == [0, 0, 0, 0, 0, 0]
    assert start["sum"].tolist() == [0, 0, 0, 0, 0, 0]
    assert stop["sum"].tolist() == [0, 0, 0, 0, 0, 0]


@pytest.mark.parametrize("normalization_method", ["raw", "window"])
@pytest.mark.parametrize("evidence_side", ["start", "stop"])
def test_one_sided_metagene_evidence_stays_finite(
    normalization_method, evidence_side, tmp_path
):
    observed = np.arange(1, 7, dtype=np.intp)
    evidence = {"chr": {30: observed.copy()}}
    start = evidence if evidence_side == "start" else {}
    stop = evidence if evidence_side == "stop" else {}

    metagene_profiling.create_metagene_figures(
        start,
        stop,
        [30],
        tmp_path,
        "fiveprime",
        normalization_method,
        2,
        4,
        [],
    )

    workbooks = {
        anchor: pd.read_excel(
            tmp_path / f"fiveprime_readcounts_{anchor}.xlsx",
            sheet_name="chr",
        )
        for anchor in ("start", "stop")
    }
    empty_side = "stop" if evidence_side == "start" else "start"
    assert np.isfinite(workbooks[empty_side]["30"]).all()
    assert workbooks[empty_side]["30"].tolist() == [0, 0, 0, 0, 0, 0]
    assert workbooks[empty_side]["sum"].tolist() == [0, 0, 0, 0, 0, 0]

    expected_observed = observed.astype(float)
    if normalization_method == "window":
        expected_observed /= observed.sum() / len(observed)
    assert np.allclose(workbooks[evidence_side]["30"], expected_observed)


@pytest.mark.parametrize(
    "spec, expected",
    [
        ("25-34", [25, 26, 27, 28, 29, 30, 31, 32, 33, 34]),
        ("30", [30]),
        ("22,23,27,34-35", [22, 23, 27, 34, 35]),
        ("34-35,22", [22, 34, 35]),
        # Reversed intervals are accepted and normalised.
        ("35-34", [34, 35]),
        # Duplicates collapse.
        ("30,30,30", [30]),
    ],
)
def test_parse_read_lengths(spec, expected):
    assert io.parse_read_lengths(spec) == expected


def test_parse_read_lengths_returns_integers():
    """Mixing str and int previously made the result unsortable."""
    assert all(isinstance(value, int) for value in io.parse_read_lengths("22,23,34-35"))


def test_equalize_dictionary_keys_balances_both_sides():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert set(start) == set(stop) == {"a", "b"}
    for coverage in (start, stop):
        for chromosome in coverage:
            assert set(coverage[chromosome]) == {30, 31}


def test_equalize_dictionary_keys_preserves_existing_data():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert start["a"][30][0] == 1
    assert stop["b"][31][0] == 1


def test_equalize_dictionary_keys_does_not_alias_windows():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)
    start["a"][31][0] = 99

    assert stop["a"][31][0] == 0, "filled windows are shared across dictionaries"
    assert start["b"][31][0] == 0, "filled windows are shared across chromosomes"
    assert start["b"][30][0] == 0, "filled windows are shared across read lengths"


def test_equalize_dictionary_keys_window_length():
    start = {"a": {30: np.ones(10, dtype=np.intp)}}
    stop = {"b": {31: np.ones(10, dtype=np.intp)}}

    start, stop = misc.equalize_dictionary_keys(start, stop, 4, 6)

    assert len(stop["a"][31]) == 4 + 6


def test_configured_read_lengths_exclude_observed_lengths_and_keep_exact_gaps():
    start = {
        "chr": {
            25: np.ones(10, dtype=np.intp),
            26: np.full(10, 2, dtype=np.intp),
            99: np.full(10, 3, dtype=np.intp),
        }
    }
    stop = {"chr": {25: np.ones(10, dtype=np.intp)}}

    start = misc.retain_read_lengths(start, [25, 27])
    stop = misc.retain_read_lengths(stop, [25, 27])
    start, stop = misc.equalize_dictionary_keys(
        start, stop, 4, 6, read_lengths=[25, 27]
    )

    assert set(start["chr"]) == {25, 27}
    assert set(stop["chr"]) == {25, 27}
    assert start["chr"][25].tolist() == [1] * 10
    assert not start["chr"][27].any()
    assert not stop["chr"][27].any()
