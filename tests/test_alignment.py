"""Regression tests for alignment interval coordinates."""

from pathlib import Path

import pysam
import pytest

from lib.alignment import IntervalReader
from lib import misc
from lib import library


def _write_bam(path: Path, cigar: str, *, is_reverse: bool = False) -> None:
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": 1000}],
    }
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        read = pysam.AlignedSegment()
        read.query_name = "read"
        read.query_sequence = "A" * 10
        read.query_qualities = pysam.qualitystring_to_array("I" * 10)
        read.flag = 16 if is_reverse else 0
        read.reference_id = 0
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigarstring = cigar
        read.set_tag("NH", 1)
        output.write(read)
    pysam.index(str(path))


@pytest.mark.parametrize(
    ("cigar", "is_reverse", "expected_stop", "expected_blocks"),
    [
        pytest.param("10M", False, 109, ((100, 109),), id="all-matched"),
        pytest.param("2S8M", False, 107, ((100, 107),), id="soft-clipped"),
        pytest.param(
            "4M2D6M",
            False,
            111,
            ((100, 103), (106, 111)),
            id="deletion",
        ),
        pytest.param("8M2S", True, 107, ((100, 107),), id="reverse-soft-clipped"),
    ],
)
def test_interval_reader_uses_cigar_aware_reference_end_and_query_length(
    tmp_path: Path,
    cigar: str,
    is_reverse: bool,
    expected_stop: int,
    expected_blocks: tuple[tuple[int, int], ...],
) -> None:
    bam = tmp_path / "reads.bam"
    _write_bam(bam, cigar, is_reverse=is_reverse)

    intervals, accepted_reads = IntervalReader(bam).output()

    strand = "-" if is_reverse else "+"
    assert list(intervals[("chr1", strand)]) == [
        (100, expected_stop, 10, expected_blocks)
    ]
    assert accepted_reads == {"chr1": 1}


@pytest.mark.parametrize(
    ("cigar", "expected"),
    [
        ("2S8M", [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]),
        ("4M2D6M", [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0]),
    ],
)
def test_global_metagene_counts_only_cigar_aligned_blocks(
    tmp_path: Path, cigar: str, expected: list[int]
) -> None:
    from lib import metagene

    bam = tmp_path / "reads.bam"
    _write_bam(bam, cigar)
    intervals, _ = IntervalReader(bam).output()

    coverage = metagene.metagene_mapping_start(
        {"+": {"chr1": [(102, 104)]}},
        intervals,
        positions_out_ORF=2,
        positions_in_ORF=12,
        mapping_method="global",
    )

    assert coverage["chr1"][10].tolist() == expected


def test_global_read_count_does_not_treat_a_cigar_gap_as_overlap(
    tmp_path: Path,
) -> None:
    bam = tmp_path / "reads.bam"
    _write_bam(bam, "4M2D6M")
    intervals, _ = IntervalReader(bam).output()

    assert misc.count_reads(intervals, "chr1", "+", 104, 105, "global") == 0
    assert misc.count_reads(intervals, "chr1", "+", 103, 104, "global") == 1


def test_global_track_coverage_excludes_cigar_gaps(tmp_path: Path) -> None:
    bam = tmp_path / "reads.bam"
    _write_bam(bam, "4M2D6M")

    [(reference, mappings)] = library.compute_seqid_mapping(
        bam,
        read_count_splitting=True,
        mapping_add_function=library.add_aln_mapping,
        strand_swap=False,
        clip_length=0,
    )

    assert reference == "chr1"
    assert mappings["forward"].sum() == 10
    assert mappings["forward"][100:112].tolist() == [
        1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1
    ]


def test_centered_track_normalizes_only_cigar_aligned_bases(tmp_path: Path) -> None:
    bam = tmp_path / "reads.bam"
    _write_bam(bam, "4M2D6M")

    [(_, mappings)] = library.compute_seqid_mapping(
        bam,
        read_count_splitting=True,
        mapping_add_function=library.add_centered_mapping,
        strand_swap=False,
        clip_length=1,
    )

    # One reference base is clipped from each outer end; the deletion remains
    # excluded, and the one read is distributed over eight aligned positions.
    assert mappings["forward"].sum() == pytest.approx(1.0)
    assert mappings["forward"][100:112].tolist() == pytest.approx(
        [0, 0.125, 0.125, 0.125, 0, 0, 0.125, 0.125, 0.125, 0.125, 0.125, 0]
    )
