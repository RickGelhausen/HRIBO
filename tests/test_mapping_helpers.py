"""Focused regressions for legacy mapping helper scripts."""

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pysam
import pytest

from lib import library
from lib import misc


REPO = Path(__file__).resolve().parent.parent
TOTAL_MAPPED_READS = REPO / "workflow" / "scripts" / "total_mapped_reads.py"
SAM_STRAND_INVERTER = REPO / "workflow" / "scripts" / "sam_strand_inverter.py"


def _write_alignment(handle, name, reference_id, start, sequence, flag=0, nh=None):
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = sequence
    alignment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    alignment.flag = flag
    alignment.reference_id = reference_id
    alignment.reference_start = start
    alignment.mapping_quality = 255 if reference_id >= 0 else 0
    alignment.cigartuples = [(0, len(sequence))] if reference_id >= 0 else None
    if nh is not None:
        alignment.set_tag("NH", nh)
    handle.write(alignment)


def test_total_mapped_reads_writes_exact_counts_and_average_lengths(tmp_path):
    bam = tmp_path / "sample.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrA", "LN": 1000}, {"SN": "chrB", "LN": 500}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        _write_alignment(handle, "a-20", 0, 10, "A" * 20)
        _write_alignment(handle, "a-30-reverse", 0, 100, "C" * 30, flag=16)
        _write_alignment(handle, "b-15", 1, 20, "G" * 15)
        _write_alignment(handle, "unmapped", -1, -1, "T" * 10, flag=4)
    pysam.index(str(bam))

    mapped = tmp_path / "mapped.txt"
    average_lengths = tmp_path / "average_lengths.txt"
    subprocess.run(
        [
            sys.executable,
            str(TOTAL_MAPPED_READS),
            "-b",
            str(bam),
            "-m",
            str(mapped),
            "-l",
            str(average_lengths),
        ],
        check=True,
    )

    assert mapped.read_text() == "sample\tchrA\t2\nsample\tchrB\t1\n"
    assert average_lengths.read_text() == "sample\tchrA\t25.0\nsample\tchrB\t15.0\n"


def test_total_mapped_reads_weights_multimappers_and_lengths_by_inverse_nh(
    tmp_path,
):
    bam = tmp_path / "fractional.bam"
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrA", "LN": 1000}, {"SN": "chrB", "LN": 500}],
    }
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        _write_alignment(handle, "unique", 0, 10, "A" * 20, nh=1)
        _write_alignment(handle, "multi", 0, 100, "C" * 30, nh=2)
        _write_alignment(handle, "multi", 1, 20, "C" * 30, nh=2)
        # Older aligners may omit NH for unique reads; that means NH=1.
        _write_alignment(handle, "implicit-unique", 1, 100, "G" * 10)
    pysam.index(str(bam))

    mapped = tmp_path / "mapped.txt"
    average_lengths = tmp_path / "average_lengths.txt"
    subprocess.run(
        [
            sys.executable,
            str(TOTAL_MAPPED_READS),
            "-b",
            str(bam),
            "-m",
            str(mapped),
            "-l",
            str(average_lengths),
        ],
        check=True,
    )

    assert mapped.read_text() == (
        "fractional\tchrA\t1.5\n"
        "fractional\tchrB\t1.5\n"
    )
    average_rows = {
        contig: float(value)
        for _, contig, value in csv.reader(
            average_lengths.open(newline=""), delimiter="\t"
        )
    }
    assert average_rows == pytest.approx(
        {"chrA": 70 / 3, "chrB": 50 / 3}
    )


@pytest.mark.parametrize("nh", [0, -1, 1.5])
def test_total_mapped_reads_rejects_invalid_nh(nh, tmp_path):
    bam = tmp_path / f"invalid-{nh}.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr", "LN": 1000}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as handle:
        _write_alignment(handle, "invalid", 0, 10, "A" * 20, nh=nh)
    pysam.index(str(bam))

    result = subprocess.run(
        [
            sys.executable,
            str(TOTAL_MAPPED_READS),
            "-b",
            str(bam),
            "-m",
            str(tmp_path / "mapped.txt"),
            "-l",
            str(tmp_path / "lengths.txt"),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "invalid NH tag" in result.stderr
    assert "positive integer" in result.stderr


def test_mapped_read_summary_aggregates_fractional_counts_by_library(tmp_path):
    summary = tmp_path / "mapped.txt"
    summary.write_text(
        "sample\tchrA\t1.5\n"
        "sample\tchrB\t2.25\n"
        "other\tchrA\t1\n"
        "other\tchrB\t1\n"
    )

    by_contig, totals = misc.read_mapped_read_summary(summary)

    assert by_contig[("sample", "chrA")] == 1.5
    assert totals == {"sample": 3.75, "other": 2.0}
    assert library.get_read_count_dict(summary, "sample") == (3.75, 2.0)

    with pytest.raises(ValueError, match="no mapped-read rows for library 'missing'"):
        library.get_read_count_dict(summary, "missing")


@pytest.mark.parametrize(
    "contents, message",
    [
        ("sample\tchr\n", "expected three tab-separated fields"),
        ("sample\tchr\tnan\n", "must be finite"),
        ("sample\tchr\t-1\n", "must not be negative"),
        ("sample\tchr\t0\n", "zero mapped reads"),
        ("sample\tchr\t1\nsample\tchr\t2\n", "duplicate mapped-read row"),
    ],
)
def test_mapped_read_summary_rejects_invalid_normalization_totals(
    contents, message, tmp_path
):
    summary = tmp_path / "invalid.txt"
    summary.write_text(contents)

    with pytest.raises(ValueError, match=message):
        misc.read_mapped_read_summary(summary)


def test_tracks_use_one_library_wide_factor_for_every_contig(monkeypatch, tmp_path):
    summary = tmp_path / "mapped.txt"
    summary.write_text(
        "sample\tchrA\t3\n"
        "sample\tchrB\t1\n"
        "other\tchrA\t1\n"
        "other\tchrB\t1\n"
    )
    library_total, minimum_total = library.get_read_count_dict(summary, "sample")

    monkeypatch.setattr(
        library,
        "compute_seqid_mapping",
        lambda *args: [
            ("chrA", {"forward": np.array([1.0]), "reverse": np.array([0.0])}),
            ("chrB", {"forward": np.array([1.0]), "reverse": np.array([0.0])}),
        ],
    )
    for normalization in ("raw", "min", "mil"):
        (tmp_path / normalization).mkdir()

    library.compute_wig(
        "unused.bam",
        f"{tmp_path}/",
        "sample",
        library_total,
        minimum_total,
        mapping_style="global",
    )

    minimum_track = (tmp_path / "min" / "sample.min.forward.wig").read_text()
    per_million_track = (
        tmp_path / "mil" / "sample.mil.forward.wig"
    ).read_text()
    assert minimum_track.count("1 0.5") == 2
    assert per_million_track.count("1 250000.0") == 2


def test_sam_strand_inverter_toggles_only_reverse_strand_flag(tmp_path):
    source = tmp_path / "source.sam"
    destination = tmp_path / "inverted.sam"
    rows = [
        ["@HD", "VN:1.6", "SO:unknown"],
        ["forward", "0", "chr1", "1", "255", "4M", "*", "0", "0", "ACGT", "IIII"],
        ["reverse", "16", "chr1", "2", "255", "4M", "*", "0", "0", "ACGT", "IIII"],
        ["paired-forward", "99", "chr1", "3", "255", "4M", "=", "7", "8", "ACGT", "IIII"],
        ["paired-reverse", "147", "chr1", "7", "255", "4M", "=", "3", "-8", "ACGT", "IIII"],
    ]
    with source.open("w", newline="") as handle:
        csv.writer(handle, dialect="excel-tab").writerows(rows)

    subprocess.run(
        [
            sys.executable,
            str(SAM_STRAND_INVERTER),
            "--sam_in_filepath",
            str(source),
            "--sam_out_filepath",
            str(destination),
        ],
        check=True,
    )

    with destination.open(newline="") as handle:
        inverted = list(csv.reader(handle, dialect="excel-tab"))

    expected = [row.copy() for row in rows]
    for row, flag in zip(expected[1:], (16, 0, 115, 131), strict=True):
        row[1] = str(flag)
    assert inverted == expected
