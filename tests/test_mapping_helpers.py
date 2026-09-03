"""Focused regressions for legacy mapping helper scripts."""

import csv
import subprocess
import sys
from pathlib import Path

import pysam


REPO = Path(__file__).resolve().parent.parent
TOTAL_MAPPED_READS = REPO / "workflow" / "scripts" / "total_mapped_reads.py"
SAM_STRAND_INVERTER = REPO / "workflow" / "scripts" / "sam_strand_inverter.py"


def _write_alignment(handle, name, reference_id, start, sequence, flag=0):
    alignment = pysam.AlignedSegment()
    alignment.query_name = name
    alignment.query_sequence = sequence
    alignment.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
    alignment.flag = flag
    alignment.reference_id = reference_id
    alignment.reference_start = start
    alignment.mapping_quality = 255 if reference_id >= 0 else 0
    alignment.cigartuples = [(0, len(sequence))] if reference_id >= 0 else None
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
