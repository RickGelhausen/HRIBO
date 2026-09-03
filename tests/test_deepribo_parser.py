"""Tests for HRIBO's strict DeepRibo bedGraph compatibility loader."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from deepribo_data_parser import (
    BedGraphError,
    load_image_parser,
    load_signal,
    read_bedgraph,
    validate_tracks,
)


def write_track(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_first_and_only_interval_is_preserved_and_empty_track_is_zero(tmp_path):
    coverage = write_track(tmp_path / "coverage.bedgraph", "chr1\t2\t5\t7\n")
    elongating = write_track(tmp_path / "empty.bedgraph", "")

    signal, elongating_signal = load_signal(
        str(coverage), str(elongating), "chr1", 8
    )

    np.testing.assert_array_equal(signal, [0, 0, 7, 7, 7, 0, 0, 0])
    np.testing.assert_array_equal(elongating_signal, np.zeros(8))


def test_absent_contig_is_zero_without_discarding_other_contigs(tmp_path):
    coverage = write_track(
        tmp_path / "coverage.bedgraph",
        "chr2\t0\t2\t3\nchr2\t4\t6\t5\n",
    )
    elongating = write_track(
        tmp_path / "elongating.bedgraph",
        "# generated fixture\ntrack type=bedGraph\nchr2\t1\t3\t2\n",
    )

    missing, missing_elongating = load_signal(
        str(coverage), str(elongating), "chr1", 6
    )
    present, present_elongating = load_signal(
        str(coverage), str(elongating), "chr2", 6
    )

    np.testing.assert_array_equal(missing, np.zeros(6))
    np.testing.assert_array_equal(missing_elongating, np.zeros(6))
    np.testing.assert_array_equal(present, [3, 3, 0, 0, 5, 5])
    np.testing.assert_array_equal(present_elongating, [0, 2, 2, 0, 0, 0])


def test_malformed_or_incompatible_tracks_fail_clearly(tmp_path):
    malformed_cases = (
        ("short", "chr1\t0\t2\n", "expected four"),
        ("coordinate", "chr1\t2\t2\t1\n", "invalid half-open interval"),
        ("count", "chr1\t0\t2\tNaN\n", "invalid count"),
        (
            "overlap",
            "chr1\t0\t3\t1\nchr1\t2\t4\t1\n",
            "overlaps or precedes",
        ),
    )

    for name, content, message in malformed_cases:
        path = write_track(tmp_path / f"{name}.bedgraph", content)
        with pytest.raises(BedGraphError, match=message):
            read_bedgraph(str(path.resolve()))

    fasta = write_track(tmp_path / "genome.fa", ">chr1 description\nAAAAAA\n")
    empty = write_track(tmp_path / "empty.bedgraph", "")
    unknown = write_track(tmp_path / "unknown.bedgraph", "chrX\t0\t2\t1\n")
    beyond = write_track(tmp_path / "beyond.bedgraph", "chr1\t5\t7\t1\n")
    with pytest.raises(BedGraphError, match="unknown contig 'chrX'"):
        validate_tracks([str(unknown), str(empty), str(empty), str(empty)], str(fasta))
    with pytest.raises(BedGraphError, match="beyond chr1's length 6"):
        validate_tracks([str(beyond), str(empty), str(empty), str(empty)], str(fasta))


def test_image_parser_checksum_guard_rejects_unexpected_source(tmp_path):
    parser = tmp_path / "DataParser.py"
    parser.write_text("raise RuntimeError('must never be imported')\n", encoding="utf-8")

    with pytest.raises(RuntimeError, match="unexpected DeepRibo parser"):
        load_image_parser(parser)
