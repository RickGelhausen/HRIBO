"""End-to-end test of the TIS advisor against a simulated library.

The simulation plants a known P-site offset in a known set of read lengths, so
this checks that the whole path -- BAM, annotation filtering, metagene profiling,
offset estimation, recommendation -- recovers ground truth, not just that it
completes.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pysam")
pytest.importorskip("interlap")

import simulate_riboseq as sim

REPO = Path(__file__).resolve().parent.parent
ADVISOR = REPO / "workflow" / "scripts" / "tis_advisor.py"


def run_advisor(reference_dir, bam, out_dir):
    result = subprocess.run(
        [
            sys.executable, str(ADVISOR),
            "-b", str(bam),
            "-a", str(reference_dir / "annotation.gff"),
            "-g", str(reference_dir / "genome.fa"),
            "-o", str(out_dir),
            "-r", "24-35",
            "--rpkm_threshold", "0",
            "--filtering_methods", "length",
            "--include_plotly_js", "online",
        ],
        capture_output=True,
        text=True,
        cwd=str(REPO / "workflow" / "scripts"),
    )
    assert result.returncode == 0, result.stderr
    return json.loads((out_dir / "tis_recommendation.json").read_text())


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    directory = tmp_path_factory.mktemp("reference")
    sim.write_genome(directory / "genome.fa")
    sim.write_annotation(directory / "annotation.gff")
    return directory


def test_recovers_planted_offset_with_periodicity(reference, tmp_path):
    bam = tmp_path / "RIBO-A-1.bam"
    sim.write_bam(bam, periodic=True, signal=True, seed=1)
    payload = run_advisor(reference, bam, tmp_path / "out")

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"], "no recommendation was made for a clear signal"
    assert set(recommendation["read_lengths"]) <= sim.GOOD_LENGTHS
    assert set(recommendation["offsets"].values()) == {sim.PLANTED_OFFSET}
    assert recommendation["confidence"] in {"high", "medium"}


def test_recovers_planted_offset_without_periodicity(reference, tmp_path):
    """The bacterial case: a clear initiation peak but no 3-nt periodicity."""
    bam = tmp_path / "RIBO-A-2.bam"
    sim.write_bam(bam, periodic=False, signal=True, seed=2)
    payload = run_advisor(reference, bam, tmp_path / "out")

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"]
    assert set(recommendation["offsets"].values()) == {sim.PLANTED_OFFSET}


def test_no_recommendation_without_signal(reference, tmp_path):
    """An RNA-seq-like library must not produce a confident setup."""
    bam = tmp_path / "RNA-A-1.bam"
    sim.write_bam(bam, periodic=False, signal=False, seed=3)
    payload = run_advisor(reference, bam, tmp_path / "out")

    recommendation = payload["recommendation"]
    assert recommendation["read_lengths"] == []
    assert recommendation["confidence"] == "none"
    assert recommendation["warnings"]


def test_outputs_are_written(reference, tmp_path):
    bam = tmp_path / "RIBO-A-1.bam"
    sim.write_bam(bam, periodic=True, signal=True, seed=1)
    out = tmp_path / "out"
    run_advisor(reference, bam, out)

    for name in ("tis_recommendation.html", "tis_recommendation.json", "read_length_evidence.tsv"):
        assert (out / name).is_file(), f"{name} was not written"

    evidence = (out / "read_length_evidence.tsv").read_text().splitlines()
    assert len(evidence) > 1
    assert evidence[0].startswith("read_length\t")


def test_report_contains_a_pasteable_config(reference, tmp_path):
    bam = tmp_path / "RIBO-A-1.bam"
    sim.write_bam(bam, periodic=True, signal=True, seed=1)
    out = tmp_path / "out"
    run_advisor(reference, bam, out)

    html = (out / "tis_recommendation.html").read_text()
    assert "psiteOffsets" in html
    assert "readLengths" in html
