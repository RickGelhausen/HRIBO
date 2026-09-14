"""Regression tests for cross-condition genome-browser exports."""

import json
import shutil
import subprocess
from urllib.parse import unquote

import pytest

from lib.diffex_browser import export_tracks


def _fixture():
    rows = [
        {
            "feature_id": "g/2;odd",
            "conditions": {
                "WT": {
                    "RNA": {"state": "detected", "normalized_count": 12.5},
                    "RIBO": {"state": "not_detected"},
                },
                "Stress / Δ": {
                    "RNA": {"state": "uncertain"},
                    "RIBO": {"state": "detected", "replicate_fraction": 1.0},
                },
            },
            "contrasts": {
                "Stress / Δ-WT": {
                    "RNA": {"state": "up", "log2fc": 2.5, "padj": 0.001},
                    "RIBO": {"state": "not_significant", "log2fc": 0.2},
                    "TE": {"state": "down", "log2fc": -1.3, "padj": 0.04},
                },
            },
        },
        {
            "feature_id": "g1",
            "conditions": {
                "WT": {
                    "RNA": {"state": "detected"},
                    "RIBO": {"state": "detected"},
                },
            },
            "contrasts": {
                "Stress / Δ-WT": {
                    "RNA": {"state": "down", "log2fc": -2.0, "padj": 0.02},
                    "TE": {"state": "not_tested"},
                },
            },
        },
        {
            "feature_id": "missing-coords",
            "conditions": {"WT": {"RNA": {"state": "detected"}}},
            "contrasts": {},
        },
    ]
    coords = {
        "g/2;odd": {
            "seqid": "chr;2", "start": 100, "end": 199, "strand": "-",
            "type": "CDS", "name": "alpha;beta=1%",
        },
        "g1": {
            "seqid": "chr1", "start": 5, "end": 50, "strand": "+",
            "type": "gene", "name": "gene one",
        },
    }
    return rows, coords


def _track(manifest, output_dir, kind, label, assay, state):
    entry = next(
        track for track in manifest["tracks"]
        if (track["kind"], track["label"], track["assay"], track["state"])
        == (kind, label, assay, state)
    )
    return entry, output_dir / entry["file"]


def test_exports_only_positive_states_with_attributes_and_manifest(tmp_path):
    rows, coords = _fixture()
    outputs = export_tracks(
        rows, tmp_path, ["WT", "Stress / Δ"], ["Stress / Δ-WT"], coords
    )
    manifest = json.loads((tmp_path / "tracks_manifest.json").read_text())

    # Two assays per condition and up/down tracks for three contrast assays.
    assert len(outputs) == 11
    assert len(manifest["tracks"]) == 10
    assert manifest["unmapped_feature_ids"] == ["missing-coords"]
    assert all(path.exists() for path in outputs)
    assert all("/" not in track["file"] for track in manifest["tracks"])

    entry, path = _track(manifest, tmp_path, "condition", "WT", "RNA", "detected")
    assert entry["feature_count"] == 2
    records = path.read_text().splitlines()[1:]
    assert len(records) == 2
    assert records[0].split("\t")[3:5] == ["5", "50"]
    fields = records[1].split("\t")
    assert len(fields) == 9
    assert fields[:8] == ["chr%3B2", "HRIBO", "CDS", "100", "199", ".", "-", "."]
    assert "ID=g%2F2%3Bodd" in fields[8]
    assert "Name=alpha%3Bbeta%3D1%25" in fields[8]
    assert "gene_id=g%2F2%3Bodd" in fields[8]
    assert "condition=WT;assay=RNA;state=detected" in fields[8]
    assert "normalized_count=12.5" in fields[8]

    _, path = _track(manifest, tmp_path, "condition", "WT", "RIBO", "detected")
    assert "g%2F2%3Bodd" not in path.read_text()
    _, path = _track(manifest, tmp_path, "condition", "Stress / Δ", "RIBO", "detected")
    assert "replicate_fraction=1.0" in path.read_text()

    _, path = _track(manifest, tmp_path, "contrast", "Stress / Δ-WT", "RNA", "up")
    assert "state=up;log2fc=2.5;padj=0.001" in path.read_text()
    assert "contrast=Stress%20%2F%20%CE%94-WT" in path.read_text()
    _, path = _track(manifest, tmp_path, "contrast", "Stress / Δ-WT", "RNA", "down")
    assert "ID=g1" in path.read_text()
    _, path = _track(manifest, tmp_path, "contrast", "Stress / Δ-WT", "TE", "down")
    assert "log2fc=-1.3;padj=0.04" in path.read_text()


def test_empty_tracks_are_valid_header_only_files(tmp_path):
    rows, coords = _fixture()
    export_tracks(rows, tmp_path, ["WT"], ["Stress / Δ-WT"], coords)
    manifest = json.loads((tmp_path / "tracks_manifest.json").read_text())
    entry, path = _track(manifest, tmp_path, "contrast", "Stress / Δ-WT", "RIBO", "up")
    assert entry["feature_count"] == 0
    assert path.read_text() == "##gff-version 3\n"


def test_deterministic_order_and_safe_distinct_filenames(tmp_path):
    rows, coords = _fixture()
    rows[0]["conditions"]["A/B"] = {"RNA": {"state": "detected"}}
    rows[0]["conditions"]["A B"] = {"RNA": {"state": "detected"}}
    labels = ["A/B", "A B", "WT"]
    first = tmp_path / "first"
    second = tmp_path / "second"
    export_tracks(rows, first, labels, [], coords)
    export_tracks(list(reversed(rows)), second, list(reversed(labels)), [], coords)
    files = sorted(path.name for path in first.iterdir())
    assert files == sorted(path.name for path in second.iterdir())
    assert len(files) == len(set(files))
    for name in files:
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_numeric_metadata_omits_nonfinite_values(tmp_path):
    rows, coords = _fixture()
    rows[0]["contrasts"]["Stress / Δ-WT"]["RNA"]["padj"] = float("nan")
    export_tracks(rows, tmp_path, [], ["Stress / Δ-WT"], coords)
    manifest = json.loads((tmp_path / "tracks_manifest.json").read_text())
    _, path = _track(manifest, tmp_path, "contrast", "Stress / Δ-WT", "RNA", "up")
    assert "log2fc=2.5" in path.read_text()
    assert "padj=" not in path.read_text()


def test_duplicate_feature_ids_are_rejected(tmp_path):
    rows, coords = _fixture()
    with pytest.raises(ValueError, match="Duplicate feature ID"):
        export_tracks([rows[0], rows[0]], tmp_path, ["WT"], [], coords)


def test_invalid_coordinates_are_rejected_for_exported_feature(tmp_path):
    rows, coords = _fixture()
    coords["g1"]["start"] = 0
    with pytest.raises(ValueError, match="Invalid coordinates"):
        export_tracks(rows, tmp_path, ["WT"], [], coords)


def test_gff_attribute_values_round_trip(tmp_path):
    rows, coords = _fixture()
    export_tracks(rows, tmp_path, ["WT"], [], coords)
    manifest = json.loads((tmp_path / "tracks_manifest.json").read_text())
    _, path = _track(manifest, tmp_path, "condition", "WT", "RNA", "detected")
    special = next(line for line in path.read_text().splitlines() if "g%2F2%3Bodd" in line)
    attributes = dict(item.split("=", 1) for item in special.split("\t")[8].split(";"))
    assert unquote(attributes["ID"]) == "g/2;odd"
    assert unquote(attributes["Name"]) == "alpha;beta=1%"


@pytest.mark.skipif(shutil.which("gt") is None, reason="GenomeTools is unavailable")
def test_tracks_pass_genometools_gff3_validator(tmp_path):
    rows, coords = _fixture()
    paths = export_tracks(rows, tmp_path, ["WT", "Stress / Δ"], ["Stress / Δ-WT"], coords)
    for path in paths:
        if path.suffix != ".gff3":
            continue
        result = subprocess.run(
            [shutil.which("gt"), "gff3validator", str(path)],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"{path.name} is invalid GFF3:\n{result.stderr}"
