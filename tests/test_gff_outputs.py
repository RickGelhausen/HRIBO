"""Golden-output tests for the annotation-transforming scripts.

These rewrite GFF files in the middle of the workflow and had no tests. The
snapshots in tests/golden_gff were captured before the attribute parsing was
consolidated.

Regenerate intentional changes with:
    python3 tests/regenerate_golden_gff.py
"""

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("pandas")

import make_gff_fixture

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "workflow" / "scripts"
GOLDEN = Path(__file__).resolve().parent / "golden_gff"


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    return make_gff_fixture.build(tmp_path_factory.mktemp("gff_inputs"))


def command(name, inputs, output):
    annotation = str(inputs / "annotation_gff3.gff")
    commands = {
        "gtf2gff3": ["gtf2gff3.py", "-a", str(inputs / "annotation_gtf2.gtf")],
        "gtf2gff3_passthrough": ["gtf2gff3.py", "-a", annotation],
        "enriched": ["enrich_annotation.py", "-a", annotation],
        "united": ["annotation_unite.py", "-a", str(inputs / "total_annotation.gff")],
        "reannotated": ["reannotate_orfs.py", "-a", annotation,
                        "-c", str(inputs / "reparation_tracks_typed.gff")],
        "reparation_merged": ["merge_duplicates_reparation.py",
                              "-i", str(inputs / "reparation_tracks_typed.gff")],
        "deepribo_merged": ["merge_duplicates_deepribo.py",
                            "-i", str(inputs / "deepribo_tracks.gff"), "-a", annotation],
    }
    script, *rest = commands[name]
    return [sys.executable, str(SCRIPTS / script), *rest, "-o", str(output)]


def run(name, inputs, tmp_path):
    output = tmp_path / f"{name}.gff"
    result = subprocess.run(
        command(name, inputs, output), capture_output=True, text=True, cwd=str(SCRIPTS)
    )
    assert result.returncode == 0, f"{name} failed:\n{result.stderr}"
    assert output.is_file(), f"{name} produced no output"
    return output


SCRIPT_NAMES = ["gtf2gff3", "gtf2gff3_passthrough", "enriched", "united",
                "reannotated", "reparation_merged", "deepribo_merged"]


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_output_matches_golden(name, inputs, tmp_path):
    output = run(name, inputs, tmp_path)
    actual = output.read_text()
    expected = (GOLDEN / f"{name}.gff").read_text()
    assert actual == expected, f"{name}.gff no longer matches tests/golden_gff/{name}.gff"


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_output_is_not_empty(name, inputs, tmp_path):
    output = run(name, inputs, tmp_path)
    assert output.stat().st_size > 0, f"{name} produced an empty file"


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_output_has_nine_columns(name, inputs, tmp_path):
    """Every row must still be a valid GFF record."""
    output = run(name, inputs, tmp_path)
    for number, line in enumerate(output.read_text().splitlines(), start=1):
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        assert len(fields) == 9, f"{name}.gff line {number} has {len(fields)} columns, not 9"


def test_reannotate_handles_empty_orf_type(inputs, tmp_path):
    """Reparation writes "ORF_type=;" when it has no type for an ORF.

    Splitting the attributes on both ";" and "=" and dropping empty fields left
    an odd number of items, and rebuilding them pairwise then ran off the end of
    the list. That crashed the whole workflow at reannotatedORFs, since
    tracks/reparation_annotated.gff feeds read counting and the overview table.
    """
    output = tmp_path / "empty_orf_type.gff"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / "reannotate_orfs.py"),
         "-a", str(inputs / "annotation_gff3.gff"),
         "-c", str(inputs / "reparation_tracks.gff"),
         "-o", str(output)],
        capture_output=True, text=True, cwd=str(SCRIPTS),
    )
    assert result.returncode == 0, (
        "reannotate_orfs.py still fails on an empty ORF_type:\n" + result.stderr
    )
    assert output.is_file() and output.stat().st_size > 0


def test_empty_orf_type_keeps_the_same_records(inputs, tmp_path):
    """The empty type must not change which ORFs are reported."""
    typed = run("reannotated", inputs, tmp_path)
    untyped = tmp_path / "untyped.gff"
    subprocess.run(
        [sys.executable, str(SCRIPTS / "reannotate_orfs.py"),
         "-a", str(inputs / "annotation_gff3.gff"),
         "-c", str(inputs / "reparation_tracks.gff"),
         "-o", str(untyped)],
        capture_output=True, text=True, cwd=str(SCRIPTS), check=True,
    )

    def coordinates(path):
        return [
            tuple(line.split("\t")[:5])
            for line in Path(path).read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ]

    assert coordinates(typed) == coordinates(untyped)


# --------------------------------------------------------------------------
# Regressions
# --------------------------------------------------------------------------


def test_merge_outputs_are_deterministic(inputs, tmp_path):
    """The evidence field was joined from a set, so it varied between runs.

    Python randomises string hashing per process, so the same input produced a
    different file every time, which defeats the point of a reproducible
    workflow.
    """
    for name in ("reparation_merged", "deepribo_merged"):
        outputs = []
        for attempt in range(3):
            output = tmp_path / f"{name}_{attempt}.gff"
            subprocess.run(
                command(name, inputs, output), capture_output=True, text=True,
                cwd=str(SCRIPTS), check=True,
            )
            outputs.append(output.read_text())
        assert outputs[0] == outputs[1] == outputs[2], f"{name} is not deterministic"


def test_enrich_annotation_survives_an_unresolvable_parent(inputs, tmp_path):
    """A Parent naming a feature absent from the file used to raise KeyError.

    enrich_annotation feeds auxiliary/enriched_annotation.gff, which every read
    counting rule depends on, so the crash took the whole workflow down.
    """
    output = tmp_path / "enriched_orphan.gff"
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / "enrich_annotation.py"),
         "-a", str(inputs / "annotation_gff3.gff"), "-o", str(output)],
        capture_output=True, text=True, cwd=str(SCRIPTS),
    )
    assert result.returncode == 0, (
        "enrich_annotation.py still fails on an unresolvable Parent:\n" + result.stderr
    )
    assert output.stat().st_size > 0


def test_deepribo_old_locus_tag_is_not_the_locus_tag(inputs, tmp_path):
    """old_locus_tag was read from the "locus_tag" attribute by mistake."""
    output = run("deepribo_merged", inputs, tmp_path)

    checked = 0
    for line in output.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 9:
            continue
        attributes = fields[8]
        pairs = dict(
            field.split("=", 1) for field in attributes.split(";") if "=" in field
        )
        if "old_locus_tag" in pairs and "locus_tag" in pairs:
            assert pairs["old_locus_tag"] != pairs["locus_tag"], (
                f"old_locus_tag duplicates locus_tag: {attributes}"
            )
            checked += 1
    assert checked > 0, "no rows carried both tags, so the check proved nothing"
