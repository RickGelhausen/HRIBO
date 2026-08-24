"""Golden-output tests for the excel generating scripts.

These scripts had no tests and produce spreadsheets people read directly, so the
refactor that deduplicated them had to preserve their output exactly. The
snapshots in tests/golden were captured from the implementation as it stood
before that refactor.

If a change here is intentional, regenerate with:
    python3 tests/regenerate_golden.py
and review the diff, which is plain CSV text.
"""

import subprocess
import sys
from pathlib import Path

import pytest

pytest.importorskip("Bio")
pytest.importorskip("xlsxwriter")
pytest.importorskip("openpyxl")

import excel_snapshot
import make_excel_fixture

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "workflow" / "scripts"
GOLDEN = Path(__file__).resolve().parent / "golden"


@pytest.fixture(scope="module")
def inputs(tmp_path_factory):
    return make_excel_fixture.build(tmp_path_factory.mktemp("excel_inputs"))


def command(name, inputs, output):
    """The command line for one script, mirroring how the rules invoke it."""
    genome = ["-g", str(inputs / "genome.fa")]
    totals = ["-t", str(inputs / "total_mapped_reads.txt")]
    annotation = ["-a", str(inputs / "annotation.gff")]

    commands = {
        "annotation": ["generate_excel.py", *genome, *totals,
                       "-r", str(inputs / "total_annotation.gtf")],
        "reparation": ["generate_excel_reparation.py", *genome, *totals,
                       "-r", str(inputs / "reparation_annotation.gff")],
        "deepribo": ["generate_excel_deepribo.py", *genome, *totals,
                     "-r", str(inputs / "deepribo_annotation.gff")],
        "readtable": ["generate_read_table.py",
                      "-r", str(inputs / "total_annotation.gtf"), *totals],
        "riborex": ["generate_excel_riborex.py", *annotation, *genome,
                    "-i", str(inputs / "riborex_all.csv")],
        "xtail": ["generate_excel_xtail.py", *annotation, *genome,
                  "-i", str(inputs / "xtail_all.csv")],
        "deltate": ["generate_excel_deltate.py", *annotation, *genome,
                    "-i", str(inputs / "deltaRibo.txt"),
                    "-r", str(inputs / "deltaRNA.txt"),
                    "-t", str(inputs / "deltaTE.txt")],
    }
    script, *rest = commands[name]
    return [sys.executable, str(SCRIPTS / script), *rest, "-o", str(output)]


def run(name, inputs, tmp_path):
    output = tmp_path / f"{name}.xlsx"
    result = subprocess.run(
        command(name, inputs, output), capture_output=True, text=True, cwd=str(SCRIPTS)
    )
    assert result.returncode == 0, f"{name} failed:\n{result.stderr}"
    assert output.is_file(), f"{name} produced no output"
    return output


SCRIPT_NAMES = ["annotation", "reparation", "deepribo", "readtable",
                "riborex", "xtail", "deltate"]


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_output_matches_golden(name, inputs, tmp_path):
    """Byte-for-byte equal sheets, headers and cell values."""
    output = run(name, inputs, tmp_path)
    actual = excel_snapshot.to_snapshot(output)
    expected = (GOLDEN / f"{name}.csv").read_text()

    assert actual == expected, (
        f"{name}.xlsx no longer matches tests/golden/{name}.csv\n"
        + excel_snapshot.diff_summary(expected, actual)
    )


@pytest.mark.parametrize("name", SCRIPT_NAMES)
def test_output_is_not_empty(name, inputs, tmp_path):
    """A silently empty spreadsheet would still match a stale golden."""
    output = run(name, inputs, tmp_path)
    sheets = excel_snapshot.read_workbook(output)
    assert sheets, f"{name} produced no sheets"
    assert any(len(frame) > 0 for frame in sheets.values()), (
        f"{name} produced only empty sheets"
    )


def test_translational_efficiency_columns_are_present(inputs, tmp_path):
    """The fixture pairs RIBO with RNA, so TE columns must be generated."""
    output = run("reparation", inputs, tmp_path)
    columns = list(excel_snapshot.read_workbook(output)["CDS"].columns)
    assert any(column.endswith("_TE") for column in columns)
    assert any(column.endswith("_rpkm") for column in columns)
    # A replicate average is added only when a condition has more than one.
    assert any("avg" in column for column in columns)


def test_annotation_splits_features_into_sheets(inputs, tmp_path):
    output = run("annotation", inputs, tmp_path)
    sheets = excel_snapshot.read_workbook(output)
    assert {"CDS", "gene", "rRNA", "tRNA", "sRNA"} <= set(sheets)
