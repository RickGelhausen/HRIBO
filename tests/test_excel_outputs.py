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

import pandas as pd
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


def command(name, inputs, output, overview_contrasts=("B-A",)):
    """The command line for one script, mirroring how the rules invoke it."""
    genome = ["-g", str(inputs / "genome.fa")]
    totals = ["-t", str(inputs / "total_mapped_reads.txt")]
    annotation = ["-a", str(inputs / "annotation.gff")]
    overview_contrast_args = (
        [] if overview_contrasts is None else ["-c", *overview_contrasts]
    )

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
        # The overview pulls in most of excel_utils, including the annotation,
        # prediction and pooled differential expression readers that nothing
        # else exercises.
        "overview": ["generate_excel_overview.py", *genome, *totals,
                     "-a", str(inputs / "total_annotation.gtf"),
                     "--mapped_reads_reparation", str(inputs / "reparation_annotation.gff"),
                     "--mapped_reads_deepribo", str(inputs / "deepribo_annotation.gff"),
                     "--riborex", str(inputs / "riborex_pooled.csv"),
                     "--xtail", str(inputs / "xtail_pooled.csv"),
                     "--deltate", str(inputs / "deltate_pooled.csv"),
                     *overview_contrast_args],
    }
    script, *rest = commands[name]
    return [sys.executable, str(SCRIPTS / script), *rest, "-o", str(output)]


def run(name, inputs, tmp_path, overview_contrasts=("B-A",)):
    output = tmp_path / f"{name}.xlsx"
    result = subprocess.run(
        command(name, inputs, output, overview_contrasts),
        capture_output=True,
        text=True,
        cwd=str(SCRIPTS),
    )
    assert result.returncode == 0, f"{name} failed:\n{result.stderr}"
    assert output.is_file(), f"{name} produced no output"
    return output


SCRIPT_NAMES = ["annotation", "reparation", "deepribo", "readtable",
                "riborex", "xtail", "deltate", "overview"]

POOLED_DIFFEX_COLUMNS = {
    "riborex": [
        "gene_id", "baseMean", "log2FC", "log2FC_SE", "stat", "pvalue",
        "pvalue_adjusted", "contrast",
    ],
    "xtail": [
        "gene_id", "mRNA_log2FC", "RPF_log2FC", "log2FC_TE_v1",
        "pvalue_v1", "log2FC_TE_v2", "pvalue_v2", "log2FC_TE_final",
        "pvalue_final", "pvalue_adjusted", "contrast",
    ],
    "deltate": [
        "gene_id", "RIBO_baseMean", "RIBO_log2FC", "RIBO_log2FC_SE",
        "RIBO_pvalue", "RIBO_pvalue_adjusted", "RNA_baseMean", "RNA_log2FC",
        "RNA_log2FC_SE", "RNA_pvalue", "RNA_pvalue_adjusted", "TE_baseMean",
        "TE_log2FC", "TE_log2FC_SE", "TE_stat", "TE_pvalue",
        "TE_pvalue_adjusted", "contrast",
    ],
}


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


@pytest.mark.parametrize("name", ["riborex", "xtail", "deltate"])
def test_generated_diffex_workbook_can_be_pooled(name, inputs, tmp_path):
    """The pooler must consume the exact standardized workbook schema."""
    workbook = tmp_path / "B-A_sorted.xlsx"
    generated = subprocess.run(
        command(name, inputs, workbook),
        capture_output=True,
        text=True,
        cwd=str(SCRIPTS),
    )
    assert generated.returncode == 0, generated.stderr

    pooled = tmp_path / f"{name}_all.csv"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "merge_differential_expression.py"),
            str(workbook),
            "--output_csv",
            str(pooled),
            "--tool",
            name,
        ],
        capture_output=True,
        text=True,
        cwd=str(SCRIPTS),
    )
    assert result.returncode == 0, result.stderr

    actual = pd.read_csv(pooled)
    source = excel_snapshot.read_workbook(workbook)["all"]
    assert list(actual.columns) == POOLED_DIFFEX_COLUMNS[name]
    assert actual["gene_id"].tolist() == source["Identifier"].tolist()
    assert set(actual["contrast"]) == {f"{name}_B-A"}


def test_translational_efficiency_columns_are_present(inputs, tmp_path):
    """The fixture pairs RIBO with RNA, so TE columns must be generated."""
    output = run("reparation", inputs, tmp_path)
    columns = list(excel_snapshot.read_workbook(output)["CDS"].columns)
    assert any(column.endswith("_TE") for column in columns)
    assert any(column.endswith("_rpkm") for column in columns)
    # A replicate average is added only when a condition has more than one.
    assert any("avg" in column for column in columns)


def test_annotation_rpkm_uses_each_two_contig_library_total(inputs, tmp_path):
    output = run("annotation", inputs, tmp_path)
    actual = excel_snapshot.read_workbook(output)["all"]
    source = pd.read_csv(inputs / "total_annotation.gtf", sep="\t", header=None)
    totals = pd.read_csv(
        inputs / "total_mapped_reads.txt",
        sep="\t",
        header=None,
        names=["library", "contig", "count"],
    ).groupby("library")["count"].sum()

    for contig in make_excel_fixture.CONTIGS:
        source_row = source[source[0] == contig].iloc[0]
        identifier = f"{source_row[0]}:{source_row[3]}-{source_row[4]}:{source_row[6]}"
        actual_row = actual[actual["Identifier"] == identifier].iloc[0]
        for index, wildcard in enumerate(make_excel_fixture.WILDCARDS):
            expected = round(
                source_row[9 + index]
                * 1_000_000_000
                / (totals[wildcard] * (source_row[4] - source_row[3] + 1)),
                2,
            )
            assert actual_row[f"{wildcard}_rpkm"] == expected


@pytest.mark.parametrize(
    "name, reads_flag",
    [("annotation", "-r"), ("readtable", "-r"), ("overview", "-a")],
)
def test_feature_tables_reject_a_summary_count_column_mismatch(
    name, reads_flag, inputs, tmp_path
):
    malformed_reads = tmp_path / "missing-library-count.gff"
    rows = []
    for line in (inputs / "total_annotation.gtf").read_text().splitlines():
        rows.append("\t".join(line.split("\t")[:-1]))
    malformed_reads.write_text("\n".join(rows) + "\n")

    output = tmp_path / f"{name}.xlsx"
    argv = command(name, inputs, output)
    reads_argument = argv.index(reads_flag) + 1
    argv[reads_argument] = str(malformed_reads)
    result = subprocess.run(
        argv, capture_output=True, text=True, cwd=str(SCRIPTS)
    )

    assert result.returncode != 0
    assert "found 7 library read-count columns" in result.stderr
    assert "8 libraries are present in the mapped-read summary" in result.stderr


def test_annotation_splits_features_into_sheets(inputs, tmp_path):
    output = run("annotation", inputs, tmp_path)
    sheets = excel_snapshot.read_workbook(output)
    assert {"CDS", "gene", "rRNA", "tRNA", "sRNA"} <= set(sheets)


def test_overview_honors_explicit_contrast_orientation(inputs, tmp_path):
    """A requested B-A contrast must neither be renamed nor replaced by A-B."""
    output = run("overview", inputs, tmp_path, overview_contrasts=("B-A",))
    sheets = excel_snapshot.read_workbook(output)

    for sheet_name, frame in sheets.items():
        contrast_columns = [
            column
            for column in frame.columns
            if column.startswith(("xtail_", "riborex_", "deltaTE_"))
        ]
        assert contrast_columns, f"{sheet_name} has no differential-expression columns"
        assert all("_B-A_" in column for column in contrast_columns)
        assert not any("_A-B_" in column for column in contrast_columns)

    # The fixture contains pooled B-A rows. Checking a populated column makes
    # this a value-propagation regression, not merely a header spelling test.
    assert sheets["all"]["xtail_B-A_TE_log2FC"].notna().any()


def test_overview_infers_pairwise_contrasts_when_none_are_passed(inputs, tmp_path):
    """Standalone use without -c retains the deterministic pairwise fallback."""
    output = run("overview", inputs, tmp_path, overview_contrasts=None)
    columns = excel_snapshot.read_workbook(output)["all"].columns

    assert "xtail_A-B_TE_log2FC" in columns
    assert "xtail_B-A_TE_log2FC" not in columns


def test_overview_writes_every_declared_side_output(inputs, tmp_path):
    output = run("overview", inputs, tmp_path)

    tsv = output.with_suffix(".tsv")
    gff_outputs = (
        output.with_suffix(".gff"),
        output.with_name(f"{output.stem}_misc.gff"),
    )
    for expected in (tsv, *gff_outputs):
        assert expected.is_file()
        assert expected.stat().st_size > 0

    pd.testing.assert_frame_equal(
        pd.read_csv(tsv, sep="\t"),
        excel_snapshot.read_workbook(output)["all"],
        check_dtype=False,
    )
    for gff in gff_outputs:
        lines = gff.read_text().splitlines()
        assert lines[0] == "##gff-version 3"
        for line in lines[1:]:
            fields = line.split("\t")
            assert len(fields) == 9
            assert 1 <= int(fields[3]) <= int(fields[4])
