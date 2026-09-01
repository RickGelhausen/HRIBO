"""Regression tests for empty and repeat prediction aggregation runs."""

import subprocess
import sys
from pathlib import Path

import pytest


pd = pytest.importorskip("pandas")

REPO = Path(__file__).resolve().parent.parent
SCRIPTS = REPO / "workflow" / "scripts"
GFF3_HEADER = "##gff-version 3\n"


def gff_record(
    seqid,
    start,
    end,
    attributes,
    *,
    source="test",
    score=".",
    strand="+",
    phase="0",
):
    return (
        "\t".join(
            [
                seqid,
                source,
                "CDS",
                str(start),
                str(end),
                score,
                strand,
                phase,
                attributes,
            ]
        )
        + "\n"
    )


def run_script(name, *arguments):
    result = subprocess.run(
        [sys.executable, str(SCRIPTS / name), *map(str, arguments)],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def records(path):
    return [
        line.split("\t")
        for line in Path(path).read_text().splitlines()
        if line and not line.startswith("#")
    ]


def attributes(row):
    return {
        key: value
        for field in row[8].split(";")
        if "=" in field
        for key, value in [field.split("=", 1)]
    }


def test_concatenate_is_sorted_and_nonempty_to_empty_rerun_truncates(tmp_path):
    earlier = tmp_path / "earlier.gff"
    later = tmp_path / "later.gff"
    output_a = tmp_path / "combined-a.gff"
    output_b = tmp_path / "combined-b.gff"

    earlier.write_text(GFF3_HEADER + gff_record("chr", 10, 20, "ID=early;"))
    later.write_text(GFF3_HEADER + gff_record("chr", 30, 40, "ID=late;"))

    run_script("concatenate_gff.py", later, earlier, "-o", output_a)
    run_script("concatenate_gff.py", earlier, later, "-o", output_b)

    assert output_a.read_text() == output_b.read_text()
    assert [int(row[3]) for row in records(output_a)] == [10, 30]
    assert output_a.read_text().count(GFF3_HEADER) == 1

    # Reuse the first output path after predictions disappear. Both common
    # empty representations must replace, not preserve, the earlier records.
    earlier.write_text(GFF3_HEADER)
    later.write_text("")
    run_script("concatenate_gff.py", earlier, later, "-o", output_a)
    assert output_a.read_text() == GFF3_HEADER


def test_concatenate_failure_does_not_replace_existing_output(tmp_path):
    malformed = tmp_path / "malformed.gff"
    output = tmp_path / "existing.gff"
    previous = GFF3_HEADER + gff_record("chr", 1, 3, "ID=previous;")
    malformed.write_text("chr\ttest\tCDS\t1\t3\t.\t+\t.\n")
    output.write_text(previous)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "concatenate_gff.py"),
            str(malformed),
            "-o",
            str(output),
        ],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert "expected 9 tab-separated GFF3 columns" in result.stderr
    assert output.read_text() == previous


@pytest.mark.parametrize(
    "empty_content", ["", GFF3_HEADER], ids=["zero-byte", "header-only"]
)
def test_deepribo_empty_rerun_writes_and_truncates_both_outputs(
    tmp_path, empty_content
):
    predictions = tmp_path / "deepribo-all.gff"
    annotation = tmp_path / "annotation.gff"
    merged = tmp_path / "explicit-main-name.gff"
    plus = tmp_path / "explicit-declared-plus-name.gff"
    reversed_predictions = tmp_path / "deepribo-reversed.gff"
    reversed_merged = tmp_path / "deepribo-reversed-main.gff"
    reversed_plus = tmp_path / "deepribo-reversed-plus.gff"

    rows = [
        gff_record(
            "chr",
            10,
            30,
            "ID=p1;deepribo_distance=-1;condition=A;method=deepribo;replicate=1;",
            source="deepribo",
            score="0.9",
        ),
        gff_record(
            "chr",
            10,
            30,
            "ID=p1;Condition=A;Method=deepribo;Replicate=2;",
            source="deepribo",
            score="0.8",
            # Legacy intermediates stored distance in this column. The merger
            # must consume it but emit a required complete-ORF CDS phase of 0.
            phase="0",
        ),
        gff_record(
            "chr",
            2,
            6,
            "ID=p2;deepribo_distance=0;condition=B;method=deepribo;replicate=1;",
            source="deepribo",
            score="0.5",
        ),
        gff_record(
            "chr",
            100,
            120,
            "ID=p3;deepribo_distance=0;condition=B;method=deepribo;replicate=2;",
            source="deepribo",
            score="0.5",
        ),
    ]
    predictions.write_text(GFF3_HEADER + "".join(rows))
    reversed_predictions.write_text(GFF3_HEADER + "".join(reversed(rows)))
    annotation.write_text(
        GFF3_HEADER + "chr\ttest\tgene\t100\t200\t.\t+\t.\tID=gene1;Name=gene1;\n"
    )

    command = (
        "merge_duplicates_deepribo.py",
        "-i",
        predictions,
        "-o",
        merged,
        "--plus-output",
        plus,
        "-a",
        annotation,
    )
    run_script(*command)
    merged_rows = records(merged)
    plus_rows = records(plus)
    assert [int(row[3]) for row in plus_rows] == [2, 10, 100]
    assert all(row[7] == "0" for row in merged_rows + plus_rows)
    assert [row[5] for row in merged_rows] == ["1", "2", "3"]
    assert [int(row[3]) for row in merged_rows[1:]] == [2, 100]

    best = attributes(merged_rows[0])
    assert best["pred_value"] == "0.9"
    assert best["evidence"] == "A-1 A-2"
    assert best["deepribo_distance"] == "-1"
    assert best["novel_rank"] == "1"
    assert not any(key[0].isupper() for key in best if key not in {"ID", "Name"})

    run_script(
        "merge_duplicates_deepribo.py",
        "-i",
        reversed_predictions,
        "-o",
        reversed_merged,
        "--plus-output",
        reversed_plus,
        "-a",
        annotation,
    )
    assert merged.read_text() == reversed_merged.read_text()
    assert plus.read_text() == reversed_plus.read_text()

    predictions.write_text(empty_content)
    run_script(*command)
    assert merged.read_text() == GFF3_HEADER
    assert plus.read_text() == GFF3_HEADER


def test_reparation_merge_is_sorted_and_empty_rerun_truncates(tmp_path):
    predictions = tmp_path / "reparation-all.gff"
    reversed_predictions = tmp_path / "reparation-reversed.gff"
    output = tmp_path / "reparation.gff"
    reversed_output = tmp_path / "reparation-reversed-output.gff"

    rows = [
        gff_record(
            "chr",
            30,
            40,
            "ID=late;ORF_type=sORF;Length=11;Ribo_count=60;Prob=0.6;Condition=B;Replicate=1;Method=reparation;",
            source="reparation",
            score="0.6",
        ),
        gff_record(
            "chr",
            10,
            20,
            "ID=early;ORF_type=sORF;Length=11;Ribo_count=70;Prob=0.7;Condition=A;Replicate=2;Method=reparation;",
            source="reparation",
            score="0.7",
        ),
        gff_record(
            "chr",
            10,
            20,
            "ID=early;ORF_type=sORF;Length=11;Ribo_count=80;Prob=0.8;Condition=A;Replicate=1;Method=reparation;",
            source="reparation",
            score="0.8",
        ),
    ]
    predictions.write_text(GFF3_HEADER + "".join(rows))
    reversed_predictions.write_text(GFF3_HEADER + "".join(reversed(rows)))

    run_script("merge_duplicates_reparation.py", "-i", predictions, "-o", output)
    run_script(
        "merge_duplicates_reparation.py",
        "-i",
        reversed_predictions,
        "-o",
        reversed_output,
    )

    assert output.read_text() == reversed_output.read_text()
    output_rows = records(output)
    assert [int(row[3]) for row in output_rows] == [10, 30]
    assert output.read_text().startswith(GFF3_HEADER)
    early = output_rows[0]
    early_attributes = attributes(early)
    assert early[5] == "0.8"
    assert early_attributes["prob"] == "0.8"
    assert early_attributes["ribo_count"] == "80"
    assert early_attributes["replicate"] == "1"
    assert early_attributes["evidence"] == "reparation-A-1 reparation-A-2"
    assert not any(
        key[0].isupper()
        for key in early_attributes
        if key not in {"ID", "Name"}
    )

    predictions.write_text("")
    run_script("merge_duplicates_reparation.py", "-i", predictions, "-o", output)
    assert output.read_text() == GFF3_HEADER

    output.write_text(GFF3_HEADER + gff_record("chr", 1, 3, "ID=stale;"))
    predictions.write_text(GFF3_HEADER)
    run_script("merge_duplicates_reparation.py", "-i", predictions, "-o", output)
    assert output.read_text() == GFF3_HEADER

    annotation = tmp_path / "annotation.gff"
    reannotated = tmp_path / "reparation-annotated.gff"
    annotation.write_text(
        GFF3_HEADER + "chr\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=gene1;\n"
    )
    reannotated.write_text(GFF3_HEADER + gff_record("chr", 1, 3, "ID=stale;"))
    run_script(
        "reannotate_orfs.py",
        "-a",
        annotation,
        "-c",
        output,
        "-o",
        reannotated,
    )
    assert reannotated.read_text() == GFF3_HEADER


def test_zero_predictions_reach_schema_bearing_excel_workbooks(tmp_path):
    pytest.importorskip("xlsxwriter")
    pytest.importorskip("openpyxl")

    empty_annotation = tmp_path / "empty-predictions.gff"
    raw_counts = tmp_path / "prediction-counts.raw"
    mapped_counts = tmp_path / "prediction-counts.gff"
    ribo_bam = tmp_path / "RIBO-A-1.bam"
    rna_bam = tmp_path / "RNA-A-1.bam"
    empty_annotation.write_text(GFF3_HEADER)

    # No BAM or featureCounts executable is needed: an empty annotation is a
    # complete result and returns before the external counter is invoked.
    run_script(
        "call_featurecounts.py",
        "-b",
        ribo_bam,
        rna_bam,
        "-a",
        empty_annotation,
        "-o",
        raw_counts,
    )
    raw_text = raw_counts.read_text()
    assert raw_text.startswith("#hribo-read-counts-v1\t")
    assert "RIBO-A-1" in raw_text and "RNA-A-1" in raw_text

    run_script(
        "map_reads_to_annotation.py",
        "-i",
        raw_counts,
        "-a",
        empty_annotation,
        "-o",
        mapped_counts,
    )
    mapped_text = mapped_counts.read_text()
    assert mapped_text.startswith("#hribo-gff-read-counts-v1\t")
    assert records(mapped_counts) == []

    genome = tmp_path / "genome.fa"
    totals = tmp_path / "totals.txt"
    genome.write_text(">chr\nATGATGATG\n")
    totals.write_text("RIBO-A-1\tchr\t10\nRNA-A-1\tchr\t10\n")

    expected_columns = {
        "generate_excel_reparation.py": "Reparation_probability",
        "generate_excel_deepribo.py": "Novel_rank",
    }
    for script, expected_column in expected_columns.items():
        workbook = tmp_path / f"{script}.xlsx"
        run_script(
            script,
            "-t",
            totals,
            "-r",
            mapped_counts,
            "-g",
            genome,
            "-o",
            workbook,
        )
        frame = pd.read_excel(workbook, sheet_name="CDS", engine="openpyxl")
        assert frame.empty
        assert expected_column in frame.columns
        assert "RIBO-A-1_rpkm" in frame.columns


def test_reannotate_failure_preserves_existing_output(tmp_path):
    annotation = tmp_path / "annotation.gff"
    malformed = tmp_path / "malformed.gff"
    output = tmp_path / "existing.gff"
    previous = GFF3_HEADER + gff_record("chr", 1, 3, "ID=previous;")
    annotation.write_text(
        GFF3_HEADER + "chr\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=gene1;\n"
    )
    malformed.write_text("chr\ttest\tCDS\t1\t3\t.\t+\t.\n")
    output.write_text(previous)

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "reannotate_orfs.py"),
            "-a",
            str(annotation),
            "-c",
            str(malformed),
            "-o",
            str(output),
        ],
        cwd=SCRIPTS,
        capture_output=True,
        text=True,
    )

    assert result.returncode != 0
    assert output.read_text() == previous


def test_deepribo_output_pair_rolls_back_if_second_replace_fails(
    tmp_path, monkeypatch
):
    import merge_duplicates_deepribo as deepribo_merge

    main = tmp_path / "main.gff"
    plus = tmp_path / "plus.gff"
    main.write_text("old main\n")
    plus.write_text("old plus\n")

    real_replace = deepribo_merge.os.replace
    calls = 0

    def fail_once_on_second_install(source, destination):
        nonlocal calls
        calls += 1
        # Two old files move to backups first, then the staged main is installed.
        if calls == 4:
            raise OSError("simulated second-output replacement failure")
        return real_replace(source, destination)

    monkeypatch.setattr(deepribo_merge.os, "replace", fail_once_on_second_install)

    with pytest.raises(OSError, match="simulated"):
        deepribo_merge.atomic_write_outputs(
            [(main, "new main\n"), (plus, "new plus\n")]
        )

    assert main.read_text() == "old main\n"
    assert plus.read_text() == "old plus\n"


def test_prediction_rules_have_no_hidden_append_state_and_declare_plus_output():
    deepribo_rule = (REPO / "workflow" / "rules" / "deepribo.smk").read_text()
    merge_rule = (REPO / "workflow" / "rules" / "merge.smk").read_text()

    assert "condition=set(" not in deepribo_rule
    assert "condition=set(" not in merge_rule
    assert "condition=conditions" in deepribo_rule
    assert "condition=conditions" in merge_rule
    assert ".unsorted" not in merge_rule
    assert ">>" not in merge_rule
    assert "--plus-output {output.plus:q}" in deepribo_rule
