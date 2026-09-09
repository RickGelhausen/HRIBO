"""Regression tests for the reusable real-data release comparator."""

import importlib.util
import json
import struct
import sys
from pathlib import Path

import openpyxl
import pysam
import pytest


REPO = Path(__file__).resolve().parent.parent
COMPARATOR_PATH = REPO / ".github" / "scripts" / "compare_real_data_outputs.py"
SPEC = importlib.util.spec_from_file_location(
    "compare_real_data_outputs", COMPARATOR_PATH
)
assert SPEC is not None and SPEC.loader is not None
COMPARATOR = importlib.util.module_from_spec(SPEC)
sys.modules["compare_real_data_outputs"] = COMPARATOR
SPEC.loader.exec_module(COMPARATOR)


def write_text_result(root: Path, relative: str, text: str) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


def write_workbook(
    root: Path, relative: str, rows: list[list[object]], sheet: str = "all"
) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook = openpyxl.Workbook()
    worksheet = workbook.active
    worksheet.title = sheet
    for row in rows:
        worksheet.append(row)
    workbook.save(path)
    workbook.close()
    return path


def write_bam(root: Path, relative: str, starts: list[int]) -> Path:
    return write_named_bam(
        root,
        relative,
        [(f"read-{index}", start) for index, start in enumerate(starts)],
    )


def write_named_bam(
    root: Path,
    relative: str,
    records: list[tuple[str, int]],
    *,
    program: dict[str, str] | None = None,
) -> Path:
    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr", "LN": 200}]}
    if program is not None:
        header["PG"] = [program]
    with pysam.AlignmentFile(path, "wb", header=header) as output:
        for name, start in records:
            alignment = pysam.AlignedSegment()
            alignment.query_name = name
            alignment.query_sequence = "A" * 10
            alignment.flag = 0
            alignment.reference_id = 0
            alignment.reference_start = start
            alignment.mapping_quality = 60
            alignment.cigar = [(0, 10)]
            alignment.query_qualities = pysam.qualitystring_to_array("I" * 10)
            output.write(alignment)
    pysam.index(str(path))
    return path


def write_bigwig(
    root: Path,
    relative: str,
    *,
    total: float = 10.0,
    chromosome: str = "chr",
) -> Path:
    """Write the minimal BigWig structures decoded by the comparator."""

    path = root / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    endian = "<"
    key_size = 8
    chromosome_tree_offset = 64
    root_offset = chromosome_tree_offset + 32
    root_size = 4 + key_size + 8
    total_summary_offset = root_offset + root_size
    full_data_offset = total_summary_offset + 40
    full_index_offset = full_data_offset + 1
    header = struct.pack(
        f"{endian}IHHQQQHHQQIQ",
        COMPARATOR.BIGWIG_MAGIC,
        4,
        0,
        chromosome_tree_offset,
        full_data_offset,
        full_index_offset,
        0,
        0,
        0,
        total_summary_offset,
        0,
        0,
    )
    tree_header = struct.pack(
        f"{endian}IIIIQQ",
        COMPARATOR.BIGWIG_CHROM_TREE_MAGIC,
        256,
        key_size,
        8,
        1,
        0,
    )
    key = chromosome.encode().ljust(key_size, b"\0")
    leaf = struct.pack(f"{endian}BBH", 1, 0, 1) + key + struct.pack(
        f"{endian}II", 0, 200
    )
    summary = struct.pack(f"{endian}Qdddd", 10, 0.0, 2.0, total, 20.0)
    path.write_bytes(header + tree_header + leaf + summary + b"\0\0")
    return path


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ("TE log2FC", "te_log2fc"),
        ("pvalue-adjusted", "pvalue_adjusted"),
        ("RIBO_log2FoldChange", "ribo_log2fc"),
        ("RNA_lfcSE", "rna_log2fc_se"),
        ("TE_padj", "te_pvalue_adjusted"),
        ("identifer", "identifier"),
        ("Pred_value", "deepribo_score"),
        ("Pred_rank", "deepribo_rank"),
        ("Pred_probability", "reparation_probability"),
        ("15nt upstream", "upstream_15nt"),
        ("  Genome / Start  ", "genome_start"),
    ],
)
def test_header_normalization_handles_historical_punctuation(header, expected):
    assert COMPARATOR.normalize_header(header) == expected


def test_blank_workbook_heading_is_rejected_instead_of_becoming_text_none():
    with pytest.raises(COMPARATOR.ComparisonError, match="empty column heading"):
        COMPARATOR.build_table(["Identifier", None], [["gene-a", 1]])


def test_release_18_prediction_headers_match_current_workbook_headers(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "auxiliary/predictions.xlsx"
    write_workbook(
        baseline,
        relative,
        [
            [
                "identifer",
                "Pred_value",
                "Pred_rank",
                "Pred_probability",
                "15nt upstream",
            ],
            ["gene-a", 2.5, 1, 0.8, "ATGATGATGATGATG"],
        ],
    )
    write_workbook(
        candidate,
        relative,
        [
            [
                "Identifier",
                "Deepribo_score",
                "Deepribo_rank",
                "Reparation_probability",
                "Upstream_15nt",
            ],
            ["gene-a", 2.5, 1, 0.8, "ATGATGATGATGATG"],
        ],
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is False


def test_identical_table_runs_are_unchanged(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    contents = "Identifier,count\ngene-a,10\ngene-b,20\n"
    write_text_result(baseline, "readcounts/counts.csv", contents)
    write_text_result(candidate, "readcounts/counts.csv", contents)

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"] == {
        "baseline_artifacts": 1,
        "candidate_artifacts": 1,
        "comparable_artifacts": 1,
        "unchanged_artifacts": 1,
        "changed_artifacts": 0,
        "candidate_only_artifacts": 0,
        "baseline_only_artifacts": 0,
        "invalid_artifacts": 0,
        "errors": 0,
        "warnings": 0,
        "review_required": False,
    }


def test_r_row_names_and_legacy_diffex_headers_are_canonicalized(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "riborex/B-A_deseq2.csv",
        '"",baseMean,log2FoldChange,lfcSE,padj\n'
        "gene-b,20,2,0.2,0.02\n"
        "gene-a,10,-1,0.1,0.01\n",
    )
    write_text_result(
        candidate,
        "riborex/B-A_deseq2.csv",
        "Identifier,baseMean,log2FC,log2FC_SE,pvalue_adjusted\n"
        "gene-a,10,-1,0.1,0.01\n"
        "gene-b,20,2,0.2,0.02\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    table = artifact["comparison"]["tables"]["table"]

    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is False
    assert table["candidate"]["normalized_columns"] == [
        "identifier",
        "basemean",
        "log2fc",
        "log2fc_se",
        "pvalue_adjusted",
    ]


def test_pooled_diffex_uses_identifier_and_contrast_as_composite_key(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    header = "gene_id,log2FoldChange,padj,contrast\n"
    write_text_result(
        baseline,
        "riborex/riborex_all.csv",
        header + "gene-a,1,0.01,riborex_B-A\ngene-a,-2,0.02,riborex_C-A\n",
    )
    write_text_result(
        candidate,
        "riborex/riborex_all.csv",
        header + "gene-a,1,0.01,riborex_B-A\ngene-a,-3,0.02,riborex_C-A\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    detail = artifact["comparison"]["tables"]["table"]["comparison"]

    assert detail["baseline_key_columns"] == ["identifier", "contrast"]
    assert artifact["baseline"]["tables"]["table"]["duplicate_keys"] == 0
    assert detail["common_keys"] == 2
    assert detail["numeric_columns"]["log2fc"]["changed"] == 1


def test_annotation_rows_use_identifier_and_feature_as_order_stable_key(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "auxiliary/annotation_total.xlsx"
    header = ["Identifier", "Feature", "Genome", "Start", "Stop", "value"]
    gene = ["chr:1-9:+", "gene", "chr", 1, 9, 1]
    coding_sequence = ["chr:1-9:+", "CDS", "chr", 1, 9, 2]
    unrelated = ["chr:20-28:+", "gene", "chr", 20, 28, 9]
    write_workbook(baseline, relative, [header, gene, unrelated, coding_sequence])
    write_workbook(candidate, relative, [header, coding_sequence, gene, unrelated])

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    detail = artifact["comparison"]["tables"]["all"]["comparison"]

    assert detail["baseline_key_columns"] == ["identifier", "feature"]
    assert artifact["baseline"]["tables"]["all"]["duplicate_keys"] == 0
    assert artifact["status"] == "unchanged"


def test_sample_rows_use_method_condition_replicate_composite_key(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    content = (
        "method,condition,replicate,value\n"
        "RIBO,A,1,10\n"
        "RIBO,A,2,11\n"
    )
    write_text_result(baseline, "readcounts/samples.csv", content)
    write_text_result(
        candidate,
        "readcounts/samples.csv",
        "method,condition,replicate,value\nRIBO,A,2,11\nRIBO,A,1,10\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    detail = artifact["comparison"]["tables"]["table"]["comparison"]

    assert detail["baseline_key_columns"] == ["method", "condition", "replicate"]
    assert artifact["status"] == "unchanged"


def test_pca_special_formats_have_explicit_headers_and_keys(tmp_path):
    root = tmp_path / "run"
    meta = write_text_result(
        root,
        "pca/meta.csv",
        "\tsampletype\texpr\nRIBO_A_1\tRIBO_A\tA\n",
    )
    correlation = write_text_result(
        root,
        "pca/rld_cor.tsv",
        "\tRIBO_A_1\tRIBO_A_2\n"
        "RIBO_A_1\t1\t0.8\nRIBO_A_2\t0.8\t1\n",
    )
    variance = write_text_result(
        root, "pca/variance_percentages.tsv", "0.7\n0.2\n0.1\n"
    )

    meta_table = COMPARATOR.inspect_artifact(meta, "pca/meta.csv")["tables"]["table"]
    correlation_table = COMPARATOR.inspect_artifact(
        correlation, "pca/rld_cor.tsv"
    )["tables"]["table"]
    variance_table = COMPARATOR.inspect_artifact(
        variance, "pca/variance_percentages.tsv"
    )["tables"]["table"]

    assert meta_table["key_columns"] == ("sample",)
    assert correlation_table["key_columns"] == ("sample",)
    assert variance_table["columns"] == ("component", "variance_percentage")
    assert list(variance_table["keyed_rows"]) == ['["PC1"]', '["PC2"]', '["PC3"]']


def test_pca_scores_use_sample_name_instead_of_pc1_as_key(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    header = "PC1\tPC2\tgroup\tsampletype\tname\n"
    write_text_result(
        baseline,
        "pca/rld.tsv",
        header + "-1\t2\tA\tRIBO_A\tRIBO_A_1\n1\t-2\tA\tRIBO_A\tRIBO_A_2\n",
    )
    write_text_result(
        candidate,
        "pca/rld.tsv",
        header + "1\t-2\tA\tRIBO_A\tRIBO_A_2\n-1\t2\tA\tRIBO_A\tRIBO_A_1\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    detail = artifact["comparison"]["tables"]["table"]["comparison"]

    assert detail["baseline_key_columns"] == ["sampletype", "name"]
    assert artifact["status"] == "unchanged"


def test_changed_table_reports_key_overlap_numeric_deltas_and_direction(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "deltate/result.csv",
        "Identifier,value,TE log2FC\na,1,-2\nb,2,1\nc,3,2\n",
    )
    write_text_result(
        candidate,
        "deltate/result.csv",
        "Identifier,value,TE-log2FC,extra\na,1.1,-1.8,x\nb,2.2,1.2,x\nc,3.3,2.2,x\nd,4,1,x\n",
    )

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = report["artifacts"][0]
    detail = artifact["comparison"]["tables"]["table"]["comparison"]

    assert artifact["status"] == "changed"
    assert artifact["review_hint"].startswith("Differential results may change")
    assert detail["key_jaccard"] == 0.75
    assert detail["candidate_only_keys"] == 1
    assert detail["candidate_only_columns"] == ["extra"]
    assert detail["numeric_columns"]["value"]["spearman"] == 1.0
    assert detail["numeric_columns"]["te_log2fc"]["sign_agreement"] == 1.0


def test_optional_key_overlap_threshold_turns_large_drift_into_error(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/counts.csv", "Identifier,count\na,1\nb,2\n")
    write_text_result(candidate, "readcounts/counts.csv", "Identifier,count\na,1\nc,2\n")

    report = COMPARATOR.compare_runs(
        baseline, candidate, minimum_key_overlap=0.5
    )

    assert report["summary"]["errors"] == 1
    assert "key Jaccard" in report["issues"][0]["message"]


def test_missing_candidate_artifact_is_a_structural_error(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/one.csv", "Identifier,count\na,1\n")
    write_text_result(baseline, "readcounts/two.csv", "Identifier,count\na,2\n")
    write_text_result(candidate, "readcounts/one.csv", "Identifier,count\na,1\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["baseline_only_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert report["issues"][0]["path"] == "readcounts/two.csv"


def test_intentionally_removed_output_can_be_allowlisted(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/one.csv", "Identifier,count\na,1\n")
    write_text_result(baseline, "readcounts/legacy.csv", "Identifier,count\na,2\n")
    write_text_result(candidate, "readcounts/one.csv", "Identifier,count\na,1\n")

    report = COMPARATOR.compare_runs(
        baseline, candidate, allow_missing=("readcounts/legacy.*",)
    )

    assert report["summary"]["errors"] == 0
    legacy = next(item for item in report["artifacts"] if "legacy" in item["path"])
    assert legacy["allowed"] is True


def test_new_candidate_artifact_is_validated_but_not_an_error(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    common = "Identifier,count\na,1\n"
    write_text_result(baseline, "readcounts/counts.csv", common)
    write_text_result(candidate, "readcounts/counts.csv", common)
    write_text_result(
        candidate,
        "tis_advice/RIBO-A-1/tis_recommendation.json",
        '{"mapping_method": "fiveprime", "offset": 12}\n',
    )

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["candidate_only_artifacts"] == 1
    assert report["summary"]["errors"] == 0
    added = next(item for item in report["artifacts"] if item["kind"] == "json")
    assert added["candidate"]["top_level_keys"] == ["mapping_method", "offset"]


def test_malformed_candidate_gff_fails_validation(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    valid = "##gff-version 3\nchr\tHRIBO\tCDS\t1\t9\t.\t+\t0\tID=one;\n"
    write_text_result(baseline, "tracks/result.gff", valid)
    write_text_result(candidate, "tracks/result.gff", "chr\ttoo\tfew\tcolumns\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert "instead of 9" in report["artifacts"][0]["error"]


def test_malformed_delimited_row_shape_is_reported_as_invalid(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/counts.csv", "Identifier,count\na,1\n")
    write_text_result(
        candidate, "readcounts/counts.csv", "Identifier,count\na,1,unexpected\n"
    )

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert "3 columns instead of 2" in report["artifacts"][0]["error"]


def test_invalid_utf8_is_reported_as_an_invalid_artifact(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/counts.csv", "Identifier,count\na,1\n")
    invalid = candidate / "readcounts" / "counts.csv"
    invalid.parent.mkdir(parents=True)
    invalid.write_bytes(b"Identifier,count\na,\xff\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert "cannot parse delimited table" in report["artifacts"][0]["error"]


@pytest.mark.parametrize(
    ("fields", "message"),
    [
        (["chr", "x", "CDS", "1", "9", ".", "+", "0", "ID"], "key=value"),
        (
            ["chr", "x", "CDS", "1", "9", ".", "+", "0", "ID=a;ID=b"],
            "occurs more than once",
        ),
        (["chr", "x", "CDS", "1", "9", ".", "+", "0", "=a"], "empty key"),
        (["chr", "x", "CDS", "1", "9", "NaN", "+", "0", "ID=a"], "non-finite score"),
        (["chr", "x", "CDS", "1", "9", "high", "+", "0", "ID=a"], "non-numeric score"),
        (["chr", "x", "CDS", "1", "9", ".", "+", "3", "ID=a"], "invalid phase"),
        (["", "x", "CDS", "1", "9", ".", "+", "0", "ID=a"], "seqid"),
        (["chr", "", "CDS", "1", "9", ".", "+", "0", "ID=a"], "empty source"),
        (["chr", "x", ".", "1", "9", ".", "+", "0", "ID=a"], "feature"),
        (["chr", "x", "CDS", "1", "9", ".", "+", "0", ""], "attributes column"),
    ],
)
def test_strict_gff_field_validation_reports_malformed_records(
    tmp_path, fields, message
):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    valid = "chr\tx\tCDS\t1\t9\t.\t+\t0\tID=a\n"
    write_text_result(baseline, "tracks/result.gff", valid)
    write_text_result(candidate, "tracks/result.gff", "\t".join(fields) + "\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert message in report["artifacts"][0]["error"]


@pytest.mark.parametrize(
    "relative",
    [
        "tracks/A.deepribo.gff",
        "tracks/deepribo_all.gff",
        "tracks/deepribo_merged.gff",
        "tracks/deepribo_merged_plus.gff",
        "tracks/totalAnnotation.gff",
        "tracks/updated_annotation.gff",
    ],
)
def test_known_release_18_deepribo_phases_are_normalized_with_warning(
    tmp_path, relative
):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    baseline_records = (
        "chr\tdeepribo\tCDS\t1\t9\t0.8\t+\t999999\tID=one;Pred_value=0.8;\n"
        "chr\tdeepribo\tCDS\t20\t28\t0.7\t-\t1\tID=two;Pred_value=0.7;\n"
    )
    candidate_records = baseline_records.replace("\t999999\t", "\t0\t").replace(
        "\t1\tID=two", "\t0\tID=two"
    )
    write_text_result(baseline, relative, baseline_records)
    write_text_result(candidate, relative, candidate_records)

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = report["artifacts"][0]
    compatibility = artifact["baseline_compatibility"]

    assert artifact["status"] == "unchanged"
    assert report["summary"]["errors"] == 0
    assert report["summary"]["warnings"] == 1
    assert report["summary"]["review_required"] is True
    assert compatibility == {
        "records_normalized": 2,
        "original_phase_counts": {"1": 1, "999999": 1},
        "normalized_phase": "0",
    }
    assert artifact["baseline"]["legacy_deepribo_phase_compatibility"] == (
        compatibility
    )
    assert "known HRIBO 1.8 DeepRibo CDS phase misuse" in report["issues"][0][
        "message"
    ]


def test_candidate_deepribo_phase_remains_strict_on_legacy_named_path(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tracks/deepribo_merged.gff"
    write_text_result(
        baseline,
        relative,
        "chr\tdeepribo\tCDS\t1\t9\t0.8\t+\t0\tID=one;deepribo_distance=-1;\n",
    )
    write_text_result(
        candidate,
        relative,
        "chr\tdeepribo\tCDS\t1\t9\t0.8\t+\t999999\tID=one;\n",
    )

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert "invalid phase '999999'" in report["artifacts"][0]["error"]


@pytest.mark.parametrize(
    ("relative", "source", "feature", "phase"),
    [
        ("tracks/unrelated.gff", "deepribo", "CDS", "999999"),
        ("tracks/deepribo_all.gff", "another-tool", "CDS", "999999"),
        ("tracks/deepribo_all.gff", "deepribo", "gene", "999999"),
        ("tracks/deepribo_all.gff", "deepribo", "CDS", "legacy-rank"),
        ("tracks/deepribo_all.gff", "deepribo", "CDS", "९"),
    ],
)
def test_legacy_phase_allowance_does_not_relax_other_baseline_gff(
    tmp_path, relative, source, feature, phase
):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        relative,
        f"chr\t{source}\t{feature}\t1\t9\t0.8\t+\t{phase}\tID=one;\n",
    )
    write_text_result(
        candidate,
        relative,
        f"chr\t{source}\t{feature}\t1\t9\t0.8\t+\t0\tID=one;\n",
    )

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert f"invalid phase {phase!r}" in report["artifacts"][0]["error"]


@pytest.mark.parametrize(
    "relative",
    [
        "tracks/A.deepribo.gff",
        "tracks/A.reparation.gff",
        "tracks/A.merged.gff",
    ],
)
def test_release_18_zero_byte_condition_gff_is_explicit_empty_compatibility(
    tmp_path, relative
):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, relative, "")
    write_text_result(candidate, relative, "##gff-version 3\n")

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = report["artifacts"][0]
    compatibility = {
        "records_normalized": 0,
        "interpretation": "empty feature set",
    }

    assert artifact["status"] == "unchanged"
    assert artifact["baseline"]["bytes"] == 0
    assert artifact["baseline"]["legacy_zero_byte_gff_compatibility"] == (
        compatibility
    )
    assert artifact["baseline_zero_result_compatibility"] == compatibility
    assert report["summary"]["errors"] == 0
    assert report["summary"]["warnings"] == 1
    assert report["summary"]["review_required"] is True
    assert "zero-byte per-condition GFF" in report["issues"][0]["message"]


def test_zero_byte_candidate_condition_gff_remains_invalid(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tracks/A.deepribo.gff"
    write_text_result(baseline, relative, "##gff-version 3\n")
    write_text_result(candidate, relative, "")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert report["artifacts"][0]["error"] == "file is empty"


def test_header_only_candidate_gff_is_a_valid_empty_feature_set(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tracks/result.gff"
    write_text_result(baseline, relative, "##gff-version 3\n")
    write_text_result(
        candidate,
        relative,
        "##gff-version 3\n# A zero-prediction result has no feature lines.\n",
    )

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["errors"] == 0
    assert report["artifacts"][0]["status"] == "unchanged"
    assert report["artifacts"][0]["candidate"]["records"] == 0


@pytest.mark.parametrize("contents", ["# no predictions\n", "##gff-version 2\n"])
def test_comment_only_candidate_without_gff3_version_is_invalid(tmp_path, contents):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tracks/result.gff"
    write_text_result(baseline, relative, "##gff-version 3\n")
    write_text_result(candidate, relative, contents)

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert "does not declare '##gff-version 3'" in report["artifacts"][0]["error"]


@pytest.mark.parametrize(
    "relative",
    [
        "tracks/deepribo_all.gff",
        "tracks/deepribo_merged.gff",
        "tracks/reparation.gff",
        "tracks/all.gff",
        "tracks/arbitrary.gff",
    ],
)
def test_zero_byte_global_or_arbitrary_baseline_gff_remains_invalid(
    tmp_path, relative
):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, relative, "")
    write_text_result(candidate, relative, "##gff-version 3\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert report["artifacts"][0]["error"] == "file is empty"


def test_nonempty_malformed_legacy_condition_gff_remains_invalid(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tracks/A.reparation.gff"
    write_text_result(baseline, relative, "not a GFF record\n")
    write_text_result(candidate, relative, "##gff-version 3\n")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert report["summary"]["errors"] == 1
    assert "instead of 9" in report["artifacts"][0]["error"]


def test_gff_comparison_reports_exact_and_three_nt_matches(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "tracks/predictions.gff",
        "##gff-version 3\nchr\tdeepribo\tCDS\t10\t30\t.\t+\t0\tID=one;\n"
        "chr\tdeepribo\tCDS\t50\t70\t.\t-\t0\tID=two;\n",
    )
    write_text_result(
        candidate,
        "tracks/predictions.gff",
        "##gff-version 3\nchr\tdeepribo\tCDS\t10\t30\t.\t+\t0\tID=one;\n"
        "chr\tdeepribo\tCDS\t53\t67\t.\t-\t0\tID=two-shifted;\n",
    )

    report = COMPARATOR.compare_runs(baseline, candidate)
    comparison = report["artifacts"][0]["comparison"]

    assert comparison["exact_coordinate_matches"] == 1
    assert comparison["within_3nt_matches"] == 2
    assert comparison["within_3nt_baseline_fraction"] == 1.0
    assert comparison["within_3nt_candidate_fraction"] == 1.0


def test_gff_tolerance_matching_is_maximum_cardinality_not_greedy(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    # Taking the exact 10-20 pair first leaves 13-23 unmatched.  The maximum
    # matching instead pairs 10-20 with 7-17 and 13-23 with 10-20.
    write_text_result(
        baseline,
        "tracks/predictions.gff",
        "chr\tx\tCDS\t10\t20\t.\t+\t0\tID=one;\n"
        "chr\tx\tCDS\t13\t23\t.\t+\t0\tID=two;\n",
    )
    write_text_result(
        candidate,
        "tracks/predictions.gff",
        "chr\tx\tCDS\t10\t20\t.\t+\t0\tID=one;\n"
        "chr\tx\tCDS\t7\t17\t.\t+\t0\tID=shifted;\n",
    )

    comparison = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0][
        "comparison"
    ]

    assert comparison["exact_coordinate_matches"] == 1
    assert comparison["within_3nt_matches"] == 2


def test_gff_status_is_semantic_across_line_and_attribute_order(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "tracks/result.gff",
        "chr\tx\tCDS\t1\t9\t.\t+\t0\tID=one;Name=a;\n"
        "chr\tx\tCDS\t20\t28\t.\t-\t0\tID=two;Name=b;\n",
    )
    write_text_result(
        candidate,
        "tracks/result.gff",
        "chr\tx\tCDS\t20\t28\t.\t-\t0\tName=b;ID=two;\n"
        "chr\tx\tCDS\t1\t9\t.\t+\t0\tName=a;ID=one;\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is False
    assert artifact["comparison"]["semantic_equal"] is True


def test_header_only_zero_result_workbooks_are_valid(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    rows = [["Identifier", "Reparation probability"]]
    baseline_book = write_workbook(
        baseline, "auxiliary/predictions_reparation.xlsx", rows, "CDS"
    )
    candidate_book = candidate / "auxiliary" / "predictions_reparation.xlsx"
    candidate_book.parent.mkdir(parents=True)
    candidate_book.write_bytes(baseline_book.read_bytes())

    report = COMPARATOR.compare_runs(baseline, candidate)
    summary = report["artifacts"][0]["candidate"]["tables"]["CDS"]

    assert report["summary"]["errors"] == 0
    assert summary["rows"] == 0
    assert summary["key_columns"] == ["identifier"]


def test_indexed_bams_are_compared_by_bounded_reference_statistics(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_bam(baseline, "maplink/RIBO-A-1.bam", [10, 30])
    write_bam(candidate, "maplink/RIBO-A-1.bam", [10, 30, 50])

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = report["artifacts"][0]
    metric = artifact["comparison"]["tables"]["index_statistics"]["comparison"]

    assert artifact["status"] == "changed"
    assert artifact["baseline"]["mapped"] == 2
    assert artifact["candidate"]["mapped"] == 3
    assert metric["numeric_columns"]["mapped"]["max_absolute_delta"] == 1.0


def test_bam_alignment_digest_detects_changes_hidden_by_index_counts(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_bam(baseline, "maplink/RIBO-A-1.bam", [10, 30])
    write_bam(candidate, "maplink/RIBO-A-1.bam", [20, 40])

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["baseline"]["mapped"] == artifact["candidate"]["mapped"] == 2
    assert artifact["comparison"]["tables"]["index_statistics"]["status"] == "unchanged"
    assert artifact["comparison"]["alignment_digest_equal"] is False
    assert artifact["status"] == "changed"


def test_bam_alignment_digest_is_independent_of_equal_coordinate_order(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    records = [("read-a", 10), ("read-b", 10), ("read-c", 10)]
    write_named_bam(baseline, "maplink/RIBO-A-1.bam", records)
    write_named_bam(candidate, "maplink/RIBO-A-1.bam", list(reversed(records)))

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["comparison"]["alignment_digest_equal"] is True
    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is None


def test_bam_program_provenance_does_not_change_semantic_status(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    records = [("read-a", 10), ("read-b", 20)]
    write_named_bam(
        baseline,
        "maplink/RIBO-A-1.bam",
        records,
        program={"ID": "mapper", "PN": "mapper", "VN": "1.0"},
    )
    write_named_bam(
        candidate,
        "maplink/RIBO-A-1.bam",
        records,
        program={"ID": "mapper", "PN": "mapper", "VN": "2.0"},
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["comparison"]["header_equal"] is True
    assert artifact["comparison"]["program_provenance_equal"] is False
    assert artifact["comparison"]["alignment_digest_equal"] is True
    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is None


def test_maplink_internal_bam_symlink_uses_index_beside_logical_path(tmp_path):
    for run_name in ("baseline", "candidate"):
        root = tmp_path / run_name
        target = write_bam(root, "stored/RIBO-A-1.bam", [10, 30])
        target_index = Path(f"{target}.bai")
        logical = root / "maplink" / "RIBO-A-1.bam"
        logical.parent.mkdir(parents=True)
        logical.symlink_to(target, target_is_directory=False)
        Path(f"{logical}.bai").symlink_to(target_index)

    report = COMPARATOR.compare_runs(
        tmp_path / "baseline", tmp_path / "candidate"
    )

    assert report["artifacts"][0]["status"] == "unchanged"
    assert report["artifacts"][0]["baseline"]["mapped"] == 2


def test_bigwig_discovery_and_decoded_summary_comparison(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    directories = (
        "globaltracks/raw",
        "centeredtracks/raw",
        "fiveprimetracks/raw",
        "threeprimetracks/raw",
    )
    for directory in directories:
        relative = f"{directory}/RIBO-A-1.raw.forward.bw"
        write_bigwig(baseline, relative, total=10.0)
        write_bigwig(candidate, relative, total=10.0)
    changed = "globaltracks/raw/RIBO-A-1.raw.forward.bw"
    write_bigwig(candidate, changed, total=11.0)

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = next(item for item in report["artifacts"] if item["path"] == changed)

    assert report["summary"]["comparable_artifacts"] == 4
    assert report["summary"]["changed_artifacts"] == 1
    assert artifact["kind"] == "bigwig"
    assert artifact["candidate"]["chromosomes"] == {"chr": {"id": 0, "length": 200}}
    assert artifact["candidate"]["total_summary"]["sum"] == 11.0
    assert artifact["comparison"]["decoded_summary_equal"] is False
    assert "compressed-byte equality" in artifact["comparison"]["comparison_mode"]


def test_invalid_bigwig_magic_is_a_structural_error(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "globaltracks/raw/RIBO-A-1.raw.forward.bw"
    write_bigwig(baseline, relative)
    write_text_result(candidate, relative, "not a BigWig")

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["invalid_artifacts"] == 1
    assert "BigWig magic" in report["artifacts"][0]["error"]


def test_identical_indexed_bam_summaries_are_not_forced_changed(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_bam(baseline, "maplink/RIBO-A-1.bam", [10, 30])
    write_bam(candidate, "maplink/RIBO-A-1.bam", [10, 30])

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["artifacts"][0]["status"] == "unchanged"
    assert report["summary"]["review_required"] is False


def test_numeric_missing_nonfinite_and_type_mismatches_are_reported(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "readcounts/counts.csv",
        "Identifier,value\na,1\nb,NaN\nc,3\nd,Infinity\n",
    )
    write_text_result(
        candidate,
        "readcounts/counts.csv",
        "Identifier,value\na,\nb,Infinity\nc,broken\nd,4\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    metric = artifact["comparison"]["tables"]["table"]["comparison"][
        "numeric_columns"
    ]["value"]

    assert metric["missing_mismatches"] == 1
    assert metric["nonfinite_value_mismatches"] == 1
    assert metric["type_mismatches"] == 1
    assert metric["nonfinite_mismatches"] == 1
    assert len(metric["mismatch_examples"]) == 4
    assert artifact["status"] == "changed"


def test_established_na_tokens_remain_missing_in_numeric_columns(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "readcounts/counts.csv",
        "Identifier,value\na,1\nb,NA\nc,2\nd,3\n",
    )
    write_text_result(
        candidate,
        "readcounts/counts.csv",
        "Identifier,value\na,1\nb,NA_character_\nc,2\nd,3\n",
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]
    metric = artifact["comparison"]["tables"]["table"]["comparison"][
        "numeric_columns"
    ]["value"]

    assert metric["pairs"] == 3
    assert metric["missing_pairs"] == 1
    assert artifact["status"] == "unchanged"


def test_correlation_threshold_accepts_unchanged_constant_column(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    contents = "Identifier,value\na,5\nb,5\nc,5\n"
    write_text_result(baseline, "readcounts/counts.csv", contents)
    write_text_result(candidate, "readcounts/counts.csv", contents)

    report = COMPARATOR.compare_runs(
        baseline, candidate, minimum_correlation=0.5
    )

    assert report["summary"]["errors"] == 0
    assert report["summary"]["review_required"] is False
    assert report["artifacts"][0]["status"] == "unchanged"


def test_correlation_threshold_rejects_undefined_constant_candidate(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(
        baseline,
        "readcounts/counts.csv",
        "Identifier,value\na,1\nb,2\nc,3\n",
    )
    write_text_result(
        candidate,
        "readcounts/counts.csv",
        "Identifier,value\na,5\nb,5\nc,5\n",
    )

    report = COMPARATOR.compare_runs(
        baseline, candidate, minimum_correlation=0.5
    )

    assert report["summary"]["errors"] == 1
    assert "Spearman is undefined for 3 finite pairs" in report["issues"][0][
        "message"
    ]


def test_json_status_is_semantic_and_retains_byte_identity(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tis_advice/RIBO-A-1/tis_recommendation.json"
    write_text_result(baseline, relative, '{"offset":12,"method":"fiveprime"}\n')
    write_text_result(
        candidate, relative, '{\n  "method": "fiveprime",\n  "offset": 12\n}\n'
    )

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["status"] == "unchanged"
    assert artifact["byte_identical"] is False
    assert artifact["comparison"]["semantic_equal"] is True


def test_json_non_boolean_number_spellings_remain_semantically_equal(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tis_advice/RIBO-A-1/tis_recommendation.json"
    write_text_result(baseline, relative, '{"offset":12}\n')
    write_text_result(candidate, relative, '{"offset":12.0}\n')

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["status"] == "unchanged"
    assert artifact["comparison"]["semantic_equal"] is True


def test_json_boolean_and_number_are_not_semantically_equal(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "tis_advice/RIBO-A-1/tis_recommendation.json"
    write_text_result(baseline, relative, '{"usable":true}\n')
    write_text_result(candidate, relative, '{"usable":1}\n')

    artifact = COMPARATOR.compare_runs(baseline, candidate)["artifacts"][0]

    assert artifact["status"] == "changed"
    assert artifact["comparison"]["semantic_equal"] is False


def test_global_stop_metagene_change_carries_scientific_review_hint(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = (
        "metageneprofiling/RIBO-A-1/raw/global_readcounts_stop.xlsx"
    )
    write_workbook(baseline, relative, [["Position", "28"], [-1, 0], [0, 4]], "chr")
    write_workbook(candidate, relative, [["Position", "28"], [-1, 4], [0, 0]], "chr")

    report = COMPARATOR.compare_runs(baseline, candidate)
    artifact = report["artifacts"][0]

    assert "dropped minus/start and plus/stop windows" in artifact["review_hint"]


def test_cli_writes_machine_readable_report_and_returns_structural_status(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    report_path = tmp_path / "reports" / "comparison.json"
    write_text_result(baseline, "readcounts/one.csv", "Identifier,count\na,1\n")
    write_text_result(baseline, "readcounts/two.csv", "Identifier,count\na,2\n")
    write_text_result(candidate, "readcounts/one.csv", "Identifier,count\na,1\n")

    return_code = COMPARATOR.main(
        [str(baseline), str(candidate), "--report", str(report_path)]
    )

    assert return_code == 1
    saved = json.loads(report_path.read_text())
    assert saved["schema_version"] == 1
    assert saved["summary"]["errors"] == 1
    assert not list(report_path.parent.glob("*.tmp"))


def test_discovery_includes_sorted_diffex_correlation_and_overview_outputs(tmp_path):
    root = tmp_path / "run"
    write_workbook(
        root,
        "xtail/B-A_sorted.xlsx",
        [["Identifier", "log2FC_TE_final"], ["gene-a", 1.0]],
    )
    write_text_result(
        root,
        "figures/SpearmanCorr_readCounts.tab",
        "#plotCorrelation --outFileCorMatrix\n"
        "\tRIBO_A_1\tRIBO_A_2\n"
        "RIBO_A_1\t1\t0.8\nRIBO_A_2\t0.8\t1\n",
    )
    write_text_result(
        root, "auxiliary/overview.tsv", "Identifier\tvalue\ngene-a\t1\n"
    )
    gff = "chr\tHRIBO\tCDS\t1\t9\t.\t+\t0\tID=one;\n"
    write_text_result(root, "auxiliary/overview.gff", gff)
    write_text_result(root, "auxiliary/overview_misc.gff", gff)

    discovered = COMPARATOR.discover_outputs(root.resolve())

    assert set(discovered) == {
        "auxiliary/overview.gff",
        "auxiliary/overview.tsv",
        "auxiliary/overview_misc.gff",
        "figures/SpearmanCorr_readCounts.tab",
        "xtail/B-A_sorted.xlsx",
    }
    matrix = COMPARATOR.inspect_artifact(
        discovered["figures/SpearmanCorr_readCounts.tab"],
        "figures/SpearmanCorr_readCounts.tab",
    )["tables"]["table"]
    assert matrix["key_columns"] == ("sample",)


def test_historical_root_correlation_matrix_maps_to_current_figures_path(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    matrix = "\tA\tB\nA\t1\t0.8\nB\t0.8\t1\n"
    write_text_result(baseline, "SpearmanCorr_readCounts.tab", matrix)
    write_text_result(candidate, "figures/SpearmanCorr_readCounts.tab", matrix)

    report = COMPARATOR.compare_runs(baseline, candidate)

    assert report["summary"]["comparable_artifacts"] == 1
    assert report["summary"]["baseline_only_artifacts"] == 0
    assert report["summary"]["candidate_only_artifacts"] == 0
    assert report["artifacts"][0]["path"] == "figures/SpearmanCorr_readCounts.tab"
    assert report["artifacts"][0]["status"] == "unchanged"


def test_comparator_rejects_equal_resolved_roots(tmp_path):
    root = tmp_path / "run"
    write_text_result(root, "readcounts/counts.csv", "Identifier,count\na,1\n")

    with pytest.raises(COMPARATOR.ComparisonError, match="same directory"):
        COMPARATOR.compare_runs(root, root / ".")


def test_discovery_rejects_output_symlink_that_escapes_run_root(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    write_text_result(baseline, "readcounts/counts.csv", "Identifier,count\na,1\n")
    external = write_text_result(
        tmp_path, "outside.csv", "Identifier,count\na,1\n"
    )
    link = candidate / "readcounts" / "counts.csv"
    link.parent.mkdir(parents=True)
    link.symlink_to(external)

    with pytest.raises(COMPARATOR.ComparisonError, match="escapes run root"):
        COMPARATOR.compare_runs(baseline, candidate)


def test_discovery_rejects_bam_index_symlink_that_escapes_run_root(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    relative = "maplink/RIBO-A-1.bam"
    write_bam(baseline, relative, [10])
    candidate_bam = write_bam(candidate, relative, [10])
    candidate_index = Path(f"{candidate_bam}.bai")
    external_index = tmp_path / "outside.bai"
    candidate_index.replace(external_index)
    candidate_index.symlink_to(external_index)

    with pytest.raises(COMPARATOR.ComparisonError, match="BAM index symlink escapes"):
        COMPARATOR.compare_runs(baseline, candidate)


def test_output_symlink_loop_is_reported_as_a_clean_setup_error(tmp_path, capsys):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    report_path = tmp_path / "report.json"
    write_text_result(baseline, "readcounts/counts.csv", "Identifier,count\na,1\n")
    loop = candidate / "tracks" / "loop.gff"
    loop.parent.mkdir(parents=True)
    loop.symlink_to(loop.name)

    return_code = COMPARATOR.main(
        [str(baseline), str(candidate), "--report", str(report_path)]
    )

    assert return_code == 2
    assert "cannot resolve output 'tracks/loop.gff'" in capsys.readouterr().err
    assert not report_path.exists()


def test_comparator_rejects_roots_without_supported_outputs(tmp_path):
    baseline = tmp_path / "baseline"
    candidate = tmp_path / "candidate"
    baseline.mkdir()
    candidate.mkdir()

    with pytest.raises(COMPARATOR.ComparisonError, match="no supported outputs"):
        COMPARATOR.compare_runs(baseline, candidate)
