"""Direct regression tests for differential-expression input and merge helpers."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
import yaml
from pandas.testing import assert_frame_equal


REPO = Path(__file__).resolve().parent.parent
PREPARE_DIFFEX = REPO / "workflow" / "scripts" / "prepare_diffex_input.py"
PREPARE_DELTATE = REPO / "workflow" / "scripts" / "prepare_deltate_input.py"
MERGE_DIFFEX = REPO / "workflow" / "scripts" / "merge_differential_expression.py"
DELTATE_RULE = (REPO / "workflow" / "rules" / "diffex_deltate.smk").read_text()
XTAIL_RULE = (REPO / "workflow" / "rules" / "diffex_xtail.smk").read_text()
XTAIL_SCRIPT = (REPO / "workflow" / "scripts" / "xtail.R").read_text()
XTAIL_ENV = (REPO / "workflow" / "envs" / "xtail.yaml").read_text()
XTAIL_DEPLOY = (
    REPO / "workflow" / "envs" / "xtail.post-deploy.sh"
).read_text()
RIBOREX_SCRIPT = (REPO / "workflow" / "scripts" / "riborex.R").read_text()


def run_script(script: Path, *arguments: object) -> subprocess.CompletedProcess[str]:
    """Run a repository helper exactly as its Snakemake rule does."""
    return subprocess.run(
        [sys.executable, str(script), *(str(argument) for argument in arguments)],
        capture_output=True,
        text=True,
    )


def test_xtail_resources_are_validated_and_forwarded():
    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    schema = yaml.safe_load(
        (REPO / "workflow" / "schemas" / "config.schema.yaml").read_text()
    )
    settings_schema = schema["properties"]["differentialExpressionSettings"]
    bins_schema = settings_schema["properties"]["xtailBins"]
    min_mean_schema = settings_schema["properties"]["xtailMinMeanCount"]

    assert config["differentialExpressionSettings"]["xtailBins"] == 10000
    assert bins_schema["type"] == "integer"
    assert bins_schema["minimum"] == 1
    assert bins_schema["default"] == 10000
    assert "xtailBins" not in settings_schema["required"]
    assert config["differentialExpressionSettings"]["xtailMinMeanCount"] == 1
    assert min_mean_schema["type"] == "number"
    assert min_mean_schema["minimum"] == 1
    assert min_mean_schema["default"] == 1
    assert "xtailMinMeanCount" not in settings_schema["required"]

    assert 'get("xtailBins", 10000)' in XTAIL_RULE
    assert 'get(\n            "xtailMinMeanCount", 1\n        )' in XTAIL_RULE
    assert "--threads {threads}" in XTAIL_RULE
    assert "--bins {params.bins:q}" in XTAIL_RULE
    assert "--min_mean_count {params.min_mean_count:q}" in XTAIL_RULE

    assert 'c("--threads"), type = "integer", default = 1L' in XTAIL_SCRIPT
    assert 'c("--bins"), type = "integer", default = 10000L' in XTAIL_SCRIPT
    assert 'c("--min_mean_count"), type = "double", default = 1' in XTAIL_SCRIPT
    assert "options$threads < 1L" in XTAIL_SCRIPT
    assert "options$bins < 1L" in XTAIL_SCRIPT
    assert "options$min_mean_count < 1" in XTAIL_SCRIPT
    assert "threads = options$threads" in XTAIL_SCRIPT
    assert "bins = options$bins" in XTAIL_SCRIPT
    assert "minMeanCount = options$min_mean_count" in XTAIL_SCRIPT
    assert "test.tab[, result.columns, drop = FALSE]" in XTAIL_SCRIPT
    assert 'default = "NULL"' not in XTAIL_SCRIPT

    assert "official 1.2.0 source release" in XTAIL_ENV
    assert "xtail_1.2.0-source.tar.gz" in XTAIL_DEPLOY
    assert (
        "5975e7b9ea692be69ebaee85a64f710ca742acf6f2464d961c764fe5f1cccc95"
        in XTAIL_DEPLOY
    )
    assert '--library="$target_library"' in XTAIL_DEPLOY
    assert 'find.package("xtail", lib.loc = library_path)' in XTAIL_DEPLOY
    assert 'packageVersion("xtail", lib.loc = library_path)' in XTAIL_DEPLOY


def test_riborex_cli_only_requires_the_result_it_writes():
    assert 'default = "NULL"' not in RIBOREX_SCRIPT
    assert "--riborexdeseq_result_path" in RIBOREX_SCRIPT
    assert "--riborexedgeR_result_path" not in RIBOREX_SCRIPT
    assert "--riborexvoom_result_path" not in RIBOREX_SCRIPT
    assert "(-r, -m, -c, -x)" in RIBOREX_SCRIPT


@pytest.mark.parametrize(
    ("tool", "expected_ribo_columns", "expected_rna_columns", "condition_vector"),
    [
        (
            "riborex",
            ["Identifier", "RIBO-A-2", "RIBO-A-1", "RIBO-B-2", "RIBO-B-1"],
            ["Identifier", "RNA-A-2", "RNA-A-1", "RNA-B-1", "RNA-B-2"],
            "treated,treated,control,control\n",
        ),
        (
            "xtail",
            ["Identifier", "RIBO-B-2", "RIBO-B-1", "RIBO-A-2", "RIBO-A-1"],
            ["Identifier", "RNA-B-1", "RNA-B-2", "RNA-A-2", "RNA-A-1"],
            "control,control,treated,treated\n",
        ),
    ],
)
def test_prepare_diffex_writes_ordered_count_tables_and_condition_vector(
    tmp_path, tool, expected_ribo_columns, expected_rna_columns, condition_vector
):
    read_counts = pd.DataFrame(
        {
            "Identifier": ["gene-1", "gene-2"],
            "RIBO-B-2": [42, 43],
            "RNA-A-2": [12, 13],
            "RIBO-A-2": [22, 23],
            "RNA-B-1": [31, 32],
            "RIBO-A-1": [21, 24],
            "RNA-A-1": [11, 14],
            "RIBO-B-1": [41, 44],
            "RNA-B-2": [33, 34],
        }
    )
    source = tmp_path / "all read counts.csv"
    output = tmp_path / f"{tool} input"
    read_counts.to_csv(source, index=False)

    result = run_script(
        PREPARE_DIFFEX,
        "--read_count_file",
        source,
        "--contrast",
        "A-B",
        "--output_folder",
        output,
        "--tool",
        tool,
    )

    assert result.returncode == 0, result.stderr
    ribo = pd.read_csv(output / "A-B_ribo_readcount_table.tsv", sep="\t")
    rna = pd.read_csv(output / "A-B_rna_readcount_table.tsv", sep="\t")
    assert_frame_equal(ribo, read_counts[expected_ribo_columns])
    assert_frame_equal(rna, read_counts[expected_rna_columns])
    assert (output / "A-B_condition_vector.csv").read_text() == condition_vector


def test_prepare_deltate_explicit_bams_ignore_stale_matching_files(tmp_path):
    bam_folder = tmp_path / "bam inputs"
    bam_folder.mkdir()
    selected_samples = [
        f"{method}-{condition}-{replicate}"
        for method in ("RIBO", "RNA")
        for condition in ("A", "B")
        for replicate in ("1", "2")
    ]
    for sample in [*selected_samples, "RIBO-B-9", "RIBO-unused-1"]:
        (bam_folder / f"{sample}.bam").touch()

    read_counts = pd.DataFrame(
        {
            "Identifier": ["gene-1", "gene-2"],
            **{
                sample: [index, index + 100]
                for index, sample in enumerate(selected_samples, start=1)
            },
            "RIBO-B-9": [998, 998],
            "RIBO-unused-1": [999, 999],
        }
    )
    source = tmp_path / "all counts.csv"
    output = tmp_path / "deltaTE input"
    read_counts.to_csv(source, index=False)

    result = run_script(
        PREPARE_DELTATE,
        "--bam_files",
        *(bam_folder / f"{sample}.bam" for sample in selected_samples),
        "--read_count_file",
        source,
        "--contrast",
        "A-B",
        "--output_folder",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert (output / "samples_info.txt").read_text() == (
        "SampleID\tCondition\tSeqType\tBatch\n"
        "RIBO-B-1\t1\tRIBO\t1\n"
        "RIBO-B-2\t1\tRIBO\t2\n"
        "RIBO-A-1\t2\tRIBO\t1\n"
        "RIBO-A-2\t2\tRIBO\t2\n"
        "RNA-B-1\t1\tRNA\t1\n"
        "RNA-B-2\t1\tRNA\t2\n"
        "RNA-A-1\t2\tRNA\t1\n"
        "RNA-A-2\t2\tRNA\t2\n"
    )
    assert (output / "has_replicates.txt").read_text() == "True"
    assert (output / "ribo_counts.txt").read_text() == (
        "RIBO-B-1\tRIBO-B-2\tRIBO-A-1\tRIBO-A-2\n"
        "gene-1\t3\t4\t1\t2\n"
        "gene-2\t103\t104\t101\t102\n"
    )
    assert (output / "rna_counts.txt").read_text() == (
        "RNA-B-1\tRNA-B-2\tRNA-A-1\tRNA-A-2\n"
        "gene-1\t7\t8\t5\t6\n"
        "gene-2\t107\t108\t105\t106\n"
    )


def test_prepare_deltate_rule_passes_only_declared_bams():
    assert "--bam_files {input.bam:q}" in DELTATE_RULE
    assert "-b bam/" not in DELTATE_RULE


def test_prepare_deltate_marks_an_unreplicated_design(tmp_path):
    bam_folder = tmp_path / "bam"
    bam_folder.mkdir()
    samples = ["RIBO-A-1", "RIBO-B-1", "RNA-A-1", "RNA-B-1"]
    for sample in samples:
        (bam_folder / f"{sample}.bam").touch()
    source = tmp_path / "counts.csv"
    pd.DataFrame(
        {"Identifier": ["gene"], **{sample: [1] for sample in samples}}
    ).to_csv(source, index=False)
    output = tmp_path / "output"

    result = run_script(
        PREPARE_DELTATE,
        "-b",
        bam_folder,
        "-r",
        source,
        "-c",
        "A-B",
        "-o",
        output,
    )

    assert result.returncode == 0, result.stderr
    assert (output / "has_replicates.txt").read_text() == "False"


MERGE_CASES = {
    "riborex": {
        "columns": [
            "Identifier",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
        ],
        "output_columns": [
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj",
            "contrast",
        ],
    },
    "xtail": {
        "columns": [
            "Identifier",
            "mRNA_log2FC",
            "RPF_log2FC",
            "log2FC_TE_v1",
            "pvalue_v1",
            "log2FC_TE_v2",
            "pvalue_v2",
            "log2FC_TE_final",
            "pvalue_final",
            "pvalue_adjusted",
        ],
        "output_columns": [
            "gene_id",
            "mRNA_log2FC",
            "RPF_log2FC",
            "log2FC_TE_v1",
            "pvalue_v1",
            "log2FC_TE_v2",
            "pvalue_v2",
            "log2FC_TE_final",
            "pvalue_final",
            "pvalue_adjust",
            "contrast",
        ],
    },
    "deltate": {
        "columns": [
            "Identifier",
            "RIBO_baseMean",
            "RIBO_log2FoldChange",
            "RIBO_lfcSE",
            "RIBO_pvalue",
            "RIBO_padj",
            "RNA_baseMean",
            "RNA_log2FoldChange",
            "RNA_lfcSE",
            "RNA_pvalue",
            "RNA_padj",
            "TE_baseMean",
            "TE_log2FoldChange",
            "TE_lfcSE",
            "TE_stat",
            "TE_pvalue",
            "TE_padj",
        ],
        "output_columns": [
            "gene_id",
            "RIBO_baseMean",
            "RIBO_log2FoldChange",
            "RIBO_lfcSE",
            "RIBO_pvalue",
            "RIBO_padj",
            "RNA_baseMean",
            "RNA_log2FoldChange",
            "RNA_lfcSE",
            "RNA_pvalue",
            "RNA_padj",
            "TE_baseMean",
            "TE_log2FoldChange",
            "TE_lfcSE",
            "TE_stat",
            "TE_pvalue",
            "TE_padj",
            "contrast",
        ],
    },
}


@pytest.mark.parametrize("tool", ["riborex", "xtail", "deltate"])
def test_merge_diffex_pools_every_contrast_with_its_origin(tmp_path, tool):
    case = MERGE_CASES[tool]
    value_columns = case["columns"][1:]
    first = pd.DataFrame(
        [
            ["gene-1", *range(10, 10 + len(value_columns))],
            ["shared", *range(20, 20 + len(value_columns))],
        ],
        columns=case["columns"],
    )
    second = pd.DataFrame(
        [
            ["shared", *range(30, 30 + len(value_columns))],
            ["gene-2", *range(40, 40 + len(value_columns))],
        ],
        columns=case["columns"],
    )
    first_path = tmp_path / "A-B_sorted.xlsx"
    second_path = tmp_path / "C-D_sorted.xlsx"
    first.to_excel(first_path, index=False)
    second.to_excel(second_path, index=False)
    output = tmp_path / f"{tool} pooled.csv"

    result = run_script(
        MERGE_DIFFEX,
        first_path,
        second_path,
        "--output_csv",
        output,
        "--tool",
        tool,
    )

    assert result.returncode == 0, result.stderr
    actual = pd.read_csv(output)
    expected = pd.DataFrame(
        [
            ["gene-1", *range(10, 10 + len(value_columns)), f"{tool}_A-B"],
            ["shared", *range(20, 20 + len(value_columns)), f"{tool}_A-B"],
            ["shared", *range(30, 30 + len(value_columns)), f"{tool}_C-D"],
            ["gene-2", *range(40, 40 + len(value_columns)), f"{tool}_C-D"],
        ],
        columns=case["output_columns"],
    )
    assert_frame_equal(actual, expected)
