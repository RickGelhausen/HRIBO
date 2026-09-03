"""Execute the production xTail and Riborex rules on a bounded count fixture."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml


REPO = Path(__file__).resolve().parent.parent
SNAKEFILE = REPO / "workflow" / "Snakefile"
CONTRAST = "A-B"
NEUTRAL_IDS = tuple(f"neutral_{index:03d}" for index in range(48))
UP_IDS = tuple(f"up_{index:03d}" for index in range(16))
DOWN_IDS = tuple(f"down_{index:03d}" for index in range(16))
LOW_COUNT_IDS = ("low_count_000",)
EFFECT_IDS = (*NEUTRAL_IDS, *UP_IDS, *DOWN_IDS)
ALL_IDS = (*EFFECT_IDS, *LOW_COUNT_IDS)

XTAIL_COLUMNS = (
    "mRNA_log2FC",
    "RPF_log2FC",
    "log2FC_TE_v1",
    "pvalue_v1",
    "log2FC_TE_v2",
    "pvalue_v2",
    "log2FC_TE_final",
    "pvalue_final",
    "pvalue.adjust",
)
RIBOREX_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)


def _negative_binomial_count(rng: np.random.Generator, mean: float) -> int:
    """Draw a reproducible overdispersed count while keeping filtered genes."""

    shape = 80.0
    probability = shape / (shape + mean)
    return max(2, int(rng.negative_binomial(shape, probability)))


def _write_count_fixture(path: Path) -> None:
    """Write balanced positive/negative TE effects without a library-size shift."""

    rng = np.random.default_rng(20260903)
    rows = []
    replicate_scales = {"1": 0.97, "2": 1.03}

    for gene_index, identifier in enumerate(EFFECT_IDS):
        if identifier.startswith("neutral_"):
            paired_index = gene_index
            a_ribo_effect = b_ribo_effect = 1.0
        elif identifier.startswith("up_"):
            paired_index = gene_index - len(NEUTRAL_IDS)
            a_ribo_effect, b_ribo_effect = 4.0, 1.0
        else:
            paired_index = gene_index - len(NEUTRAL_IDS) - len(UP_IDS)
            a_ribo_effect, b_ribo_effect = 1.0, 4.0

        # Up/down genes at the same paired index have identical expected base
        # abundance, keeping the condition-level library totals balanced.
        rna_mean = 180.0 + 11.0 * (paired_index % 17)
        ribo_mean = 150.0 + 13.0 * (paired_index % 19)
        row: dict[str, int | str] = {"Identifier": identifier}

        for condition, effect in (("A", a_ribo_effect), ("B", b_ribo_effect)):
            for replicate, scale in replicate_scales.items():
                row[f"RIBO-{condition}-{replicate}"] = _negative_binomial_count(
                    rng, ribo_mean * effect * scale
                )
                row[f"RNA-{condition}-{replicate}"] = _negative_binomial_count(
                    rng, rna_mean * scale
                )
        rows.append(row)

    # Both mean counts are below xTail 1.2.0's new implicit cutoff of 10 but
    # above HRIBO's historical cutoff of 1. This row makes preserving that
    # analysis boundary an executed (rather than merely static) assertion.
    rows.append(
        {
            "Identifier": LOW_COUNT_IDS[0],
            "RIBO-A-1": 4,
            "RNA-A-1": 5,
            "RIBO-A-2": 6,
            "RNA-A-2": 7,
            "RIBO-B-1": 4,
            "RNA-B-1": 5,
            "RIBO-B-2": 6,
            "RNA-B-2": 7,
        }
    )

    path.parent.mkdir(parents=True)
    pd.DataFrame(rows).to_csv(path, index=False)


def _write_diffex_config(
    path: Path,
    genome_file: Path,
    annotation_file: Path,
    samples: pd.DataFrame,
) -> None:
    sample_path = path.parent / "samples.tsv"
    samples.to_csv(sample_path, sep="\t", index=False, na_rep="")

    config = yaml.safe_load((REPO / "config" / "config.yaml").read_text())
    config["biologySettings"].update(
        {
            "genome": str(genome_file),
            "annotation": str(annotation_file),
            "samples": str(sample_path),
        }
    )
    config["differentialExpressionSettings"].update(
        {"contrasts": [CONTRAST], "xtailBins": 1000}
    )
    config["predictionSettings"]["deepribo"] = "off"
    config["workflowSettings"]["stages"] = ["differential_expression"]
    path.write_text(yaml.safe_dump(config, sort_keys=False))


def _assert_probability_columns(table: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce")
        assert values.notna().all(), f"{column} contained missing/non-numeric values"
        assert np.isfinite(values).all(), f"{column} contained non-finite values"
        assert values.between(0.0, 1.0, inclusive="both").all()


def test_diffex_rules_execute_real_xtail_and_riborex(
    snakemake_command,
    genome_file,
    annotation_file,
    samples,
    tmp_path,
):
    workdir = tmp_path / "differential expression run with spaces"
    workdir.mkdir()
    config_path = workdir / "fixture config.yaml"
    _write_diffex_config(config_path, genome_file, annotation_file, samples)
    _write_count_fixture(
        workdir / "readcounts" / "differential_expression_read_counts.csv"
    )

    configured_prefix = os.environ.get("HRIBO_TEST_CONDA_PREFIX")
    conda_prefix = (
        Path(configured_prefix)
        if configured_prefix
        else REPO / ".snakemake" / "conda"
    ).resolve()
    environment = {
        **os.environ,
        "PATH": os.pathsep.join(
            [str(Path(sys.executable).parent), os.environ.get("PATH", "")]
        ),
        "XDG_CACHE_HOME": str(workdir / ".cache"),
        "OPENBLAS_NUM_THREADS": "1",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }

    xtail_path = workdir / "xtail" / f"{CONTRAST}.csv"
    xtail_fc_plot = workdir / "xtail" / f"fc_{CONTRAST}.pdf"
    xtail_r_plot = workdir / "xtail" / f"r_{CONTRAST}.pdf"
    riborex_path = workdir / "riborex" / f"{CONTRAST}_deseq2.csv"
    xtail_conditions = (
        workdir
        / "diffex_input"
        / "xtail"
        / f"{CONTRAST}_condition_vector.csv"
    )
    command = [
        *snakemake_command,
        str(xtail_path.relative_to(workdir)),
        str(xtail_fc_plot.relative_to(workdir)),
        str(xtail_r_plot.relative_to(workdir)),
        str(riborex_path.relative_to(workdir)),
        "--cores",
        "1",
        "--printshellcmds",
        "--show-failed-logs",
        "--rerun-incomplete",
        "--notemp",
        "--software-deployment-method",
        "conda",
        "--conda-prefix",
        str(conda_prefix),
        "--allowed-rules",
        "contrastInput",
        "prepareXtailInput",
        "prepareRiborexInput",
        "xtail",
        "riborex",
        "--snakefile",
        str(SNAKEFILE),
        "--directory",
        str(workdir),
        "--configfile",
        str(config_path),
    ]
    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=environment,
        timeout=1200,
    )
    rendered = result.stdout + result.stderr
    assert result.returncode == 0, rendered
    assert "rule xtail" in rendered
    assert "rule riborex" in rendered
    assert not (workdir / "bam").exists(), "mapping escaped the fixture boundary"
    assert xtail_conditions.read_text() == "control,control,treated,treated\n"

    xtail = pd.read_csv(xtail_path, index_col=0)
    riborex = pd.read_csv(riborex_path, index_col=0)
    assert tuple(xtail.columns) == XTAIL_COLUMNS
    assert tuple(riborex.columns) == RIBOREX_COLUMNS
    assert len(xtail.index) == len(ALL_IDS)
    assert len(riborex.index) == len(ALL_IDS)
    assert set(xtail.index.astype(str)) == set(ALL_IDS)
    assert set(riborex.index.astype(str)) == set(ALL_IDS)

    _assert_probability_columns(
        xtail.loc[list(EFFECT_IDS)],
        ("pvalue_v1", "pvalue_v2", "pvalue_final", "pvalue.adjust"),
    )
    _assert_probability_columns(riborex.loc[list(EFFECT_IDS)], ("pvalue", "padj"))
    assert pd.to_numeric(
        xtail.loc[list(LOW_COUNT_IDS), "pvalue_final"], errors="coerce"
    ).notna().all()

    assert xtail.loc[list(UP_IDS), "log2FC_TE_final"].median() > 0
    assert xtail.loc[list(DOWN_IDS), "log2FC_TE_final"].median() < 0
    assert abs(xtail.loc[list(NEUTRAL_IDS), "log2FC_TE_final"].median()) < 0.5
    assert riborex.loc[list(UP_IDS), "log2FoldChange"].median() > 0
    assert riborex.loc[list(DOWN_IDS), "log2FoldChange"].median() < 0
    assert abs(riborex.loc[list(NEUTRAL_IDS), "log2FoldChange"].median()) < 0.5

    assert xtail_fc_plot.read_bytes().startswith(b"%PDF")
    assert xtail_r_plot.read_bytes().startswith(b"%PDF")

    outputs = (xtail_path, xtail_fc_plot, xtail_r_plot, riborex_path)
    before = {path: (path.read_bytes(), path.stat().st_mtime_ns) for path in outputs}
    rerun = subprocess.run(
        command,
        capture_output=True,
        text=True,
        env=environment,
        timeout=1200,
    )
    rerun_rendered = rerun.stdout + rerun.stderr
    assert rerun.returncode == 0, rerun_rendered
    assert "Nothing to be done" in rerun_rendered
    assert {
        path: (path.read_bytes(), path.stat().st_mtime_ns) for path in outputs
    } == before
