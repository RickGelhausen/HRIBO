"""Regression tests for component-aware PCA reporting."""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from plot_PCA import plot_scatter_2D, plot_scatter_3D


pytest.importorskip("plotly")

REPO = Path(__file__).resolve().parent.parent
SCRIPT = REPO / "workflow" / "scripts" / "plot_PCA.py"
R_SCRIPT = REPO / "workflow" / "scripts" / "analyse_variance.R"


@pytest.mark.parametrize("component_count", [1, 2, 3])
def test_pca_report_uses_every_available_component(tmp_path, component_count):
    sample_count = max(2, component_count + 1)
    scores = {
        f"PC{index}": [
            float((sample - (sample_count - 1) / 2) * index)
            for sample in range(sample_count)
        ]
        for index in range(1, component_count + 1)
    }
    scores.update(
        {
            "group": [f"RIBO_{chr(65 + sample)}" for sample in range(sample_count)],
            "sampletype": [
                f"RIBO_{chr(65 + sample)}" for sample in range(sample_count)
            ],
            "name": [
                f"RIBO_{chr(65 + sample)}_1" for sample in range(sample_count)
            ],
        }
    )
    score_path = tmp_path / "rld.tsv"
    pd.DataFrame(scores).to_csv(score_path, sep="\t", index=False)

    variance_path = tmp_path / "variance.tsv"
    variance_path.write_text(
        "".join(f"{1 / component_count}\n" for _ in range(component_count))
    )
    correlation_path = tmp_path / "correlation.tsv"
    sample_names = scores["name"]
    correlation = [
        [1.0 if row == column else 0.8 for column in range(sample_count)]
        for row in range(sample_count)
    ]
    pd.DataFrame(
        correlation,
        index=sample_names,
        columns=sample_names,
    ).to_csv(correlation_path, sep="\t")

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "-r",
            str(score_path),
            "-p",
            str(variance_path),
            "-c",
            str(correlation_path),
            "-o",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr
    assert (tmp_path / "diffex_QC.html").stat().st_size > 0
    pca_output = (tmp_path / "PCA_3D.html").read_text()
    if component_count < 3:
        assert "showing the lower-dimensional PCA" in pca_output
    else:
        assert "showing the lower-dimensional PCA" not in pca_output
    assert "RIBO-A-1" in (tmp_path / "diffex_QC.html").read_text()


def test_pca_trace_dimension_matches_available_components():
    two_dimensional = pd.DataFrame(
        {
            "PC1": [-1.0, 1.0],
            "PC2": [0.0, 0.0],
            "group": ["RIBO_A", "RIBO_B"],
            "name": ["RIBO_A_1", "RIBO_B_1"],
        }
    )
    three_dimensional = two_dimensional.assign(PC3=[-0.5, 0.5])

    figure_2d = plot_scatter_2D(two_dimensional, [0.8, 0.2])
    figure_3d = plot_scatter_3D(three_dimensional, [0.7, 0.2, 0.1])

    assert {trace.type for trace in figure_2d.data} == {"scatter"}
    assert {trace.type for trace in figure_3d.data} == {"scatter3d"}


def test_r_preprocessing_selects_only_available_components():
    source = R_SCRIPT.read_text()

    assert "effective_rank <- sum(pca$sdev > tolerance)" in source
    assert "component_count <- min(3, max(1, effective_rank))" in source
    assert "pca$x[, 3]" not in source
