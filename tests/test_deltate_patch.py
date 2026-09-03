"""Regression tests for the checksum-guarded deltaTE compatibility patch."""

from __future__ import annotations

import hashlib
import stat
import subprocess

import pytest

from patch_deltate import PatchError, materialize_patched_script, patch_source


PLOT_BLOCK = """### Examples for each class of genes
par(mfrow=c(2,2))
goi = forwarded[1]
plot(res[goi,2])
goi = exclusive[1]
plot(res[goi,2])
goi = buffered[1]
plot(res[goi,2])
goi = intensified[1]
plot(res[goi,2])
dev.off()
"""


def synthetic_source() -> str:
    return "\n".join(
        (
            "#!/usr/bin/env Rscript",
            "vsd <- vst(ddsMat_ribo)",
            "vsd <- vst(ddsMat_ribo)",
            "vsd <- vst(ddsMat_rna)",
            "vsd <- vst(ddsMat_rna)",
            "max_val = max(res_ribo[,2],res_rna[,2],na.rm = T)",
            PLOT_BLOCK,
        )
    )


def test_patch_supports_small_inputs_and_empty_regulatory_classes():
    patched = patch_source(synthetic_source())

    assert "vst(ddsMat_" not in patched
    assert patched.count("varianceStabilizingTransformation(ddsMat_ribo") == 2
    assert patched.count("varianceStabilizingTransformation(ddsMat_rna") == 2
    assert "max(abs(c(res_ribo[,2],res_rna[,2]))" in patched
    assert "if (!is.finite(max_val) || max_val == 0)" in patched
    assert 'if (length(genes) == 0L)' in patched
    assert '"No genes in this class"' in patched
    assert "goi = forwarded[1]" not in patched
    assert 'plot_class_example(forwarded, "Forwarded")' in patched
    assert patched.rstrip().endswith("dev.off()")


def test_materializer_is_atomic_checksum_guarded_and_executable(tmp_path):
    source = tmp_path / "DTEG.R"
    output = tmp_path / "nested" / "patched DTEG.R"
    source.write_text(synthetic_source(), encoding="utf-8")
    expected_hash = hashlib.sha256(source.read_bytes()).hexdigest()

    materialize_patched_script(source, output, expected_hash)

    assert "varianceStabilizingTransformation" in output.read_text(encoding="utf-8")
    assert output.stat().st_mode & stat.S_IXUSR
    assert list(output.parent.glob(f".{output.name}.*")) == []
    syntax = subprocess.run(
        [
            "Rscript",
            "-e",
            "parse(file=commandArgs(trailingOnly=TRUE)[1])",
            str(output),
        ],
        capture_output=True,
        text=True,
    )
    assert syntax.returncode == 0, syntax.stdout + syntax.stderr

    output.write_text("previous verified output\n", encoding="utf-8")
    with pytest.raises(PatchError, match="unexpected deltaTE source"):
        materialize_patched_script(source, output, "0" * 64)
    assert output.read_text(encoding="utf-8") == "previous verified output\n"


def test_patch_rejects_structural_source_drift():
    missing_vst = synthetic_source().replace(
        "vsd <- vst(ddsMat_ribo)", "vsd <- changed(ddsMat_ribo)", 1
    )
    with pytest.raises(PatchError, match="expected 2 occurrences"):
        patch_source(missing_vst)

    missing_marker = synthetic_source().replace(
        "### Examples for each class of genes", "### changed marker"
    )
    with pytest.raises(PatchError, match="example-plot marker"):
        patch_source(missing_marker)
