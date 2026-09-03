#!/usr/bin/env python3
"""Apply HRIBO's audited compatibility fixes to the pinned deltaTE script.

The deltaTE container predates support for small gene sets and assumes that
every regulatory class contains at least one gene.  Rather than redistribute
the third-party script, this tool verifies the exact installed source and
materializes a narrowly patched copy for one workflow run.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import stat
import tempfile
from pathlib import Path


EXPECTED_SHA256 = (
    "e214c9f59de11aab71299fee512ebff36ac0dd4542bacdff4818263498df02b6"
)
PLOT_MARKER = "### Examples for each class of genes"
SCATTER_MAXIMUM = "max_val = max(res_ribo[,2],res_rna[,2],na.rm = T)"

SAFE_EXAMPLE_PLOTS = r'''### Examples for each class of genes
plot_class_example <- function(genes, class_name) {
  if (length(genes) == 0L) {
    plot.new()
    title(main=paste(class_name, "gene"))
    text(0.5, 0.5, "No genes in this class")
    return(invisible(NULL))
  }

  goi <- genes[[1L]]
  values <- c(
    RibOnly=as.numeric(res[goi,2]),
    Ribo=as.numeric(res_ribo[goi,2]),
    RNA=as.numeric(res_rna[goi,2])
  )
  if (any(!is.finite(values))) {
    plot.new()
    title(main=paste(class_name, "gene"))
    text(0.5, 0.5, "No finite fold changes")
    return(invisible(NULL))
  }

  limits <- range(c(0, values))
  if (limits[1] == limits[2]) {
    padding <- max(1, abs(limits[1]) * 0.1)
    limits <- limits + c(-padding, padding)
  }
  plot(c(1,2), c(0,values[["RibOnly"]]), type="l", col="red", xaxt="n",
       xlab="Conditions", ylim=limits, ylab="Log2 Fold Change",
       main=paste(class_name, "gene"))
  lines(c(1,2), c(0,values[["Ribo"]]), col="gray")
  lines(c(1,2), c(0,values[["RNA"]]), col="blue")
  axis(1, at=c(1,2), labels=c(1,2), las=1)
  legend("bottomleft", c("RNA", "Ribo", "RibOnly"),
         fill=c("blue", "gray", "red"), cex=1, border=NA, bty="n")
  invisible(NULL)
}

par(mfrow=c(2,2))
plot_class_example(forwarded, "Forwarded")
plot_class_example(exclusive, "Exclusive")
plot_class_example(buffered, "Buffered")
plot_class_example(intensified, "Intensified")

dev.off()'''


class PatchError(RuntimeError):
    """The installed deltaTE source is unexpected or cannot be patched safely."""


def replace_exact(text: str, old: str, new: str, expected_count: int) -> str:
    """Replace an audited expression only when its occurrence count is exact."""

    actual_count = text.count(old)
    if actual_count != expected_count:
        raise PatchError(
            f"expected {expected_count} occurrences of {old!r}, found {actual_count}"
        )
    return text.replace(old, new)


def patch_source(text: str) -> str:
    """Return the compatibility-patched source after structural checks."""

    text = replace_exact(
        text,
        "vsd <- vst(ddsMat_ribo)",
        "vsd <- varianceStabilizingTransformation(ddsMat_ribo, blind=TRUE)",
        2,
    )
    text = replace_exact(
        text,
        SCATTER_MAXIMUM,
        "\n".join(
            (
                "max_val = max(abs(c(res_ribo[,2],res_rna[,2])),na.rm = T)",
                "if (!is.finite(max_val) || max_val == 0) {",
                "  max_val <- 1",
                "}",
            )
        ),
        1,
    )
    text = replace_exact(
        text,
        "vsd <- vst(ddsMat_rna)",
        "vsd <- varianceStabilizingTransformation(ddsMat_rna, blind=TRUE)",
        2,
    )

    if text.count(PLOT_MARKER) != 1:
        raise PatchError(
            f"expected one deltaTE example-plot marker, found {text.count(PLOT_MARKER)}"
        )
    marker_start = text.index(PLOT_MARKER)
    terminator = "dev.off()"
    terminator_start = text.find(terminator, marker_start)
    if terminator_start < 0:
        raise PatchError("deltaTE example-plot block has no dev.off() terminator")
    block_end = terminator_start + len(terminator)
    if text[block_end:].strip():
        raise PatchError("unexpected source follows deltaTE's example-plot block")

    return text[:marker_start] + SAFE_EXAMPLE_PLOTS + "\n"


def materialize_patched_script(
    source: Path,
    output: Path,
    expected_sha256: str = EXPECTED_SHA256,
) -> None:
    """Verify, patch, and atomically publish an executable deltaTE script."""

    try:
        source_bytes = source.read_bytes()
    except OSError as exc:
        raise PatchError(f"cannot read deltaTE source {source}: {exc}") from exc
    actual_sha256 = hashlib.sha256(source_bytes).hexdigest()
    if actual_sha256 != expected_sha256:
        raise PatchError(
            f"unexpected deltaTE source at {source}: SHA-256 {actual_sha256}, "
            f"expected {expected_sha256}"
        )
    try:
        text = source_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PatchError(f"deltaTE source {source} is not UTF-8") from exc
    patched = patch_source(text).encode("utf-8")

    temporary_path: Path | None = None
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=f".{output.name}.",
            dir=output.parent,
            delete=False,
        ) as handle:
            temporary_path = Path(handle.name)
            handle.write(patched)
            handle.flush()
            os.fsync(handle.fileno())
        source_mode = stat.S_IMODE(source.stat().st_mode)
        temporary_path.chmod(source_mode | 0o111)
        temporary_path.replace(output)
        temporary_path = None
    except OSError as exc:
        raise PatchError(f"cannot publish patched deltaTE script {output}: {exc}") from exc
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Verify and patch the deltaTE script installed in its container."
    )
    parser.add_argument("source", type=Path, help="installed DTEG.R path")
    parser.add_argument("output", type=Path, help="patched output path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        materialize_patched_script(args.source, args.output)
    except PatchError as exc:
        raise SystemExit(f"patch_deltate: error: {exc}") from exc


if __name__ == "__main__":
    main()
