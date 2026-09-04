#!/usr/bin/env python3
"""Apply audited safety guards to DeepRibo's pinned S-curve estimator.

The estimator bundled with DeepRibo keeps removing observations until it finds
a lower bend, but it does not stop before too few observations remain.  This
launcher verifies the installed source byte-for-byte and creates a narrowly
patched copy for HRIBO's parameter-estimation rule.
"""

import argparse
import hashlib
import os
import stat
import tempfile
from pathlib import Path


EXPECTED_SHA256 = (
    "00bcb14bc71a55e79176614eb772e5811b2f3d35fbbe6e4c33dbe8e97c0f9750"
)

LOOP_START = "  while (MINCOV>0.60){"
MINIMUM_GUARD = """  while (MINCOV>0.60){
    if (fit_idx < 5L) {
      stop(paste(
        "DeepRibo S-curve did not reach coverage <= 0.60 before",
        "fewer than five observations remained"
      ))
    }"""

BEND_RESULT = "    MINCOV <- round(predict(bent_curve, MINRPKM), digits = 6)"
FINITE_GUARD = """    MINCOV <- round(predict(bent_curve, MINRPKM), digits = 6)
    if (!is.finite(MINRPKM) || !is.finite(MINCOV)) {
      stop("DeepRibo S-curve produced a non-finite bend")
    }"""

OLD_AXIS = 'plot(df$rpk_elo, df$coverage_elo, xlab="RPKM",ylab="Coverage")'
NEW_AXIS = (
    'plot(df$rpk_elo, df$coverage_elo, '
    'xlab="Log mean A-site occupancy per nucleotide",ylab="Coverage")'
)


class PatchError(RuntimeError):
    """The installed estimator is unexpected or cannot be patched safely."""


def replace_exact(text, old, new, expected_count=1):
    """Replace one audited fragment only when its occurrence count is exact."""

    actual_count = text.count(old)
    if actual_count != expected_count:
        raise PatchError(
            "expected {} occurrences of {!r}, found {}".format(
                expected_count, old, actual_count
            )
        )
    return text.replace(old, new)


def patch_source(text):
    """Return source with the bounded-loop, finite-result, and label fixes."""

    text = replace_exact(text, LOOP_START, MINIMUM_GUARD)
    text = replace_exact(text, BEND_RESULT, FINITE_GUARD)
    return replace_exact(text, OLD_AXIS, NEW_AXIS)


def materialize_patched_script(
    source, output, expected_sha256=EXPECTED_SHA256
):
    """Verify, patch, and atomically publish an executable R script."""

    source = Path(source)
    output = Path(output)
    try:
        source_bytes = source.read_bytes()
    except OSError as exc:
        raise PatchError("cannot read DeepRibo estimator {}: {}".format(source, exc))

    actual_sha256 = hashlib.sha256(source_bytes).hexdigest()
    if actual_sha256 != expected_sha256:
        raise PatchError(
            "unexpected DeepRibo estimator at {}: SHA-256 {}, expected {}".format(
                source, actual_sha256, expected_sha256
            )
        )
    try:
        text = source_bytes.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise PatchError("DeepRibo estimator {} is not UTF-8".format(source)) from exc

    patched = patch_source(text).encode("utf-8")
    temporary_path = None
    try:
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(
            mode="wb",
            prefix=".{}.".format(output.name),
            dir=str(output.parent),
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
        raise PatchError(
            "cannot publish patched DeepRibo estimator {}: {}".format(output, exc)
        ) from exc
    finally:
        if temporary_path is not None:
            try:
                temporary_path.unlink()
            except FileNotFoundError:
                pass


def parse_args():
    parser = argparse.ArgumentParser(
        description="Verify and patch DeepRibo's installed S-curve estimator."
    )
    parser.add_argument("source", type=Path, help="installed estimator path")
    parser.add_argument("output", type=Path, help="patched output path")
    return parser.parse_args()


def main():
    args = parse_args()
    try:
        materialize_patched_script(args.source, args.output)
    except PatchError as exc:
        raise SystemExit("patch_deepribo_scurve: error: {}".format(exc))


if __name__ == "__main__":
    main()
