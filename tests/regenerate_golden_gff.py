#!/usr/bin/env python3
"""Regenerate the GFF golden snapshots from the current scripts.

Run only when an output change is intentional, and review the resulting diff.
"""

import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import make_gff_fixture
from test_gff_outputs import SCRIPT_NAMES, command

GOLDEN = Path(__file__).resolve().parent / "golden_gff"
SCRIPTS = Path(__file__).resolve().parent.parent / "workflow" / "scripts"


def main():
    GOLDEN.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as workspace:
        workspace = Path(workspace)
        inputs = make_gff_fixture.build(workspace / "inputs")
        for name in SCRIPT_NAMES:
            output = workspace / f"{name}.gff"
            result = subprocess.run(
                command(name, inputs, output), capture_output=True, text=True, cwd=str(SCRIPTS)
            )
            if result.returncode != 0:
                sys.exit(f"{name} failed:\n{result.stderr}")
            shutil.copy(output, GOLDEN / f"{name}.gff")
            print(f"updated tests/golden_gff/{name}.gff")


if __name__ == "__main__":
    main()
