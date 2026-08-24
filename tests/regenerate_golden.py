#!/usr/bin/env python3
"""Regenerate the excel golden snapshots from the current scripts.

Run this only when an output change is intentional, and review the resulting
diff: the snapshots are plain CSV, so the change is readable.
"""

import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import excel_snapshot
import make_excel_fixture
from test_excel_outputs import SCRIPT_NAMES, command

GOLDEN = Path(__file__).resolve().parent / "golden"
SCRIPTS = Path(__file__).resolve().parent.parent / "workflow" / "scripts"


def main():
    GOLDEN.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as workspace:
        workspace = Path(workspace)
        inputs = make_excel_fixture.build(workspace / "inputs")
        for name in SCRIPT_NAMES:
            output = workspace / f"{name}.xlsx"
            result = subprocess.run(
                command(name, inputs, output), capture_output=True, text=True, cwd=str(SCRIPTS)
            )
            if result.returncode != 0:
                sys.exit(f"{name} failed:\n{result.stderr}")
            excel_snapshot.write_snapshot(output, GOLDEN / f"{name}.csv")
            print(f"updated tests/golden/{name}.csv")


if __name__ == "__main__":
    main()
