"""Read an .xlsx into comparable text, for golden-output testing.

The excel scripts are consumed by people reading spreadsheets, so a refactor has
to preserve them exactly. Comparing the binary files is useless (they embed
timestamps), so each sheet is rendered to CSV and compared as text.
"""

from pathlib import Path

import pandas as pd


def read_workbook(path):
    """{sheet name: DataFrame}, with values as written."""
    return pd.read_excel(path, sheet_name=None, engine="openpyxl")


def to_snapshot(path):
    """Render a workbook as one deterministic text blob."""
    sheets = read_workbook(path)
    parts = []
    for name in sheets:
        frame = sheets[name]
        parts.append(f"### sheet: {name} ({len(frame)} rows)")
        parts.append(frame.to_csv(index=False, lineterminator="\n").rstrip("\n"))
    return "\n".join(parts) + "\n"


def write_snapshot(xlsx_path, snapshot_path):
    Path(snapshot_path).parent.mkdir(parents=True, exist_ok=True)
    Path(snapshot_path).write_text(to_snapshot(xlsx_path))


def diff_summary(expected, actual, limit=15):
    """First differing lines, for a readable assertion message."""
    expected_lines = expected.splitlines()
    actual_lines = actual.splitlines()
    if len(expected_lines) != len(actual_lines):
        header = f"line count differs: expected {len(expected_lines)}, got {len(actual_lines)}"
    else:
        header = "content differs"

    differences = []
    for number, (left, right) in enumerate(zip(expected_lines, actual_lines), start=1):
        if left != right:
            differences.append(f"  line {number}:\n    expected: {left}\n    actual:   {right}")
        if len(differences) >= limit:
            differences.append("  ...")
            break
    return header + "\n" + "\n".join(differences)


if __name__ == "__main__":
    import sys

    write_snapshot(sys.argv[1], sys.argv[2])
    print(f"wrote {sys.argv[2]}")
