#!/usr/bin/env python3
"""Fail CI when pytest reports a skipped test in its JUnit XML."""

from __future__ import annotations

import argparse
import sys
import xml.etree.ElementTree as ET
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("junit_xml", type=Path)
    parser.add_argument("--expected-tests", type=int, required=True)
    args = parser.parse_args()

    root = ET.parse(args.junit_xml).getroot()
    cases = root.findall(".//testcase")
    if len(cases) != args.expected_tests:
        print(
            f"Expected {args.expected_tests} tests, but JUnit contains {len(cases)}. "
            "Update the CI baseline only after confirming the collection change.",
            file=sys.stderr,
        )
        return 1

    skipped = []
    for case in cases:
        marker = case.find("skipped")
        if marker is None:
            continue
        nodeid = "::".join(
            part
            for part in (case.get("classname"), case.get("name"))
            if part
        )
        reason = marker.get("message") or (marker.text or "").strip() or "no reason"
        skipped.append((nodeid, reason))

    if not skipped:
        print(f"pytest ran all {len(cases)} expected tests with no skips.")
        return 0

    print("CI requires every collected Python test to run; skipped tests:", file=sys.stderr)
    for nodeid, reason in skipped:
        print(f"- {nodeid}: {reason}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
