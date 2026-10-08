#!/usr/bin/env python3
"""Convert existing HRIBO advisor JSONs to ORFBounder read-length/offset JSONs."""

import argparse
import json
from pathlib import Path

from lib.orfbounder import export_orfbounder_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--recommendations", nargs="+", type=Path, required=True,
                        help="Existing tis_recommendation.json files, one per library.")
    parser.add_argument("-o", "--output-dir", type=Path, required=True,
                        help="Destination for end-specific JSON pairs and manifest.json.")
    args = parser.parse_args()
    try:
        payloads = [json.loads(path.read_text(encoding="utf-8")) for path in args.recommendations]
        manifest = export_orfbounder_inputs(payloads, args.output_dir)
    except (OSError, ValueError, TypeError) as error:
        parser.error(str(error))
    for end, result in manifest["read_ends"].items():
        print(f"{end}: {len(result['samples'])} calibrated libraries exported")


if __name__ == "__main__":
    main()
