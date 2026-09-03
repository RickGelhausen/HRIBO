#!/usr/bin/env python3
"""Parse shipped YAML and validate the example config/sample sheet schemas."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml
from snakemake.utils import validate


ROOT = Path(__file__).resolve().parents[2]
YAML_PATTERNS = (
    "*.yaml",
    "*.yml",
    "config/*.yaml",
    "config/*.yml",
    "workflow/envs/*.yaml",
    "workflow/envs/*.yml",
    "workflow/schemas/*.yaml",
    "workflow/schemas/*.yml",
    ".github/workflows/*.yaml",
    ".github/workflows/*.yml",
)


def yaml_files() -> list[Path]:
    """Return each maintained YAML file once, in stable display order."""
    return sorted({path for pattern in YAML_PATTERNS for path in ROOT.glob(pattern)})


def main() -> None:
    for path in yaml_files():
        with path.open(encoding="utf-8") as handle:
            yaml.safe_load(handle)

    config = yaml.safe_load((ROOT / "config/config.yaml").read_text(encoding="utf-8"))
    validate(config, str(ROOT / "workflow/schemas/config.schema.yaml"))

    samples = pd.read_csv(
        ROOT / "config/samples.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    validate(samples, str(ROOT / "workflow/schemas/samples.schema.yaml"))

    print(f"Parsed {len(yaml_files())} YAML files; config and samples match their schemas.")


if __name__ == "__main__":
    main()
