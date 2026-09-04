#!/usr/bin/env python3
"""Regenerate every HRIBO Linux Conda pin with the approved lock tool."""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Callable
from pathlib import Path

from check_conda_pins import PinError, conda_lock_hashes, validate_pin


PLATFORM = "linux-64"
CONDA_LOCK_VERSION = "4.0.2"
ROOT_ENVIRONMENTS = ("environment.yaml", "environment-dev.yaml")


def environment_files(repository: Path) -> list[Path]:
    """Return the stable set of launcher, development, and rule manifests."""
    environments = [repository / name for name in ROOT_ENVIRONMENTS]
    environments.extend(sorted((repository / "workflow" / "envs").glob("*.yaml")))
    missing = [path for path in environments if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing environment manifests: {missing}")
    return environments


def solver_arguments(solver: Path) -> list[str]:
    """Select conda-lock's matching backend for an explicit solver binary."""
    name = solver.name.lower()
    if "micromamba" in name:
        return ["--conda", str(solver), "--no-mamba", "--micromamba"]
    if "mamba" in name:
        return ["--conda", str(solver), "--mamba", "--no-micromamba"]
    return ["--conda", str(solver), "--no-mamba", "--no-micromamba"]


def regenerate(
    repository: Path,
    conda_lock: Path,
    solver: Path,
    hash_provider: Callable[[Path, Path], set[str]] = conda_lock_hashes,
) -> int:
    version = subprocess.run(
        [str(conda_lock), "--version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if version != f"conda-lock, version {CONDA_LOCK_VERSION}":
        raise RuntimeError(
            f"expected conda-lock {CONDA_LOCK_VERSION}, received {version!r}"
        )

    environments = environment_files(repository)
    virtual_packages = repository / "workflow" / "conda-lock" / "virtual-packages.yml"
    if not virtual_packages.is_file():
        raise RuntimeError(f"missing virtual-package specification: {virtual_packages}")

    process_environment = os.environ.copy()
    process_environment["CONDA_CHANNEL_PRIORITY"] = "strict"
    staging_parent = repository / ".snakemake"
    staging_parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(
        prefix="conda-pins-", dir=staging_parent
    ) as temporary:
        temporary_directory = Path(temporary)
        generated: list[tuple[Path, Path]] = []
        for manifest in environments:
            destination = manifest.with_suffix(f".{PLATFORM}.pin.txt")
            template = temporary_directory / f"{manifest.stem}.{{platform}}.pin.txt"
            command = [
                str(conda_lock),
                "lock",
                "--file",
                str(manifest),
                "--platform",
                PLATFORM,
                "--kind",
                "explicit",
                "--filename-template",
                str(template),
                "--virtual-package-spec",
                str(virtual_packages),
                "--strip-auth",
                *solver_arguments(solver),
            ]
            subprocess.run(
                command,
                cwd=repository,
                env=process_environment,
                check=True,
            )
            rendered = temporary_directory / f"{manifest.stem}.{PLATFORM}.pin.txt"
            if not rendered.is_file() or rendered.stat().st_size == 0:
                raise RuntimeError(f"conda-lock did not render {rendered.name}")
            validate_pin(manifest, rendered, virtual_packages, hash_provider)
            generated.append((rendered, destination))

        for rendered, destination in generated:
            rendered.replace(destination)
    return len(environments)


def executable(value: str) -> Path:
    resolved = shutil.which(value)
    if resolved is None:
        raise argparse.ArgumentTypeError(f"executable not found: {value}")
    return Path(resolved).resolve()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--repository",
        type=Path,
        default=Path(__file__).resolve().parents[2],
    )
    parser.add_argument("--conda-lock", type=executable, default="conda-lock")
    parser.add_argument("--solver", type=executable, default="mamba")
    args = parser.parse_args()
    try:
        count = regenerate(
            args.repository.resolve(),
            executable(args.conda_lock) if isinstance(args.conda_lock, str) else args.conda_lock,
            executable(args.solver) if isinstance(args.solver, str) else args.solver,
        )
    except (OSError, PinError, RuntimeError, subprocess.CalledProcessError) as error:
        print(f"Conda pin regeneration failed: {error}", file=sys.stderr)
        return 1
    print(f"Regenerated {count} {PLATFORM} Conda pins.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
