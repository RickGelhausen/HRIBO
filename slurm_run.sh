#!/usr/bin/env bash
# Submit HRIBO to SLURM from a project directory. Activate the HRIBO launcher
# environment first; cluster-specific module setup belongs outside this script.
set -euo pipefail

hribo_dir=$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

exec snakemake \
    -s "${hribo_dir}/workflow/Snakefile" \
    --directory "${PWD}" \
    --profile "${hribo_dir}/workflow/profiles/slurm" \
    "$@"
