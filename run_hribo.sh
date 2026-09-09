#!/usr/bin/env bash
# Run HRIBO locally from a project directory. The checkout may live elsewhere.
set -euo pipefail

hribo_dir=$(CDPATH='' cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)

exec snakemake -p -k \
    --sdm conda apptainer \
    --resources reparation_instances=1 \
    --scheduler-greediness 0 \
    -s "${hribo_dir}/workflow/Snakefile" \
    --directory "${PWD}" \
    -j 5 \
    --latency-wait 60 \
    "$@"
