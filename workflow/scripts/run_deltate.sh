#!/usr/bin/env bash

# Execute deltaTE while keeping "the tool failed" distinct from "Snakemake saw
# an output pathname". The caller validates replication before scheduling this
# script, so every scheduled run must produce real, non-empty result artifacts.
set -euo pipefail

if [[ $# -ne 9 ]]; then
    echo "Usage: run_deltate.sh RIBO_COUNTS RNA_COUNTS SAMPLES RESULT_DIR RIBO_OUT RNA_OUT TE_OUT FIGURE_SOURCE FIGURE_OUT" >&2
    exit 64
fi

ribo_counts=$1
rna_counts=$2
samples=$3
result_dir=${4%/}
ribo_output=$5
rna_output=$6
te_output=$7
figure_source=$8
figure_output=$9

# A retry must not be able to satisfy the job with artifacts left by an older
# successful run. Snakemake removes declared outputs after failures; the source
# figure is internal to deltaTE, so it is cleared explicitly as well.
rm -f -- "$ribo_output" "$rna_output" "$te_output" "$figure_source" "$figure_output"

DTEG.R "$ribo_counts" "$rna_counts" "$samples" 0 "$result_dir/"

for required in "$ribo_output" "$rna_output" "$te_output" "$figure_source"; do
    if [[ ! -s "$required" ]]; then
        echo "deltaTE completed without required non-empty output: $required" >&2
        exit 1
    fi
done

cp -- "$figure_source" "$figure_output"
