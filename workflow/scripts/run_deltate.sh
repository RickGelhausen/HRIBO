#!/usr/bin/env bash

# Execute deltaTE while keeping "the tool failed" distinct from "Snakemake saw
# an output pathname". The caller validates replication before scheduling this
# script, so every scheduled run must produce real, non-empty result artifacts.
set -euo pipefail

if [[ $# -ne 10 ]]; then
    echo "Usage: run_deltate.sh RIBO_COUNTS RNA_COUNTS SAMPLES RESULT_DIR RIBO_OUT RNA_OUT TE_OUT FIGURE_SOURCE FIGURE_OUT DTEG_SCRIPT" >&2
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
engine=${10}

# A retry must not be able to satisfy the job with artifacts left by an older
# successful run. Snakemake removes declared outputs after failures; the source
# figure is internal to deltaTE, so it is cleared explicitly as well.
rm -f -- "$ribo_output" "$rna_output" "$te_output" "$figure_source" "$figure_output"

"$engine" "$ribo_counts" "$rna_counts" "$samples" 0 "$result_dir/"

for required in "$ribo_output" "$rna_output" "$te_output" "$figure_source"; do
    if [[ ! -s "$required" ]]; then
        echo "deltaTE completed without required non-empty output: $required" >&2
        exit 1
    fi
done

validate_table() {
    local path=$1
    local expected_header=$2
    local expected_fields=$3
    local header
    IFS= read -r header < "$path"
    if [[ "$header" != "$expected_header" ]]; then
        echo "deltaTE output has an unexpected header: $path" >&2
        exit 1
    fi
    if ! awk -F '\t' -v fields="$expected_fields" \
        'NR > 1 && (NF != fields || $1 == "") { bad = 1 } END { exit (NR < 2 || bad) }' \
        "$path"; then
        echo "deltaTE output has no valid data rows: $path" >&2
        exit 1
    fi
}

validate_table "$ribo_output" $'baseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj' 6
validate_table "$rna_output" $'baseMean\tlog2FoldChange\tlfcSE\tpvalue\tpadj' 6
validate_table "$te_output" $'baseMean\tlog2FoldChange\tlfcSE\tstat\tpvalue\tpadj' 7

pdf_magic=$(head -c 5 -- "$figure_source")
if [[ "$pdf_magic" != '%PDF-' ]]; then
    echo "deltaTE figure is not a PDF: $figure_source" >&2
    exit 1
fi

cp -- "$figure_source" "$figure_output"
