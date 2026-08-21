#!/bin/bash
# Run HRIBO locally. Execute from your project directory, with HRIBO cloned into it.
snakemake -p -k \
    --sdm conda apptainer \
    --resources reparation_instances=1 \
    --scheduler-greediness 0 \
    -s HRIBO/workflow/Snakefile \
    --directory "${PWD}" \
    -j 5 \
    --latency-wait 60
