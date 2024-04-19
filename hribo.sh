#!/bin/bash

snakemake -p -k --use-conda --resources reparation_instances=1 --use-singularity --singularity-args "\"-c\"" --greediness 0 -s HRIBO/Snakefile --directory ${PWD} -j 5 --latency-wait 60