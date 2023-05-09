#!/bin/bash
module load gcc12-env/12.1.0
module load miniconda3/4.12.0

conda activate snakemake

snakemake --profile slurm_profile