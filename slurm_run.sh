#!/bin/bash
# Submit HRIBO to SLURM. Run from your project directory, with HRIBO cloned into it.
#
# The SLURM executor plugin does not use a custom job script, so anything the
# compute nodes need must be in the environment at submission time; sbatch
# propagates it with its default --export=ALL.
module load gcc12-env/12.1.0
module load miniconda3/4.12.0
module load singularity/3.8.7

export SINGULARITY_BIND="${SINGULARITY_BIND},${TMPDIR}"

conda activate snakemake

snakemake \
    -s HRIBO/workflow/Snakefile \
    --directory "${PWD}" \
    --profile HRIBO/workflow/profiles/slurm
