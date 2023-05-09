#!/bin/bash
#properties = {properties}
module load gcc12-env/12.1.0
module load singularity/3.8.7
export SINGULARITY_BIND=$SINGULARITY_BIND,$TMPDIR
{exec_job}
