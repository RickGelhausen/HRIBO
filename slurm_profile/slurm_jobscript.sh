#!/bin/bash
#properties = {properties}
module load singularity/3.8.7
export SINGULARITY_BIND=$SINGULARITY_BIND,$TMPDIR
{exec_job}
