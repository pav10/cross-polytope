#!/bin/bash
# Stage A for all 39 strata as a job array (Pavlo's cluster).  Submit from the campaign root:
#   sbatch --array=0-38 scripts/slurm_stageA.sbatch
# NOTE (PLAN §2): with (A1)-(A5) only, mid-size strata give ~1e4 sets and run for hours; the output is
# written only at the end, so a 12 h kill loses the stratum.  Strengthen Stage A first.
#SBATCH --partition=ffa
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=cp4A
#SBATCH --output=results/stageA_%a.log
python3 -u scripts/stageA_run.py ${SLURM_ARRAY_TASK_ID}
