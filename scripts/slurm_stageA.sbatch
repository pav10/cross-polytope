#!/bin/bash
# Stage A as a job array (Pavlo's cluster), for the 14 strata not excluded by proof (cp4.strata.REMAINING).
# Submit from the campaign root:
#   sbatch --array=0-5,8-11,14,15,17,30 scripts/slurm_stageA.sbatch
# (stageA_run.py refuses an excluded stratum, so a task outside this list exits at once.)
# NOTE (PLAN §2): with (A1)-(A5) only, mid-size strata give ~1e4 sets and run for hours; the output is
# written only at the end, so a 12 h kill loses the stratum.  Strengthen Stage A first.
#SBATCH --partition=ffa
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=cp4A
#SBATCH --output=results/stageA_%a.log
python3 -u scripts/stageA_run.py ${SLURM_ARRAY_TASK_ID}
