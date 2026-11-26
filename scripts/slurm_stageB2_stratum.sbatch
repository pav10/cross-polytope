#!/bin/bash
# Stage B2 over the B1 leaves of one stratum, as a job array.  Prerequisite: B1 complete
# (python3 scripts/stageB_stratum.py <SID> --summary shows "partial": 0, "untouched": 0).  Submit from the campaign
# root, e.g. for stratum 11:
#   SID=11 sbatch --array=0-199 --export=ALL scripts/slurm_stageB2_stratum.sbatch
# Task i handles the B1 sets k = i (mod array size); per-set results results/stageB2_<SID>/set_<k>.json are rewritten
# every 25 leaves and resumed, so RESUBMIT THE SAME COMMAND until every set is complete, then merge:
#   python3 scripts/stageB2_stratum.py <SID> --merge        # -> results/stageB2_<SID>.json (sets_done / sets_total)
# Leaves of dim 0-3 are decided exactly (sympy); dim >= 4 are only recorded as 'pending' (msolve, PLAN §3) unless
# MAXDIM is raised (dim 4 with sympy is slow; never change MAXDIM between resubmissions - the script refuses).
# Expected cost, stratum 11: 47 830 leaves of dim <= 3 (41 176 of dim 3 at ~20-50 s), i.e. a few hundred CPU-h.
#SBATCH --partition=ffa
#SBATCH --mem=6G
#SBATCH --time=12:00:00
#SBATCH --job-name=cp4B2
#SBATCH --output=results/stageB2_%A_%a.log
SID=${SID:-11}
NT=${NTASKS:-$SLURM_ARRAY_TASK_COUNT}
python3 -u scripts/stageB2_stratum.py $SID --task $SLURM_ARRAY_TASK_ID --ntasks $NT --seconds 41400 --max-dim ${MAXDIM:-3}
