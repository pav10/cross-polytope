#!/bin/bash
# Stage B1 (supporting facets, STATUS §3 Thm B.1', Lemmas B.3-B.6) over ALL Stage-A contact sets of one stratum, as one
# job array.  Prerequisite: results/stageA_<SID>.json (python3 scripts/stageA_run.py <SID>; stratum 11: 1 469 sets,
# ~2 min).  Submit from the campaign root, e.g. for stratum 11:
#   SID=11 sbatch --array=0-299 --export=ALL scripts/slurm_stageB_stratum.sbatch
# Task i searches the sets k = i (mod array size).  Every set has its own checkpoint results/stageB_<SID>/set_<k>.pkl
# (rewritten atomically at every heartbeat); a task stops itself after 11.5 h.  RESUBMIT THE SAME COMMAND (same array
# size, after the previous array has finished - never two arrays at once) until
#   python3 scripts/stageB_stratum.py <SID> --summary
# reports "untouched": 0, "partial": 0.  Then run B2: python3 scripts/stageB2_stratum.py <SID>.
# Expected cost, stratum 11 (Knuth through the real code path, 24 of 1 469 sets; PLAN §3): ~1e9 nodes at ~50-70
# nodes/s, i.e. roughly 3-6e3 CPU-h before the filter order, less with it - about 10-20 resubmissions of a 300-task
# array.  Heavy-tailed: a few sets take far longer than the median; their tasks simply continue next time.
# Needs python3 with numpy, scipy, sympy >= 1.12 and, strongly recommended (~3x faster LPs), highspy:
#   python3 -m pip install --user highspy
#SBATCH --partition=ffa
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --job-name=cp4B1s
#SBATCH --output=results/stageB1_%A_%a.log
SID=${SID:-11}
NT=${NTASKS:-$SLURM_ARRAY_TASK_COUNT}
python3 -c 'import highspy' 2>/dev/null || echo "WARNING: highspy missing - scipy fallback, about 3x slower"
# optional: SLACK_TOL=1e-6 (2x fewer nodes on zeta's S, same leaves; STATUS §4) - must stay fixed across resubmissions
python3 -u scripts/stageB_stratum.py $SID --task $SLURM_ARRAY_TASK_ID --ntasks $NT --seconds 41400 --heartbeat 5000 \
    ${SLACK_TOL:+--slack-tol $SLACK_TOL}
