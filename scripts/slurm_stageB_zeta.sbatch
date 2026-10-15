#!/bin/bash
# Stage B1 FULL run on zeta's contact set (gate G6), resumable.  Submit once from the campaign root:
#   sbatch scripts/slurm_stageB_zeta.sbatch
# A chain of up to 6 sequential 12 h tasks (array %1).  Each task resumes from the checkpoint
# results/stageB_faces_zeta.pkl (rewritten atomically at every heartbeat), stops itself after 11.5 h, and
# exits at once if the run is already complete.  When done, `python3 verify/run_all.py` reports G6.
# Expected cost: Knuth estimate ~5.6e5 nodes (results/stageB_faces_zeta_estimate.json); measured ~12 nodes/s
# at the top of the tree; the facet-only pass took ~4.7x its estimate, so plan for roughly 10-40 h.
# Needs python3 with numpy, scipy, sympy >= 1.12 and, strongly recommended (~3x faster LPs), highspy:
#   python3 -m pip install --user highspy
#SBATCH --partition=ffa
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --array=0-5%1
#SBATCH --job-name=cp4B1zeta
#SBATCH --output=results/stageB_faces_zeta_full_%a.log
python3 -c 'import highspy' 2>/dev/null || echo "WARNING: highspy missing - scipy fallback, about 3x slower"
python3 -u scripts/stageB_faces_zeta.py 41400 --static --determine --resume --heartbeat 2000
