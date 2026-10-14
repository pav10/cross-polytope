# CLAUDE.md — brief for Claude Code sessions on the C₄ lattice-packing campaign

## Mission
Prove **Theorem Z**: the densest lattice packing of the 4-dimensional cross-polytope
`C₄ = {‖x‖₁ ≤ 1}` has determinant `D_pack = D* = (16−9ζ+6ζ²)/16 = 0.80822000090353537…`
(`ζ` = root of `z³−6z+3` in `(1/2,2/3)`), attained only by signed-permutation images of `Λ_ζ = ℤ⁴B(ζ)`.
Owner: Pavlo Yatsyna (Charles University). Everything except the finite computation is already proved.

## Read order (do not skip; do not re-prove)
1. `docs/STATUS.md` — every established fact with its tag ([P], [P+C], [Cite], [Num], [Conj]) and §0 conventions.
2. `docs/PLAN.md` — the reduction tree, the gates, Stages A/B/C, measured costs, decision rules.
3. `docs/DEAD_ENDS.md` — false statements and failed computations. **Never revive one without new math.**
4. As needed: `docs/THEORY.md` (full proofs), `docs/audit_note.tex` (external audit), `docs/SLICING.md`, `docs/REFERENCES.md`.

## Hard rules
* All thresholds are **determinants** (`D*`, window `[2p/3, p·D*)` in frame coordinates), never densities.
* Accept/reject decisions in **exact arithmetic** (`cp4.exact`, `fractions`, sympy rationals). Floats only to
  search or to prune conservatively (the LP prunes in B1 use a strict-slack test; closure LPs are conservative).
* Admissibility via `cp4.lattice.exact_short_vectors` (rigorous coefficient bound). **Never** fixed coefficient boxes.
* Validate every new tool on a known answer first (gates in PLAN §1; `python3 verify/run_all.py`).
* Every long job prints a heartbeat and has a time limit or checkpoint; launch detached:
  `setsid nohup python3 -u script.py args > results/x.log 2>&1 < /dev/null &` and poll the log.
* Separate [P] / [P+C] / [Num] in everything written. Record new results in `docs/STATUS.md` (tagged) and new
  failures in `docs/DEAD_ENDS.md` with the measured reason.
* Do not edit `data/audit_certificate.json`, `docs/audit_note.tex`, `verify/verify_audit.py` (read-only, external).

## Claude Code workflow
* Use Plan Mode before new features, significant design tradeoffs, and bugs that survive the first attempted fix.

## Layout
```
src/cp4/exact.py      exact det/inv/rref/nullspace/solve_affine
src/cp4/lattice.py    conventions, M(z), B(z), ZETA, DSTAR, exact & fast short vectors
src/cp4/strata.py     HNFs, Lee criterion, 39 strata, stabilisers, candidate points
src/cp4/contacts.py   Stage A (contact sets S; Lemma A.1 / Cor A.2)
src/cp4/faces.py      Stage B1 (face data, LP pruning, ambient-orbit reduction, affine families)
src/cp4/critical.py   Stage B2 (eliminant of critical values of det on an affine family)
scripts/              make_strata, zeta_frame, stageA_run, stageB_faces_zeta, stageB_critical_zeta, classify_catalog
data/                 strata.json, zeta.json, audit_certificate.json, wolff test sets, catalog
results/              outputs + logs of runs
verify/               verify_audit.py (external), run_all.py (regression gates)
```

## Commands
```
python3 verify/run_all.py              # all cheap gates, ~2 min; must print ALL PASS
python3 scripts/stageA_run.py <sid>    # Stage A for one stratum -> results/stageA_<sid>.json
python3 scripts/stageB_faces_zeta.py [seconds]   # Stage B1 gate on zeta's S
python3 scripts/stageB_critical_zeta.py          # Stage B2 gate (zeta cubic, 21 s)
```
Dependencies: python3, numpy, scipy (HiGHS via linprog), sympy ≥ 1.12. Recommended to add: `highspy`
(warm-started LPs), `msolve` (B2 for dim A = 5, 6).

## Environment
* Pavlo's cluster: SLURM, partition `ffa`, 5–10 GB RAM, 12 h wall time. Pavlo runs cluster jobs; Claude writes
  code and math. Give him a single sbatch script per job array, with checkpointing and heartbeat.
* The claude.ai Project "Cross-polytopes" holds the literature PDFs and older scripts (several superseded —
  see DEAD_ENDS A1–A4 before using any of them).

## Where we are (14 Oct 2026)
Gates G0–G3, G5 pass. Stage A is too weak on mid-size strata (≥ 7 000 contact sets in stratum 11 and growing);
Stage B1 costs ~5 nodes/s with an unknown, large tree per S. **Next work = PLAN §2 decision rule and §3
remedies (merged search §4), measured on stratum 11 before any full run.**
