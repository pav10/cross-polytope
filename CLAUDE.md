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
src/cp4/contacts.py   Stage A (contact sets S; Lemma A.1 / Cor A.2; swap= frame-swap rule A6)
src/cp4/frameswap.py  exact frame-swap classification (FRAME_SWAP.md): FrameSwapRule(strata, sid).ok / pool_excluded
src/cp4/faces.py      Stage B1 (face data, LP pruning, ambient-orbit reduction, affine families; supporting=True,
                      ball=True: half-open facet cells + ||cY||<=2 rows, STATUS §3 Thm B.1', Lemmas B.4-B.5;
                      eutaxy=True: eutaxy leaf prune, Lemma B.6)
src/cp4/covectors.py  certified covectors of the contact arrangement; LP-free column pre-filter of B1 (Lemma B.3)
src/cp4/critical.py   Stage B2 (eliminant of critical values of det on an affine family)
src/cp4/merged.py     merged Stage A + B1 search prototype (PLAN §4; infeasible frame-first, DEAD_ENDS B13)
scripts/              make_strata, zeta_frame, stageA_run, stageB_faces_zeta, stageB_critical_zeta, classify_catalog,
                      probe_A6_frame_facets, merged_run (--estimate K = Knuth tree-size estimate), merged_witness_zeta,
                      check_covectors_zeta (G8), compare_prefilter_zeta (pre-filter vs none),
                      stageB_witness_supporting (G6s: zeta and D4 paths in supporting mode),
                      stageB_stratum (stratum-wide B1 driver, per-set checkpoints; --summary),
                      stageB2_stratum (B2 window check over B1 leaves), slurm_stageB_stratum.sbatch (job array)
data/                 strata.json, zeta.json, audit_certificate.json, wolff test sets, catalog
results/              outputs + logs of runs
verify/               verify_audit.py (external), run_all.py (regression gates), verify_support_reduction /
                      verify_structural_exclusions / verify_coordinate_sections (exact stratum exclusions, gate G9),
                      verify_frame_swap (FRAME_SWAP counts), verify_stratum11 (11.A counterexample) - both in G9
```

## Commands
```
python3 verify/run_all.py              # all cheap gates, ~2 min; must print ALL PASS
python3 scripts/stageA_run.py <sid>    # Stage A for one stratum (frame-swap rule on) -> results/stageA_<sid>.json
                                       # stratum 11: 1 469 sets in ~2 min
python3 scripts/stageB_faces_zeta.py [seconds] [--facet-only] [--static] [--determine] [--prefilter]
                                     [--frame-symmetry] [--supporting] [--ball] [--estimate K] [--resume] [--out F]
                                                 # Stage B1 on zeta's S (G6 / G6f / G6S); the complete default is
                                                 # --static --determine --prefilter --frame-symmetry --supporting --ball
python3 scripts/stageB_critical_zeta.py          # Stage B2 gate (zeta cubic, 21 s)
SID=11 sbatch --array=0-299 --export=ALL scripts/slurm_stageB_stratum.sbatch   # B1 over all sets (resubmit until done)
python3 scripts/stageB_stratum.py 11 --summary   # B1 progress; --only k runs one set locally
python3 scripts/stageB2_stratum.py 11            # B2 on the complete sets' leaves -> results/stageB2_11.json
```
Dependencies: python3, numpy, scipy (HiGHS via linprog), sympy ≥ 1.12, and `highspy` (used by `cp4.lp` when
installed: ~3× faster LPs, same verdicts — `scripts/check_lp_backend.py`). Recommended to add: `msolve`
(B2 for dim A = 5, 6).

## Environment
* Pavlo's cluster: SLURM, partition `ffa`, 5–10 GB RAM, 12 h wall time. Pavlo runs cluster jobs; Claude writes
  code and math. Give him a single sbatch script per job array, with checkpointing and heartbeat.
* The claude.ai Project "Cross-polytopes" holds the literature PDFs and older scripts (several superseded —
  see DEAD_ENDS A1–A4 before using any of them).

## Where we are (13 Nov 2026)
**Theory first:** 26 of the 39 strata are excluded by proofs with exact certificates — THEORY_REVIEW (15),
BLOCK_COUPLING + FACET_LAYERS (7), COORDINATE_SECTIONS (6, 7, 25), COMMON_FACET §3c (30, det ≥ 196/225).
**13 remain: 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17** (STATUS §3). Search rejections by float LPs are not
yet certified (THEORY_REVIEW §4).

**Stratum 11 (ζ's)** is now a search problem of cluster size:
* **Conjecture 11.A is false.** An exact frame of det 115/144 satisfies the 7-point and cube constraints; its short
  vector lies outside the cube (STRATUM11 §3, DEAD_ENDS A14, `verify/verify_stratum11.py`). So there is no
  cube-only relaxation proof.
* **Frame-swap lemma** (docs/FRAME_SWAP.md, `verify/verify_frame_swap.py`): unimodular quadruples of minimal
  vectors give other max-det frames, which prune Stage A. It is now Stage-A rule (A6) (`cp4.frameswap`,
  `enumerate_contact_sets(swap=...)`, on by default in `stageA_run.py`). Stratum 11: **1 469 sets in ~2 min**
  (`results/stageA_11.json`, gate G4).
* **Supporting-facet B1** (STATUS §3 Thm B.1′, Lemmas B.4–B.5; PLAN §3): non-frame contacts take only the 16
  facets in half-open cells, and every LP carries `‖cY‖₁ ≤ 2` for the undecided contacts. Complete; gates G6s/G6S.
  ζ's S is **complete in 55 613 nodes (38 min, 15 leaves, ζ found)**; the original B1 was incomplete at 1.7·10⁶.
  This is the new default (`supporting=True, ball=True`, static order; the LP-free filter order gave no gain,
  DEAD_ENDS B20).
* **Eutaxy leaf prune** (Lemma B.6, `eutaxy=True`): an RLT relaxation of THEORY Thm 5.2 drops 97 % of the
  (numerous, dim 5–6) leaves of small contact sets. It is useless at internal nodes (B19).
* **The pipeline is ready for the cluster** (PLAN §3):
  * `scripts/slurm_stageB_stratum.sbatch`: a job array over all sets, with per-set checkpoints; resubmit until
    `stageB_stratum.py 11 --summary` is complete;
  * then `scripts/stageB2_stratum.py 11` checks the window [2p/3, p·D*) exactly for the leaves of dim ≤ 3 (both
    signs of det). Leaves of dim ≥ 4 are *pending (msolve)*.
  * **On ζ's own S, B2 finds no critical value in the window except ζ itself** (both orientations; gate G5S).
  * Estimated B1 cost ≈ 1.1·10⁹ nodes, ≈ 3–6·10³ CPU-h (Knuth, 24 sets).
  * `--slack-tol 1e-6` (`SLACK_TOL=1e-6` for the sbatch) halves that, with the same leaves on ζ's S. Half the
    accepted LP cells are float noise around the 10⁻⁹ threshold (STATUS §4). This is Pavlo's call; the exact
    alternative is an exact LP for borderline cells.

**Next work:**
* Pavlo submits the stratum-11 B1 array.
* Install msolve on the cluster for dim 4–6 B2 leaves.
* Stage C for any window root.
* Certificates for the float LP prunes.
* Faster LPs (warm start) and prefix sharing across sets (PLAN §3).
* Theory: stratum 17 (COMMON_FACET §6.3) and the fiber-offset argument for strata 1, 8, 9, 10, 15.

The original full B1 on ζ's S (G6, 46 leaves, incomplete) is superseded by G6S; its cluster job
(`scripts/slurm_stageB_zeta.sbatch`) is no longer needed.

Earlier B1 history (PLAN §3): highspy, static order, determine-then-check, resume, the covector pre-filter (G8,
Lemma B.3), frame-symmetry isomorph rejection (`frame_group`; stratum 0 needs only the star A4 set,
≈ 4·10⁴ nodes). Rejected: option 4(a) (DEAD_ENDS B12), frame-first merged search (B13), column-wise B1 (B15),
closed facet cells (B18). Knuth-estimate any new search before building it.
