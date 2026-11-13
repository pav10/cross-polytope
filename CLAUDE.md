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
src/cp4/covectors.py  certified covectors of the contact arrangement; LP-free column pre-filter of B1 (Lemma B.3)
src/cp4/critical.py   Stage B2 (eliminant of critical values of det on an affine family)
src/cp4/merged.py     merged Stage A + B1 search prototype (PLAN §4; infeasible frame-first, DEAD_ENDS B13)
scripts/              make_strata, zeta_frame, stageA_run, stageB_faces_zeta, stageB_critical_zeta, classify_catalog,
                      probe_A6_frame_facets, merged_run (--estimate K = Knuth tree-size estimate), merged_witness_zeta,
                      check_covectors_zeta (G8), compare_prefilter_zeta (pre-filter vs none)
data/                 strata.json, zeta.json, audit_certificate.json, wolff test sets, catalog
results/              outputs + logs of runs
verify/               verify_audit.py (external), run_all.py (regression gates), verify_support_reduction /
                      verify_structural_exclusions / verify_coordinate_sections (exact stratum exclusions, gate G9)
```

## Commands
```
python3 verify/run_all.py              # all cheap gates, ~2 min; must print ALL PASS
python3 scripts/stageA_run.py <sid>    # Stage A for one stratum -> results/stageA_<sid>.json
python3 scripts/stageB_faces_zeta.py [seconds] [--facet-only] [--static] [--determine] [--prefilter]
                                     [--frame-symmetry] [--estimate K] [--resume] [--out F]
                                                 # Stage B1 on zeta's S (G6 / G6f); facet-only static --prefilter ~6 min
sbatch scripts/slurm_stageB_zeta.sbatch          # full B1 on zeta's S on the cluster (resumable chain)
python3 scripts/stageB_critical_zeta.py          # Stage B2 gate (zeta cubic, 21 s)
```
Dependencies: python3, numpy, scipy (HiGHS via linprog), sympy ≥ 1.12, and `highspy` (used by `cp4.lp` when
installed: ~3× faster LPs, same verdicts — `scripts/check_lp_backend.py`). Recommended to add: `msolve`
(B2 for dim A = 5, 6).

## Environment
* Pavlo's cluster: SLURM, partition `ffa`, 5–10 GB RAM, 12 h wall time. Pavlo runs cluster jobs; Claude writes
  code and math. Give him a single sbatch script per job array, with checkpointing and heartbeat.
* The claude.ai Project "Cross-polytopes" holds the literature PDFs and older scripts (several superseded —
  see DEAD_ENDS A1–A4 before using any of them).

## Where we are (10 Nov 2026)
**Theory first:** 26 of the 39 strata are excluded by proofs with exact certificates — THEORY_REVIEW (15),
BLOCK_COUPLING + FACET_LAYERS (7), COORDINATE_SECTIONS (6, 7, 25), COMMON_FACET §3c (30, det ≥ 196/225).
**13 remain: 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17** (STATUS §3). Search rejections by float LPs are not
yet certified (THEORY_REVIEW §4). Priority target: stratum 11 (ζ's stratum) by theory, which would remove its search entirely (docs/STRATUM11.md:
7-point reformulation, Conjecture 11.A). Next: stratum 17 by the stratum-30 programme
(COMMON_FACET §6.3); extend the fiber-offset argument (COORDINATE_SECTIONS §6) to one axis row
(strata 1, 8, 9, 10, 15).
Gates G0–G3, G5, G6a, G6f, G7a, G8, G9 pass. Stage A is too weak on mid-size strata (≥ 9 800 contact sets in stratum 11
and growing). The PLAN §2 decision rule has been run: option 4(a) kills nothing (DEAD_ENDS B12) and the
frame-first merged search is ~60–1 000× Stage A's tree on stratum 0 (DEAD_ENDS B13). B1 remedies (PLAN §3):
highspy LPs, static order, facet-only pass, determine-then-check, resume. The facet-only pass on ζ's S is complete
and contains ζ's datum (G6f). The column-wise (covector) alternative is ~10²–10³× worse (no prune before the
4th column; DEAD_ENDS B15), but the covector list as an LP-free pre-filter inside row-wise B1 is now built in
(`prefilter=True` / `--prefilter`, exact certificate, STATUS §3 Lemma B.3, gate G8): identical trees and leaf
sets, 4.2× (facet-only) and 6–10× (full, static order). The full run (G6) has found ζ's datum (46 leaves after
≈ 1.73 M nodes, 288 of 1 633 depth-3 subtrees) but is not complete; the Knuth estimate was ≥ 3× low (heavy tail).
Finish it on the cluster from the committed checkpoint (`scripts/slurm_stageB_zeta.sbatch`; long runs die in the
cloud container, DEAD_ENDS C).
B1 has frame-symmetry isomorph rejection (`frame_group`, PLAN §3): ζ's S ~9×; stratum 0 needs only the star
A4 set, ≈ 4·10⁴ nodes. Stage A refuses excluded strata (`cp4.strata.EXCLUDED`).
**Next work: finish G6 (cluster); B1 is still far too slow campaign-wide
(~10⁴ S per stratum) — a structurally better enumeration is needed (PLAN §3, last paragraph); Knuth-estimate
any new search before building it.**
