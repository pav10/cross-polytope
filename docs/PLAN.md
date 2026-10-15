# PLAN — the computation that proves Theorem Z

Read `STATUS.md` first (what is proved), then this file (what to compute), then `DEAD_ENDS.md`
(what not to try). Every number below marked *measured* was produced by the scripts named.

## 0. The proof skeleton (all reductions already [P])

```
critical lattice Λ (exists, Thm 4.1)
  └─ max-det frame X, Γ = ΛX⁻¹ ⊇ ℤ⁴, index p            (Thm 6.1)
       └─ Γ is one of 39 strata                          (Thm 6.2)            → data/strata.json
            └─ contact set S ⊆ candidates(Γ) satisfying (A1)–(A5)   (Lemma A.1, Cor A.2)   → Stage A
                 └─ face datum (s_c)_{c∈S}, affine space A, dim A ≤ 6   (Thm B.1)          → Stage B1
                      └─ |det X| = |critical value of det|_A|          (Thm B.1)          → Stage B2
                           └─ Cor B.2: no critical value in [2p/3, p·D*)   ⇒ Theorem Z
```
Theorem Z ⟸ for every leaf, the eliminant R(D) has no real root in `[2p/3, p·D*)`, and the root `p·D*`
occurs only on leaves of Λ_ζ (stratum 11). Leaves violating this go to Stage C.

## 1. Gates (each already passed unless marked) — rerun with `python3 verify/run_all.py`

| gate | expected | status |
|---|---|---|
| G0 audit checker | PASS | passed |
| G1 strata | 39 orbits / 1040 superlattices / p-set; strict Lee: 6 / 54 | passed |
| G2 ζ located | stratum 11, p = 4, 15 contacts | passed |
| G3 Stage A, p = 1 | exactly 3 sets, all = A₄ (size 10) | passed (13 s) |
| G4 Stage A, stratum 11 contains ζ's S | ζ's canonical S ∈ output | see STATUS §4 |
| G5 Stage B2 on ζ's face datum | eliminant contains `64D³−1920D²+13584D−26011` with root `4·D*` | **passed** (21 s, sympy grevlex+FGLM); other factors `D−6`, `D³+62D²+8D+200` — no root in `[8/3, 4D*)` |
| G6a B1 witness | ζ's datum survives every prune contact by contact | **passed** |
| G6f Stage B1 facet-only pass on ζ's S | ζ's face datum among the facet-only leaves | **passed** (65 148 nodes, 23 min; 6 leaves, all dim 3) |
| G6 Stage B1 on ζ's S (complete run) | ζ's face datum among leaves (up to ambient group) | not run yet: ~10–40 h, `scripts/slurm_stageB_zeta.sbatch` (resumable) |
| G7a merged-search witness | ζ's labels (15 in, 53 out) survive every prune of `src/cp4/merged.py`; leaf dim A = 3 | **passed** (4 s) |

## 2. Stage A — contact sets (`src/cp4/contacts.py`, `scripts/stageA_run.py <sid>`)

Include/exclude branching over the candidate list with hull-closure forcing; output reduced mod the
stratum stabiliser; heartbeat every 20 000 nodes; output `results/stageA_<sid>.json`
(`contact_sets`: list of sets of numerator 4-tuples; `stats`).

*Measured:* stratum 0 — 54 726 nodes, 3 sets. Stratum 11 — ~1 000 nodes/s, ~10 sets per 1 000 nodes,
> 6 000 sets after 560 000 nodes and still growing. Pruning profile: `dead_det` ≈ 82 % of kills,
`dead_interior` ≈ 4 %, **volume (A4) never fires**, `dead_excluded_forced` ≈ 5 %.
**Conclusion: (A1)–(A5) alone leave 10⁴–10⁵ sets per mid-size stratum.** Stage A must be strengthened
before running all 39 strata. Options, in order of expected value:

1. **Realisability needs face data.** The convex relaxation "∃X: ‖cX‖₁ ≤ 2 on S" is trivially feasible
   (X = 0) and `‖eᵢX‖₁ = 2`, `F(c) ≥ 2` off S are non-convex (DEAD_ENDS A8). So the only strong
   realisability test is B1 itself ⇒ the **merged search of §4** is the main remedy, not a new Stage-A rule.
2. Already implied (do not re-add): the parity rule (Lemma 1.6 ⇒ for `c, c′ ∈ S` with `(c−c′)/2 ∈ Γ`,
   both `(c±c′)/2` are minimal) is a special case of hull closure (A3b) — the midpoint lies in P.
3. Rank ≥ 10 / eutaxy (Thm 5.1/5.2) need face data → B1.
4. Possible cheap additions to test: (a) the frame rows' sign vectors (≤ 15⁴/384 choices up to the
   ambient group) + LP "F(c) ≤ 2 on S with equality on the frame" — non-trivial and convex once the frame
   signs are fixed; (b) volume (A4) is useless as implemented — drop it or strengthen via Minkowski's
   second theorem.
5. Symmetry: results are canonical under the stabiliser of Γ; the ambient group acts in B1.

Decision rule: try 4(a) (≤ 1 day); re-measure stratum 11 (target < 2 000 sets). Otherwise implement
the merged search of §4.

*Outcome (14 Oct 2026):* **4(a) fails.** In its sound form (one closed facet `uᵢ ∈ {±1}⁴` per frame row,
330 patterns mod the ambient group) it kills **0** of the 200 Stage-A sets of stratum 0 and **0** of 931
of stratum 11 (`scripts/probe_A6_frame_facets.py`; DEAD_ENDS B12). The merged search of §4 was prototyped
next; in its frame-first form it is infeasible too (§4; DEAD_ENDS B13). Next step: undecided.

## 3. Stage B — face data and critical values

**B1** (`src/cp4/faces.py`, `scripts/stageB_faces_zeta.py`): for fixed S, branch over sign vectors
`s_c ∈ {−1,0,1}⁴∖0`; dynamic order = contact with fewest coordinatewise-allowed signs (8 closure LPs
per contact, conservative); orbit reduction under the subgroup of the 384 ambient signed permutations
fixing the chosen signs (stabiliser chain — complete and non-redundant; contacts are not moved by the
ambient group, so the dynamic order is equivariant); prune by strict LP feasibility (max common slack);
leaves: exact `solve_affine`, discard `dim A > 6`.
*Measured on ζ's S (15 contacts):* ~5 nodes/s; 7 596 nodes in 1 500 s, all dead by depth 14, no leaf
yet, still inside the first few depth-3 branches (8 at depth 3) ⇒ a full B1 run on one S is hours. With
~10⁴ sets S per stratum, **B1 as implemented is infeasible campaign-wide; this is the bottleneck.**
Structural reason: after ~13 facet equations the system has rank 16 (X determined), so most work is
enumerating sign vectors of 13 contacts that are later killed by the overdetermined remainder.
Remedies:

* **Determine-then-check.** Branch only until the equation rank reaches 16 − d for the target d (or
  until X is unique), then compute X (or the small affine family) exactly and *read off* the remaining
  sign vectors instead of branching on them; and choose the first contacts to maximise rank growth
  (independent equations first). Expected gain: the tree is cut at depth ≈ rank, not 15.
* **Eutaxy prune (Thm 5.2).** At any node where the contacts with fixed signs already have rank ≥ 10,
  the KKT condition `I ∈ cone{σ_cᵀ(cX)}` over *all* contacts is necessary; a Farkas certificate on the
  partial data does not prune soundly (more contacts can be added), but at a leaf it is a cheap exact check.

* **Cache by A.** The certificate depends only on the affine space A (Thm B.1): store eliminants keyed by
  the RREF of A's equations; different (S, face datum) with equal A cost nothing. Also: A is determined
  by any subset of contacts whose equations span; enumerate face data contact-by-contact and stop
  branching (go straight to B2) once `rank = 16 − dim A` stops growing? **No** — remaining contacts
  still constrain the cell (sign conditions); they matter for Stage C only. So: once rank reaches its
  final value is unknowable; keep full branching but cache B2.
* **Cheaper LPs:** replace 8 range LPs per undecided contact by one LP for the current point + a
  sign-pattern guess; or use a warm-started HiGHS model (`highspy`) instead of `scipy.linprog`.
* **Facet-only first pass:** at ζ every contact is interior to a facet (all `s_c` full). Leaves with
  some `s_cj = 0` have more equations ⇒ smaller A; enumerate them too (they are needed), but a
  facet-only pass gives the fastest complete picture for the dominant case.

*Measured on ζ's S (15 Oct 2026; Knuth estimates, 300 probes, `results/stageB_faces_zeta_estimate.json`):*
the original dynamic order spends ~116 LPs per node on closure ranges; a **static order** needs ~35 and
gives a tree of similar size (5.6·10⁵ vs 6.2·10⁵ nodes, medians ~3·10³), so it is ~2× cheaper overall.
**highspy** (`src/cp4/lp.py`) gives identical verdicts on 3 011 recorded LPs at 0.66 vs 2.06 ms/LP
(`scripts/check_lp_backend.py`). **Determine-then-check** is implemented (`determine=True`) but never fired
on ζ's S (no branch reached rank 16). A norm-range lookahead (min/max of `‖cX‖₁` over the cell for undecided
contacts) cut the tree only 2–3× at 1.5–2.5× the cost per node — no net gain (DEAD_ENDS B14). The
**facet-only pass is complete**: 65 148 nodes, 1 404 s, 6 leaves (all dim A = 3), ζ's datum among them (G6f);
an interrupted and resumed run gives the identical leaf set. The **full run** (zeros allowed) is ≈ 12 nodes/s at
the top of the tree; expected 10–40 h ⇒ cluster (`scripts/slurm_stageB_zeta.sbatch`, chained `--resume`).
Campaign-wide (10⁴ sets S per stratum) B1 remains far too slow; a structurally better enumeration is needed
(e.g. column sign patterns = covectors of the contact arrangement in ℝ⁴, instead of row-by-row branching).

**B2** (`src/cp4/critical.py`, `scripts/stageB_critical_zeta.py`): `f = det(x0 + Σ t_k N_k)`
(degree 4, `d = dim A ≤ 6` variables); eliminant of `{∂f/∂t_k = 0, D − f = 0}` by sympy grevlex +
FGLM → lex; exact real roots in the window. *Measured:* d = 3 in 21 s. For d = 5, 6 sympy will be too
slow (Bézout 3^d critical points): use **msolve** (`apt`/conda, or build from source — it is a single C
library; `msolve -P 1` gives the eliminating polynomial with exact real-root isolation) or Singular/M2
on the cluster. Positive-dimensional critical sets (non-zero-dimensional ideal) are reported and must be
handled by restricting to the component (they arise e.g. when A contains a direction along which f is
constant).

**Stage C** (only for leaves with a root in the band): locate real critical points (msolve), check the
strict sign conditions of the face datum (cell membership), the frame condition (Thm 6.1: all `|cᵢ|≤1`,
quadruple dets ≤ 1), then exact admissibility of `ΓX` via `cp4.lattice.is_admissible_exact`
(Fincke–Pohst bound, never fixed coefficient boxes). A surviving critical point below `p·D*` would be a
counterexample to Theorem Z — check it three ways before believing it.

## 4. Merged search (fallback if Stage A stays large)

Branch directly on pairs (candidate c, sign vector s_c) with the LP in X-space as the only state; a
candidate is "excluded" ⇔ the constraint `F(c) > 2` is imposed later via its own face datum (i.e. its
sign vector with `Σ s_j(cX)_j > 2`). Hull closure (A3) and quadruple dets (A2) still prune combinatorially.
This removes the Stage A / B1 duplication (every S is re-realised in B1 anyway).

*Prototype (14 Oct 2026, `src/cp4/merged.py`, `scripts/merged_run.py`):* labels `in(s)` / lazy `out` per
candidate; root frame pattern reduced mod (stratum stabiliser × ambient group) — 927 roots (stratum 0),
10 511 (stratum 11); exact frame-sign filter; A2, A3, A5, strict-slack LP; determine-then-check at rank 16;
exclusion test at full leaves; Knuth estimator. ζ witness passes (G7a). *Measured* (Knuth, 200 probes,
`results/merged_{0,11}_estimate.json`): stratum 0 ≈ 5.3·10⁷ nodes (median 3.4·10⁶) against Stage A's
54 726; stratum 11 ≈ 4.2·10⁸ (median 1.8·10⁷), ≈ 4.8·10⁵ leaves; ~240 nodes/s, 85 % of the time in LPs.
**Infeasible as designed.** Reason: branching on the frame signs first replays Stage A's combinatorial
tree once per frame pattern, and partial face data are rarely LP-infeasible at shallow depth (the B1
depth-14 effect), so the LP prune cannot pay for the replication. A revival must defer sign branching
until S is (nearly) fixed and be Knuth-estimated before it is built. The full stratum-0 gate (G7) was not
added: at this cost it would take more than a day.

## 5. Order of work (recommended)

1. `python3 verify/run_all.py` (≤ 2 min). Read the result lines in STATUS §4.
2. Finish G6 (B1 on ζ's S): if it does not complete in 1 h, implement the B1 remedies first.
3. ~~Implement (A6) in Stage A~~ (A6 = §2 option 4(a): kills nothing, DEAD_ENDS B12); strengthen
   Stage A by other means; re-measure stratum 11; then run all strata in increasing candidate-pool
   size (`data/strata.json`, field `candidates`), on the cluster (one SLURM job per stratum,
   partition `ffa`, ≤ 10 GB, 12 h; checkpoint the include/exclude path; heartbeat).
4. B1 + B2 for every S; cache eliminants by A; record per-leaf certificates as JSON
   (`{stratum, S, signs, A_rref, eliminant_factors, window_roots: []}`) — this JSON *is* the proof.
5. Independent checker: a separate script that re-derives A from (S, signs), recomputes det(A) and the
   eliminant (or verifies it via a Groebner certificate), and confirms no root in the window.

## 6. Cost model (to be filled in as measurements arrive)

`total ≈ Σ_strata #S(stratum) × #leaves(S) × cost(B2)`. Known: #S(0) = 3; #S(11) ≥ 6 000 (A1–A5 only);
cost(B2, d=3) ≈ 20 s (sympy). Unknown: #leaves(S). Do not launch the full campaign before
#leaves(ζ's S) and #S(11) with a working Stage-A strengthening are measured ((A6) does not work).
Merged search, frame-first (§4): ≈ 5·10⁷ nodes (stratum 0) and ≈ 4·10⁸ (stratum 11) at ~240 nodes/s.
B1 per S (ζ's S, §3): facet-only 23 min; full ≈ 5.6·10⁵ nodes (Knuth), 10–40 h single-core.
