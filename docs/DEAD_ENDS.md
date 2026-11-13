# DEAD ENDS — do not retry (each cost real time; reasons recorded)

## A. False or unproved mathematics (never use in a proof)

1. **"L₄ (Miller's rational lattice) is optimal."** FALSE. It is the point z = 11/21 of an admissible family and is not even stationary (d′(11/21) = −2/9261). The old documents `PROOF-COMPLETE.md`, `crosspolytope.tex`, `WRITEUP-SKELETON.md`, `dim4-*.tex` all assert this and are **superseded** (kept only in `legacy/` of the outputs folder; do not import anything from them without checking).
2. **"Exclude the Wolff types 2–6."** FALSE target. Types are *configurations*, not disjoint lattice classes: L₄ (and Λ_ζ) contain clean minimal frames of index 1, 2, 3 and 4. Any exclusion of a type must use an extra minimality (least or greatest determinant frame). The month-long `typebnb.py` cluster campaign chased this impossible statement.
3. **"Descent Theorem / Corollary A-i: every admissible lattice with a minimal basis has a {−1,0,1} good basis."** UNPROVED (4 gaps: populated case may return an already-good vector; one swap does not control Φ over all minimal vectors; finiteness of ρ assumed; flat family not covered). The "6,246,480-sector certificate" is not reconstructible. Replaced by the max-det frame (STATUS §2, Thm 6.1).
4. **"Saturation: one minimal pair per parity class."** Not available. Two minimal pairs can share a class (then their ambient supports are disjoint). There are 40 ternary pairs, not 15. The 131 072-configuration enumeration and the catalog `results.csv` assume it.
5. **Parity floor Φ₀ = 32** is basis dependent (GL₄(𝔽₂) does not preserve Hamming weight).
6. **Swap formula** `c′ⱼ = sgn(bⱼ)cⱼ` is wrong; correct: `c′ⱼ = cⱼ`, `c′ᵢ = cᵢ − cⱼb′ᵢ` with `x′ⱼ = sgn(bⱼ)y`.
7. **Chain lemma, dependent case** is false (x_r = −x_s, k = 2).
8. **"Fixing the signs of 4 basis vectors makes all norm inequalities linear."** False: `‖qB‖₁ ≥ 2` is a *disjunction* of halfspaces.
9. **Determinant dichotomy "singular feasible points exist"**: an artifact of an incorrect constraint encoding. Correct finite test (THEORY Thm 1.7) forces nonsingularity.
10. **Lemma "det Λ ≥ Δ₁/p" (superlattice reduction)** is true but weaker by the factor p; it removes nothing.
11. **Bantegnie's paper does not exclude types 2–6**; it confirms max clean index 5 and bounds the index of admissible superlattices of ℤ⁴ by 16 (cyclic) / 18 (non-cyclic) — consistent with our strata (max p = 18).
12. **"The determinant minimum over a polytope is at a vertex."** False (ζ is a cubic irrational). Kills every vertex-LP approach.
13. **"The three near-shell vectors (norm 10/9) of stratum 30 suffice."** FALSE (10 Nov 2026, [Num]). With only them, |det X| drops to ≈ 7.31 (det ≈ 0.812, at the threshold). In that frame the facet rows leak 0.24, 0.48, 0.44 onto X_1's main coordinate. The three norm-11/9 vectors v_s − e_j cut this off and are needed (COMMON_FACET §3b). Likewise, "a single facet-row leak onto t keeps |det X| = 8" holds only for the three-vector system. The summed six-vector constraint alone is far weaker still (|det X| ≈ 2).
14. **"Conjecture 11.A: the 19 constraints of the 7-point system force det ≥ D* on stratum 11."** FALSE (13 Nov 2026, [P]).
    The rational frame X = [[−23/12, 1/12, 0, 0], [1/2, −1/4, −5/4, 0], [0, 2, 0, 0], [−1/3, −1/4, −3/4, −2/3]]
    satisfies them, and even ‖cX‖₁ ≥ 2 for all 68 cube candidates of Γ₁₁, with det = 115/144 < D*. Its short
    vector (½, ¾, −¼, −5/4), of norm 4/3, lies outside the cube (`verify/verify_stratum11.py`). The 10 Nov
    "evidence" came from penalty Nelder–Mead, which never left 4D*. Sequential polyhedral descent finds this point
    in 1–3 of every 25–60 random starts (STRATUM11 §3). Any relaxation proof for stratum 11 must include vectors
    outside [−1,1]⁴.

## B. Computational approaches that failed (measured)

1. **Interval branch-and-bound over shape space** (`bnb.py`, `typebnb.py`): non-convergent. Frontier flat or growing for months; cost ~ w⁻⁹ to resolve the 9-dimensional boundary; targets sat 8×10⁻⁵ above the optimum so det-pruning never bit. p12 finished; p22/p13 never would.
2. **Lazy sign branching over test vectors** (`critdet.py`): at the root only 0.3 % of test vectors have forced signs; ~16^|T| leaves.
3. **Single-contact LP pruning over the Wolff test set** (`contactbf.py`): 98 % of extensions feasible for the first 4–5 levels (18200 → 16038 → 12477 → 6129 → 1622 → … → 0 at depth 11 for ζ's chain). Needs combinatorial pruning first (→ Stage A).
4. **Fixed coefficient boxes for admissibility**: near-singular lattices hide short vectors with large coefficients; produced det = 0 "solutions" three separate times. Always use the rigorous inverse-matrix coefficient bound or Fincke–Pohst (`cp4.lattice`).
5. **Unconstrained or penalty optimisation followed by Newton** for the p = 1 stratum: collapses toward singular matrices, or converges to spurious det = 0 critical points (adjugate vanishes at rank ≤ 2). Only use constrained optimisation inside an exact family with admissibility enforced, and validate on ζ's family first.
6. **SLSQP / differential evolution "infeasible" = exclusion**: never. An optimiser failing to find a point proves nothing.
7. **Uniform slicing** (every lattice has a good dual direction): refuted by a volume estimate (short by ×2–4). Slicing is a *prune*; measured threshold det ≈ 1.0 with trivial Δ₃ (≈ 0.94 with sharp Δ₃) — it leaves exactly the dangerous band.
8. **Primal (projection) slicing**: ratio 0.855 at Λ_ζ, can never be tight there; do not invest.
9. **Cohn–Elkies LP bound**: not sharp (≈0.984 vs 18/19 in dim 3; ≈0.913 vs 0.8249 in dim 4). Worth a separate paper as an upper bound for all packings, not a route to Theorem Z. Grid LP is unstable beyond degree 5; needs SOS/SDP.
10. **Treating cluster silence as a stall**: scripts that print only at the end look dead. Every long job must emit a heartbeat (all current scripts do).
11. **Checkpoint races**: never run two writers on one checkpoint; gzip patch must import gzip (the old header lacked the substring `import pickle`).
12. **Frame-facet LP as a Stage-A prune ("A6", PLAN §2 option 4(a))**: fix a closed facet `uᵢ ∈ {±1}⁴` of each frame row (330 patterns mod the ambient group) and require `uᵢ·Xᵢ = 2`, `‖cX‖₁ ≤ 2` on S. Sound and hereditary, but it kills **0 of 200** Stage-A sets (stratum 0) and **0 of 931** (stratum 11, 30 000 nodes); sampled sets keep 187–270 of 330 patterns feasible, ζ's S 295. Reason: `‖cX‖₁ ≤ 2` is only an upper bound, met by shrinking `cX`; the binding conditions (`F(c) = 2` on S, `F ≥ 2` off S) are non-convex (A8). `scripts/probe_A6_frame_facets.py`.
13. **Merged search with frame signs branched first (PLAN §4 as written)**: `src/cp4/merged.py` is sound (ζ witness G7a passes) but the Knuth estimate (200 probes) is ≈ 5·10⁷ nodes for stratum 0 (median 3.4·10⁶; Stage A needs 54 726) and ≈ 4·10⁸ for stratum 11 (median 1.8·10⁷), at ~240 nodes/s with 85 % of the time in LPs. Reason: the root splits into 927 / 10 511 frame patterns, each replaying most of Stage A's tree, while partial face data are rarely LP-infeasible before ~13 contacts. Do not retry without deferring sign branching until S is (nearly) fixed, and Knuth-estimate first (`scripts/merged_run.py <sid> --estimate K`).
14. **Norm-range lookahead in B1** (prune a node when, for some undecided contact c, `min ‖cX‖₁ > 2` or `max ‖cX‖₁ < 2` over the cell closure): sound, but on ζ's S it cuts the Knuth-estimated tree only 2–3× (to 2.9–4.0·10⁵ nodes) while raising the cost per node 1.5–2.5× (130–240 LPs/node) — no net gain; not implemented (PLAN §3).
15. **Column-wise ("covector") B1 enumeration**: choose each column of the face-datum sign matrix as a covector of the arrangement `{c^⊥ : c ∈ S}` in ℝ⁴ (ζ's S: N = 852 classes mod ±, T = 186 topes; ζ's columns and all 6 facet-only leaves' columns are among them), as a multiset mod the ambient group, LP at the leaf. Sound, but **no convex prune exists before the 4th column**: the partial LP (strict cells, `Σ_{j≤k} |c·Yⱼ| ≤ 2`) was feasible for 1 000/1 000 random 1-, 2- and 3-tuples, because a partial assignment scales towards 0. So every one of C(N+3,4) = **2.2·10¹⁰** multisets (5.2·10⁷ for topes) reaches the leaf LP (0 of 5 000 random leaves feasible; ~5 ms/LP ⇒ ~3 years, facet-only ~3 days), against row-wise B1's ~5.6·10⁵ nodes (facet-only 65 148, done in 23 min). `scripts/estimate_covector_B1.py`, `results/covector_B1_estimate.json`. A revival needs a non-convex partial prune. (The covector list is useful the other way round: as an LP-free pre-filter inside row-wise B1 it saves 4.5–10× — PLAN §3, not a dead end.)

16. **Facet-class layer prune inside B1** (discard a node when the decided contacts whose sign vectors put cX on
    one closed facet ±σ generate Γ, so σ·Λ ⊆ 2ℤ and det ≥ 5/6 by FACET_LAYERS §1): sound, exact and LP-free, but
    nearly inert on ζ's S. Facet-only complete pass 65 148 → 64 094 nodes (−1.6 %, same 6 leaves). Full static
    tree, 900 s each way: the same 25 depth-3 subtrees finished with 90 708 vs 90 819 nodes, although 2 322
    children were cut. None of the 46 + 6 known leaves is killed (ζ's leaves are spread over at most 3
    contacts per facet class). Reason: B1's dead subtrees die by LP infeasibility around depth 13–14,
    before any facet class contains ≥ 4 contacts generating the index-4 group Γ. Not worth building in.
    `scripts/measure_layer_prune.py`, `results/layer_prune_zeta.json`.

17. **Realization-wide "pair lemma"** (any admissible Λ containing minimal u, v with disjoint supports and
    (u+v)/2 ∈ Λ has det ≥ 81/100; it would exclude strata 1, 8, 9, 10, 15 at once): numerically almost false.
    With the other two frame rows of norm 2, constrained searches reach det 0.853 (supports 2+2), 0.867 (1+3),
    0.922 (1+2) [Num], about 5 % above the target. The coordinate-section bound gives only det ≥ ¼‖u‖∞‖v‖∞
    here. Per-stratum residue structure, or the max-det condition, is needed instead.

18. **Supporting-facet B1 with closed cells or a fixed tie side** (13 Nov 2026, scratchpad Knuth on stratum-11 sets,
    60–400 probes). Theorem B.1′ (STATUS §3) lets non-frame contacts take facets only. But with **closed** cells
    (σ_j(cX)_j ≥ 0 everywhere) every forced zero coordinate is explored on both sides, with duplicate
    boundaries. Trees: ζ's S 1.3·10⁶ (8.6·10⁶ with ball rows); set #275 6.2·10⁸ (max 3.7·10¹⁰), against 2.5·10⁴
    for the strict (incomplete) facet pass. With the **fixed** tie side ρ_j = sign of the first nonzero frame
    entry of column j, cells pinned to (cX)_j = 0 survive whenever ρ_j is the infeasible side: 3.2·10⁵ (ζ's S),
    1.6·10⁶ (#1165), 3.3·10⁶ (#129). Use the half-open cells with a per-node tie side (Lemma B.4): with ball rows,
    5.4·10⁴ (ζ's S), 3.4·10⁵ (#1165), 2.1·10⁵ (#275), against 1.2·10⁶, 8.3·10⁶, 3.2·10⁶ for the original B1.

## C. Process rules learned the hard way

* Validate every new tool on a **known answer** before any long run (p = 1 ⇒ A₄; stratum 11 ⇒ contains ζ's S; Stage B on ζ's S ⇒ recovers the cubic and D*).
* Exact arithmetic for every accept/reject decision; floats only to search or to prune *conservatively*.
* Separate [P] from [Num] in every write-up; an external audit already caught a false main theorem once.
* Do not quote δ(L₄) or the density where a determinant threshold is meant.
* **Long runs in the Claude Code cloud container are not reliable**: the container is reclaimed a few minutes after the
  session goes idle (16 Oct 2026: G6 legs of 2 h 10 min, 2 h 9 min, then 6.5 min and 5 min; each later leg lost its
  unfinished depth-3 subtree). Run anything longer than a working session on the cluster from a committed checkpoint.
