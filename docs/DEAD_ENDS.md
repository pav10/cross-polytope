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

## C. Process rules learned the hard way

* Validate every new tool on a **known answer** before any long run (p = 1 ⇒ A₄; stratum 11 ⇒ contains ζ's S; Stage B on ζ's S ⇒ recovers the cubic and D*).
* Exact arithmetic for every accept/reject decision; floats only to search or to prune *conservatively*.
* Separate [P] from [Num] in every write-up; an external audit already caught a false main theorem once.
* Do not quote δ(L₄) or the density where a determinant threshold is meant.
