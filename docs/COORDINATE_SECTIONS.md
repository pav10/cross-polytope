# Coordinate sections and fiber offsets: strata 6, 7 and 25 excluded

9 November 2026. The lemmas are theoretical [P]; the finite facts they use about the residue groups, and one
inequality about a single piecewise-linear function (§3), are checked exactly by
`verify/verify_coordinate_sections.py` [P+C]. Conventions as in `THEORY_REVIEW.md`: Γ ⊇ ℤ⁴ of index p, X
nonsingular with frame rows X_i of ℓ¹ norm 2, Λ = ΓX admissible (nonzero vectors have ℓ¹ norm ≥ 2).
β_d is a valid lower bound for the determinant of an admissible lattice in (ℝ^d, ℓ¹): β_1 = β_2 = 2
(Minkowski). Nothing here uses local optimality or max-det properties of the frame: every admissible
realization of these strata has determinant ≥ 81/100 > D*. Together with THEORY_REVIEW, BLOCK_COUPLING and
FACET_LAYERS this leaves **14 strata: 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17, 30.**

## 1. The coordinate-section lemma [P]

**Lemma 1.** Let T ⊆ {1,2,3,4}, |T| = r, and suppose Λ_T := Λ ∩ ℝ^T has rank r. Let μ_T be the ℓ¹ covering
radius of Λ_T in ℝ^T and π the coordinate projection onto ℝ^{T^c}. Then π(Λ) is a lattice of rank 4−r, every
nonzero vector of π(Λ) has ℓ¹ norm ≥ 2−μ_T, and det Λ = det Λ_T · det π(Λ). Hence, if μ_T < 2,

    det Λ ≥ det Λ_T · β_{4−r} · ((2 − μ_T)/2)^{4−r}.                                   (1)

*Proof.* Λ_T is the kernel of π on Λ and has full rank in ker π = ℝ^T, so π(Λ) is a lattice and the
determinant factorises (π is an orthogonal projection). If v ∈ Λ has π(v) ≠ 0, choose w ∈ Λ_T with
‖v_T − w‖₁ ≤ μ_T. Then v − w ∈ Λ∖0, so 2 ≤ ‖v − w‖₁ = ‖v_T − w‖₁ + ‖π(v)‖₁ ≤ μ_T + ‖π(v)‖₁. Scaling π(Λ) by
2/(2−μ_T) gives an admissible lattice of dimension 4−r. ∎

The block bound of THEORY_REVIEW §2 is the case where Λ splits; (1) does not need Λ to split.

**Axis rows.** Call X_i an *axis row* if X_i = ±2e_k. Let R be a set of axis rows occupying the coordinates T
(|T| = |R|), and Γ_R = Γ ∩ ℝ^R (vectors of Γ supported on the rows R). Then **Λ_T = Γ_R X_R exactly**: if cX
is supported in T, the rows outside R restricted to T^c form a nonsingular square block (X is nonsingular and
the rows of R vanish on T^c), so c vanishes outside R. Up to a signed permutation of coordinates (an ℓ¹
isometry), Λ_T = 2Γ_R.

**Corollary 2 (an axis pair).** If rows i, j are axis rows and Γ_{ij} = ℤ² + ℤ(½,½), then **det Λ ≥ 1**.
Indeed Λ_T ≅ 2Γ_{ij} = D₂ = {(u,v) ∈ ℤ² : u+v even}, with det 2. Its ℓ¹ covering radius is 1: in the
coordinates u+v, u−v the lattice D₂ is 2ℤ², and ‖(u,v)‖₁ = max(|u+v|,|u−v|). By (1),
det Λ ≥ 2 · 2 · (1/2)² = 1.

**Forcing axis rows.** If three rows are pairwise support-disjoint (THEORY_REVIEW §2, from the forced unit
shell), their supports are three disjoint nonempty subsets of a 4-element set, of sizes (1,1,1) or (1,1,2) in
some order; so at least two of them are axis rows.

## 2. Stratum 6 [P]

p = 4, Γ = ℤ⁴ + ⟨(0,0,½,½), (0,½,0,½)⟩ (residues (0,0,2,2), (0,2,0,2), (0,2,2,0) over 4). The unit shell
contains (0,0,½,±½), (0,½,±½,0), (0,½,0,±½), so rows 2, 3, 4 are pairwise support-disjoint and two of them,
say i and j, are axis rows. For **every** pair {i,j} ⊆ {2,3,4}, Γ_{ij} = ℤ² + ℤ(½,½). Corollary 2: det Λ ≥ 1.

## 3. Stratum 7 [P+C]

p = 4, Γ = ℤ⁴ + ℤ(0,½,¼,¼); every residue has first coordinate 0, so Γ = ℤe_1 ⊕ Γ_{234}. The unit shell
contains (0,0,½,±½) and (0,½,±¼,±¼) (equal signs in the last two places), so rows 2, 3, 4 are pairwise
support-disjoint. Γ_{34} = ℤ² + ℤ(½,½), while Γ_{23} = Γ_{24} = ℤ².

*Patterns with X_3, X_4 both axis rows* (this includes (1,1,1)): det Λ ≥ 1 by Corollary 2.

*Residual pattern.* Otherwise X_2 is an axis row and exactly one of X_3, X_4 is; Γ is invariant under swapping
coordinates 3 and 4, so take X_3 axis and X_4 with a two-element support. After a signed permutation of the
ambient coordinates (a, b, c, d):

    X_2 = 2e_a, X_3 = 2e_b, X_4 = (0, 0, y_c, y_d) with y_c, y_d > 0, y_c + y_d = 2, X_1 = (α, β, x_c, x_d).

Then |det X| = 4D with D := |x_c y_d − x_d y_c|, so **det Λ = D**. Put e := max(y_c, y_d) − 1 ∈ [0,1).

*Fiber offsets.* The vectors of Λ with frame coefficient c_1 = 1 are X_1 + 2c_2e_a + 2c_3e_b + c_4X_4,
(c_2,c_3,c_4) ∈ Γ_{234} = ℤ³ + ℤ(½,¼,¼). Their (c,d)-part is x + c_4y, with c_4 = k/4, k ∈ ℤ; for fixed k the
pair 2(c_2,c_3) runs over a coset o_k + 2ℤ², with o_0 = (0,0), o_1 = (1,½), o_2 = (0,1), o_3 = (1,3/2)
(k mod 4). Writing δ₂(t) for the distance from t to 2ℤ and f_k = δ₂(α + o_{k,1}) + δ₂(β + o_{k,2}),
admissibility gives

    f_{k mod 4} + g(k) ≥ 2 for every k ∈ ℤ, where g(k) := ‖x + (k/4) y‖₁.                 (2)

Since δ₂(t) + δ₂(t+1) = 1 for every real t, **f_0 + f_1 + f_2 + f_3 = 4** identically, so (2) implies

    g(k) + g(k+1) + g(k+2) + g(k+3) ≥ 4 for every k ∈ ℤ.                                  (3)

*Closed form of g.* With |a|+|b| = max(|a+b|, |a−b|), S = x_c + x_d and Δ = x_c − x_d,
g(k) = max(|S + k/2|, |Δ + (y_c − y_d)k/4|). Substituting s = k + 2S, and then reflecting s ↦ −s and changing
the sign inside the second absolute value as needed (both preserve the set of 4-windows), gives

    g = max(|s|/2, |D + e s/2|) on the points s ∈ δ + ℤ, for some δ ∈ [0,1),                (4)

because Δ − σeS = x_c y_d − x_d y_c for σ = sign(y_c − y_d).

**Lemma 3 [P+C].** If D ∈ [0, 81/100], e ∈ [0,1] and δ ∈ [0,1), then some window
{δ+j, δ+j+1, δ+j+2, δ+j+3}, j ∈ {−4,…,1}, has sum of (4) below 4. *Certificate:* on the rational grid of step
1/40 in (D, e, δ) the smallest window sum is at most 5903/1600 = 3.69 < 4 − 1/5. Each window sum is
8-Lipschitz in the sup-distance to the nearest grid point: |∂g/∂D| ≤ 1, |∂g/∂e| ≤ |s|/2 ≤ 5/2 and
|∂g/∂δ| ≤ 1/2 on these windows, summed over four terms with half a step in each parameter.
Without the grid, the bound found numerically is about D ≥ 0.959 [Num]; it is about 1 at e = 0.

By (3) and Lemma 3, D ≥ 81/100 in the residual pattern. **Every admissible realization of stratum 7 has
det Λ ≥ 81/100 > D*.** [Num] A direct search of the residual pattern (`s7probe` in the session log; 400
Nelder–Mead runs with full admissibility) found no determinant below 1.0036; the infimum is approached as X_4
becomes an axis row.

## 4. Stratum 25 [P+C]

p = 8. Its unit shell is that of stratum 7, (0,0,½,±½) and (0,½,±¼,±¼), so rows 2, 3, 4 are pairwise
support-disjoint, and Γ_{34} = ℤ² + ℤ(½,½): patterns with X_3, X_4 axis give det Λ ≥ 1. Γ_25 contains the
group of stratum 7 with index 2; every residue has first coordinate 0 or 4 (over 8), and Γ is again invariant
under swapping coordinates 3 and 4.
In the residual pattern (same normalisation as §3) now **det Λ = 4D/8 = D/2**. Use the coset line c_1 = ½. The
vectors ½X_1 + 2c_2e_a + 2c_3e_b + c_4X_4 have c_4 ∈ ⅛ + ¼ℤ, and the four classes
c_4 ≡ ⅛, ⅜, ⅝, ⅞ have offsets (½,5/4), (3/2,7/4), (½,¼), (3/2,¾). Their distance sum to 2ℤ² is again
identically 4. The (c,d)-parts ½x + c_4y lie on a line at ℓ¹ distance (D/2)/(1+e) from 0. The same
computation as for (4) gives g = max(|s|/2, |D/2 + e s/2|) on s ∈ δ' + ℤ. Lemma 3 with D replaced by D/2
gives det Λ = D/2 ≥ 81/100.

## 5. Partial restrictions for other strata [P]

Corollary 2 also excludes, in strata whose rows 3, 4 are support-disjoint and Γ_{34} ∋ (½,½) — strata 1, 8,
9, 10, 15 — every realization in which **X_3 and X_4 are both axis rows**; in stratum 8 also every realization
in which X_1 and X_2 are both axis rows. These are not exclusions of whole strata, since the rows need not be
axis rows. A search of these strata may assume them.

## 6. Scope and next steps

The fiber-offset argument of §3 needs two axis rows whose coordinate section is 2ℤ², so that the fibers over
the remaining coordinates are cosets of 2ℤ² with a computable offset sum. Strata 1, 8, 9, 10 and 15 force only
pairwise disjointness of one or two row pairs. Strata 4, 17 and 30 force one positive unit-shell vector on rows
2–4. Strata 0, 2, 3, 5, 11 and 14 have an empty unit shell beyond the frame. Extending the fiber-offset method
to one axis row, where the fibers are cosets of a one-dimensional lattice over a three-dimensional projection,
is the natural next target.
