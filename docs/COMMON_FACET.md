# Strata 17 and 30: structure, necessary conditions, and a route to exclusion

9 November 2026. Work in progress. §§1–4 are proved [P]; §5 is numerical [Num]; §6 is the plan. No stratum is
excluded here yet. Conventions as in `THEORY_REVIEW.md` (Γ ⊇ ℤ⁴ of index p, frame rows X_i of ℓ¹ norm 2,
Λ = ΓX admissible).

## 1. Common facet and the coset structure [P]

Both strata contain the unit-shell vector (0,⅓,⅓,⅓). By equality in the triangle inequality, rows 2, 3 and 4
lie on one closed facet σ (FACET_LAYERS §2 argument). After coordinate sign changes σ = (1,1,1,1), so
**X_2, X_3, X_4 ≥ 0 entrywise, each with coordinate sum 2.** Row X_1 is arbitrary with ‖X_1‖₁ = 2. Put
c = (X_2+X_3+X_4)/3 (a contact on the facet), a_i = X_i − c, and

    Λ'' = Λ ∩ V,  V = span(X_2, X_3, X_4),  Λ'' = ℤc ⊕ Λ''_0,  Λ''_0 = ℤa_2 + ℤa_3  (σ·Λ''_0 = 0),

because Γ ∩ {c_1 = 0} = ℤ³ + ℤ(⅓,⅓,⅓) on rows 2–4 in both strata. Then Λ = Λ'' + ℤh with

* stratum 17 (p = 6, Γ = ℤ⁴ + ℤ(0,⅓,⅓,⅓) + ℤ(½,0,½,½)): h ≡ ½X_1 + ½a_2 (mod Λ''; h = (½,0,½,½)X), and X_1 ≡ 2h;
* stratum 30 (p = 9, Γ = ℤ⁴ + ℤ(⅓,1/9,4/9,7/9)): h = (⅓,1/9,4/9,7/9)X = ⅓X_1 + ⅓(a_4 − a_2) + (4/3)c, and
  X_1 ≡ 3h (mod Λ'').

If N ∈ ℝ⁴ is the vector of signed 3×3 minors of the rows X_2, X_3, X_4 (a normal of V), then

    det Λ = det_V(Λ'') · dist₂(h, V) = |N·h|/3 = d_V(h)·‖N‖_∞/3,                      (1)

where d_V is the ℓ¹ distance to V. This uses det_V(Λ'') = ‖N‖₂/3 and the fact that the ℓ¹ distance to the
hyperplane N^⊥ is |N·x|/‖N‖_∞.

## 2. The cancellation-budget lemma [P]

For c ∈ Γ put C(c) = Σ_i |c_i|·‖X_i‖₁ − ‖cX‖₁ ≥ 0 (the cancellation of the triangle inequality). Coordinatewise,
C(c) = 2 Σ_k min(P_k, N_k), where P_k and N_k are the totals of the positive and negative terms c_i X_ik.
**Admissibility gives C(c) ≤ 2‖c‖₁ − 2.** For ‖c‖₁ = 1 this is the forced-shell lemma of THEORY_REVIEW §2
(no cancellation at all). For a *near-shell* vector, ‖c‖₁ = 1 + ε, the total cancellation is at most 2ε. The
inequality min(a+b, c+d) ≥ min(a,d) + min(b,c) (a, b, c, d ≥ 0) splits it into pairwise overlap budgets.

## 3. Consequences for stratum 30 [P]

Γ₃₀ has three near-shell cosets (up to sign) with unique minimal representatives of ℓ¹ norm 10/9:

    v = ⅓X_1 + (1/9)X_i + (4/9)X_j − (2/9)X_l,  (i,j,l) ∈ {(2,3,4), (4,2,3), (3,4,2)}.

Each has cancellation budget 2/9. With R_k = X_2k + X_3k + X_4k (so Σ_k R_k = 6), summing the three budgets and
splitting X_1 = A − B into positive and negative parts gives

    Σ_k min(3|X_1k|, 2R_k) ≤ 3,                                                            (2)
    Σ_k min(X_ik + 4X_jk, 2X_lk) ≤ 1   for each cyclic (i,j,l).                            (3)

In particular (3) gives pairwise overlaps O_ab = Σ_k min(X_ak, X_bk) ≤ ½ for all pairs of rows among 2–4
(min(4x, 2y) ≥ 2 min(x, y)). Since a+b+c ≤ max + (sum of the pairwise minima) for nonnegative a, b, c, rows 2–4
then put total mass ≥ 6 − 3/2 = 4.5 on their dominant coordinates (Σ_k max_i X_ik ≥ 4.5). That alone allows
leakage up to ½ per row, which is too much for a diagonal-dominance determinant bound (1.5⁴ ≈ 5 < 9·0.81). So
the X_1 coupling (2) and the slicing family (5) must carry the rest.
So X_1 must largely avoid the coordinates carrying rows 2–4: each cyclic shift demands the opposite sign of X_1
on one row's support, so the three demands conflict. Rows 2–4 must also be nearly disjoint in a weighted sense.
(Layer special case: if σ·X_1 = −2 then σ·Λ = 2ℤ and det ≥ 5/6 by FACET_LAYERS.)

## 4. Consequences for stratum 17 [P]

The near-shell vectors are ½X_1 ± (⅓X_2 − ⅙X_3 − ⅙X_4), of frame norm 7/6 and budget ⅓ each. They give

    Σ_k min(3|X_1k|, R_k) ≤ 2,      Σ_k min(2X_2k, X_3k + X_4k) ≤ 1.                       (4)

The vectors ½(X_1 ± X_3 ∓ X_4) (frame norm 3/2, budget ½) add the weaker Σ_k min(X_3k, X_4k) ≤ 1.
(Layer special case: if σ·X_1 = 0 then σ·Λ = 2ℤ and det ≥ 5/6.)

## 5. Calibration and the slicing family [Num] / [P]

*Slicing family [P].* For every m ∈ Γ* (m·Γ ⊆ ℤ), w = X⁻¹m lies in Λ*, so by THEORY 8.1

    det Λ ≥ vol₃(C_4 ∩ w^⊥) / ‖w‖₂.                                                        (5)

FACET_LAYERS is the case w = σ/2. At the axis frame X = 2I of stratum 30, m = (−1,−1,1,0) ∈ Γ₃₀* gives
w = m/2 and exactly 1.0 (s(1110) = √3/2, ‖w‖₂ = √3/2).

*[Num].* Constrained searches (`scripts/probe_common_facet.py`; basin hopping from the best points) give det
≈ 0.968 for stratum 17 (best point with 7 minimal pairs, not yet a local minimum) and ≈ 1.03 for stratum 30
(7 pairs). At the best stratum-30 frame, which is near-axis (diagonal 1.68–1.84, leakage ≤ 0.32 per row), (5) with
m = (−1,−1,1,0) gives ≈ 0.851 > 81/100 (Monte Carlo section volume).

## 6. Plan to finish the two exclusions

1. From (2)–(4), prove that every admissible frame of the stratum is *near-axis*: after a permutation of
   coordinates, X_ii ≥ 2 − ℓ with an explicit leakage bound ℓ. Rows 2–4 are nearly disjoint by (3)/(4), and
   then X_1 is pushed onto the remaining coordinate by (2)/(4). Handle the residual cases where rows 2–4 share
   coordinates by the same budgets, ending either in a contradiction or in a determinant bound.
2. In the near-axis region, bound ‖X⁻¹m‖₂ for a fixed short m ∈ Γ* (for stratum 30, m = (−1,−1,1,0) and its
   images under the cyclic symmetry of rows 2–4), and the section volume from below, so that (5) gives ≥ 81/100.
   The margin at the axis is 1.0, and at the best-known point 0.851.
3. Same for stratum 17 with its short dual vectors.
Step 1 is pure inequality work. Step 2 needs explicit perturbation bounds for X⁻¹ and for section volumes; THEORY 8.1
has the exact (Hermite–Genocchi) section formula.
