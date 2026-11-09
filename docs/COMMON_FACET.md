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

### 3a. Stratum 30 with rows 2–4 on the axes: det ≥ 8/9, sharp [P]

Suppose X_2, X_3, X_4 are axis rows 2e_{k_2}, 2e_{k_3}, 2e_{k_4} (distinct coordinates, positive by §1); let t be
the remaining coordinate, x = X_1, A_r = (x_{k_r})₊ and B_r = (x_{k_r})₋. For the shift (i,j,l),
vX = ⅓x + (2/9)e_{k_i} + (8/9)e_{k_j} − (4/9)e_{k_l}. Its total is 20/9, and ‖vX‖₁ ≥ 2 caps the cancellation:
min(B_i, 2/3) + min(B_j, 8/3) + min(A_l, 4/3) ≤ ⅓. Each cap exceeds ⅓, so this reads

    B_2 + B_3 + A_4 ≤ ⅓,   B_4 + B_2 + A_3 ≤ ⅓,   B_3 + B_4 + A_2 ≤ ⅓.

Summing: (A_2+A_3+A_4) + 2(B_2+B_3+B_4) ≤ 1. So X_1 carries mass ≤ 1 on the row coordinates, and |x_t| ≥ 1.
Since |det X| = 8|x_t|, **det Λ ≥ 8/9 > D***. Equality forces B_r = 0, A_r = ⅓ and |x_t| = 1, i.e.
X_1 = (⅓ on the row coordinates, ±1 at t). That lattice *is* admissible: with X_2 = 2e_3, X_3 = 2e_1,
X_4 = 2e_4, X_1 = (⅓, ∓1, ⅓, ⅓) the exact short-vector test gives λ₁ = 2 with 16 minimal pairs, and
det = 8/9 [P+C]. So in this case 8/9 is the exact minimum. The relaxed search (14 constraint vectors,
random starts) found nothing lower anywhere (best 0.912, at a frame with rows 2–4 within 0.012 of the axes)
[Num]. **Conjecture: every admissible realization of stratum 30 has det ≥ 8/9.** What remains is the
non-axis case: show that leakage of rows 2–4 off their main coordinates does not lower the bound. That is the
perturbation step of §6.

## 4. Consequences for stratum 17 [P]

The near-shell vectors are ½X_1 ± (⅓X_2 − ⅙X_3 − ⅙X_4), of frame norm 7/6 and budget ⅓ each. They give

    Σ_k min(3|X_1k|, R_k) ≤ 2,      Σ_k min(2X_2k, X_3k + X_4k) ≤ 1.                       (4)

The vectors ½(X_1 ± X_3 ∓ X_4) (frame norm 3/2, budget ½) add the weaker Σ_k min(X_3k, X_4k) ≤ 1.
(Layer special case: if σ·X_1 = 0 then σ·Λ = 2ℤ and det ≥ 5/6.)
*Axis case.* The same computation with rows 2–4 on the axes gives B_{k_2} + A_{k_3} + A_{k_4} ≤ ⅓ and
A_{k_2} + B_{k_3} + B_{k_4} ≤ ⅓, so |x_t| ≥ 4/3 and det Λ = 8|x_t|/6 ≥ 16/9. But the searches reach 0.968 [Num], so
stratum 17's low configurations are **not** near-axis. A relaxed random-start search (16 vectors, ‖c‖₁ ≤ 3/2)
found det 0.958 at a frame where one facet row is split, X_4 ≈ (1.12, 0.68, 0.20, 0), overlapping the main
coordinate of X_1 ≈ (0, 1.79, 0, 0.20); the other rows are near-axis [Num]. Stratum 17 needs its own case
analysis, organized by which facet row splits and where.

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
