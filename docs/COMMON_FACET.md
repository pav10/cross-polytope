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
[Num]. **Conjecture: every admissible realization of stratum 30 has det ≥ 8/9.** §3b proves it in the whole
near-axis region (the "main configuration"). What remains is the configuration lemma of §3b.

### 3b. Six short vectors; the main configuration has det ≥ 8/9 [P]

*The six vectors.* For each cyclic shift s = (i,j,l), Γ₃₀ contains v_s = ⅓e_1 + (1/9)e_i + (4/9)e_j − (2/9)e_l
(norm 10/9) and v'_s = v_s − e_j (norm 11/9). With x = X_1, P_r = X_r (r = 2, 3, 4) and
u_s = P_i + 4P_j − 2P_l, u'_s = P_i − 5P_j − 2P_l, admissibility gives

    (N_s)  ‖3x + u_s‖₁ ≥ 18,        (M_s)  ‖3x + u'_s‖₁ ≥ 18.                                   (6)

[Num] The three N_s alone are **not** enough. With only them, the minimum of |det X| drops to ≈ 7.31 (det ≈ 0.812,
at the threshold). The minimizer is a frame whose three facet rows leak 0.24, 0.48 and 0.44 onto X_1's main
coordinate, and it violates all three M_s. With all six, take the exact inner minimum over X_1 (a MILP for each
frame) and run Nelder–Mead over rows 2–4 from near-axis starts. Four runs end at 8.006–8.018, next to the axis.
The minimum also grows linearly along leakage directions: 9.44 for a 0.2 leak onto t, 9.12 for a 0.1 cross
leak, and over 6·10⁴ random leak patterns min (|det X| − 8)/L ≈ 1.6.

*Main configuration.* Label the coordinates k_2, k_3, k_4, t (distinct). Put d_r = P_r(k_r),
e_{r→r'} = P_r(k_{r'}), f_r = P_r(t), ℓ_r = 2 − d_r, L = Σℓ_r = E + F (E the cross leaks, F the leaks onto
t), y_r = x(k_r). Let M be the 3×3 block of rows 2–4 on columns k_2, k_3, k_4, and
C_r = det(M with column r replaced by (f_2, f_3, f_4)ᵀ). The frame is in the *main configuration with sign τ*
if, for every s, the sign vector of 3x + u_s on (k_i, k_j, k_l, t) is (+, +, −, τ) and that of 3x + u'_s is
(+, −, −, τ). Zero entries count as either sign.

If ε is the sign vector of z, then ‖z‖ = ε·z, and ε·x ≤ 2 − 2·(the mass of x where ε_k x_k < 0). So in the
main configuration with τ = + the six constraints (6) are equivalent to linear inequalities:

    q_l := 2d_l − e_{i→l} − 4e_{j→l} − 3   ≥ 3(y_l⁺ + y_i⁻ + y_j⁻ + x_t⁻),
    m_s := 5S_j + 2S_l − S_i − 12          ≥ 3(y_j⁺ + y_l⁺ + y_i⁻ + x_t⁻),   S_r = P_r(k_j) + P_r(k_l).

**Budget identity.** Σ_l q_l + Σ_s m_s = 9(1 − L) identically, for any rows of sum 2. Every unit of leakage
costs exactly 9: each cross leak has N-weight 1 or 4 and M-weight 6 or 3 (sum 7), plus 2 from its own row.
Summing the six inequalities gives 9Σ|y_r| + 18x_t⁻ ≤ 9(1 − L). With ‖x‖₁ = 2 this forces

    x_t ≥ 1 + L   and   L ≤ 1.                                                                 (7)

The axis case of §3a is L = 0.

**Determinant.** Schur's complement and Cramer give det X = −(x_t det M − Σ_r y_r C_r). Hence
|det X| ≥ x_t det M − (2 − x_t) max_r |C_r|, which increases in x_t, so by (7)

    |det X| ≥ (1 + L) det M − (1 − L) max_r |C_r|.

**Lemma.** If P ≥ 0, the rows have sum 2 and L ≤ 1, then (1 + L) det M − (1 − L)|C_r| ≥ 8 for each r.
*Proof.* Write the other two rows as b, c, with u = ℓ_b + ℓ_c and E_r = e_{r→b} + e_{r→c}.
(a) C_r ≥ 0. *Column trick:* if x ≥ 0 and Σx = 2, then det X = −2 det(M − 𝟙yᵀ). To see this, add all
columns to column t (it becomes (2, 2, 2, 2)ᵀ) and subtract row 1 from the others. Take
x = (1 − L)e_{k_r} + (1 + L)e_t. Then (1 + L) det M − (1 − L)C_r = 2 det M′ with M′ = M − (1 − L)𝟙e_rᵀ, whose
column r is (1 + u, −α, −β) with α = 1 − L − e_{b→r} and β = 1 − L − e_{c→r}. Expanding along column r,

    det M′ = (1 + u)(d_b d_c − e_{b→c}e_{c→b}) + e_{r→b}(α d_c − β e_{b→c}) + e_{r→c}(β d_b − α e_{c→b}).

Since e_{b→c} ≤ ℓ_b and e_{c→b} ≤ ℓ_c, the first term is ≥ (1 + u)(4 − 2u). Using 0 ≤ d_c ≤ 2 and β ≤ 1 − L,
αd_c − βe_{b→c} ≥ −2(e_{b→r} − (1 − L))⁺ − (1 − L)e_{b→c} ≥ −2ℓ_b; symmetrically, the other bracket is
≥ −2ℓ_c. So det M′ ≥ 4 + 2u(1 − u − E_r) ≥ 4, because E_r + u ≤ L ≤ 1.
(b) C_r < 0. Let D_r = det M^{(r)}, where M^{(r)} is M with the t-column added to column r; its rows sum to 2.
Since det M = D_r − C_r, the expression equals (1 + L)D_r − 2L C_r ≥ (1 + L)D_r. Write M^{(r)} = 2I − ℒ with ℒ a
Laplacian (zero row sums, off-diagonals ≤ 0). Then D_r = 8 − 4 tr ℒ + 2σ₂(ℒ) ≥ 8 − 4 tr ℒ, since the
principal 2×2 minors of ℒ are ≥ 0 and det ℒ = 0. Also tr ℒ = E_r + u ≤ L, so
(1 + L)D_r ≥ (1 + L)(8 − 4L) ≥ 8. ∎

So **in the main configuration with τ = +, |det X| ≥ 8, i.e. det Λ ≥ 8/9.** It is sharp at the axis
(§3a). The lemma needs only L ≤ 1, and the six constraints enter only through their sum (7).

*τ = −.* The linear forms become q⁻_l = −u_s(k_l) − u_s(t) − 3 ≥ 3(y_l⁺ + y_i⁻ + y_j⁻ + x_t⁺) and
m⁻_s = u′_s(k_i) = d_i − 5e_{j→i} − 2e_{l→i} ≥ 3(y_j⁺ + y_l⁺ + y_i⁻ + x_t⁺). They satisfy, identically,

    Σq⁻ + Σm⁻ = 9 − 9E − 6F,   2Σq⁻ − Σm⁻ = −9E₃ − 9F,   2Σm⁻ − Σq⁻ = 9 − 9E₆ + 3F,

where E₃ (E₆) is the cross leakage of N-weight 4 (1). The plain sum forces x_t < 0. Now |det X| =
|x_t| det M + Σ y_r C_r ≥ 2 det M − (det M + c⁻)ΣA − (det M + c⁺)ΣB, where A, B are the positive and negative
parts of y and c^± = max_r (C_r)^±. Weight the N's by α = (det M + 2c⁺ − c⁻)/9 and the M's by
β = (det M + 2c⁻ − c⁺)/9; then

    |det X| ≥ det M (1 + E + ⅔F) + c⁺(E₃ + F) − c⁻(1 − E₆ + ⅓F).

Crude bounds finish it:
* det M ≥ 8 − 4L (the same Laplacian expansion with diagonal 2 − f_a);
* c⁻ ≤ 2 f_max E (only the off-diagonal cofactors can make C_r negative);
* f_j ≤ ¼ (from q⁻_l ≥ 0) and 3E + 5F ≤ 3 (from Σq⁻ ≥ 0, since Σn ≥ E);
* β ≥ 0, since det M ≥ 4 ≥ 4F ≥ c⁺.

Together, |det X| ≥ min over that triangle of (8 − 4E − 4F)(1 + E + ⅔F) − 2 min(¼, F)E(1 + ⅓F) = 7.84 > 7.29,
attained at E = 0, F = 3/5. Numerically the displayed bound itself is ≥ 8 (68 000 feasible leak samples).
**So in the main configuration with τ = −, det Λ ≥ 0.871 > 81/100.**

*The remaining step: the configuration lemma.* It is not proved that every admissible realization of
stratum 30 is in a main configuration for some labelling and some τ. Ways out of it, and what each needs:
* a row without a coordinate of mass > 1;
* |y_r| beyond the flat range of its column;
* |u(t)| > 3|x_t|.
The summed constraint alone, (S): Σ_k Φ(x_k, p_k) ≥ 108 with Φ(ξ, p) = Σ over the six vectors of |3ξ + u(p)|,
is far too weak: [Num] it admits |det X| ≈ 2. So the configuration lemma must use the six constraints
individually. Their exact form is a disjunction: (N_s) holds iff, for some T ⊆ {1,2,3,4} (the coordinates
where 3x + u_s < 0), 2P_l(T) − P_i(T) − 4P_j(T) − 3 ≥ 3(x⁺(T) + x⁻(Tᶜ)); likewise (M_s) with
5P_j + 2P_l − P_i − 12. Facts used so far:
* a row l has P_l(T_l) ≥ 3/2 on its N-set;
* if d_r > ½ then k_r ∈ T_r, so y_r⁺ ≤ ⅓;
* if d_r > 1 then k_r lies in no other N-set, so y_r⁻ ≤ ⅓.

Exact checks of the identities above are in `verify/verify_structural_exclusions.py`.

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

1. **Stratum 30: the configuration lemma (§3b) is the only open step.** Show that every admissible frame is in
   a main configuration for some labelling and some τ, or else has |det X| ≥ 7.29. The tools are the exact
   ∃T forms of the six constraints (6), the row-only bounds (3) (pairwise overlaps ≤ ½), the 4/3-vectors, and
   ‖X_1 ± X_r‖ ≥ 2. A small slack is enough: the main-configuration bounds hold up to the total "mismatch"
   δ (the summed |entries| of the six vectors 3x + u whose sign differs from the main pattern). Their sum
   weakens by 2δ, so x_t ≥ 1 + L − 2δ/9, and since the target is 7.29 and not 8, δ up to about 0.35 can be
   absorbed. If a hand proof stalls, the fallback is an exact enumeration of the sign configurations (the
   membership type of each coordinate in the six T-sets, with the sign conditions as linear constraints). That
   needs certified (rational Farkas) LP infeasibility, so it is a [P+C] route.
2. The slicing route of the earlier plan is superseded for stratum 30: at the extremal lattice the best
   slicing bound is only ≈ 0.745.
3. Stratum 17: same programme. Find its short vectors (norm ≤ 3/2: ½X_1 ± (⅓X_2 − ⅙X_3 − ⅙X_4) and
   ½(X_1 ± X_3 ∓ X_4), §4), check numerically which subset already forces |det X| ≥ 6·0.81 = 4.86, look for the
   analogous budget identity in its low configuration (one facet row split), then do the determinant algebra.
