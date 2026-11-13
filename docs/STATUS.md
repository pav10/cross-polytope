# STATUS — what is established (read this before doing anything)

Date: 14 October 2026.  Tags: **[P]** proved (proof here or in `THEORY.md`, section given);
**[P+C]** proved with an exact finite computation (script named, reproducible);
**[Cite]** literature; **[Num]** numerical evidence only; **[Conj]** open.

Do **not** re-derive anything tagged [P] or [P+C]. Do **not** use anything not tagged [P]/[P+C]/[Cite] in a proof.

---

## §0 Conventions (load-bearing — two earlier campaigns died on normalisation slips)

* `‖x‖ = ‖x‖₁` on ℝ⁴; `C₄ = {‖x‖ ≤ 1}`, `vol C₄ = 2/3`; body `K = 2C₄`.
* Lattice `Λ = ℤ⁴Y`, **rows** of `Y` form a basis; coefficient row `q` ↦ `v = qY`.
* Admissible ⟺ `‖v‖ ≥ 2` for all `v ∈ Λ∖{0}` (Lemma 1.1: ⟺ translates of `C₄` pack). Normalise `λ₁ = 2`.
* `D_pack = inf det Λ` over admissible `Λ` ( = Δ(2C₄) = 16·Δ(C₄); do not call it "the critical determinant of C₄").
* Density `δ = (2/3)/det`. **Every threshold in code is a DETERMINANT: `D* = 0.80822000090353537…`, never the density `0.82485791730…`.**
* `ζ` = root of `z³−6z+3` in `(1/2,2/3)` ≈ 0.52397639708186596553; isolating interval `[5239/10000, 131/250]`.
* `M(z) = [[2−2z,2+z,−2−z,−2],[2z−4,−2z,2,−2],[−2z,−2z,2z−4,4−2z],[2−2z,2+z,2,2+z]]`, `B(z)=M(z)/4`, `Λ_ζ = ℤ⁴B(ζ)`.
* `D* = det Λ_ζ = (16−9ζ+6ζ²)/16 = 0.80822000090353537300824717153…`; density `0.82485791730145054648…`.
* **Frame coordinates.** For a frame `X` (4×4, rows = 4 independent minimal vectors), `Γ = ΛX⁻¹ ⊇ ℤ⁴`, `p = [Γ:ℤ⁴] = |det X|/det Λ`, and `F(c) = ‖cX‖₁` (so `F(eᵢ) = 2`). In data files a frame-coordinate point is stored as its **integer numerator** `n = p·c`.
* Symmetries: signed permutations (order 384) act (a) on ambient coordinates = columns of `X`, (b) on frame coordinates = rows of `X` / coordinates of `c`. Both preserve everything.

**Target Theorem Z [Conj].** `D_pack = D*`, attained exactly by the signed-permutation images of `Λ_ζ`.

---

## §1 Settled facts about the candidate (all exact; `verify/verify_audit.py` re-checks them)

| # | statement | tag | where |
|---|---|---|---|
| 1.1 | Miller's `L₄ = ℤ⁴A/84` (det `2514935/3111696`) is **not optimal**: `B(131/250)` is admissible with det `50513750079/62500000000` < det L₄. | [P+C] | audit §2 |
| 1.2 | `ℤ⁴B(z)` admissible for all `z∈[1/2,2/3]`; `det B(z) = −(16−12z+12z²−z⁴)/16`; unique minimum at `ζ`. Hence `D_pack ≤ D*`. | [P+C] | THEORY Thm 2.1 |
| 1.3 | `L₄ = B(11/21)` and `d′(11/21) = −2/9261 < 0`: L₄ is not even stationary. | [P] | THEORY Rem 2.2 |
| 1.4 | `Λ_ζ` has exactly 15 minimal pairs (on the cell `[5239/10000,131/250]`), contact rank 13, kernel `{M′, U, V}`; identity `−det(M+uU+vV) = 16(16−12z+12z²−z⁴)+16(z²−2z+2)u²+16z²uv+8(z²−2z+4)v²`; 15 strictly positive KKT multipliers. | [P+C] | audit §3, `data/zeta.json` |
| 1.5 | `Λ_ζ` is a **strict local minimum** of det among admissible lattices. | [P+C] | THEORY Thm 3.1 |
| 1.6 | Among matrices with ζ's contact signs (the whole 3-parameter cell), `|det|` is uniquely minimised at `M(ζ)`. | [P+C] | audit Prop 3.2 |
| 1.7 | Minor distribution of the 15 pairs: `{0:420,1:831,2:104,3:9,4:1}`; every independent quadruple is clean; unique max-det frame, index **p = 4**; Λ_ζ lies in **stratum 11** (§3). | [P+C] | `scripts/zeta_frame.py` |

Status: `2/3 ≤ D_pack ≤ D* < det L₄`.

---

## §2 General theory (all [P], full proofs in `THEORY.md`)

| result | content | THEORY § |
|---|---|---|
| Lemma 1.1 | packing ⟺ admissibility | §1 |
| Cor 1.4 | `det Λ ≥ 2/3` (Minkowski) | §1 |
| Lemma 1.5/1.6 | minimal vectors primitive; same-parity minimal pairs have disjoint supports | §1 |
| **Thm 1.7** | finite admissibility test for an **arbitrary** frame lattice `Γ ⊇ ℤ⁴` (not only clean frames); finiteness via van der Corput, `|cᵢ| ≤ (nᵢ+2)·4!/p`; Cor 1.8: passing the test ⟹ nonsingular, `|det X| ≥ 2p/3` | §1.7 |
| Thm 4.1/4.2 | critical lattice exists [Cite Mahler]; its minimal vectors span | §4 |
| **Thm 5.1** | at a local minimum, for any choice of supporting facets the contact functionals have rank ≥ 10 (`E ∩ B·Sym = 0`); ≥ 10 minimal pairs | §5 |
| **Thm 5.2** | ℓ¹-eutaxy: `Σ λ_q σ_qᵀ v_q = I`, `λ ≥ 0`, `Σλ = 2`, for **every** facet choice | §5 |
| Thm 5.3 | [Cite Miller arXiv:2501.09245] at most 20 minimal pairs; 20 only for D₄⁺-equivalents | §5 |
| **Thm 6.1** | max-det frame: all minimal vectors have frame coords in `Γ∩[−1,1]⁴`; all quadruple dets ≤ 1 in frame coords; `Γ ∩ {‖c‖₁<1} = {0}`; `p ≤ 24`; minimal vectors are p-modular in a basis | §6 |
| **Thm 6.2** | exactly **39 strata** (Γ mod signed perms), `p ∈ {1,…,10,12,16,18}` | §6, `scripts/make_strata.py` [P+C] |
| **Thm 7.1** | stratum p = 1: exactly 10 minimal pairs, forming `A₄` (`{±Xᵢ, ±(Xᵢ−Xⱼ)}`, an ℓ¹-equilateral 5-point set) | §7 [P+C] |
| Thm 8.1 | slicing: `det Λ ≥ Δ₃(K∩ŵ^⊥)/‖w‖₂` for primitive dual `w` (a pruning tool, not a proof route; see DEAD_ENDS) | §8, `SLICING.md` |
| **Thm 9.1** | Reduction: `D_pack` is the minimum of finitely many explicit critical values; algebraic; Theorem Z ⟺ finite Claim C | §9 |
| six types | [Cite Brunngraber/Wolff/Mordell] + independent exact HNF certificate (54 clean superlattices, 6 orbits) | audit §5 [P+C] |

---

## §3 New results of this campaign set-up (14 Oct 2026) — proofs here

**Stratum 30 excluded, 10 November 2026 [P] (exact checks: G9). 26 of 39 strata excluded; 13 remain:
0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17.** `docs/COMMON_FACET.md` §3–3c. Every admissible realization of
stratum 30 has **det Λ ≥ 196/225 = 0.871 > 81/100 > D***; 8/9 is attained (sharp axis lattice, 16 minimal pairs).
The proof uses only six vectors of Γ₃₀: three of norm 10/9 and three of norm 11/9. [Num] The three 10/9 vectors
alone allow |det X| ≈ 7.31, so the 11/9 ones are essential. The steps:
* *Per-coordinate form.* Each constraint reads Σ_k min(Pos_k + 3x_k⁺, Neg_k + 3x_k⁻) ≤ 1 or ≤ 2.
* *Structure.* Facet-row overlaps are ≤ ½, so every facet row has a main entry d_r > ½ on its own coordinate.
  Then |X_1(k_r)| ≤ ⅓, so |x_t| ≥ 1, and then every d_r > 1.
* *Linearity.* Every cost then takes its linear value, except M_s at k_{i(s)} ("failure").
* *Budget.* A budget identity (each unit of leakage costs exactly 9), and exact rational certificates for the
  8 failure patterns, give x_t ≥ 1 + L when x_t > 0.
* *Determinant.* A 3×3 determinant lemma (column trick plus Laplacian bound) then gives |det X| ≥ 8.
  For x_t < 0 a weighted sum gives ≥ 196/25 (≥ 8 numerically), or a failure forces X_1 ≤ 0 on the row
  coordinates, which gives ≥ 8.

**Frame-swap lemma, 13 November 2026 [P] / [P+C].** `docs/FRAME_SWAP.md`, exact check
`verify/verify_frame_swap.py`. If four minimal vectors Q ⊆ S have |det Q| = 1 in frame coordinates, then QX is
another max-det frame, with frame lattice ΓQ⁻¹. Consequences:
* (a) Γ′ has no point of frame norm < 1. This is already implied by Stage A's hull-interior rule.
* (b) If Γ′ is an excluded stratum, then det ≥ 0.81.
* (c) Assign each critical lattice to the first remaining stratum, in a fixed order, among those of its max-det
  frames. A stratum's search may then discard Q whose Γ′ lies in an earlier stratum.
The rules are hereditary. Single-vector swaps shrink every remaining candidate pool except stratum 0's (table in
the note); for stratum 11, 24 of 68 candidates go by (b)/(c). **[Num] Stage A for stratum 11 now completes:**
* pool rules alone: 2 477 contact sets in 126 s;
* pool rules plus the dynamic quadruple rule: **1 469 sets in 113 s**.
Before, it stood at ≥ 9 800 sets after 25 min and was still growing. ζ's S survives, and every set passes the
exact re-verification. This was a scratchpad prototype; it is not yet integrated into `src/cp4/contacts.py`.

**Stratum 11, a 7-point reformulation, 10 November 2026 [P] / [Num] / [Conj].** `docs/STRATUM11.md`.
Put A = (½,¼,¼,¼)X, B = (−½,¼,¼,¼)X and C = A + B. Then Λ = ℤ⟨A, B, X_2, X_3, X_4⟩ with the single relation
X_2 + X_3 + X_4 = 2C, and det Λ = |det[A, B, X_2, X_3]|. ζ's 15 contacts are exactly the pairs at distance 2 in
{O, A, B, C, X_2, X_3, X_4}: every pair except O–C and the X_r–X_s (exact check) [P].
Conjecture 11.A (these 19 pair constraints alone imply det ≥ D*) is **false** (13 Nov 2026) [P]: an exact rational
frame satisfies them, and ‖cX‖₁ ≥ 2 for all 68 candidates of the cube, with det = 115/144 < D*. The lattice is not
admissible, because the vector (½, ¾, −¼, −5/4), outside the cube, has norm 4/3 (`verify/verify_stratum11.py`,
STRATUM11 §3, DEAD_ENDS A14). Programme to settle stratum 11 without search (it must use vectors outside the cube):
* a local configuration theorem on ζ's face-datum region, extending audit Prop 3.2;
* a budget configuration lemma, as for stratum 30.

**Strata 17 and 30: structure and necessary conditions, 9 November 2026 [P].**
Rows 2–4 lie on one facet, Λ = (ℤc ⊕ Λ''_0) + ℤh, and det Λ = d_V(h)·‖N‖_∞/3.
*Cancellation-budget lemma:* a vector c ∈ Γ with ‖c‖₁ = 1 + ε allows total cancellation ≤ 2ε. (The case ε = 0
is the forced-shell lemma.) Consequences:
* stratum 30: Σ_k min(3|X_1k|, 2R_k) ≤ 3 and pairwise row overlaps ≤ ½;
* stratum 17: Σ_k min(3|X_1k|, R_k) ≤ 2 and Σ_k min(2X_2k, X_3k+X_4k) ≤ 1.
For stratum 17 the axis case gives 16/9, but its low configurations are not near-axis [Num] (relaxed minimum
≈ 0.958, with one facet row split). The plan for 17 is COMMON_FACET §6.3: the same programme as stratum 30.

**Stratum 0 needs one contact set, 9 November 2026 [P].** A local minimum in stratum 0 has the A4 contact
set {±(x_i − x_j)} (THEORY 7.1). In stratum 0 every four independent minimal vectors form a basis (max-det
index 1), so for any vertex v the star X_i = x_i − x_v is a max-det frame. Hence B1 on the **star embedding** of
A4 (one of the three Stage A sets, the one with frame stabilizer of order 48) covers all of stratum 0. The path
and fork embeddings are redundant. With frame symmetry its Knuth estimate is ≈ 3.7·10⁴ nodes (PLAN §3).

**Numerical triage of the next targets, 9 November 2026 [Num].** Constrained Nelder–Mead from feasible
starts, with full admissibility via short vectors.
* *Pair problem*: admissible Λ containing minimal u, v with disjoint supports and (u+v)/2 ∈ Λ, which covers
  strata 1, 8, 9, 10, 15. With the two other frame rows of norm 2 (stratum 1's form) the search found
  det 0.853 (supports 2+2), 0.867 (1+3) and 0.922 (1+2) (`scripts/probe_pair_stratum1.py`). A realization-wide "pair lemma" at 81/100 would
  therefore have only ~5 % slack; not pursued in that form (DEAD_ENDS B17).
* *Common-facet strata* (shell vector (0,⅓,⅓,⅓), rows 2–4 on one facet; the search keeps rows 2–4 ≥ 0):
  `scripts/probe_common_facet.py`; best found stratum 4: 0.891 (40 restarts); stratum 17: 1.030; stratum 30: 1.028 (partial runs). So 17 and
  30 may have up to ~27 % slack over 81/100; stratum 4 up to ~10 %. *Caveat:* these are upper bounds on each
  stratum's infimum from non-converged searches. The best stratum-30 point (det 1.035) has only 7 minimal
  pairs, and a local minimum needs ≥ 10 (Thm 5.1), so the true infima may be lower. The same caveat applies to
  the pair-problem values.
  Old catalog (`data/catalog_local_optima.json`): every local optimum found lies in stratum 11, except one in
  stratum 14 with det exactly 1.0.

**Coordinate sections, 9 November 2026 [P] / [P+C]: strata 6, 7, 25 excluded; 14 remain.**
`docs/COORDINATE_SECTIONS.md`. *Lemma [P]:* if Λ ∩ ℝ^T has full rank r with ℓ¹ covering radius μ_T, then the
projection of Λ to ℝ^{T^c} has ℓ¹ minimum ≥ 2 − μ_T, so det Λ ≥ det(Λ∩ℝ^T)·β_{4−r}((2−μ_T)/2)^{4−r}. Two axis rows
(X_i = ±2e_k) with Γ_{ij} = ℤ²+ℤ(½,½) give det ≥ 1. Three pairwise support-disjoint rows force two axis rows.
This excludes **stratum 6** (det ≥ 1). In **strata 7 and 25** it leaves one support pattern. There a
fiber-offset argument applies: the distances of the four coset offsets to 2ℤ² sum to 4 identically, so
every 4-window of g(k) = ‖x + (k/4)y‖₁ sums to ≥ 4. With a certified one-dimensional inequality
(`verify/verify_coordinate_sections.py`, exact rational grid with a Lipschitz margin) this gives
**det ≥ 81/100 > D*** [P+C]. [Num] The residual pattern of stratum 7 has inf det ≈ 1. Also, in strata 1, 8, 9,
10, 15 every realization with X_3, X_4 both axis rows (in stratum 8 also X_1, X_2) has det ≥ 1.
**Remaining 14: 0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 14, 15, 17, 30**, with indices {1,2,3,4,5,6,9}.

**Review of the 8–9 November notes (9 November 2026).** Every data-dependent claim of BLOCK_COUPLING and
FACET_LAYERS is now checked exactly by `verify/verify_structural_exclusions.py`: cyclic generators and
residue groups, shell vectors, the identities for aX and bX, the planar identity, the section volume 20/3 by
exact triangulation, and D* < D(13/25) = 8089/10000 < 81/100. All proofs check, including COMMON_DESCENT
against THEORY 1.7 and INDEX_ONE_GEOMETRY.
*Remarks.*
(i) The layer theorem is THEORY 8.1 with w = σ/2: s(1111) = 5/6 there, ×8 = 20/3. The new ingredient is the
arithmetic proof that σ/2 ∈ Λ*.
(ii) The generators of strata 28, (1,1,3,3)/8, and 31, (2,4,1,3)/10, are also positive with numerators
summing to p, so FACET_LAYERS alone excludes them too (BLOCK_COUPLING's bounds 9/8 and 9/10 are stronger).
(iii) THEORY_REVIEW §4.1 is right that B1/Stage A rejections by float LPs (slack ≤ 1e-9) are not certified.
That also covers the G6 run, and closing it needs exact Farkas certificates for every rejection.
(iv) The layer theorem as an in-search prune of B1 is inert (DEAD_ENDS B16).
(v) For ζ's S the frame stabilizer has order 24. B1's orbit reduction uses only the ambient column group,
so symmetry breaking on the frame rows could save up to 12× on the G6 gate (not on generic S).

**Common descent, 9 November 2026 [P].** For an equilateral A4 configuration,
the intersection of the finite admissibility test sets over all spanning-tree
bases is sufficient. The proof uses the common integer-valued affine mass
mu(b)=sum |b_i|/2; see `docs/COMMON_DESCENT.md`. At mass two, only the 15
disjoint-edge tests can fail. Higher masses have not been eliminated, and the
intersection size is not yet known. No additional stratum is excluded.

**Index-one geometric restrictions, 9 November 2026 [P].**
`docs/INDEX_ONE_GEOMETRY.md` records the general facet-generation obstruction,
its consequence that no facet can contain four independent contacts in an
index-one candidate below 5/6, and the 15 strict separation inequalities for
midpoints of disjoint edges of its equilateral simplex. These are necessary
conditions, not a new determinant bound. The remaining count stays **17**.

**Facet-layer obstruction, 8 November 2026 [P].** If an admissible lattice has
sigma.Lambda=2Z for a sign vector sigma, the central section of 2C4 has volume
20/3 and layer spacing 1, giving **det Lambda >= 5/6** by three-dimensional
Minkowski. A strictly positive unit-shell cyclic generator forces this layer
structure. This excludes strata **12, 13, 16, 18, 19**, without enumeration.
Full proof, including the section-volume calculation: `docs/FACET_LAYERS.md`.
Combined with the reductions below, **22 of 39 strata are excluded; 17 remain**:
**0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 25, 30**.
The possible indices are **{1,2,3,4,5,6,8,9}**. Stratum 11 remains. These are
theoretical exclusions; existing search code and its certification status are
unchanged.

**Further theoretical reduction, 8 November 2026 [P].** The explicit residue
coupling in strata **28** and **31** gives determinant lower bounds **9/8** and
**9/10**, respectively. These exceed D*, so both strata are excluded. The proofs
use forced-contact sign compatibility and two planar row-width inequalities;
no new enumeration or numerical optimization is needed. See
`docs/BLOCK_COUPLING.md`. Together with the next paragraph, **17 of 39 strata
are excluded, leaving 22**, with indices in **{1,2,3,4,5,6,7,8,9}**. Stratum 11
remains. The existing support-reduction checker covers the earlier 15 exclusions
only; the two new exclusions are proved directly in the note.

**Update, 8 November 2026 — forced-support reduction [P+C].** Equality in the
triangle inequality for forced unit-shell contacts forces selected frame rows
to have disjoint supports. Block determinant bounds exclude strata **20, 21,
22, 23, 24, 26, 27, 29, 32, 33, 34, 35, 36, 37, 38** before Stage A, leaving
**24** strata (including 11). All indices 12, 16, 18 are excluded. Proof and
limitations: `docs/THEORY_REVIEW.md` §§2–3. Exact integer/rational checker:
`python3 verify/verify_support_reduction.py` (PASS). The search code has not yet
been changed to apply these exclusions.

**Certification qualification, 8 November 2026.** Existing floating LP/hull
branch rejections do not yet carry exact certificates or proved numerical error
bounds; rechecking surviving leaves does not establish completeness. Read
`docs/THEORY_REVIEW.md` §4 before treating Stage A/B enumeration outputs as a
global proof certificate. This does not retract the separately checked exact
candidate identities, strata regeneration, or covector certificate.

Notation: `Λ` critical, `X` a max-det frame, `Γ`, `p`, `F` as in §0, and
`S` = set of **all** minimal vectors of Λ in frame coordinates (one per ± pair).

**Lemma A.1 (hull closure) [P].** Let `P = conv(±S)`. Then
(a) `Γ ∩ int P = {0}`; (b) every nonzero point of `Γ ∩ P` lies in `±S`.

*Proof.* `F` is convex with `F = 2` on `±S`, so `F ≤ 2` on `P`.
(b) For `c ∈ Γ∩P`, `c≠0`: admissibility gives `F(c) ≥ 2`, so `F(c) = 2`, i.e. `c` is minimal.
(a) Let `c ∈ int P`, `c ≠ 0`. If `F(c) ≥ 2` then `F(c) = 2`; for small `ε>0`, `(1+ε)c ∈ P` but
`F((1+ε)c) = 2(1+ε) > 2`, contradicting `F ≤ 2` on `P`. So `F(c) < 2`, contradicting admissibility. ∎

**Corollary A.2 (volume) [P].** `vol(P) ≤ 16/p`. *(Minkowski for Γ, `det Γ = 1/p`, applied to (a).)*

**Consequences used by Stage A.** (A1) `S ⊇ {c ∈ Γ : ‖c‖₁ = 1}` (such `c` has `F(c) ≤ 2`); (A2) quadruple dets ≤ 1 (Thm 6.1); (A3) Lemma A.1; (A4) Cor A.2; (A5) `10 ≤ |S| ≤ 20` (Thm 5.1, Miller). (A2),(A3a),(A4) and `|S| ≤ 20` are hereditary; (A3b) is a forcing rule. Candidate points: `S ⊆ Γ∩[−1,1]⁴∩{‖c‖₁ ≥ 1}` (stored per stratum in `data/strata.json`).

**Theorem B.1 (critical-value reduction) [P].** Let `s_c ∈ {−1,0,1}⁴` be the exact sign vector of `cX` for `c ∈ S` (the *face datum*), and `A` the affine space of matrices `Y` with
`Σⱼ s_cj (cY)ⱼ = 2` and `(cY)ⱼ = 0` whenever `s_cj = 0` (all `c ∈ S`). Then
(i) `X ∈ A`, `dim A ≤ 6`; (ii) `X` is a critical point of `det|_A`; (iii) hence `p·D_pack = |det X|` is the
absolute value of a **critical value** of the polynomial `det|_A`.

*Proof.* For `Y ∈ A` close to `X`, nonzero coordinates of `cY` keep their signs and zero ones stay zero, so
`‖cY‖₁ = Σ s_cj(cY)ⱼ = 2` for `c ∈ S`; non-contacts keep a uniform gap above 2 (THEORY Lemma 5.0).
So `ΓY` is admissible and `|det Y|/p ≥ D_pack = |det X|/p`: `X` is a local minimum of `|det|` on `A`,
and `det X ≠ 0`, giving (ii). All directions of `A` are two-sided feasible, so the argument of Thm 5.1 gives
`(dir A) ∩ X·Sym₄ = 0`, i.e. `dim A ≤ 6`. ∎

**Corollary B.2 [P].** Theorem Z follows if, for every stratum `Γ`, every `S` passing (A1)–(A5), and every
face datum whose system is consistent with `dim A ≤ 6` and whose strict sign conditions are feasible,
no critical value `v` of `det|_A` satisfies `2p/3 ≤ |v| < p·D*`, and `|v| = p·D*` occurs only for
the data of `Λ_ζ`. *(Lower cut `2p/3` from Cor 1.4.)* **The certificate per leaf is a univariate polynomial
(elimination of the critical-value variable) with no real root in the band — no need to locate critical
points, verify cells, or test admissibility, except for leaves with roots in the band.**

**Lemma B.3 (covector column pre-filter) [P].** Let `X` be nonsingular with face datum `(s_c)_{c∈S}` and let
`Yⱼ` be the j-th column of `X`. Then each column `(s_cj)_{c∈S}` is the covector `(sign(c·Yⱼ))_{c∈S}` of the
central arrangement `{c^⊥ : c ∈ S}` in ℝ⁴ at a point `Yⱼ ≠ 0`. Hence a B1 node (decided contacts `D`) with a
partial column that is not the restriction to `D` of a covector realised by a nonzero vector has no
descendant cell containing a nonsingular `X`, and can be discarded without an LP; Cor B.2 only concerns
cells containing a critical (so nonsingular) `X`. The test commutes with the orbit reduction of B1.
*Proof.* `(cX)ⱼ = c·Yⱼ`, so `s_cj = sign(c·Yⱼ)` by the definition of the face datum, and `det X ≠ 0` gives
`Yⱼ ≠ 0`. An ambient element `g = (π, ε)` fixing the decided sign vectors maps column `k` of the partial
matrix to `ε_k` times column `π(k)`, and the covector set is closed under negation. ∎

**Covector certificate [P+C]** (`cp4.covectors.certified_covectors`). The float LP list of covectors is
accepted only if every listed cell contains an exact rational point (exact sign check), every zero set is a
flat (exact ranks), and for every flat `F` the number of listed covectors with zero set exactly `F` equals
`r_F = Σ_{G ⊇ F} |μ(F,G)|`, the number of regions of the arrangement restricted to `⋂_{i∈F} c_i^⊥`
(Zaslavsky 1975 [Cite]). Distinct genuine covectors with the right count on every flat are all of them.
ζ's S: 117 flats (by rank 1, 15, 55, 45, 1), **1 705 faces, 372 topes**, 13 s (gate G8,
`scripts/check_covectors_zeta.py`).


**Theorem B.1′ (supporting facets) [P]** (13 Nov 2026). Let `X` be a critical frame with contact set `S` and face
datum `(s_c)`. Fix any `T ⊆ S` (in B1: the four frame vectors) and keep `T` on its exact faces. For `c ∈ S∖T`
choose any **supporting facet** `σ_c ∈ {±1}⁴`, i.e. `σ_cj = s_cj` whenever `s_cj ≠ 0` (free where `(cX)_j = 0`).
Let `A′` be cut out by the equations of Theorem B.1 for `c ∈ T` and by `σ_c·(cY) = 2` for `c ∈ S∖T`. Then `X ∈ A′`,
`X` is a local minimum of `|det|` on `A′`, and (ii), (iii) of B.1, Cor B.2 and `dim A′ ≤ 6` hold with `A′` in
place of `A`.
*Proof.* For `Y ∈ A′` near `X`: `c ∈ T` as in B.1; `c ∈ S∖T`: `‖cY‖₁ ≥ σ_c·(cY) = 2` (Hölder, `‖σ_c‖_∞ = 1`);
the non-contacts by THEORY Lemma 5.0. So `ΓY` is admissible for all `Y ∈ A′` near `X`, in both directions of
every line; the rest is the proof of B.1 and THEORY Thm 5.1 (which is already stated for every choice of
supporting facets). ∎
*Use.* Non-frame contacts need only the 16 facets, never the 64 lower-dimensional faces with zero coordinates.
That was the dominant cost of B1: on ζ's S the facet-only pass is complete in 14 304 nodes (frame symmetry) while
the full one was not, at ≥ 1.7·10⁶. The leaf families have one equation per non-frame contact (zero equations
only for frame rows), so they can be larger than the B.1 families, but `dim ≤ 6` still holds.

**Lemma B.4 (half-open facet cells) [P].** For a contact `c` and any `ρ ∈ {±1}⁴` (the *tie side*), the 16 cells
`H_ρ(σ) = {Y : σ·(cY) = 2; σ_j(cY)_j ≥ 0 if σ_j = ρ_j, σ_j(cY)_j > 0 if σ_j ≠ ρ_j}`, `σ ∈ {±1}⁴`, partition
`{Y : ‖cY‖₁ = 2}`. `Y` lies in the cell with `σ_j = sign((cY)_j)`, and `ρ_j` where `(cY)_j = 0`; that `σ` is a
supporting facet of `cY`. Hence B1 may branch each non-frame contact over the 16 cells with any `ρ`, chosen
afresh at each node, and every critical `X` reaches exactly one leaf. That leaf's affine space is an `A′` of
Theorem B.1′.
* *Symmetry.* Below the frame no orbit reduction is applied (trivial group), so no compatibility of `ρ` with the
  ambient group is needed; the frame rows keep the exact faces, the stabiliser chain and the frame-symmetry
  isomorph rejection unchanged.
* *Pre-filter.* Lemma B.3 extends: a decided entry `σ_j = ρ_j` allows the true sign `σ_j` **or 0**
  (`ColumnFilter.allowed(..., zero_ok)`).
* *Choice of ρ* (efficiency only): `ρ_j = +1` iff `(cY)_j > 0` is possible on the node's cell (covector filter,
  else two closure LPs). This avoids spurious cells pinned to `(cY)_j = 0`. With a fixed `ρ` (frame-derived), or
  with closed cells and no tie rule, the trees are 2–20× (fixed `ρ`) and up to 10²–10⁴× (closed, duplicated
  boundaries, heavy tails) larger (scratchpad Knuth, DEAD_ENDS B18).
* *Determine-then-check.* At rank 16, decided non-frame contacts are checked against their half-open cells, and
  an undecided one only needs `‖cX‖₁ = 2` exactly.

**Lemma B.5 (ball rows) [P].** At a critical `X`, `‖cX‖₁ = 2` for every `c ∈ S`. So the rows
`−2 ≤ τ·(cY) ≤ 2` (`τ ∈ {1}×{±1}³`), i.e. `‖cY‖₁ ≤ 2`, may be added for **every** `c ∈ S` to every LP of B1
(decided or not) without losing `X`. *Not* DEAD_ENDS B14 (look-ahead LPs per undecided contact) or B12 (the same
inequality as a Stage-A prune, inert there). Inside the B1 LPs it is the convex half of the undecided contacts'
`F = 2`. Measured, B1 (Knuth, 400 probes): 2–20× fewer nodes at 1.4–2× the LP cost.

Implementation: `enumerate_face_data` / `estimate_face_data(..., supporting=True, ball=True)` (static order, frame
first); witness gate G6s (`scripts/stageB_witness_supporting.py`: ζ, and D4 with the frame (e₁±e₂, e₃±e₄), whose
non-frame contacts have 28 zero coordinates); complete run on ζ's S, gate G6S. Defaults are unchanged: the
facet-only frame-symmetry pass reproduces 14 304 nodes and the same LP-dead, filter-dead and symmetry-dead counts.
---

## §4 Computations already done, with their exact outcomes

| computation | outcome | script / data |
|---|---|---|
| audit checker | PASS (all items of §1, 54 clean superlattices, test-set counts 1156/182/196/65/29/16) | `verify/verify_audit.py` |
| strata | 39 orbits, 1040 superlattices, p ∈ {1..10,12,16,18}; candidate pools 40–220, forced points 4–20 | `data/strata.json` |
| ζ's stratum | stratum **11** (p=4, orbit 16, |stab| 24, 68 candidates, 4 forced) | `data/zeta.json` |
| Stage A, stratum 0 (p=1) | **exactly 3** contact sets, all of size 10 = the A₄ configuration in its 3 frame embeddings; 54 726 nodes, 13 s | `results/stageA_0.json` |
| Stage A, stratum 11 | still running at handover: ≥ 9 500 contact sets after 1.46 M nodes (~25 min, ~1 000 nodes/s); kills: det 82 %, interior 4 %, excluded-forced 5 %, volume 0 — **(A1)–(A5) leave ~10⁴ sets per mid-size stratum** (PLAN §2). Output file written only at the end; G4 (ζ's S present) checked by `verify/run_all.py` once it exists | `results/stageA_11.log/.json` |
| Stage B1 witness, ζ's S | ζ's face datum (all 15 contacts interior to facets) survives every prune when added contact by contact; allowed-sign counts 80,80,80,35,53,53,80,53,35,35,35,12,2,6,4 | `scripts/stageB_witness_zeta.py` |
| Stage B1 full, ζ's S | 1 500 s: 7 596 nodes, 0 leaves, incomplete; per-depth nodes 1,1,1,8,17,54,45,2455,1721,1550,891,514,274,64 — **every branch explored so far dies at depth 14**; ~5 nodes/s. 6 h rerun launched at handover | `results/stageB_faces_zeta*.log/.pkl` |
| Stage B2, ζ's face datum | dim A = 3; eliminant `(D−6)(D³+62D²+8D+200)(64D³−1920D²+13584D−26011)`, sympy grevlex+FGLM, 21 s; the last factor has the root `4·D*` (= `p·D*`, p = 4); no other root in `[8/3, 4D*)` — the leaf certificate of Cor B.2 works as designed | `scripts/stageB_critical_zeta.py`, `src/cp4/critical.py` |
| catalog `results.csv` (2169 searched lattices) | 26 distinct density values; **all 25 non-trivial ones have ζ's exact contact set** (stratum 11, 15 pairs); several are points of ζ's own family (z ≈ 0.5037, 0.548, 0.6254, 2/3). It is the output of the earlier ternary-configuration search, **not** an independent catalog of local optima. | `data/catalog_local_optima.json` |
| PLAN §2 option 4(a) ("A6", frame-facet LP) | kills **0** of 200 Stage-A sets (stratum 0, complete) and **0** of 931 (stratum 11, 30 000 nodes); sampled sets keep 187–270 of 330 facet patterns, ζ's S keeps 295 ⇒ fails the decision rule (DEAD_ENDS B12) | `scripts/probe_A6_frame_facets.py` |
| merged search prototype (PLAN §4, frame signs first) | G7a: ζ's labels survive every prune (15 in / 53 out, leaf dim A = 3, rank 13). Knuth estimate, 200 probes: stratum 0 ≈ 5.3·10⁷ nodes (median 3.4·10⁶), stratum 11 ≈ 4.2·10⁸ (median 1.8·10⁷); ~240 nodes/s ⇒ infeasible as designed (DEAD_ENDS B13) | `src/cp4/merged.py`, `scripts/merged_run.py --estimate`, `results/merged_{0,11}_estimate.json` |
| LP backend for float prunes | highspy and scipy give **identical verdicts on 3 011** LPs recorded from B1 walks on ζ's S; 0.66 vs 2.06 ms/LP | `src/cp4/lp.py`, `scripts/check_lp_backend.py` |
| Stage B1 on ζ's S, Knuth estimates (300 probes) | facet-only: 3.5·10⁴ (static) / 1.1·10⁴ (dynamic) nodes; full: 5.6·10⁵ (static, ~35 LPs/node) / 6.2·10⁵ (dynamic, ~116 LPs/node) | `results/stageB_faces_zeta_estimate.json` |
| **Stage B1 facet-only pass on ζ's S (G6f)** | **complete**: 65 148 nodes, 1 404 s, 6 leaves, all dim A = 3; **ζ's face datum is among them**; an interrupted + resumed run gives the identical leaf set. Facet-only is a validation pass, not complete for the proof (zero coordinates excluded) | `scripts/stageB_faces_zeta.py --facet-only --static --determine`, `results/stageB_faces_zeta_facet.{pkl,log}` |
| Facet-class layer prune inside B1 (ζ's S) | facet-only 65 148 → 64 094 nodes (−1.6 %), same 6 leaves; full static tree: same 25 depth-3 subtrees in 90 708 vs 90 819 nodes over 900 s; 0 of 46 + 6 known leaves killed — inert (DEAD_ENDS B16) | `scripts/measure_layer_prune.py`, `results/layer_prune_zeta.json` |
| B1 frame symmetry (ζ's S, stratum 0) | facet-only exact: 14 304 vs 82 069 nodes, leaves agree up to G_S×AMB, ζ found; Knuth (full): ζ 3.5·10⁶ → 3.8·10⁵, star-A4 1.8·10⁶ → 3.7·10⁴ | `scripts/check_frame_symmetry.py`, `results/frame_symmetry_*.json` |
| Structural verifiers | `verify_support_reduction.py` (15 strata), `verify_structural_exclusions.py` (BLOCK_COUPLING, FACET_LAYERS: 7 strata), `verify_coordinate_sections.py` (strata 6, 7, 25): all PASS, exact | `verify/` |
| Stage B1 full run on ζ's S (G6) | **incomplete; to be finished on the cluster** (`scripts/slurm_stageB_zeta.sbatch`, resumes from the committed checkpoint `results/stageB_faces_zeta.pkl`). Run here 16 Oct 2026 (`--static --determine --prefilter --resume`, ≈ 100–120 nodes/s): ≈ 1.73 M nodes in two long legs (2 h 10 min, 2 h 9 min), **46 leaves, ζ's face datum among them**; 288 of the **1 633 depth-3 subtrees** finished, all in the first depth-1 branch (contact 0 in a facet, 381 subtrees); every rank-16 read-off dead. The session container is reclaimed a few minutes after the session goes idle, so later legs (≈ 5 min) could not finish subtree 289 (> 50 000 nodes); local resumption stopped. **The Knuth estimate (5.6·10⁵) is ≥ 3× too low** (heavy tail: probe maxima 10⁷–10⁸). Per-branch Knuth (150 probes; contact 0 = facet / 1 zero / 2 zeros / 3 zeros): 1.1·10⁶ / 9.1·10⁴ / 1.0·10⁵ / 6.6·10³; the first branch dominates and already exceeds its estimate. Remaining work, roughly 1–3 M nodes (3–8 h single-core), is a guess | `scripts/knuth_G6_by_branch.py`, `scripts/stageB_faces_zeta.py`, `scripts/slurm_stageB_zeta.sbatch`, `results/stageB_faces_zeta_full.log` |
| Column-wise ("covector") B1 on ζ's S, cost | arrangement `{c^⊥}`: 1 705 covectors, N = 852 classes mod ±, T = 186 topes (ζ's columns among them); partial LP never prunes (1 000/1 000 at levels 1–3) ⇒ C(N+3,4) = 2.2·10¹⁰ leaf LPs (5.2·10⁷ facet-only), ~10²–10³× row-wise B1 (DEAD_ENDS B15) | `scripts/estimate_covector_B1.py`, `results/covector_B1_estimate.json` |
| Covector pre-filter in row-wise B1 (ζ's S, static order) | facet-only complete pass **319 s vs 1 420 s** (4.5×), identical tree (65 148 nodes) and 6 leaves, 213 807 vs 1 030 938 LPs; full run **93 vs 9.6 nodes/s** (9.7×) over 600 s; LP-subsumed: 0 of 20 031 rejected candidates LP-feasible | `scripts/measure_covector_prefilter.py`, `results/covector_prefilter_*.{json,log}` |
| **Covector pre-filter built into B1** (`enumerate_face_data(prefilter=True)`, Lemma B.3 + certificate) | validated against the search without it on ζ's S: facet-only complete pass through the real code path **337 s vs 1 404 s**, identical per-depth counts (65 148 nodes) and **identical leaf set** (6 leaves, ζ hit); static full order, first 2 000 nodes: identical tree, 33 s vs 209 s (incl. the 13 s certificate); dynamic order, first 1 000 nodes: identical tree, no speed-up (closure LPs dominate). In all three, LP-dead without the filter = LP-dead with it + filter rejections **exactly** (965 791 = 148 660 + 817 131; 150 015 = 11 395 + 138 620; 5 953 = 2 392 + 3 561) | `scripts/compare_prefilter_zeta.py`, `results/covector_prefilter_validation.{json,log}`, `results/stageB_faces_zeta_facet_prefilter.{pkl,log}` |

## §5 Numerical evidence (never cite as proof)

* [Num] Stratum p = 1 (A₄): best determinant found ≈ 0.8363 (+3.5 % over D*), 263 restarts with exact shortest-vector checks; other local values 0.8674, 0.9075 seen.
* [Num] Slicing is tight at Λ_ζ: `Δ₃(K∩ŵ^⊥)/‖w‖ = 0.808220001` for `w = (−1,0,0,1)B*`.
* [Num] Cohn–Elkies grid LP (degree 5): ≈ 0.984 (octahedron, truth 18/19) and ≈ 0.913 (C₄) — not rigorous, not sharp.
* [Num] No numerical search to date has found an admissible lattice with det < D*.
