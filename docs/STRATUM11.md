# Stratum 11 (ζ's stratum): a 7-point reformulation and a programme to settle it without search

10 November 2026; §3 revised 13 November 2026. §§1–3 are proved [P] (the finite facts are checked exactly by the
scripts quoted there); §4 is the plan. Conventions as in `THEORY_REVIEW.md`.

Stratum 11 cannot be excluded, because Λ_ζ lies in it (STATUS 1.7). The goal is to prove that **every admissible
realization of stratum 11 has det ≥ D*, with equality only at Λ_ζ**. If that is done by theory, stratum 11 needs
no Stage A or B search. That search is the campaign's real bottleneck: Stage A has found ≥ 9 800 contact sets
in stratum 11, and B1 costs 10⁵–10⁶ nodes per set.

## 1. The 7-point picture [P]

Γ₁₁ = ℤ⁴ + ℤg with g = (½, ¼, ¼, ¼), cyclic of order 4 (2g ≡ (0, ½, ½, ½)). For a frame X, put

    A = gX = ½X_1 + ¼(X_2 + X_3 + X_4),   B = (−½, ¼, ¼, ¼)X,   C = A + B = ½(X_2 + X_3 + X_4).

Then A − B = X_1, so **Λ = ℤ⟨A, B, X_2, X_3, X_4⟩ with the one relation X_2 + X_3 + X_4 = 2C**. The four vectors
A, B, X_2, X_3 form a basis, and det Λ = |det[A, B, X_2, X_3]| = |det X|/4.

Take the seven lattice points O = 0, A, B, C, X_2, X_3, X_4. Admissibility puts every pair at ℓ¹-distance ≥ 2.
**At Λ_ζ, exactly 17 of the 21 pairs are at distance 2.** Their differences are the 15 contact vectors of ζ
(A−C ≡ O−B and B−C ≡ O−A as vectors). The four longer pairs are O–C and X_r–X_s.
* So O, A, C, B is a rhombus of side 2 whose diagonal A–B also has length 2, while O–C is longer.
* Each X_r is at distance 2 from all four vertices of the rhombus.
* The X_r are pairwise more than 2 apart, and they sum to 2C.

(Exact check, frame-coordinate arithmetic on `data/zeta.json` and `data/strata.json`: all seven points lie in
Γ₁₁; X_2 + X_3 + X_4 = 2C and A − B = e_1; the 17 pair differences are contacts of ζ, covering all 15 up to
sign; and the 4 others are not.)

Stated without coordinates, ζ's configuration is a **7-point set in an admissible lattice whose "distance-2
graph" is K₇ minus a triangle (on the X_r) and minus one disjoint edge (O–C), with one linear relation.**

## 2. What is already proved about stratum 11

* Λ_ζ is a strict local minimum of det among admissible lattices (STATUS 1.5, THEORY Thm 3.1) [P+C].
* On ζ's 3-parameter cell (ζ's contact signs, all 15 contacts tight), |det| is uniquely minimised at M(ζ)
  (STATUS 1.6, audit Prop 3.2) [P+C]. The identity behind it is
  −det(M + uU + vV) = 16(16 − 12z + 12z² − z⁴) + 16(z² − 2z + 2)u² + 16z²uv + 8(z² − 2z + 4)v².
* Along ζ's family B(z), det(z) = (16 − 12z + 12z² − z⁴)/16 stays below 81/100 for z up to about 0.6. So the
  sublevel set {det < 0.81} of stratum 11 is not a small neighbourhood of ζ: it contains that whole arc of the
  family, its 3-parameter cell, and the interior points near them. A proof must work with the cell, not with
  a neighbourhood of one point.

## 3. Is the 7-point system already enough? No [P]

Relaxation: keep only the 19 constraints ‖(p − q)‖₁ ≥ 2 for the 21 pairs of §1 (the 15 contacts, O–C, and the
three X_r–X_s), plus ‖X_i‖₁ = 2. Then minimise |det X|.

**Conjecture 11.A (10 Nov 2026) is false (13 Nov 2026) [P].** The rational frame

    X = [[−23/12, 1/12, 0, 0], [1/2, −1/4, −5/4, 0], [0, 2, 0, 0], [−1/3, −1/4, −3/4, −2/3]]

satisfies all 19 constraints and ‖X_i‖₁ = 2, and has det = |det X|/4 = **115/144 = 0.7986 < D***. It even
satisfies ‖cX‖₁ ≥ 2 for **all 68 candidates** of stratum 11 (every Γ₁₁ point of [−1,1]⁴ with ‖c‖₁ ≥ 1). It is
not a lattice counterexample: c = (½, ¾, −¼, −5/4) = X_1 + X_2 − X_4 − A lies in Γ₁₁ *outside the cube* and has
‖cX‖₁ = 4/3. (Exact: `verify/verify_stratum11.py`.)
* So no relaxation by the vectors of the cube [−1,1]⁴ alone can prove det ≥ D* on stratum 11. Vectors outside
  the cube can only be ignored *after* admissibility is known (they are then never minimal, by the max-det
  property). A relaxation proof must include some of them, e.g. (½, ¾, −¼, −5/4).
* Method [Num]: sequential polyhedral descent. Fix the facets σ_v = sign(vX), take the inner approximation
  ⟨σ_v, vX⟩ ≥ 2 and the frame rows in the ℓ¹ ball, minimise det by SLSQP, update σ and repeat. From random
  starts, 1–3 of every 25–60 runs end below D* (at 115/144 or 0.80787), whereas penalty Nelder–Mead in the original test
  which never left 4D*. The earlier evidence for 11.A was an artefact of that weak searcher (DEAD_ENDS A14).

## 4. Programme (theory first; mirrors the exclusion of stratum 30, COMMON_FACET §3b–3c)

1. **Per-coordinate form** of the 19 constraints. Each ‖p − q‖ ≥ 2 becomes
   Σ_k min(Pos_k + …, Neg_k + …) ≤ budget, as in COMMON_FACET (8).
2. **Local configuration theorem.** Extend audit Prop 3.2 from the cell (all 15 contacts tight) to ζ's whole
   face-datum region (the same sign patterns, contacts ≥ 2). On that region the 19 constraints are linear, so
   the claim is that det ≥ D* on a polyhedron, with equality only at ζ.
   Route: an identity of the audit type, det − 4D* = Σ λ_c(‖cX‖ − 2)·(positive) + (positive definite part),
   with ζ's 15 strictly positive KKT multipliers λ_c (STATUS 1.4) as the starting point.
3. **Configuration lemma.** Show that any admissible frame of stratum 11 with det < 0.81 either has ζ's face
   datum or has det ≥ D* by a budget argument (as COMMON_FACET §3c Steps 1–3). The 7-point picture suggests
   the structure step: the rhombus O, A, C, B and the three apexes X_r.
   **Caveat (§3):** the 19 constraints are not enough; the configuration lemma must use vectors outside the cube.
4. Consequence: stratum 11 settled without search. With stratum 0 (≈ 4·10⁴ B1 nodes, PLAN §3) this would leave
   only the strata not yet excluded by theory.
