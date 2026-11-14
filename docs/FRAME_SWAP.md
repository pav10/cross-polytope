# Frame swaps: alternative max-det frames prune the contact-set search

13 November 2026. §§1–3 are proved [P]; the per-stratum counts in §4 are exact finite computations [P+C]
(script named there). Conventions as in STATUS §0: a max-det frame X has rows that are minimal vectors and
maximises |det| among such quadruples; Γ = ΛX⁻¹ ⊇ ℤ⁴ has index p; S = all minimal vectors in frame
coordinates.

## 1. The swap lemma [P]

**Lemma 1.** Let X be a max-det frame of an admissible Λ, and let Q be four vectors of S (rows, in frame
coordinates) with |det Q| = 1. Then X′ = QX is again a max-det frame of Λ, and its frame lattice is
Γ′ = ΛX′⁻¹ = ΓQ⁻¹. It has the same index p.

*Proof.* The rows of QX are minimal vectors, and |det QX| = |det X|, which is maximal. Hence X′ is a max-det
frame. Also ΛX′⁻¹ = ΓXX⁻¹Q⁻¹ = ΓQ⁻¹. ∎

The simplest case is a single minimal vector c with |c_i| = 1: replace row i of the identity by c.

**Corollary 2 (three consequences).** In the situation of Lemma 1:
(a) **Impossible.** Γ′ cannot contain a vector c′ ≠ 0 with ‖c′‖₁ < 1. (Equivalently, the residue group of
Γ′ has no nonzero element of Lee weight < p.) Indeed ‖c′X′‖₁ ≤ Σ|c′_j|·2 < 2 would contradict
admissibility. So **no admissible lattice has these four vectors as minimal vectors**.
(b) **Excluded stratum.** If Γ′ is (up to signed permutation) a stratum excluded by proof, then det Λ ≥ 81/100.
This holds because every exclusion proof applies to any admissible realization with a frame of norm-2 rows,
and X′ is a max-det frame anyway.
(c) **Canonical stratum.** Fix a total order on the remaining strata and assign each critical lattice to the
*first* stratum among those of all its max-det frames. Then, in the search for stratum s, any Q ⊆ S whose
Γ′ lies in a remaining stratum earlier than s may be discarded. Such a lattice is assigned to that earlier
stratum and is found by its search.

*Soundness of (c).* Let Λ be critical, let s be its canonical stratum, and let X be a frame realising s. Every
unimodular quadruple Q ⊆ S(Λ, X) gives a frame of Λ, so its stratum comes no earlier than s. Hence the search
of s with rule (c) still contains (Λ, X). The other strata use the same rule, so every critical lattice is
found exactly where it is assigned. Strata excluded by proof are covered by (b). ∎

Rules (a)–(c) are **hereditary**: if S fails them, so does every superset. So they can prune Stage A's
include/exclude tree at the moment a quadruple is completed, as well as filter the candidate pool at the root.

## 2. Remark on ζ

Λ_ζ has a unique max-det frame (STATUS 1.7). So the only unimodular quadruple in ζ's S is its frame, and the
rules never touch ζ. Every rule above is satisfied by ζ's contact set.

## 3. Candidate-level consequences (single-vector swaps)

For a candidate c with |c_i| = p/p = 1 (numerator ±p), Lemma 1 applies with Q = I with row i replaced by c.
The candidate is then
* *impossible* if some such swap violates (a);
* *det ≥ 0.81* if some swap lands in an excluded stratum (b);
* *cross-stratum* if a swap lands in another remaining stratum, which rule (c) handles.

## 4. Counts [P+C] (exact; `verify/verify_frame_swap.py`, gate G9)

| stratum | p | pool | impossible | swap into excluded | cross-stratum (targets) |
|---|---|---|---|---|---|
| 1 | 2 | 58 | 2 | 0 | 24 (2, 3) |
| 2 | 2 | 52 | 4 | 0 | 20 (1, 3) |
| 3 | 2 | 48 | 8 | 0 | 28 (1, 2) |
| 4 | 3 | 64 | 9 | 0 | 18 (5) |
| 5 | 3 | 56 | 9 | 0 | 21 (4) |
| 8 | 4 | 84 | 12 | 32 | 24 (9, 10) |
| 9 | 4 | 82 | 10 | 18 | 32 (8, 10, 11) |
| 10 | 4 | 74 | 10 | 22 | 16 (8, 9, 11) |
| 11 | 4 | 68 | 14 | 15 | 9 (9, 10) |
| 14 | 5 | 72 | 16 | 18 | 0 |
| 15 | 6 | 98 | 21 | 24 | 8 (17) |
| 17 | 6 | 92 | 27 | 12 | 7 (15) |

Stratum 0 (p = 1) is unaffected, since every swap stays in Γ = ℤ⁴.

**Stratum 11.** Its 68 candidates fall into 19 classes under the stabiliser (order 24). ζ's S uses 5 classes:
the frame, A and B (5/4), C − X_r (3/2) and A − X_r, B − X_r (7/4), 15 vectors in all. Of the other 53:
* 14 are *impossible*;
* 15 force det ≥ 0.81, since they swap into stratum 7 or 12;
* 9 swap into strata 9 or 10, which rule (c) handles with order 9, 10 before 11.

Fifteen non-ζ candidates remain open, in four classes, all without a ±1 coordinate (except the last):
−C (3/2), A + C − X_r (9/4), A + C (11/4) and X_1 + X_r − X_s (3; it swaps only back into 11).

## 5. Implementation (13 November 2026)

`src/cp4/frameswap.py` (`FrameSwapRule`) implements rules (a)–(c) exactly, with the canonical order
`sorted(REMAINING)`. `cp4.contacts.enumerate_contact_sets(..., swap=FrameSwapRule(strata, sid))` applies them as
Stage-A rule (A6):
* the candidates whose single swap with the frame fails are excluded at the root;
* every new quadruple of S with |det| = p⁴ (exact) must pass `ok`.

`scripts/stageA_run.py` uses it by default (`--no-swap` to disable). `verify_frame_swap.py` checks that the root
exclusions agree with its own independent classification in every remaining stratum.

**Stratum 11 [P+C]: 1 469 contact sets in 102 s** (`results/stageA_11.json`, gate G4: ζ's S present). That is
168 523 nodes, 38 root exclusions and 19 349 branches cut by the quadruple rule. The sets are identical to the
scratchpad prototype's. Stratum 0 is unchanged (no swap applies at p = 1; G3).

