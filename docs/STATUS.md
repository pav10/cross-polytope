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

## §5 Numerical evidence (never cite as proof)

* [Num] Stratum p = 1 (A₄): best determinant found ≈ 0.8363 (+3.5 % over D*), 263 restarts with exact shortest-vector checks; other local values 0.8674, 0.9075 seen.
* [Num] Slicing is tight at Λ_ζ: `Δ₃(K∩ŵ^⊥)/‖w‖ = 0.808220001` for `w = (−1,0,0,1)B*`.
* [Num] Cohn–Elkies grid LP (degree 5): ≈ 0.984 (octahedron, truth 18/19) and ≈ 0.913 (C₄) — not rigorous, not sharp.
* [Num] No numerical search to date has found an admissible lattice with det < D*.
