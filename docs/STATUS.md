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
| Stage B1 full run on ζ's S (G6) | not run here: ≈ 12 nodes/s at the top of the tree, expected 10–40 h; resumable cluster job prepared | `scripts/slurm_stageB_zeta.sbatch` |
| Column-wise ("covector") B1 on ζ's S, cost | arrangement `{c^⊥}`: 1 705 covectors, N = 852 classes mod ±, T = 186 topes (ζ's columns among them); partial LP never prunes (1 000/1 000 at levels 1–3) ⇒ C(N+3,4) = 2.2·10¹⁰ leaf LPs (5.2·10⁷ facet-only), ~10²–10³× row-wise B1 (DEAD_ENDS B15) | `scripts/estimate_covector_B1.py`, `results/covector_B1_estimate.json` |
| Covector pre-filter in row-wise B1 (ζ's S, static order) | facet-only complete pass **319 s vs 1 420 s** (4.5×), identical tree (65 148 nodes) and 6 leaves, 213 807 vs 1 030 938 LPs; full run **93 vs 9.6 nodes/s** (9.7×) over 600 s; LP-subsumed: 0 of 20 031 rejected candidates LP-feasible | `scripts/measure_covector_prefilter.py`, `results/covector_prefilter_*.{json,log}` |
| **Covector pre-filter built into B1** (`enumerate_face_data(prefilter=True)`, Lemma B.3 + certificate) | validated against the search without it on ζ's S: facet-only complete pass through the real code path **337 s vs 1 404 s**, identical per-depth counts (65 148 nodes) and **identical leaf set** (6 leaves, ζ hit); static full order, first 2 000 nodes: identical tree, 33 s vs 209 s (incl. the 13 s certificate); dynamic order, first 1 000 nodes: identical tree, no speed-up (closure LPs dominate). In all three, LP-dead without the filter = LP-dead with it + filter rejections **exactly** (965 791 = 148 660 + 817 131; 150 015 = 11 395 + 138 620; 5 953 = 2 392 + 3 561) | `scripts/compare_prefilter_zeta.py`, `results/covector_prefilter_validation.{json,log}`, `results/stageB_faces_zeta_facet_prefilter.{pkl,log}` |

## §5 Numerical evidence (never cite as proof)

* [Num] Stratum p = 1 (A₄): best determinant found ≈ 0.8363 (+3.5 % over D*), 263 restarts with exact shortest-vector checks; other local values 0.8674, 0.9075 seen.
* [Num] Slicing is tight at Λ_ζ: `Δ₃(K∩ŵ^⊥)/‖w‖ = 0.808220001` for `w = (−1,0,0,1)B*`.
* [Num] Cohn–Elkies grid LP (degree 5): ≈ 0.984 (octahedron, truth 18/19) and ≈ 0.913 (C₄) — not rigorous, not sharp.
* [Num] No numerical search to date has found an admissible lattice with det < D*.
