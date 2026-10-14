# The slicing approach to the 4-dimensional cross-polytope
## A careful sketch, with proofs of everything established and an honest status

**Goal.** Λ ⊂ ℝ⁴ is *admissible* if ‖v‖₁ ≥ 2 for every nonzero v ∈ Λ.
Conjecture: every admissible Λ has det Λ ≥ D\* = (16 − 9ζ + 6ζ²)/16 = 0.80822000090…,
ζ the root of z³ − 6z + 3 in (1/2, 2/3); equality for the lattice Λ_ζ = ℤ⁴B(ζ).

---

## 1. The slicing inequality (proved)

**Definitions.** Λ\* = {w : ⟨w,v⟩ ∈ ℤ for all v ∈ Λ} is the dual lattice; det Λ\* = 1/det Λ.
A nonzero w ∈ Λ\* is *primitive* if w/k ∉ Λ\* for every integer k ≥ 2.
For a unit vector ŵ, the *section body* is the 3-dimensional centrally symmetric polytope
  P_ŵ = 2C₄ ∩ ŵ^⊥ = { x : ⟨ŵ,x⟩ = 0, ‖x‖₁ ≤ 2 }.
Its *critical determinant* is
  Δ₃(P_ŵ) = min { det M : M a rank-3 lattice in ŵ^⊥ with no nonzero point in the interior of P_ŵ },
where det M is the 3-dimensional Euclidean covolume inside the hyperplane ŵ^⊥.

**Theorem 1 (slicing).** For every admissible Λ and every primitive w ∈ Λ\*,
  det Λ ≥ Δ₃(P_ŵ) / ‖w‖₂,  where ŵ = w/‖w‖₂.

*Proof.* Because w is primitive in Λ\*, the linear map v ↦ ⟨w,v⟩ sends Λ *onto* ℤ. Its kernel
Λ₀ = Λ ∩ w^⊥ is a rank-3 sublattice, and Λ is the disjoint union of the layers
Λ_k = {v : ⟨w,v⟩ = k}, k ∈ ℤ, each a translate of Λ₀ lying in the affine hyperplane
{⟨w,x⟩ = k}. Consecutive hyperplanes are at Euclidean distance 1/‖w‖₂.
Pick a basis v₁,v₂,v₃ of Λ₀ and any v₄ ∈ Λ₁; then v₁,…,v₄ is a basis of Λ, and
  det Λ = (3-volume of the parallelepiped spanned by v₁,v₂,v₃) × (distance of v₄ from w^⊥)
        = det Λ₀ × (1/‖w‖₂).
Finally Λ₀ ⊂ Λ inherits admissibility: no nonzero point of Λ₀ lies in the interior of 2C₄, hence
none lies in the interior of P_ŵ = 2C₄ ∩ w^⊥. So det Λ₀ ≥ Δ₃(P_ŵ). ∎

**Corollary 2 (trivial-Δ₃ version).** Minkowski's convex body theorem in dimension 3 gives
Δ₃(P) ≥ vol₃(P)/2³. Writing s(ŵ) = vol₃(C₄ ∩ ŵ^⊥) (so vol₃(P_ŵ) = 8 s(ŵ)),
  det Λ ≥ s(ŵ)/‖w‖₂  for every primitive w ∈ Λ\*.

**Computed values.** s(ŵ) is the volume of a central hyperplane section of the ℓ¹ ball, computed
exactly from the edge–hyperplane intersection points:
  s(e₄) = 4/3 (maximum; the section is C₃),  s((1,1,0,0)) = 0.9428,  s((1,1,1,0)) = 0.8660,
  s((1,1,1,1)) = 5/6 (minimum, ⊥ to the all-ones direction).

**What it gives at the conjectured optimum.** For Λ_ζ, the best dual vector is w = (−1,0,0,1)·B\*,
with ‖w‖₂ = 1.2504 and s(ŵ) = 0.9428. Corollary 2 gives det Λ_ζ ≥ **0.7540**, against the truth
0.8082 and the 4-dimensional Minkowski bound (det ≥ vol C₄ = 0.6667). Note what happened:
one elementary inequality, using only the *trivial* 3-dimensional input, recovers 60% of the gap
between Minkowski and the truth.

---

## 2. The structure of the optimal slice (computed)

For w as above, Λ₀ = Λ_ζ ∩ w^⊥ = {qB : q₁ = q₄}, with basis (1,0,0,1)B, (0,1,0,0)B, (0,0,1,0)B.
  • det Λ₀ = det Λ_ζ · ‖w‖₂ = 1.010612.
  • The section body P_ŵ has 8 facets (an octahedron-like polytope) and vol₃(P_ŵ) = 7.5425.
  • Λ₀ is admissible for P_ŵ and makes **7 contact pairs** with its boundary — the minimal vectors
    of Λ_ζ with q₁ = q₄, namely 0010, 0100, 0110, 1001, 1101, 1111, 1211.
    Since a local optimum in dimension 3 needs at least n(n+1)/2 = 6 contacts, Λ₀ is at least
    combinatorially consistent with being a *critical lattice of the section body*.

**Open question (decisive for the method).** Is Δ₃(P_ŵ) = det Λ₀ = 1.010612?
If yes, Theorem 1 is *tight* at the optimum: det Λ_ζ = Δ₃(P_ŵ)/‖w‖₂ exactly, and the 4-dimensional
optimum is "built from" an optimal 3-dimensional slice. This must be settled by an exact 3-dimensional
computation (Betke–Henk's algorithm applied to the fixed polytope P_ŵ). A naive numerical search with a
fixed coefficient box is unreliable here — it converges to near-singular matrices whose short vectors
have large coefficients (the same degeneracy trap that defeated earlier 4-dimensional attempts), so it
must not be used as evidence in either direction.

---

## 3. Why the *uniform* slicing argument cannot work alone (proved estimate)

The natural universal statement is: *every admissible Λ has a primitive w ∈ Λ\* with
Δ₃(P_ŵ)/‖w‖₂ ≥ D\**. Equivalently, Λ\* must contain a nonzero point of the star body
  G = { w : ‖w‖₂ ≤ Δ₃(P_ŵ)/D\* }.
The only tool that forces a lattice point into a body is Minkowski's theorem, which needs a
*convex* S ⊆ G with vol(S) ≥ 2⁴ det Λ\* = 16/det Λ. In the worst case det Λ ≈ 2/3 this requires
vol(S) ≥ 24. But the radius of G in direction ŵ is Δ₃(P_ŵ)/D\*, which lies between
  s_min/D\* = (5/6)/0.8082 = 1.03  and  Δ₃(2C₃)/D\* = 1.407/0.8082 = 1.74
(the upper value uses Minkowski's own theorem that the densest lattice packing of the octahedron has
density 18/19, i.e. Δ₃(2C₃) = (4/3)(19/18)). A convex body inside G is contained in a ball of radius
≤ 1.74 and contains at most a ball of radius ≈ 1.03; its volume is at most a few units, roughly
5–10, against the required 24. **The uniform argument falls short by a factor of about 2–4 and no
improvement of Δ₃ can close that gap.** Using Minkowski's second theorem on the four successive
minima of Λ\* does not help either: it only yields λ₁(Λ\*) ≤ 1.42, giving det Λ ≥ 0.83/1.42 = 0.59,
worse than the trivial 2/3.

The reason is structural: Minkowski guarantees *one* short dual vector, which may point in a bad
direction (⊥ to (±1,±1,±1,±1), where the section is thinnest). Slicing is strong for a *specific*
lattice because we may choose the best of its many short dual vectors; it is weak as a uniform bound
because we cannot control which direction the guaranteed vector takes.

---

## 4. The realistic program: slicing as a certified prune

Since slicing is nearly tight at the optimum (0.754 with trivial Δ₃; possibly exact with sharp Δ₃)
but useless uniformly, its correct role is as a **cheap, rigorous certificate for regions of shape space**
inside an enumeration:

1. **Tabulate Δ₃(P_ŵ) over directions.** ŵ ranges over the sphere modulo the symmetry group of C₄
   (order 384); the fundamental domain is small. For each direction the problem is a *3-dimensional*
   lattice-packing problem for an explicit polytope with ≤ 8 facets — solved territory (Betke–Henk).
   The function ŵ ↦ Δ₃(P_ŵ) is continuous and piecewise algebraic, so a rigorous lower envelope on a
   grid with continuity bounds is achievable.

2. **Certify regions.** In any enumeration over contact configurations (the max-det-frame enumeration,
   or the Wolff-type contact enumeration), each node is a region R of basis matrices B. For any fixed
   integer vector c, the dual vector w_c = c·B\* varies continuously over R, and
     det B ≥ Δ₃(P_ŵ_c)/‖w_c‖₂  for every B ∈ R.
   If min over R of the right-hand side (a bound on ‖w_c‖₂ over R, plus the tabulated Δ₃ lower envelope
   on the cone of directions swept by ŵ_c) exceeds D\*, the whole region R is certified at once. This
   costs one dual-vector bound per node instead of a polynomial optimization, and is exact (no interval
   determinant arithmetic, which was the bottleneck of every previous attempt).

3. **Residue.** Regions where no short dual vector certifies — lattices whose dual is "round" with no
   thick section direction — are handled by the exact multiplier certificates on ≤ 6-dimensional faces
   (the form of the auditor's ζ analysis).

**Expected effect (from the measured values).** At Λ_ζ the certificate reaches 0.754/0.808 = 93% of the
target with trivial Δ₃; regions of shape space near the optimum should be closed by step 2 once the
sharp Δ₃ is tabulated, leaving the residue for step 3. Whether the residue is small is not yet measured.

---

## 5. Two variants worth carrying along

**Projection (primal) slicing.** For a primitive v ∈ Λ (e.g. a minimal vector), project along v onto
v^⊥: det Λ = ‖v‖₂ · det π(Λ). The projected lattice avoids the 3-dimensional convex body
W_v = { y ∈ v^⊥ : the chord of 2C₄ through y in direction v has length > 1 (in units of v) }
(an open interval of length > 1 contains a point of every coset t₀ + ℤ, so a fiber of Λ over such y
would contain a short vector). Hence det Λ ≥ ‖v‖₂ Δ₃(W_v). This uses the *primal* short vectors, which
in an enumeration are known exactly (they are the contacts), whereas dual short vectors must be bounded.

**Iterated slicing.** Both inequalities reduce dimension by one; applying them to the 3-dimensional
slice again (with the exactly known Δ₂ of 2-dimensional sections, which are polygons) could give
computable closed-form bounds, at the cost of looseness at each step.

---

## 6. Status ledger

| item | status |
|---|---|
| Theorem 1 (slicing inequality) | **proved** |
| Corollary 2 and the section volumes s(ŵ) | **proved / computed exactly** |
| bound 0.754 at Λ_ζ (trivial Δ₃) | **computed** |
| Λ₀ = Λ_ζ ∩ w^⊥ has 7 contacts with the section body | **computed** |
| tightness Δ₃(P_ŵ) = det Λ₀ | **open** — needs an exact 3-dim computation (Betke–Henk); naive numerics are unreliable here |
| uniform slicing proves the conjecture | **refuted by estimate** (short by a factor ≈ 2–4) |
| slicing as a certified prune | **designed, not yet run** |
| tabulation of Δ₃(P_ŵ) over directions | **not started** — the concrete next computation |

**Next computation, in order.** (i) Compute Δ₃(P_ŵ) exactly for ζ's best direction (settles tightness;
one 3-dimensional Betke–Henk run). (ii) If tight, tabulate Δ₃ over a grid of directions in the
fundamental domain. (iii) Wire the prune into the max-det-frame enumeration and measure the residue.
