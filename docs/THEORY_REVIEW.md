# Theory review and a reduction before Stage A

8 November 2026. This note separates a new proved reduction from certification
gaps and possible next steps. It does not establish Theorem Z. The external audit
files and existing search algorithms have not been changed.

## 1. Where the proof stands

The candidate's admissibility, strict local optimality, and uniqueness within its
specified contact region have explicit exact certificates. In particular, the
audit's contact-region proof derives 1/2 <= z <= 2/3 from the sign inequalities;
its optimization is not merely along the one-parameter family.

The global bottleneck is classification of possible contact configurations and
their realizations. Even the full B1 gate for the candidate's contact set is
unfinished. The numerical catalog is not an independent classification of local
optima. The old runtime estimates in THEORY §§10–11 are superseded by the measured
failures in STATUS and DEAD_ENDS. Faster LPs do not resolve this structural issue.

## 2. Forced contacts split the ambient coordinates [P]

Let Gamma contain Z^4, with index p, and let X be a nonsingular frame matrix
whose rows X_i have l1 norm 2. Assume Gamma X is admissible. Define

    U = {c in Gamma : ||c||_1 = 1}.

For every c in U, admissibility and the triangle inequality give

    2 <= ||cX||_1 <= sum_i |c_i| ||X_i||_1 = 2.

Consequently equality holds in every coordinate's triangle inequality. Thus,
for every ambient coordinate k, all nonzero numbers c_i X_ik have the same sign.
In particular,

    c_i c_j X_ik X_jk >= 0.

**Disjoint-support lemma.** If U contains c and d with c_i c_j > 0 and
d_i d_j < 0, then X_ik X_jk = 0 for every k. Hence rows i and j have disjoint
supports. This conclusion uses just two forced contacts, not an LP or a choice
of the unknown remaining contacts.

Construct a graph on the four row indices, joining i and j unless the lemma
certifies that they have disjoint supports. Let its connected components be
I_1,...,I_m. Every nonzero column of X is supported inside one component:
otherwise two rows from different components would overlap in that column.
If k_a columns belong to I_a, full row rank requires k_a >= |I_a|. Since there
are four columns in total, equality holds for every a. Permuting rows and columns
therefore makes X block diagonal, with square nonsingular blocks X_a.

This does **not** assert that Gamma or Gamma X splits as a direct sum. The
residue group can couple different blocks; the following argument retains p.

Put Gamma_a = Gamma intersect R^{I_a}, and let p_a be its index over Z^{I_a}.
The lattice Gamma_a X_a is admissible in its |I_a|-dimensional coordinate space.
If b_d is a valid minimum-2 l1 determinant lower bound in dimension d, then

    |det X_a| >= p_a b_|I_a|,
    det(Gamma X) = |det X|/p >= (product_a p_a b_|I_a|)/p.       (1)

Minkowski alone gives b_1 = b_2 = 2. These are the only dimensional bounds
needed for the exclusions below.

**Axis case.** If every pair of rows is disjoint, each row occupies exactly one
coordinate, with entry of absolute value 2. Thus |det X| = 16 exactly and the
lattice determinant is 16/p. For the recorded strata p <= 18, this is at least
8/9, strictly above D*. This excludes such frames for a global optimum even if
they occur in other admissible lattices.

## 3. Consequence for the 39 strata [P+C]

The small checker `verify/verify_support_reduction.py` regenerates the unit shell
from each residue group, checks it against the stored forced contacts, and uses
only integer signs and rational bounds in (1). Completeness of the 39-stratum
list remains the independently checked prerequisite G1; this script does not
claim to replace that classification.

| Stratum IDs | p | Block sizes | Lower bound for lattice determinant |
|---|---|---|---|
| 20, 22, 24, 27 | 8 | 1+1+1+1 | 2 (exact) |
| 21, 23, 26 | 8 | 2+1+1 | 1 |
| 29 | 9 | 1+1+1+1 | 16/9 (exact) |
| 32, 33 | 12 | 1+1+1+1 | 4/3 (exact) |
| 34, 35, 36, 37 | 16 | 1+1+1+1 | 1 (exact) |
| 38 | 18 | 1+1+1+1 | 8/9 (exact) |

Every bound is greater than 81/100, whereas
D* <= d(1/2) = 207/256 < 81/100. Thus **15 strata can be removed before Stage A**,
leaving 24, including the candidate's stratum 11. In particular all p=12,16,18
strata disappear. No contact-set or face enumeration is needed for these cases.
This is a reduction in the domain of the proof, not a measured runtime saving;
some existing searches might already reject some of these cases by other means.

Run: `python3 verify/verify_support_reduction.py`.

## 4. Proofreading and certification findings

1. **Search completeness is not yet certified.** `faces.lp_feasible` rejects an
   LP when the floating solver fails or reports slack <= 1e-9. Exact strict
   feasibility requires slack > 0; no separation bound or exact dual certificate
   justifies replacing this with a fixed tolerance. `allowed_signs` likewise
   discards small positive or negative ranges. Taking closures is conservative
   mathematically, but floating endpoints and these thresholds do not inherit
   that guarantee. The certified covector filter is different: its rational
   witnesses and exact cardinality check supply a completeness argument.

2. **Stage A has the same distinction.** Hull membership, forcing, volume,
   and determinant rejection use floating arithmetic. `verify_contact_set`
   checks only A1, A2, A5 on surviving outputs; its reference to
   `verify/test_regression.py` points to a file absent from this checkout.
   Verifying survivors cannot certify that rejected branches were impossible.
   These are gaps in the proof-producing implementation, not demonstrated
   incorrect outputs. Exact rejection witnesses or proved error bounds are needed.

3. **The equality case needs a fiber argument.** A critical-value eliminant can
   establish the determinant bound, but a root equal to p D* does not identify
   every lattice attaining it. Corollary B.2 should explicitly require that every
   admissible equality point represents a signed-permutation image of Lambda_zeta.
   The audit supplies that argument for the specified candidate contact region;
   each additional equality region needs an equivalence or uniqueness proof.

4. **Positive-dimensional critical loci remain an implementation gap.**
   `critical_value_polynomial` raises when its ideal is not zero-dimensional.
   This does not invalidate the mathematical finite-critical-values reduction.
   A rational polynomial is constant on each connected component of its real
   critical locus: a semialgebraic path can be partitioned into differentiable
   arcs, along each of which the derivative is zero. There are finitely many
   components. The resulting finite value set is semialgebraic over Q and hence
   consists of algebraic numbers. This also supplies the algebraicity step that
   THEORY 9.1's short proof leaves implicit. It gives no bound on elimination cost.

5. **Check both determinant signs.** The reduction concerns |det X|. A general
   leaf certificate must check both signed windows, or prove and apply an
   orientation convention. `roots_in_window` checks only its supplied interval;
   the positive-window wording in `critical.py` is not by itself a full certificate.

6. **Normalization and wording.** THEORY 3.1 writes Y*=M(zeta) while using the
   term admissible, defined earlier at minimum 2. Its proof works at minimum 8,
   as the section introduction says. State that explicitly in the theorem's
   displayed neighborhood claim, then scale by 1/4. “The index-one stratum is
   completely classified” in §7 means its contact combinatorics, not all geometric
   realizations or their optimum. Several script references in THEORY are legacy
   names; use the current files mapped in STATUS. External literature statements
   and precise equality hypotheses have not been independently source-checked here.

## 5. Next theoretical targets

The initial reduction leaves 24 strata. The subsequent theoretical proofs in
`BLOCK_COUPLING.md` exclude two more, and `FACET_LAYERS.md` excludes another
five, leaving 17. These should be approached by structure before further B1 runs.

* Extend equality-in-the-triangle-inequality reasoning to certified convex
  dependencies among contacts. If w=sum a_i v_i is a contact, a_i>0 and
  sum a_i=1, then each ambient column has no cancellation among these v_i.
  This imposes joint sign restrictions rather than independently branching over
  every contact's face. The forced-shell lemma is the special case with frame rows.
  Its usefulness on stratum 11 is unmeasured: its forced shell is only the frame.
* The 2+2 cases (strata 28 and 31) are now excluded in `BLOCK_COUPLING.md`:
  cross-block residues give bounds 9/8 and 9/10. The bound (1) alone was too weak.
  The 1+3 case is subsequently excluded by the arithmetic layer theorem in
  `FACET_LAYERS.md`. Replacing Gamma by a direct product loses essential
  information and can yield a false bound.
* Strengthen B2's necessary conditions using the determinant Hessian. For a
  two-sided feasible direction H at X, put A=X^{-1}H. Stationarity gives tr A=0,
  and second-order necessity for |det| gives tr(A^2)<=0. The whole quadratic form
  H -> -tr((X^{-1}H)^2) must be positive semidefinite on the feasible tangent
  space. This can reject stationary saddles, but a measured reduction in the
  number or complexity of leaves has not yet been established.

This note proves the initial 15-stratum exclusion. With `BLOCK_COUPLING.md` and
`FACET_LAYERS.md`, the total is 22. The remaining 17 strata and global equality
classification are open.
