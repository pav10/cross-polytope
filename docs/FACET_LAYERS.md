# A facet-layer obstruction: five further strata excluded

8 November 2026. The arguments are theoretical [P]. The stratum applications use
their explicit cyclic generators from the existing 39-stratum classification;
no new enumeration, LP, or numerical optimization is involved.

Throughout, K={x in R^4 : ||x||_1 <= 2}, and admissibility means that nonzero
lattice vectors have l1 norm at least 2. Determinants of lattices are covolumes.

## 1. The arithmetic layer theorem [P]

**Theorem.** Suppose an admissible full-rank lattice Lambda has a sign vector
sigma in {+1,-1}^4 for which

    sigma.Lambda = 2 Z.

Then det Lambda >= 5/6. In particular Lambda cannot attain or improve D*.

**Proof.** Coordinate sign changes allow us to take sigma=(1,1,1,1). Let
H={x : sum_i x_i=0} and Lambda_0=Lambda intersect H. The map
x -> (sum_i x_i)/2 is a surjective homomorphism Lambda -> Z. Thus Lambda_0
is a rank-three lattice in H and consecutive lattice layers parallel to H
have Euclidean spacing 2/||sigma||_2=1. The usual base-times-height formula gives

    det Lambda = det_H Lambda_0.

The section K intersect H is a centrally symmetric three-dimensional convex
body of Euclidean volume 20/3; a direct calculation is provided below.
Its relative interior contains no nonzero point of Lambda_0. Minkowski's
theorem in dimension three therefore yields

    20/3 <= 8 det_H Lambda_0,

and hence det Lambda >= 5/6. Finally
D* <= d(1/2)=207/256 < 5/6. QED.

**Section volume, without a numerical calculation.** Use coordinates
(x_1,x_2,x_3) on H, with x_4=-x_1-x_2-x_3. The Euclidean volume element is
2 dx_1 dx_2 dx_3: the coordinate Gram matrix is I_3+11^T and has determinant 4.
On H, the sum t of the positive coordinates equals the absolute sum of the
negative coordinates, and ||x||_1 <= 2 is equivalent to t <= 1.

* In an orthant with one positive coordinate, the three negative magnitudes
  are nonnegative and have sum at most 1. The coordinate simplex has volume
  1/6, so the Euclidean volume is 1/3. There are four such orthants, and four
  with three positive coordinates, contributing 8/3 in total.
* In an orthant with two positive coordinates, take
  x=(a,b,-c,-a-b+c), where a,b >= 0, a+b <= 1 and 0 <= c <= a+b.
  Its coordinate volume is

      integral_{a,b>=0, a+b<=1} (a+b) da db = 1/3.

  Its Euclidean volume is 2/3. There are six such orthants, contributing 4.

The orthants meet only on sets of three-dimensional measure zero. The total
is therefore 8/3+4=20/3.

## 2. Detecting the layers from a forced contact [P]

Let Gamma contain Z^4, and let X be a nonsingular matrix with minimal frame
rows X_1,...,X_4, each of l1 norm 2. Suppose Lambda=Gamma X is admissible.

**Lemma.** If Gamma contains c with c_i>0 for every i and sum_i c_i=1, then
there is a sign vector sigma such that sigma.X_i=2 for all four rows.

**Proof.** The vector v=cX is nonzero and belongs to Lambda. Hence

    2 <= ||v||_1 <= sum_i c_i ||X_i||_1 = 2.

Choose a supporting sign vector sigma with sigma.v=||v||_1=2. Each
sigma.X_i <= ||X_i||_1=2. Since their strictly positive weighted average
equals 2, every sigma.X_i equals 2. QED.

Consequently, for every a in Gamma,

    sigma.(aX) = 2 sum_i a_i.                              (1)

**Corollary.** Under the lemma's assumptions, if

    sum_i a_i is an integer for every a in Gamma,

then det(Gamma X) >= 5/6.

Indeed, (1) puts sigma.Lambda inside 2Z, and any frame row realizes 2, so
equality sigma.Lambda=2Z holds. Apply the layer theorem.

This gives a particularly simple sufficient condition:

    Gamma = Z^4 + Z c,  c_i>0,  sum_i c_i=1.              (2)

For (2), the required integrality follows directly. The argument does not
assume local optimality or maximality of the chosen frame. Positivity on all
four coordinates is essential to this proof: a zero coefficient leaves that
row's supporting-facet value undetermined.

## 3. Five applications [P]

Each of the following frame lattices satisfies (2):

| Stratum ID | p | Cyclic generator c |
|---|---|---|
| 12 | 4 | (1,1,1,1)/4 |
| 13 | 5 | (1,1,1,2)/5 |
| 16 | 6 | (2,2,1,1)/6 |
| 18 | 6 | (3,1,1,1)/6 |
| 19 | 7 | (1,1,2,3)/7 |

Every displayed generator has order p modulo Z^4 (it has a numerator
coordinate equal to 1), belongs to the recorded residue group of order p,
and therefore generates that group. All coordinates are positive and the
numerators sum to p. Thus every admissible realization in these strata has
determinant at least 5/6>D*.

The last remaining 1+3 block case, stratum 18, is included. Its block geometry
does not need separate optimization: the arithmetic facet-layer theorem
already excludes it. Stratum 16 also shows that a full block decomposition
is not necessary for this obstruction.

Together with THEORY_REVIEW §3 and BLOCK_COUPLING §§2–3, this excludes
**22 of 39 strata**, leaving exactly the following **17** IDs:

    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 25, 30.

Their indices belong to {1,2,3,4,5,6,8,9}. The candidate's stratum 11 remains.

## 4. Why the arithmetic hypothesis matters

Merely placing four independent contacts on one common facet does not imply
the layer theorem's hypothesis for the whole lattice. Additional residue
classes can create finer layers. More generally, if sigma.Lambda=(2/q)Z,
the same calculation gives only det Lambda >= 5/(6q). The cases above work
because their cyclic generator preserves the integral sum of frame coordinates.

This is a specific arithmetic application of slicing, not a proposed universal
sharp slicing bound. It avoids the earlier failed slicing approach: no choice
of an unverified short dual vector or conjectured section optimum is needed.

For the remaining strata, a contact with fewer than four nonzero coefficients
does not yet give a fixed lattice-layer spacing by this argument. The present
proof neither excludes those cases nor establishes global optimality or its
equality classification.
