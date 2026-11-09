# Two further exclusions from coupled planar blocks

8 November 2026. All arguments below are theoretical [P], applied to the explicit
residue descriptions of strata 28 and 31 in the existing 39-stratum classification.
No new enumeration or numerical optimization is used. Indices of coordinates in
the proofs are 1-based; stratum IDs retain their repository numbering.

## 1. A planar identity

Let A and B be independent planar vectors of l1 norm 2 whose coordinates have
no opposite signs: A_k B_k >= 0 for k=1,2. By changing coordinate signs they
can both be made nonnegative. Write

    A = (1+s/2, 1-s/2), B = (1+t/2, 1-t/2),
    -2 <= s,t <= 2,
    delta = |s-t|/2, m = (s+t)/2.

Then

    0 < delta <= 2,  |m| <= 2-delta,
    ||A-B||_1 = 2 delta,  |det(A,B)| = 2 delta.             (1)

The identity ||(x,y)||_1 = max(|x+y|,|x-y|) also gives

    min(||(-5A+3B)/8||_1, ||(3A-5B)/8||_1)
      = max(1/2, |delta-|m|/4|)
      <= max(1/2, delta).                                 (2)

For the equality, the coordinate sums are both -1/2 and the coordinate
differences are -m/4 plus or minus delta. For the inequality, |m|/4 <= 1/2
and |delta-|m|/4| <= max(delta,|m|/4). These formulas allow zero coordinates;
there is no assumption that the rows lie in facet interiors.

## 2. Stratum 28: determinant at least 9/8 [P]

Here p=8 and Gamma is generated over Z^4 by

    g = (1,1,3,3)/8.

In particular its forced unit-shell contacts include

    g, h=(2,2,-2,-2)/8, j=(3,3,1,1)/8.

Membership follows from h=2g-(0,0,1,1) and
j=3g-(0,0,1,1). All three have l1 norm 1 in frame coordinates.

Let X have minimal rows of norm 2, and suppose Gamma X is admissible. Equality
in the triangle inequality for g says that all four rows are coordinatewise
sign-compatible wherever their supports overlap. Equality for h says that
each of rows 1,2 has disjoint support from each of rows 3,4: the required
sign of their coordinatewise products is reversed. Invertibility therefore
splits X into two square planar blocks, as proved in THEORY_REVIEW §2.

Write their rows as A,B and C,D. Within each block the rows are sign-compatible,
so (1) applies with widths delta and epsilon. The norms of vectors in the
two blocks add, and

    det(Gamma X) = (2 delta)(2 epsilon)/8 = delta epsilon/2.

The two lattice vectors with frame coefficients

    g-e_3 = (1,1,-5,3)/8,
    g-e_4 = (1,1,3,-5)/8

have first-block norm ||(A+B)/8||_1 = 1/2. By (2), one of their second-block
norms is at most max(1/2,epsilon). Both vectors must have norm at least 2, so

    2 <= 1/2 + max(1/2,epsilon), hence epsilon >= 3/2.

Applying the same argument to j-e_1 and j-e_2 gives delta >= 3/2. Therefore

    det(Gamma X) >= (3/2)^2/2 = 9/8 > D*.

This excludes every admissible realization of this stratum at or below D*,
without using local minimality or further properties of a max-det frame.

## 3. Stratum 31: determinant at least 9/10 [P]

Here p=10. Among the forced unit-shell contacts are

    u=(2,4,1,3)/10,  v=(4,-2,-3,1)/10.

Comparing the signs of u_i u_j and v_i v_j forces disjoint support between
the row groups {1,4} and {2,3}. Within each group the two rows are
sign-compatible, by u. Thus X has two planar blocks whose ordered rows are

    A=X_1, B=X_4;  C=X_2, D=X_3.

Let their widths be delta and epsilon as in (1). Thus delta,epsilon <= 2 and

    det(Gamma X) = 4 delta epsilon/10.

The residue group contains (2,4,6,8)/10 and (4,8,2,6)/10, so Gamma contains

    a=(2,4,-4,-2)/10,
    b=(4,-2,2,-4)/10.

These membership assertions can also be read directly from the two generators
u and h=(0,0,1/2,1/2): a=u+h-e_3-e_4 and b=2u-e_2-e_4.
The two vectors separate exactly into row differences:

    aX = (A-B)/5 + 2(C-D)/5,
    bX = 2(A-B)/5 - (C-D)/5.

The summands occupy disjoint coordinate blocks. Admissibility and (1) imply

    delta + 2 epsilon >= 5,
    2 delta + epsilon >= 5.

Together with delta,epsilon <= 2, these yield delta,epsilon >= 3/2. Hence

    det(Gamma X) >= (4/10)(3/2)^2 = 9/10 > D*.

This deliberately sufficient bound avoids solving any auxiliary minimization
problem. As in §2, it holds for every admissible realization of the frame.

## 4. Consequence and scope

**Subsequent result:** `FACET_LAYERS.md` excludes five additional strata,
including the remaining 1+3 block case, leaving 17. The count below records
the intermediate consequence of the two proofs in this note.

Combining these two proofs with the previous 15-stratum exclusion leaves
**22 of the original 39 strata**. Their possible indices are now

    p in {1,2,3,4,5,6,7,8,9}.

The comparison is rational: D* <= d(1/2)=207/256 < 81/100, whereas both new
bounds exceed 81/100. The candidate's stratum 11 remains.

The key improvement is that mixed residue classes supply short-vector
constraints coupling the planar widths. Treating the lattice as a direct
product would miss exactly these constraints. This provides a concrete
theoretical route for further cases, but does not by itself treat connected
support graphs, including stratum 11, or the 1+3 block case (subsequently
excluded by the different argument in `FACET_LAYERS.md`).

The search algorithms and their certification issues are outside this note's
scope. No computational verification work is proposed or performed here.
