# Index-one geometry: necessary restrictions, not yet an exclusion

9 November 2026. All lemmas below are theoretical [P], conditional on the
established A4 contact classification in THEORY 7.1. No new determinant lower
bound or exclusion is claimed.

Let x_0=0,x_1,...,x_4 be an affinely independent equilateral set of side 2
whose differences generate an admissible lattice Lambda. At a local minimum
with maximal-frame index one, its ten minimal pairs are exactly
±(x_i-x_j), by THEORY 7.1.

## Facet-generation obstruction

If the contacts on a facet sigma.x=2 generate Lambda as an abelian group,
then sigma.Lambda=2Z. The layer theorem in FACET_LAYERS.md gives
det Lambda >= 5/6. This applies even if the generating contacts contain no
four-element lattice basis.

For the index-one stratum, any four independent minimal vectors form a lattice
basis: their determinant ratio to Lambda is a positive integer at most the
maximal-frame index 1. Consequently, a lattice of this stratum with determinant
below 5/6 cannot have four independent contacts on one facet.

Fix sigma in {±1}^4 and put h_i=sigma.x_i. Since the distances are 2,
max h_i-min h_i <= 2. If equality holds, let L and U be the sets of points at
the minimum and maximum heights. The facet contacts are exactly
x_u-x_l with u in U and l in L: equality in the scalar-product bound proves
both directions. Their graph is the complete bipartite graph on L and U.
Affine independence of the five points implies that these differences have
linear rank |L|+|U|-1 (a spanning tree is independent and spans the differences).

Thus a candidate with determinant below 5/6 must satisfy

    |L|+|U| <= 4

in every signed-coordinate direction of width 2. In particular, at least one
point must lie strictly between the two supporting planes in each such direction.
This is a necessary condition, not a determinant bound for simplices satisfying it.

## Admissibility adds a different geometric restriction

For four distinct indices i,j,k,l, the vector

    w = x_i+x_j-x_k-x_l

belongs to Lambda and is nonzero by affine independence. It is not any of the
ten contact pairs: its unique affine coefficient vector has four nonzero
entries, whereas a difference of vertices has two. Therefore, at the local
minimum under consideration, ||w||_1>2. Equivalently, the midpoints of every
two disjoint edges obey

    ||(x_i+x_j)/2 - (x_k+x_l)/2||_1 > 1.

There are 15 such inequalities up to sign: choose the omitted vertex in five
ways and partition the other four into two unordered pairs in three ways.
For an arbitrary admissible equilateral basis the weak inequality >=1 holds;
strictness uses the complete ten-contact classification.

These are exact necessary restrictions on the simplex itself. They do not
replace the complete lattice admissibility test. In particular, the facet
restriction concerns which known contacts share a supporting plane; it does
not control the lengths of these additional lattice vectors.

## The remaining analytic problem

The desired estimate is a uniform lower bound greater than D* for

    |det(x_1-x_0,...,x_4-x_0)|.

The established hypotheses are equilateral distances, the facet restriction,
the 15 midpoint-separation inequalities, and full lattice admissibility.
There is presently no proof here that either of the two simpler geometric
restrictions suffices on its own, or that together they imply admissibility.
An argument using a compact closure may replace strict midpoint inequalities
by weak ones to seek a sufficient lower bound, but must include all resulting
boundary configurations rather than assume they cannot occur.

The number of unexcluded strata remains 17. The numerical value about 0.8363
in older notes is not used as a theorem or as evidence that this restricted
analytic problem has been solved.
