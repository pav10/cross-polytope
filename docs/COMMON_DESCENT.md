# A common descent for all equilateral spanning-tree bases

9 November 2026. The result is theoretical [P]. It supplies a sufficient finite
admissibility test, but no new determinant bound or stratum exclusion.

## Statement

Let x_0,...,x_4 be affinely independent points in R^4 with all pairwise l1
distances equal to 2. Work in the coefficient lattice

    L = {b in Z^5 : sum_i b_i=0},
    F(b) = ||sum_i b_i x_i||_1.

Its roots r_ij=e_i-e_j satisfy F(r_ij)=2. The roots corresponding to any
oriented spanning tree on the five vertices form a Z-basis R of L. This follows
by removing leaves successively, or expressing every root by its tree path.

For each such basis, transport the finite test set of THEORY 1.7 from Z^4 to L;
call the resulting symmetric set T_R. Define

    J = intersection over all spanning-tree bases R of T_R.

**Theorem.** F(b)>=2 for every nonzero b in L if and only if F(b)>=2 for every
b in J. In particular J is finite and is contained in the test set for any
one chosen basis.

## Proof

Use the same auxiliary norm for every basis:

    mu(b) = (1/2) sum_i |b_i| = sum_{b_i>0} b_i.

It takes positive integer values on L minus zero. The vectors with mu=1 are
exactly the roots, and every basis vector of every tree basis has mu=1.

Suppose F has a nonzero violating vector. Choose one, b, with F(b)<2 and
mu(b) minimal. Since all roots have F=2, mu(b)>1.

Fix an arbitrary tree basis R=(r_1,...,r_4) and an index i for which b has
nonzero r_i-coordinate. Let

    V_i(b) = conv(±r_j for j!=i, ±b).

If a lattice point w in V_i(b), distinct from ±b, has nonzero b-coefficient,
write it in the unique basis (r_j for j!=i, b) as

    w = sum_{j!=i} alpha_j r_j + beta b,
    sum |alpha_j| + |beta| <= 1.

Here 0<|beta|<1; equality |beta|=1 would give w=±b. Also w is nonzero because
its r_i-coordinate is beta times the nonzero r_i-coordinate of b. The two
triangle inequalities give

    F(w) <= 2 sum |alpha_j| + |beta| F(b) < 2,
    mu(w) <= 1-|beta| + |beta| mu(b) < mu(b).

This contradicts minimality. Thus every V_i(b) has exactly the irreducibility
property required in THEORY 1.7. The coordinate l1 norm of b in the tree basis
is greater than 1: a nonzero integer coordinate vector of norm at most 1 would
be a signed basis vector, hence a root. Therefore b belongs to T_R.

The basis R was arbitrary, so b belongs to J. If F>=2 on J this is impossible,
proving sufficiency. Necessity is immediate. Finiteness follows from
J subset T_R for any fixed R, using THEORY 1.7. QED.

The key point is the common norm mu. Applying the single-basis theorem
separately would only produce possibly different witnesses in different test
sets; that alone would not justify taking their intersection.

## The first nontrivial affine mass

At mu=2, after identifying b and -b, vectors have three possible forms:

1. 2(e_i-e_j), whose F-value is 4.
2. 2e_i-e_j-e_k with distinct indices. Put u=x_i-x_j and v=x_i-x_k.
   Since ||u||=||v||=||u-v||=2, the triangle inequality applied to
   2u=(u+v)+(u-v) gives ||u+v||>=2. These tests hold automatically.
3. e_i+e_j-e_k-e_l with four distinct indices. These are exactly the 15
   disjoint-edge midpoint tests in INDEX_ONE_GEOMETRY.md.

Thus only those 15 inequalities can fail at affine mass two. This does NOT
prove that mass-two tests suffice: potential witnesses of greater mass must
still be eliminated or checked in J.

At a local minimum having exactly the ten A4 contact pairs, the inequalities
in items 2 and 3 are strict, because neither vector type is a root. Automatic
weak admissibility for item 2 must not be confused with automatic strictness.

## Consequence for the theoretical program

The next exact question is the structure of J, rather than independent
admissibility analyses for each choice of basis. A classification or analytic
bound for J could replace much of the original single-basis test. Its size,
and whether it contains essential vectors of affine mass greater than two,
have not been determined here. No claimed reduction from 1156 to 15 is justified.

This result is independent of local optimality and of the facet obstruction;
it applies to any affinely independent equilateral five-point configuration.
It leaves the desired determinant lower bound and the count of 17 remaining
strata unchanged.
