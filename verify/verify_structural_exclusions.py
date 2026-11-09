"""Exact checks of the data-dependent claims in docs/BLOCK_COUPLING.md and docs/FACET_LAYERS.md.

Run from any directory: python3 verify/verify_structural_exclusions.py   (a few seconds; prints STRUCTURAL PASS)
The proofs themselves are in the notes; this script checks every finite fact they use:
  * the cyclic generators of strata 12, 13, 16, 18, 19, 28, 31 generate the recorded residue groups;
  * the stated shell vectors and lattice memberships of strata 28 and 31, and the identities aX, bX;
  * the planar identity of BLOCK_COUPLING (1)-(2), symbolically;
  * vol_3(K cap {sum x = 0}) = 20/3 for K = {||x||_1 <= 2}, by an exact triangulation (scipy only lists facets);
  * D* < D(13/25) = 8089/10000 < 81/100 < 5/6, from the sign of z^3 - 6z + 3 at 13/25.
Residue numerators are stored canonically in [0, p) (checked).  Exact arithmetic only; no LP."""
import json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
strata = {s['id']: s for s in json.loads((ROOT / 'data/strata.json').read_text())}
assert all(0 <= x < s['p'] for s in strata.values() for r in s['residues'] for x in r)


def residues(sid):
    return {tuple(r) for r in strata[sid]['residues']}


def cyclic(g, p):
    return {tuple((k * x) % p for x in g) for k in range(p)}


def shell(sid):
    return {tuple(f) for f in strata[sid]['forced']}


def check():
    # FACET_LAYERS §3 (and strata 28, 31): positive cyclic generators with numerators summing to p.
    for sid, g in [(12, (1, 1, 1, 1)), (13, (1, 1, 1, 2)), (16, (2, 2, 1, 1)), (18, (3, 1, 1, 1)), (19, (1, 1, 2, 3)),
                   (28, (1, 1, 3, 3)), (31, (2, 4, 1, 3))]:
        p = strata[sid]['p']
        assert cyclic(g, p) == residues(sid), sid
        assert all(x > 0 for x in g) and sum(g) == p, sid
        print(f"stratum {sid:2d}: Gamma = Z^4 + Z{g}/{p}, positive, sum = p  -> det >= 5/6 (FACET_LAYERS)")
    # BLOCK_COUPLING §2 (stratum 28): g, h, j in the unit shell; g - e3, g - e4, j - e1, j - e2 in Gamma.
    sh = shell(28)
    assert {(1, 1, 3, 3), (2, 2, -2, -2), (3, 3, 1, 1)} <= sh
    # BLOCK_COUPLING §3 (stratum 31): u, v in the shell; h = (0,0,5,5)/10, (2,4,6,8)/10, (4,8,2,6)/10 in Gamma.
    sh = shell(31)
    assert (2, 4, 1, 3) in sh and ((4, -2, -3, 1) in sh or (-4, 2, 3, -1) in sh)
    assert {(0, 0, 5, 5), (2, 4, 6, 8), (4, 8, 2, 6)} <= residues(31)
    A, B, C, D = sp.symbols('A B C D')        # rows X1, X4 (block 1) and X2, X3 (block 2)
    X = {1: A, 2: C, 3: D, 4: B}
    a = [F(2, 10), F(4, 10), F(-4, 10), F(-2, 10)]
    b = [F(4, 10), F(-2, 10), F(2, 10), F(-4, 10)]
    comb = lambda c: sp.expand(sum(sp.Rational(c[i].numerator, c[i].denominator) * X[i + 1] for i in range(4)))
    assert comb(a) == sp.expand((A - B) / 5 + 2 * (C - D) / 5)
    assert comb(b) == sp.expand(2 * (A - B) / 5 - (C - D) / 5)
    # BLOCK_COUPLING (1)-(2), symbolically.
    s, t = sp.symbols('s t', real=True)
    Av, Bv = sp.Matrix([1 + s / 2, 1 - s / 2]), sp.Matrix([1 + t / 2, 1 - t / 2])
    assert sp.expand(Av[0] * Bv[1] - Av[1] * Bv[0]) == s - t
    m, d = sp.symbols('m d', real=True)
    for v, diff in (((-5 * Av + 3 * Bv) / 8, -m / 4 - d), ((3 * Av - 5 * Bv) / 8, -m / 4 + d)):
        assert sp.simplify(v[0] + v[1]) == sp.Rational(-1, 2)
        assert sp.expand((v[0] - v[1]).subs({s: m + d, t: m - d})) == sp.expand(diff)
    # FACET_LAYERS §1: the section is conv{e_i - e_j} (a cuboctahedron); exact volume via facets.
    from scipy.spatial import ConvexHull
    import numpy as np
    pts = []
    for i in range(4):
        for j in range(4):
            if i != j:
                v = [0] * 4; v[i] += 1; v[j] -= 1
                pts.append(v[:3])                 # coordinates (x1, x2, x3) on H, x4 = -x1-x2-x3
    hull = ConvexHull(np.array(pts, float))
    vol = F(0)
    for simplex in hull.simplices:                # facet triangles; cone from the origin (an interior point)
        M = [[F(x) for x in pts[k]] for k in simplex]
        dt = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
              + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
        vol += abs(dt) / 6
    assert len(pts) == 12 and len(hull.vertices) == 12
    for v in pts:                                  # every vertex lies on the boundary ||x||_1 = 2 of K
        assert sum(abs(x) for x in v) + abs(sum(v)) == 2
    assert 2 * vol == F(20, 3), vol                # Euclidean volume element on H is 2 dx1 dx2 dx3
    # D* < D(13/25) < 81/100 < 5/6: D is decreasing on (1/2, 3/4) and 13/25 < zeta.
    z = sp.Rational(13, 25)
    assert z ** 3 - 6 * z + 3 > 0 and sp.Rational(1, 2) < z            # p(1/2) > 0 > p(11/20): zeta in (13/25, 11/20)
    assert sp.Rational(11, 20) ** 3 - 6 * sp.Rational(11, 20) + 3 < 0
    Dz = (16 - 9 * z + 6 * z ** 2) / 16
    assert Dz == sp.Rational(8089, 10000) and Dz < sp.Rational(81, 100) < sp.Rational(5, 6)
    print("BLOCK_COUPLING identities, section volume 20/3, D* < 8089/10000 < 81/100: checked")
    print("STRUCTURAL PASS")


if __name__ == '__main__':
    check()
