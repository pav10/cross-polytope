"""Exact checks of the data-dependent claims in docs/BLOCK_COUPLING.md, docs/FACET_LAYERS.md and docs/COMMON_FACET.md.

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
    # COMMON_FACET §3a: the extremal stratum-30 lattice of the axis case is admissible, with det 8/9.
    import sys
    sys.path.insert(0, str(ROOT / 'src'))
    from cp4.lattice import exact_short_vectors
    from cp4.exact import det as exact_det
    p30 = strata[30]['p']
    gens = [[p30 * int(i == j) for j in range(4)] for i in range(4)] + [list(r) for r in strata[30]['residues'] if any(r)]
    M = [list(r) for r in gens]; basis = []                    # Z-basis of p*Gamma by gcd row reduction
    for col in range(4):
        while True:
            nz = [r for r in M if r[col] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(r[col])); q0 = nz[0]
            for r in nz[1:]:
                q = r[col] // q0[col]
                for k in range(4): r[k] -= q * q0[k]
            M = [r for r in M if any(r)]
        piv = [r for r in M if r[col] != 0][0]; basis.append(piv); M = [r for r in M if r is not piv]
    Bg = [[F(x, p30) for x in r] for r in basis]
    assert abs(exact_det(Bg)) == F(1, p30)
    t = F(1, 3)
    X = [[t, F(-1), t, t], [0, 0, 2, 0], [2, 0, 0, 0], [0, 0, 0, 2]]
    Y = [[sum(Bg[i][k] * X[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
    sv = exact_short_vectors(Y, 2, cap=60)
    assert min(n for _, n in sv) == 2 and sum(1 for _, n in sv if n == 2) == 16
    assert abs(exact_det(Y)) == F(8, 9)
    print("stratum 30 axis-case extremal lattice: admissible (lambda_1 = 2, 16 minimal pairs), det = 8/9: checked")
    check_main_configuration()
    print("STRUCTURAL PASS")


def check_main_configuration():
    """COMMON_FACET §3b, symbolically: the six short vectors, the budget identities, Schur/Cramer, the column trick,
    the expansion of det M', and the Laplacian determinant formula."""
    # the six vectors: v_s (norm 10/9) and v_s - e_j (norm 11/9), all in Gamma_30
    R30 = residues(30)
    for (i, j, l) in [(1, 2, 3), (3, 1, 2), (2, 3, 1)]:
        v = [F(1, 3), 0, 0, 0]; v[i] += F(1, 9); v[j] += F(4, 9); v[l] -= F(2, 9)
        w = list(v); w[j] -= 1
        for c, n in ((v, F(10, 9)), (w, F(11, 9))):
            assert sum(abs(x) for x in c) == n and tuple(int((9 * x) % 9) for x in c) in R30
    e = {(a, b): sp.Symbol(f'e{a}{b}') for a in range(3) for b in range(3) if a != b}
    f = [sp.Symbol(f'f{a}') for a in range(3)]
    ell = [f[a] + sum(e[(a, b)] for b in range(3) if b != a) for a in range(3)]
    d = [2 - ell[a] for a in range(3)]
    M = sp.Matrix(3, 3, lambda a, b: d[a] if a == b else e[(a, b)])
    P = M.row_join(sp.Matrix(f))                  # rows X2, X3, X4 on columns (k2, k3, k4, t); row sums 2
    SH = [(0, 1, 2), (2, 0, 1), (1, 2, 0)]        # (i, j, l) as row indices of P
    L = sum(ell); E = sum(ell) - sum(f); Fs = sum(f)
    E3 = e[(1, 2)] + e[(0, 1)] + e[(2, 0)]; E6 = e[(0, 2)] + e[(2, 1)] + e[(1, 0)]
    q = m = qm = mm = 0
    for (i, j, l) in SH:
        q += 2 * P[l, l] - P[i, l] - 4 * P[j, l] - 3
        S = lambda r: P[r, j] + P[r, l]
        m += 5 * S(j) + 2 * S(l) - S(i) - 12
        qm += 2 * (P[l, l] + P[l, 3]) - (P[i, l] + P[i, 3]) - 4 * (P[j, l] + P[j, 3]) - 3
        mm += P[i, i] - 5 * P[j, i] - 2 * P[l, i]
    assert sp.expand(q + m - 9 * (1 - L)) == 0
    assert sp.expand(qm + mm - (9 - 9 * E - 6 * Fs)) == 0
    assert sp.expand(2 * qm - mm - (-9 * E3 - 9 * Fs)) == 0 and sp.expand(2 * mm - qm - (9 - 9 * E6 + 3 * Fs)) == 0
    y = sp.symbols('y0:3'); xt = sp.Symbol('xt')
    X = sp.Matrix([[y[0], y[1], y[2], xt]]).col_join(P)
    C = [sp.Matrix.hstack(*[sp.Matrix(f) if k == r else M[:, k] for k in range(3)]).det() for r in range(3)]
    assert sp.expand(X.det() + xt * M.det() - sum(y[r] * C[r] for r in range(3))) == 0
    W = M - sp.ones(3, 1) * sp.Matrix([y])
    assert sp.expand(X.det().subs(xt, 2 - sum(y)) + 2 * W.det()) == 0
    r, b, c = 0, 1, 2
    Mp = M.copy(); Mp[:, r] = M[:, r] - (1 - L) * sp.ones(3, 1)
    u = ell[b] + ell[c]; al = 1 - L - e[(b, r)]; be = 1 - L - e[(c, r)]
    assert sp.expand(Mp[r, r] - (1 + u)) == 0
    assert sp.expand(Mp.det() - ((1 + u) * (d[b] * d[c] - e[(b, c)] * e[(c, b)]) + e[(r, b)] * (al * d[c] - be * e[(b, c)])
                                 + e[(r, c)] * (be * d[b] - al * e[(c, b)]))) == 0
    w = {(a, bb): sp.Symbol(f'w{a}{bb}') for a in range(3) for bb in range(3) if a != bb}
    Lap = sp.Matrix(3, 3, lambda a, bb: sum(w[(a, k)] for k in range(3) if k != a) if a == bb else -w[(a, bb)])
    s2 = sum(Lap.extract([p_, q_], [p_, q_]).det() for p_ in range(3) for q_ in range(p_ + 1, 3))
    assert sp.expand((2 * sp.eye(3) - Lap).det() - (8 - 4 * Lap.trace() + 2 * s2)) == 0
    print("COMMON_FACET §3b: six vectors in Gamma_30, budget identities, Schur/Cramer, column trick, det M', "
          "Laplacian formula: checked")


if __name__ == '__main__':
    check()
