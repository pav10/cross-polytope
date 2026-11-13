"""Exact checks for docs/STRATUM11.md §3: Conjecture 11.A is false.

The frame X below (rows = frame vectors, ambient coordinates) satisfies
  * ||X_i||_1 = 2 for the four frame rows,
  * ||cX||_1 >= 2 for the 19 vectors of the 7-point system (STRATUM11 §1) AND for all 68 candidates of
    stratum 11 (every Gamma_11 point in [-1,1]^4 with ||c||_1 >= 1),
and has |det X| / 4 = 115/144 < D*.  It is not admissible: Gamma_11 contains c = (1/2, 3/4, -1/4, -5/4)
(outside the cube) with ||cX||_1 = 4/3.  So no finite relaxation by vectors of the cube [-1,1]^4 can prove
det >= D* on stratum 11.  Run from any directory: python3 verify/verify_stratum11.py
(prints STRATUM11 PASS).
"""
import json, os, sys
from fractions import Fraction as F
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
from cp4.exact import det
from cp4.lattice import DSTAR

X = [[F(-23, 12), F(1, 12), F(0), F(0)],
     [F(1, 2), F(-1, 4), F(-5, 4), F(0)],
     [F(0), F(2), F(0), F(0)],
     [F(-1, 3), F(-1, 4), F(-3, 4), F(-2, 3)]]


def norm(c):
    return sum(abs(sum(F(c[i]) * X[i][j] for i in range(4))) for j in range(4))


def in_gamma11(c):
    """c in Z^4 + Z(1/2, 1/4, 1/4, 1/4)?"""
    g = (F(1, 2), F(1, 4), F(1, 4), F(1, 4))
    return any(all((F(c[i]) - k * g[i]).denominator == 1 for i in range(4)) for k in range(4))


def check():
    assert all(sum(abs(x) for x in row) == 2 for row in X)
    g, b = (F(1, 2), F(1, 4), F(1, 4), F(1, 4)), (F(-1, 2), F(1, 4), F(1, 4), F(1, 4))
    e = [tuple(F(int(i == j)) for j in range(4)) for i in range(4)]
    add = lambda u, v, s=1: tuple(x + s * y for x, y in zip(u, v))
    cc = add(g, b)
    seven = [g, b, cc] + [add(cc, e[r], -1) for r in (1, 2, 3)] + [add(g, e[r], -1) for r in (1, 2, 3)] \
        + [add(b, e[r], -1) for r in (1, 2, 3)] + [add(e[r], e[s], -1) for r, s in ((1, 2), (1, 3), (2, 3))] + e
    assert len(seven) == 19 and all(norm(c) >= 2 for c in seven)
    st = json.load(open(os.path.join(ROOT, 'data', 'strata.json')))[11]
    cand = [tuple(F(x, st['p']) for x in c) for c in st['candidates']]
    assert len(cand) == 68 and all(norm(c) >= 2 for c in cand)
    d = abs(det(X)) / 4
    assert d == F(115, 144) and d < F(DSTAR) - F(1, 10 ** 6)
    short = (F(1, 2), F(3, 4), F(-1, 4), F(-5, 4))
    assert in_gamma11(short) and norm(short) == F(4, 3)
    print('Conjecture 11.A counterexample: det = 115/144 < D*, all 19 + 68 constraints hold, '
          'short vector (1/2,3/4,-1/4,-5/4) of norm 4/3 (not admissible)')
    print('STRATUM11 PASS')


if __name__ == '__main__':
    check()
