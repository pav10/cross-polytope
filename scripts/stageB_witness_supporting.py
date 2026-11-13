"""Gate G6s: supporting-facet B1 (STATUS §3 Theorem B.1', Lemmas B.4-B.5) keeps a known lattice's own path.

For a lattice X (frame rows) and its contact set S (frame coordinates), follow X through the real branching code
of cp4.faces in supporting mode, without symmetry reduction: the frame rows take their exact relative-interior
faces (cp4.faces._children with the trivial group); every other contact takes the half-open facet cell that
contains X for the tie side rho chosen by the search (_children_supporting).  At every step X's cell must be
among the children (it survives the side tests, the covector pre-filter and the strict-slack LP).  At the end
X must lie in the leaf's affine family (exactly, for rational X) or pass determine-then-check.
Lattices: zeta (all contacts in facet interiors) and D4 = {x in Z^4 : sum x even} with the frame
(e1+e2, e1-e2, e3+e4, e3-e4), whose non-frame contacts 2e_i and e_i +- e_j have zero coordinates, so the half-open
cells and the zero-tolerant column filter are exercised.  Each with and without ball rows and the pre-filter.
Prints SUPPORTING WITNESS OK.  ~1-2 min."""
import sys, os, json
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4 import faces as FC
from cp4.exact import solve_affine
from cp4.lattice import B_of_z, ZETA
root = os.path.join(os.path.dirname(__file__), '..')
ID = [((0, 1, 2, 3), (1, 1, 1, 1))]


def zeta():
    z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
    S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
    X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
    pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
    return 'zeta', S, 4, [[sg[k] * x for x in X0[pm[k]]] for k in range(4)]


def d4():
    X = [[F(1), F(1), F(0), F(0)], [F(1), F(-1), F(0), F(0)], [F(0), F(0), F(1), F(1)], [F(0), F(0), F(1), F(-1)]]
    # minimal vectors of D4 (l1 norm 2), one per +- pair, in frame coordinates c = v X^{-1} = v X^T / 2, over p = 2
    vs = [tuple(2 * (k == i) for k in range(4)) for i in range(4)]
    vs += [tuple((k == i) + s * (k == j) for k in range(4)) for i in range(4) for j in range(i + 1, 4) for s in (1, -1)]
    S = []
    for v in vs:
        c = [sum(v[k] * X[i][k] for k in range(4)) for i in range(4)]       # = 2 * (v X^T / 2) = numerator over 2
        S.append(tuple(int(x) for x in c))
    return 'D4', S, 2, X


def sign(t, exact):
    if exact:
        return (t > 0) - (t < 0)
    return int(np.sign(round(float(t), 9)))


def walk(name, S, p, X, ball, prefilter):
    exact = isinstance(X[0][0], F)
    order = FC._frame_order(S, p)
    C = [tuple(F(x, p) for x in S[i]) for i in order]
    stats = {'lp_dead': 0, 'prefilter_dead': 0}
    extra = FC.ball_rows(C) if ball else None
    colf = FC._column_filter(C, stats) if prefilter else None
    vec = lambda c: [sum(c[i] * X[i][j] for i in range(4)) for j in range(4)]
    assert all(abs(sum(abs(t) for t in vec(c)) - 2) < 1e-9 for c in C), 'S must consist of contacts of X'
    eqs, st, no, signs, rhos = [], [], [], {}, {}
    for ci, c in enumerate(C):
        v = vec(c)
        if ci < 4:
            res = FC._children(C, ID, eqs, st, signs, False, FC.SIGNVECS, stats, colf, extra=extra)
            assert res is not None and res[0] == ci
            want = tuple(sign(t, exact) for t in v)
            kid = next((k for k in res[1] if k[0] == want), None)
            assert kid is not None, (name, 'frame row', ci, want)
            eqs, st = kid[2], kid[3]
        else:
            rho, kids = FC._children_supporting(C, ci, eqs, st, no, signs, rhos, stats, colf, extra)
            assert rho is not None, (name, 'dead node before contact', ci)
            want = tuple(sign(t, exact) or rho[j] for j, t in enumerate(v))
            kid = next((k for k in kids if k[0] == want), None)
            assert kid is not None, (name, 'contact', ci, c, 'cell', want, 'rho', rho, 'children', [k[0] for k in kids])
            eqs, st, no = kid[2], kid[3], kid[4]
            rhos[ci] = rho
        signs[ci] = want
    sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
    assert sol is not None
    x0, N = sol
    flat = [X[i][j] for i in range(4) for j in range(4)]
    if exact:   # X - x0 in the span of N
        A = [list(col) for col in zip(*N)] if N else [[F(0)] for _ in range(16)]
        assert solve_affine([[F(0)] * 0 + row for row in A], [flat[k] - x0[k] for k in range(16)]) is not None
    else:
        M = np.array(N, float).T if N else np.zeros((16, 1))
        r = np.array(flat, float) - np.array(x0, float)
        assert np.linalg.norm(M @ np.linalg.lstsq(M, r, rcond=None)[0] - r) < 1e-8
    zeros = sum(1 for ci in range(4, len(C)) for t in vec(C[ci]) if abs(float(t)) < 1e-12)
    print(f'{name}: ball={ball} prefilter={prefilter}: path survives ({len(C)} contacts, {zeros} zero coordinates '
          f'on non-frame contacts, leaf dim {len(N)}); stats {stats}', flush=True)


for lat in (zeta(), d4()):
    for ball in (False, True):
        for pf in (False, True):
            walk(*lat, ball=ball, prefilter=pf)
print('SUPPORTING WITNESS OK')
