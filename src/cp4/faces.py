"""Stage B (part 1): face data for a fixed contact set S, and the resulting affine families.

Unknown: the frame matrix X (rows = the 4 frame minimal vectors in ambient coordinates),
variables x[4*i+j] = X[i][j].  For a contact c (frame coordinates) put (cX)_j = sum_i c_i X[i][j].

Face datum = for each contact c a sign vector s_c in {-1,0,1}^4 minus {0} (the face of the l1 sphere
of radius 2 containing cX in its RELATIVE INTERIOR):
    sum_j s_cj (cX)_j = 2,      s_cj (cX)_j > 0  (s_cj != 0),      (cX)_j = 0  (s_cj = 0).
A critical lattice lies in the relative interior of exactly one face datum (for its own S), and
is then a critical point of det on the affine set A cut out by the equalities  [STATUS.md §B].

Branching: contacts in a fixed order; for each, the admissible sign vectors are reduced modulo the
subgroup of the ambient signed-permutation group (acting on COLUMNS of X) that fixes all signs
chosen so far (stabilizer chain => complete and non-redundant).  Pruning: LP feasibility of the
strict system (maximise a common slack t; feasible iff t > 0).  Leaves with dim A > 6 are
discarded (perfection, Thm 5.1).
"""
import numpy as np
from itertools import product, permutations
from fractions import Fraction as F
from scipy.optimize import linprog
from .exact import solve_affine, rank

AMB = [(pm, sg) for pm in permutations(range(4)) for sg in product((1, -1), repeat=4)]
SIGNVECS = [s for s in product((-1, 0, 1), repeat=4) if any(s)]

def act_sign(g, s):
    """Column action: X -> X Q with Q signed permutation; the sign vector transforms the same way."""
    pm, sg = g
    return tuple(sg[k] * s[pm[k]] for k in range(4))

def orbit_reps(group, choices):
    seen = set(); reps = []
    for s in choices:
        if s in seen:
            continue
        reps.append(s)
        for g in group:
            seen.add(act_sign(g, s))
    return reps

def _rows_for(c, s):
    """Equalities and strict inequalities for contact c (Fractions) with sign vector s."""
    eq = []; strict = []
    # (cX)_j as a 16-vector
    cols = []
    for j in range(4):
        r = [F(0)] * 16
        for i in range(4):
            r[4 * i + j] = F(c[i])
        cols.append(r)
    eq.append(([sum(s[j] * cols[j][k] for j in range(4)) for k in range(16)], F(2)))
    for j in range(4):
        if s[j] == 0:
            eq.append((cols[j], F(0)))
        else:
            strict.append([s[j] * x for x in cols[j]])
    return eq, strict

def lp_feasible(eqs, stricts):
    """Max t s.t. A_eq x = b, S x >= t, |x_k| <= 2 (frame rows have l1 norm 2), t <= 1."""
    n = 16
    A_eq = np.array([[float(v) for v in r] + [0.0] for r, _ in eqs]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub = np.array([[-float(v) for v in r] + [1.0] for r in stricts]) if stricts else None
    b_ub = np.zeros(len(stricts)) if stricts else None
    c = np.zeros(n + 1); c[-1] = -1.0
    res = linprog(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(-2, 2)] * n + [(None, 1.0)], method='highs')
    if res.status != 0:
        return False, None
    return res.x[-1] > 1e-9, res.x[:n]

def _closure_ranges(eqs, stricts, row):
    """min and max of row.x over the closure of the current polyhedron (|x|<=2)."""
    n = 16
    A_eq = np.array([[float(v) for v in r] for r, _ in eqs]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub = np.array([[-float(v) for v in r] for r in stricts]) if stricts else None
    b_ub = np.zeros(len(stricts)) if stricts else None
    out = []
    for sgn in (1.0, -1.0):
        res = linprog(sgn * np.array([float(v) for v in row]), A_ub=A_ub, b_ub=b_ub,
                      A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * n, method='highs')
        if res.status != 0:
            return None
        out.append(sgn * res.fun)
    return out[0], out[1]          # (min, max)

def allowed_signs(c, eqs, stricts, eps=1e-9):
    """Sign vectors for contact c compatible (coordinatewise) with the current polyhedron.
    Conservative: uses the closure, so never discards a feasible sign vector."""
    per = []
    for j in range(4):
        row = [F(0)] * 16
        for i in range(4):
            row[4 * i + j] = F(c[i])
        r = _closure_ranges(eqs, stricts, row)
        if r is None:
            return []
        lo, hi = r
        a = []
        if lo < -eps: a.append(-1)
        if lo <= eps and hi >= -eps: a.append(0)
        if hi > eps: a.append(1)
        per.append(a)
    return [s for s in product(*per) if any(s)]

def enumerate_face_data(S, p, max_dim=6, log_every=None, time_limit=None, dynamic=True):
    """S: list of contacts as numerator tuples over p.  Returns (leaves, stats).
    leaves: dict(signs, x0, basis (affine family), dim).
    Dynamic order: at each node branch on the undecided contact with the fewest
    coordinatewise-allowed sign vectors (computed by closure LPs; conservative)."""
    import time
    t0 = time.time()
    C = [tuple(F(x, p) for x in c) for c in S]
    leaves = []
    stats = {'nodes': 0, 'lp_dead': 0, 'dim_dead': 0, 'leaves': 0, 'per_depth': {}}

    def rec(group, eqs, stricts, signs):
        stats['nodes'] += 1
        d = len(signs)
        stats['per_depth'][d] = stats['per_depth'].get(d, 0) + 1
        if log_every and stats['nodes'] % log_every == 0:
            print(f"  [heartbeat] {time.time()-t0:.0f}s", stats, flush=True)
        if time_limit and time.time() - t0 > time_limit:
            raise TimeoutError
        undecided = [i for i in range(len(C)) if i not in signs]
        if not undecided:
            sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
            if sol is None:
                return
            x0, N = sol
            if len(N) > max_dim:
                stats['dim_dead'] += 1; return
            stats['leaves'] += 1
            leaves.append({'signs': dict(signs), 'x0': x0, 'basis': N, 'dim': len(N)})
            return
        if dynamic:
            best = None
            for i in undecided:
                al = allowed_signs(C[i], eqs, stricts)
                if best is None or len(al) < len(best[1]):
                    best = (i, al)
                if len(al) <= 1:
                    break
            ci, choices = best
        else:
            ci = undecided[0]; choices = SIGNVECS
        # orbit reduction must be over the FULL set of sign vectors, restricted afterwards:
        reps = [s for s in orbit_reps(group, SIGNVECS)]
        chs = set(choices)
        for s in reps:
            # some image of s under group must be allowed; use that image (same orbit)
            img = next((act_sign(g, s) for g in group if act_sign(g, s) in chs), None)
            if img is None:
                continue
            e, st = _rows_for(C[ci], img)
            ok, _ = lp_feasible(eqs + e, stricts + st)
            if not ok:
                stats['lp_dead'] += 1; continue
            sub = [g for g in group if act_sign(g, img) == img]
            signs[ci] = img
            rec(sub, eqs + e, stricts + st, signs)
            del signs[ci]

    try:
        rec(AMB, [], [], {})
        stats['complete'] = True
    except TimeoutError:
        stats['complete'] = False
    stats['seconds'] = time.time() - t0
    return leaves, stats
