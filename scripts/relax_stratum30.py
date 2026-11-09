"""[Num] Relaxation for stratum 30: admissibility imposed only for c in Gamma_30 with ||c||_1 <= NMAX.
Minimise |det X|/9 over rows 2-4 >= 0 (sum 2), ||X1||_1 = 2.  If the relaxed minimum stays > 0.81, those finitely
many constraints already exclude the stratum.  Usage: python3 scripts/relax_stratum30.py [NMAX] [seed] [restarts]
(starts near 2I; docs/COMMON_FACET.md §3a records random-start results too)."""
import sys, json, itertools, numpy as np
from fractions import Fraction as F
from scipy.optimize import minimize
NMAX = F(sys.argv[1]) if len(sys.argv) > 1 else F(4, 3)
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
import os
st = json.load(open(os.path.join(os.path.dirname(__file__), '..', 'data', 'strata.json')))[30]; p = 9
R = [tuple(r) for r in st['residues']]
vecs = set()
for r in R:
    for sh in itertools.product((-2, -1, 0, 1), repeat=4):
        c = tuple(F(r[i], p) + sh[i] for i in range(4))
        if any(c) and sum(abs(x) for x in c) <= NMAX:
            if tuple(-x for x in c) not in vecs: vecs.add(c)
V = np.array([[float(x) for x in c] for c in sorted(vecs)])
print(len(V), 'constraint vectors with ||c||_1 <=', NMAX, flush=True)
def frame(v):
    X = v.reshape(4, 4).copy(); X[1:] = np.abs(X[1:]); return 2 * X / np.abs(X).sum(axis=1, keepdims=True)
def obj(v):
    X = frame(v); d = abs(np.linalg.det(X)) / 9
    viol = np.maximum(0, 2 - np.abs(V @ X).sum(axis=1)).sum()
    return d + 50 * viol
best = 9
for it in range(int(sys.argv[3]) if len(sys.argv) > 3 else 60):
    v = (2 * np.eye(4) + 0.5 * rng.normal(size=(4, 4))).ravel()
    for tol in (1e-9, 1e-12):
        r = minimize(obj, v, method='Nelder-Mead', options={'maxiter': 30000, 'xatol': tol, 'fatol': tol}); v = r.x
    if r.fun < best:
        best = r.fun; X = frame(v)
        print(f"it {it}: relaxed det {abs(np.linalg.det(X))/9:.5f} violation {np.maximum(0, 2-np.abs(V@X).sum(axis=1)).sum():.2e}", flush=True)
np.set_printoptions(precision=4, suppress=True); print(frame(v))
