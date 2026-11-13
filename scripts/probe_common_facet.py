"""[Num] Common-facet strata (4, 17, 30): Gamma from data/strata.json, rows 2-4 kept >= 0 (sigma = 1, forced by
the shell vector (0,1/3,1/3,1/3)), all rows of l1 norm 2.  Minimises det subject to lambda_1 >= 2 (penalty,
Nelder-Mead from starts near 2I).  Exploration only.
Usage: python3 scripts/probe_common_facet.py <stratum> [seed] [restarts]"""
import sys, json, numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
from cp4.lattice import fast_lambda1
from scipy.optimize import minimize
sid = int(sys.argv[1]); rng = np.random.default_rng(int(sys.argv[2])); restarts = int(sys.argv[3])
st = json.load(open('data/strata.json'))[sid]; p = st['p']
def echelon_basis(rows, n=4):
    M = [list(r) for r in rows]; basis = []
    for col in range(n):
        while True:
            nz = [r for r in M if r[col] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(r[col])); q0 = nz[0]
            for r in nz[1:]:
                q = r[col] // q0[col]
                for k in range(n): r[k] -= q * q0[k]
            M = [r for r in M if any(r)]
        piv = [r for r in M if r[col] != 0]; basis.append(piv[0]); M = [r for r in M if r is not piv[0]]
    return basis
B = np.array(echelon_basis([[p*int(i==j) for j in range(4)] for i in range(4)] + [list(r) for r in st['residues'] if any(r)]), float)/p
def frame(v):
    X = v.reshape(4, 4).copy(); X[1:] = np.abs(X[1:])
    return 2*X/np.abs(X).sum(axis=1, keepdims=True)
def obj(v):
    Y = B @ frame(v); d = abs(np.linalg.det(Y))
    if d < 1e-4: return 10.0
    l = fast_lambda1(Y, cap=20)
    return 10.0 if l is None else d + 30*max(0.0, 2-l)
best = 9
for it in range(restarts):
    v0 = (2*np.eye(4) + 0.6*rng.normal(size=(4,4))).ravel()
    r = minimize(obj, v0, method='Nelder-Mead', options={'maxiter': 12000, 'xatol': 1e-10, 'fatol': 1e-12})
    r = minimize(obj, r.x, method='Nelder-Mead', options={'maxiter': 12000, 'xatol': 1e-11, 'fatol': 1e-13})
    if r.fun < best:
        best = r.fun; Y = B @ frame(r.x)
        print(f"stratum {sid} it {it}: det {abs(np.linalg.det(Y)):.5f} lambda1 {fast_lambda1(Y, cap=20):.5f}", flush=True)
print(f"stratum {sid} done best {best:.5f}")
