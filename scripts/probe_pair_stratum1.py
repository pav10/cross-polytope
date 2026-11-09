"""[Num] Pair problem in stratum 1's form (DEAD_ENDS B17): Lambda = span((u+v)/2, (u-v)/2, X1, X2) with u, v of
l1 norm 2 and disjoint supports (pattern 22, 13 or 12), X1, X2 of l1 norm 2.  Minimises det subject to
lambda_1 >= 2 (penalty, Nelder-Mead).  Exploration only.
Usage: python3 scripts/probe_pair_stratum1.py <22|13|12> [seed] [restarts]"""
import sys, numpy as np
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
from cp4.lattice import fast_lambda1
from scipy.optimize import minimize
pat = sys.argv[1]; STRAT1 = True
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
supp = {'22': ([0, 1], [2, 3]), '13': ([0], [1, 2, 3]), '12': ([0], [1, 2])}[pat]
nu, nv = len(supp[0]), len(supp[1])
def vec(params, idx):
    x = np.zeros(4); w = np.abs(params) + 1e-9; x[idx] = 2 * w / w.sum(); return x
def build(p):
    u = vec(p[:nu], supp[0]); v = vec(p[nu:nu+nv], supp[1])
    a = p[nu+nv:nu+nv+4]; b = p[nu+nv+4:nu+nv+8]
    if STRAT1: a = 2*a/np.abs(a).sum(); b = 2*b/np.abs(b).sum()
    return np.array([(u+v)/2, (u-v)/2, a, b])
def obj(p):
    Y = build(p); d = abs(np.linalg.det(Y))
    if d < 1e-4: return 10.0
    l = fast_lambda1(Y, cap=20)
    return 10.0 if l is None else d + 30*max(0.0, 2-l)
best = 9
for it in range(int(sys.argv[3]) if len(sys.argv) > 3 else 150):
    p0 = np.r_[rng.random(nu+nv), 3*rng.normal(size=8)]
    r = minimize(obj, p0, method='Nelder-Mead', options={'maxiter': 8000, 'xatol': 1e-10, 'fatol': 1e-12})
    r = minimize(obj, r.x, method='Nelder-Mead', options={'maxiter': 8000, 'xatol': 1e-11, 'fatol': 1e-13})
    if r.fun < best:
        best = r.fun; Y = build(r.x)
        print(f"{pat} it {it}: det {abs(np.linalg.det(Y)):.5f} lambda1 {fast_lambda1(Y, cap=20):.5f}", flush=True)
print(pat, 'done, best objective', best)
