"""[Num] Stratum 7, residual support pattern of docs/COORDINATE_SECTIONS.md §3:
X2 = 2e_a, X3 = 2e_b, X4 = (0,0,yc,yd), X1 free with ||X1||_1 = 2.  Minimises det(Gamma X) subject to
lambda_1 >= 2 (penalty, Nelder-Mead from random starts, fast short vectors).  Exploration only, not a proof:
the exclusion is the certified Lemma 3 (verify/verify_coordinate_sections.py).
Usage: python3 scripts/probe_stratum7_residual.py [seed]   (~15 min; prints each improvement)"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.lattice import fast_lambda1
from scipy.optimize import minimize
G = np.array([[1,0,0,0],[0,.5,.25,.25],[0,0,.5,.5],[0,0,0,1.]])   # basis of Gamma (stratum 7)
def build(p):
    yc = 2/(1+np.exp(-p[0]))              # in (0,2)
    x = p[1:5]; x = 2*x/np.abs(x).sum()
    X = np.array([x, [2,0,0,0], [0,2,0,0], [0,0,yc,2-yc]])
    return X
def obj(p):
    X = build(p); Y = G @ X
    d = abs(np.linalg.det(Y))
    l = fast_lambda1(Y, cap=20)
    if l is None: return 10.0
    return d + 50*max(0, 2-l)
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
best = (9, None)
for it in range(400):
    p0 = rng.normal(size=5)
    r = minimize(obj, p0, method='Nelder-Mead', options={'maxiter': 2000, 'xatol':1e-9, 'fatol':1e-11})
    if r.fun < best[0]:
        best = (r.fun, r.x); X = build(r.x); Y = G @ X
        print(f"it {it}: det {abs(np.linalg.det(Y)):.5f} lambda1 {fast_lambda1(Y, cap=20):.5f} X1 {np.round(X[0],4)} y {np.round(X[3,2:],4)}", flush=True)
