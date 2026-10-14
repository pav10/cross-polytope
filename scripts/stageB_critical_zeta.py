"""Gate for Stage B part 2: zeta's own face datum must yield the cubic and the critical value p*D*.
Builds the affine family directly from zeta's face datum (no dependence on the enumerator)."""
import sys, os, json, time
import numpy as np, sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fractions import Fraction as F
from cp4.faces import _rows_for
from cp4.exact import solve_affine
from cp4.lattice import B_of_z, ZETA
from cp4.critical import critical_value_polynomial
root = os.path.join(os.path.dirname(__file__), '..')
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
eqs = []
for c in S:
    s = tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / 4) @ X)
    e, _ = _rows_for(tuple(F(x, 4) for x in c), s)
    eqs += e
x0, N = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
print('dim A =', len(N))
t = time.time()
R, f, tv = critical_value_polynomial(x0, N)
print('eliminant factors:', R, f'({time.time()-t:.1f}s)')
D = sp.Symbol('D')
target = 4 * sp.Float('0.80822000090353537301', 30)
for fac, _ in R[1]:
    for r in sp.Poly(fac, D).nroots(n=25):
        if abs(sp.im(r)) < 1e-20:
            print(f'  root {sp.re(r)}   (4*D* = {target})', '<-- zeta' if abs(sp.re(r) - target) < 1e-15 else '')
