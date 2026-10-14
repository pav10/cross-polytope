"""Gate G6a: zeta's own face datum survives every Stage-B1 prune (coordinatewise allowed signs and
strict LP feasibility) when its contacts are added one by one.  Fast (~10 s)."""
import sys, json, numpy as np
import os; sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src')); os.chdir(os.path.join(os.path.dirname(__file__), '..'))
from fractions import Fraction as F
from cp4.faces import _rows_for, lp_feasible, allowed_signs
from cp4.lattice import B_of_z, ZETA
z=json.load(open('data/zeta.json'))
S=[tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
X0=np.array(z['maxdet_frame_basis_coords'],float)@np.array(B_of_z(ZETA),float)
pm=z['signed_perm_to_canonical']['perm']; sg=z['signed_perm_to_canonical']['signs']
X=np.array([sg[k]*X0[pm[k]] for k in range(4)])

eqs=[];st=[]
for k,c in enumerate(S):
    s=tuple(int(np.sign(round(v,9))) for v in (np.array(c)/4)@X)
    cf=tuple(F(x,4) for x in c)
    al=allowed_signs(cf,eqs,st)
    e,t=_rows_for(cf,s); eqs+=e; st+=t
    ok = (s in al) and lp_feasible(eqs,st)[0]
    print(k, s, 'allowed' if s in al else 'NOT ALLOWED', len(al), 'LP ok' if ok else 'LP DEAD')
    assert ok
print('WITNESS OK')
