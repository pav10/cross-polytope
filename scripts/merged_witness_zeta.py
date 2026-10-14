"""Gate G7a: zeta's own labels survive every prune of the merged search (src/cp4/merged.py).
In the search order for zeta's frame pattern: 'in' with zeta's face datum for its 15 contacts, 'out' for
the 53 other candidates of stratum 11; the full leaf must pass the exclusion test with dim A = 3.  Fast."""
import sys, os, json
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import stabilizer
from cp4.merged import Merged
from cp4.lattice import B_of_z, ZETA
root = os.path.join(os.path.dirname(__file__), '..')
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
s = json.load(open(os.path.join(root, 'data', 'strata.json')))[z['stratum_id']]
# zeta's X in canonical frame coordinates (convention of scripts/stageB_faces_zeta.py)
X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
M = Merged(s, stabilizer([tuple(r) for r in s['residues']], s['p']))
sign = lambda c: tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / s['p']) @ X)
contacts = {M.cand.index(c) if c in M.cand else M.cand.index(tuple(-x for x in c))
            for c in map(tuple, z['S_frame_numerators_over_p_canonical_stratum_coords'])}
n = M.root_pattern(tuple(sign(M.cand[i]) for i in M.frame))
for ci in M.order(n):
    if ci in contacts:
        assert sign(M.cand[ci]) in n.root[2][ci], f'frame-sign filter rejects zeta at {M.cand[ci]}'
        m = M.include(n, ci, sign(M.cand[ci]))
    else:
        m = M.exclude(n, ci)
    assert m is not None, f'zeta path dies at {M.cand[ci]}: {M.stats}'
    assert m.rank < 16, 'zeta should never reach rank 16 (its family has dim 3)'
    n = m
assert M.full_leaf(n) == 1, f'zeta leaf rejected: {M.stats}'
L = M.leaves[-1]
assert L['dim'] == 3 and len(L['S']) == 15, L
print(f"zeta path: {len(contacts)} in, {len(M.order(n)) + 4 - len(contacts)} out, leaf dim {L['dim']}, rank {n.rank}")
print('WITNESS OK')
