"""Stage B part 1 gate: enumerate face data for zeta's contact set S (stratum 11, p=4) and check that
zeta's own face datum (all 15 contacts in facets) appears among the leaves up to the ambient group.
Usage: python3 scripts/stageB_faces_zeta.py [time_limit_seconds]"""
import sys, os, json, time, pickle
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.faces import enumerate_face_data, AMB, act_sign
from cp4.lattice import B_of_z, ZETA
from collections import Counter
root = os.path.join(os.path.dirname(__file__), '..')
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
# zeta's X in canonical frame coordinates (convention: row k = sg[k] * X0[pm[k]])
X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
zsig = {i: tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / 4) @ X) for i, c in enumerate(S)}
print('zeta face datum:', zsig, flush=True)
tl = float(sys.argv[1]) if len(sys.argv) > 1 else None
leaves, st = enumerate_face_data(S, 4, log_every=1000, time_limit=tl)
hit = [k for k, L in enumerate(leaves)
       if any(all(act_sign(g, L['signs'][i]) == zsig[i] for i in zsig) for g in AMB)]
print(f"zeta S: {len(leaves)} leaves, dims {dict(Counter(l['dim'] for l in leaves))}, stats {st}")
print('zeta datum found among leaves:', hit, flush=True)
pickle.dump({'leaves': leaves, 'stats': st, 'zeta_hits': hit},
            open(os.path.join(root, 'results', 'stageB_faces_zeta.pkl'), 'wb'))
