"""Gate G8: the covector column pre-filter of Stage B1 (cp4.covectors) on zeta's contact set S.
  1. certified_covectors: the float LP enumeration is certified exactly (exact point in every cell, zero sets
     are flats, per-flat face counts = Zaslavsky's r_F); known answer 1 705 faces, 372 topes, 117 flats.
  2. zeta's four columns are covectors realised by nonzero vectors, and zeta's face datum passes the filter at
     every step of its static B1 path and of 20 random contact orders (the dynamic order decides contacts in
     a data-dependent order); the filter is not vacuous along the path.
  3. every column of every leaf of the facet-only pass (results/stageB_faces_zeta_facet.pkl, if present) is
     an admissible covector.
Usage: python3 scripts/check_covectors_zeta.py      (~20 s; prints COVECTORS OK)"""
import sys, os, json, pickle, random, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.covectors import certified_covectors, ColumnFilter
from cp4.faces import SIGNVECS
from cp4.lattice import B_of_z, ZETA
root = os.path.join(os.path.dirname(__file__), '..')
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]

t = time.time()
cov, cert = certified_covectors(S)
print(f"certificate ({time.time() - t:.0f}s): {cert}", flush=True)
assert cert == {'covectors': 1705, 'nonzero_realised': 1704, 'topes': 372, 'flats_by_rank': [1, 15, 55, 45, 1],
                'faces_zaslavsky': 1705}, cert
covset = set(cov)

X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
zsig = [tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / 4) @ X) for c in S]
assert all(all(s) for s in zsig)
assert all(tuple(s[j] for s in zsig) in covset for j in range(4)), "a column of zeta's datum is not a covector"

flt = ColumnFilter(cov)
def path_ok(order, verbose=False):
    signs, rejected = {}, 0
    for i in order:
        al = flt.allowed(signs, i)
        if not ColumnFilter.passes(al, zsig[i]):
            return None
        rejected += sum(not ColumnFilter.passes(al, s) for s in SIGNVECS)
        if verbose:
            print(f"  contact {i}: {sum(ColumnFilter.passes(al, s) for s in SIGNVECS)}/80 sign vectors pass", flush=True)
        signs[i] = zsig[i]
    return rejected
rej = path_ok(range(len(S)), verbose=True)
assert rej is not None, "zeta's datum rejected on its static path"
assert rej > 0, 'filter vacuous along the static path'
rng = random.Random(0)
for _ in range(20):
    order = list(range(len(S))); rng.shuffle(order)
    assert path_ok(order) is not None, f"zeta's datum rejected along order {order}"
msg = f"static path: {rej} of {80 * len(S)} candidates rejected LP-free; 20 random orders pass"

f = os.path.join(root, 'results', 'stageB_faces_zeta_facet.pkl')
if os.path.exists(f):
    R = pickle.load(open(f, 'rb'))
    assert all(tuple(L['signs'][i][j] for i in range(len(S))) in covset for L in R['leaves'] for j in range(4))
    msg += f"; all columns of the {len(R['leaves'])} facet-only leaves are covectors"
print(msg)
print('COVECTORS OK')
