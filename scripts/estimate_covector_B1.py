"""Cost estimate of a column-wise ("covector") Stage-B1 enumeration on zeta's contact set (DEAD_ENDS B15).
Usage: python3 scripts/estimate_covector_B1.py [partial_samples] [leaf_samples]   (defaults 1000, 5000)

Column view of a face datum: column j of the 15x4 sign matrix (s_cj) is the sign vector of (c . Y_j)_c,
Y_j = column j of X, i.e. a covector of the central arrangement {c^perp : c in S} in R^4.  A column-wise search
picks a multiset of 4 covector classes (mod the ambient group: column permutations and signs) and tests the
face-datum LP at the leaf.  This script
  1. enumerates all covectors (incremental, strict LPs) and checks that zeta's four columns are among them,
     and that every column of the facet-only B1 leaves (results/stageB_faces_zeta_facet.pkl) is a tope;
  2. tests the partial LP for 1-3 columns (strict cells, sum_{j<=k} |c.Y_j| <= 2): it never prunes, since a
     partial assignment can be scaled towards 0 -- so the tree has C(N+k-1, k) nodes at level k;
  3. samples random leaves with the face-datum LP of faces.py and times it.
Output: results/covector_B1_estimate.json."""
import sys, os, json, time, pickle, random
from math import comb
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.faces import _rows_for, lp_feasible
from cp4.lp import solve, BACKEND
from cp4.lattice import B_of_z, ZETA
root = os.path.join(os.path.dirname(__file__), '..')
n_partial = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
n_leaf = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
A = np.array(S, float) / 4
C = [tuple(F(x, 4) for x in c) for c in S]

def cone_ok(sig):
    """Is the relatively open cone {y : sign(A_i . y) = sig_i} nonempty?  (max common slack, |y| <= 1)"""
    ub, eq = [], []
    for i, s in enumerate(sig):
        (eq.append(np.r_[A[i], 0.0]) if s == 0 else ub.append(np.r_[-s * A[i], 1.0]))
    c = np.zeros(5); c[-1] = -1
    ok, x, _ = solve(c, A_ub=np.array(ub) if ub else None, b_ub=np.zeros(len(ub)) if ub else None,
                     A_eq=np.array(eq) if eq else None, b_eq=np.zeros(len(eq)) if eq else None,
                     bounds=[(-1, 1)] * 4 + [(None, 1.0)])
    return ok and (x[-1] > 1e-9 or not ub)

t0 = time.time()
cov = [()]
for m in range(len(S)):
    cov = [s + (v,) for s in cov for v in (1, 0, -1) if cone_ok(s + (v,))]
canon = lambda s: s if next(v for v in s if v) > 0 else tuple(-v for v in s)
classes = sorted({canon(s) for s in cov if any(s)})
topes = [s for s in classes if all(s)]
N, T = len(classes), len(topes)
print(f"{len(cov)} covectors (incl. 0) in {time.time() - t0:.0f}s: N = {N} classes mod ±, T = {T} topes", flush=True)

X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
covset = set(cov)
zeta_cols_ok = all(tuple(int(np.sign(round(v, 9))) for v in A @ X[:, j]) in covset for j in range(4))
facet_ok = None
f = os.path.join(root, 'results', 'stageB_faces_zeta_facet.pkl')
if os.path.exists(f):
    R = pickle.load(open(f, 'rb'))
    facet_ok = all(all(tuple(L['signs'][i][j] for i in range(15))) and
                   tuple(L['signs'][i][j] for i in range(15)) in covset for L in R['leaves'] for j in range(4))
print(f"zeta's columns are covectors: {zeta_cols_ok}; facet-only B1 leaf columns are topes: {facet_ok}", flush=True)
assert zeta_cols_ok and facet_ok is not False

def partial_ok(cols):
    k = len(cols); n = 4 * k + 1; ub, bu, eq = [], [], []
    for i in range(len(S)):
        r = np.zeros(n)
        for j, s in enumerate(cols):
            if s[i] == 0:
                e = np.zeros(n); e[4 * j:4 * j + 4] = A[i]; eq.append(e)
            else:
                g = np.zeros(n); g[4 * j:4 * j + 4] = -s[i] * A[i]; g[-1] = 1.0; ub.append(g); bu.append(0.0)
                r[4 * j:4 * j + 4] += s[i] * A[i]
        ub.append(r); bu.append(2.0)
    c = np.zeros(n); c[-1] = -1
    ok, x, _ = solve(c, A_ub=np.array(ub), b_ub=np.array(bu), A_eq=np.array(eq) if eq else None,
                     b_eq=np.zeros(len(eq)) if eq else None, bounds=[(-2, 2)] * (4 * k) + [(None, 1.0)])
    return ok and x[-1] > 1e-9

def leaf_ok(cols):
    if any(all(s[i] == 0 for s in cols) for i in range(len(S))):
        return None
    eqs, st = [], []
    for i in range(len(S)):
        e, t = _rows_for(C[i], tuple(s[i] for s in cols)); eqs += e; st += t
    return lp_feasible(eqs, st)[0]

rng = random.Random(0)
out = {'N_classes': N, 'T_topes': T, 'covectors_incl_zero': len(cov), 'backend': BACKEND,
       'zeta_columns_are_covectors': zeta_cols_ok, 'facet_leaf_columns_are_topes': facet_ok, 'partial': {}, 'leaf': {}}
for k in (1, 2, 3):
    feas = sum(partial_ok([rng.choice(classes) for _ in range(k)]) for _ in range(n_partial))
    out['partial'][k] = {'samples': n_partial, 'feasible': feas}
    print(f"level {k}: partial LP feasible for {feas}/{n_partial} random {k}-tuples", flush=True)
for name, pool in (('topes', topes), ('all', classes)):
    t = time.time(); res = [leaf_ok([rng.choice(pool) for _ in range(4)]) for _ in range(n_leaf)]
    nlp = sum(r is not None for r in res); dt = (time.time() - t) / max(nlp, 1)
    leaves = comb(len(pool) + 3, 4)
    out['leaf'][name] = {'samples': n_leaf, 'lps': nlp, 'feasible': sum(r is True for r in res),
                         'ms_per_lp': dt * 1000, 'leaf_multisets': leaves,
                         'tree_nodes': sum(comb(len(pool) + k - 1, k) for k in range(5)),
                         'hours_at_measured_rate': leaves * dt / 3600}
    print(f"{name}: {out['leaf'][name]['feasible']} feasible of {nlp} leaf LPs, {dt * 1000:.2f} ms/LP; "
          f"C({len(pool)}+3,4) = {leaves:.3g} leaves -> {leaves * dt / 3600:.3g} h", flush=True)
json.dump(out, open(os.path.join(root, 'results', 'covector_B1_estimate.json'), 'w'), indent=1,
          default=lambda o: o.item())                      # numpy scalars
