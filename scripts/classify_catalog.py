"""Classify the numerical catalog (data/catalog_results.csv) of locally optimised lattices:
for each distinct density value: #minimal pairs, max-det frame index p, stratum id, contact set.
Output: data/catalog_local_optima.json  (validation targets for Stages A/B).  Floats only:
these are NUMERICAL facts [Num], used as regression targets, never as proof."""
import csv, json, sys, os
import numpy as np
from itertools import combinations, product
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import group, canon_group, act, SIGNED_PERMS
from cp4.contacts import canon_set
from cp4.strata import stabilizer
root = os.path.join(os.path.dirname(__file__), '..')
rows = list(csv.DictReader(open(os.path.join(root, 'data', 'catalog_results.csv'))))
strata = json.load(open(os.path.join(root, 'data', 'strata.json')))
skey = {(s['p'], tuple(map(tuple, s['residues']))): s['id'] for s in strata}
TOLR = float(os.environ.get("TOLR", "1e-5"))
by = {}
for r in rows:
    by.setdefault(round(float(r['density']), 8), []).append(r)
out = []
for dens in sorted(by, reverse=True):
    reps = by[dens]
    B = np.array([float(t) for t in reps[0]['basis_flat(16)'].split(';')]).reshape(4, 4)
    Yi = np.linalg.inv(B)
    lam = None
    b = [int(np.floor(1.05 * np.max(np.abs(Yi[:, j])))) for j in range(4)]
    Q = np.array([q for q in product(*[range(-k, k + 1) for k in b]) if any(q)], float)
    N = np.abs(Q @ B).sum(axis=1); lam = N.min()
    mins = [tuple(int(x) for x in q) for q, n in zip(Q, N) if n < lam * (1 + TOLR)]
    mins = sorted({q if next(x for x in q if x) > 0 else tuple(-x for x in q) for q in mins})
    dets = {}
    for t in combinations(mins, 4):
        d = int(round(abs(np.linalg.det(np.array(t, float)))))
        dets[d] = dets.get(d, 0) + 1
    if not dets:
        print(f"{dens:.8f}: only {len(mins)} near-minimal pairs at tol {TOLR}; skipped"); continue
    p = max(dets)
    frame = next(t for t in combinations(mins, 4) if int(round(abs(np.linalg.det(np.array(t, float))))) == p)
    Xi = np.linalg.inv(np.array(frame, float))
    gens = [tuple(int(round(p * x)) % p for x in row) for row in Xi]
    G = group(gens, p)
    key = canon_group(G, p) if G is not None else None
    sid = skey.get((p, key)) if key else None
    Snum = None
    if sid is not None:
        g = next(g for g in SIGNED_PERMS if tuple(sorted(act(g, rr, p) for rr in G)) == key)
        pm, sg = g
        Sn = []
        for q in mins:
            c = np.array(q, float) @ Xi
            n = tuple(int(round(p * sg[k] * c[pm[k]])) for k in range(4))
            Sn.append(max(n, tuple(-x for x in n)))
        stab = stabilizer([tuple(x) for x in strata[sid]['residues']], p)
        Snum = [list(x) for x in canon_set(sorted(Sn), stab, p)]
    out.append({'density': dens, 'det_at_lambda1_2': (2/3)/dens, 'hits': len(reps),
                'n_min_pairs': len(mins), 'minor_distribution': dets, 'p': p,
                'stratum_id': sid, 'S_frame_canonical': Snum,
                'example_basis_lambda1_eq_1': reps[0]['basis_flat(16)']})
    print(f"{dens:.8f} hits={len(reps):4d} pairs={len(mins):2d} p={p} stratum={sid} minors={dict(sorted(dets.items()))}")
json.dump(out, open(os.path.join(root, 'data', 'catalog_local_optima.json'), 'w'), indent=1)
