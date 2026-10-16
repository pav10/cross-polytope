"""Validation of the covector pre-filter (enumerate_face_data(prefilter=True)) against the search without it
on zeta's contact set: the trees must agree (the filter only saves LPs).
  * static full order (the G6 mode: --static --determine) and dynamic order: the first N nodes of the
    depth-first traversal, with and without the filter: identical per-depth node counts and leaves;
  * facet-only: the complete pre-filtered run (results/stageB_faces_zeta_facet_prefilter.pkl, from
    stageB_faces_zeta.py --facet-only --static --determine --prefilter --out ...) against the committed
    G6f run (results/stageB_faces_zeta_facet.pkl): identical leaf sets.
Usage: python3 scripts/compare_prefilter_zeta.py [N_static] [N_dynamic]   (defaults 2000, 1000)
Output: results/covector_prefilter_validation.json."""
import sys, os, json, pickle, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.faces import enumerate_face_data
from cp4.lp import BACKEND
root = os.path.join(os.path.dirname(__file__), '..')
n_static = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
n_dynamic = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]

def stop(leaves, st):
    raise TimeoutError                               # caught by enumerate_face_data: stops after log_every nodes

def key(L):
    return tuple(sorted(L['signs'].items()))

def first_nodes(n, **opts):
    leaves, st = enumerate_face_data(S, 4, log_every=n, on_heartbeat=stop, **opts)
    return {'nodes': st['nodes'], 'per_depth': st['per_depth'], 'leaves': sorted(map(key, leaves)),
            'lp_dead': st['lp_dead'], 'prefilter_dead': st['prefilter_dead'], 'seconds': st['seconds']}

out = {'backend': BACKEND}
ok = True
for name, n, opts in (('static_full', n_static, dict(dynamic=False, determine=True)),
                      ('dynamic_full', n_dynamic, dict(dynamic=True))):
    r0 = first_nodes(n, prefilter=False, **opts)
    r1 = first_nodes(n, prefilter=True, **opts)
    same = r0['per_depth'] == r1['per_depth'] and r0['leaves'] == r1['leaves']
    ok &= same
    out[name] = {'nodes': n, 'identical': same, 'per_depth': r0['per_depth'], 'leaves': len(r0['leaves']),
                 'off': {k: r0[k] for k in ('lp_dead', 'seconds')},
                 'on': {k: r1[k] for k in ('lp_dead', 'prefilter_dead', 'seconds')},
                 'speedup': r0['seconds'] / r1['seconds']}
    print(name, json.dumps(out[name]), flush=True)

fa, fb = (os.path.join(root, 'results', f) for f in ('stageB_faces_zeta_facet.pkl', 'stageB_faces_zeta_facet_prefilter.pkl'))
if os.path.exists(fa) and os.path.exists(fb):
    A, B = pickle.load(open(fa, 'rb')), pickle.load(open(fb, 'rb'))
    same = (A['stats']['complete'] and B['stats']['complete'] and B['options'].get('prefilter') is True
            and sorted(map(key, A['leaves'])) == sorted(map(key, B['leaves'])))
    ok &= same
    out['facet_only'] = {'identical_leaf_sets': same, 'leaves': len(B['leaves']), 'zeta_hits': B['zeta_hits'],
                         'nodes_off': A['stats']['nodes'], 'nodes_on': B['stats']['nodes'],
                         'seconds_off': A['stats']['seconds'], 'seconds_on': B['stats']['seconds'],
                         'prefilter_dead': B['stats']['prefilter_dead'], 'lp_dead_off': A['stats']['lp_dead'],
                         'lp_dead_on': B['stats']['lp_dead'], 'speedup': A['stats']['seconds'] / B['stats']['seconds']}
    print('facet_only', json.dumps(out['facet_only']), flush=True)
out['all_identical'] = bool(ok)
json.dump(out, open(os.path.join(root, 'results', 'covector_prefilter_validation.json'), 'w'), indent=1)
print('PREFILTER VALIDATION OK' if ok else 'PREFILTER VALIDATION FAILED')
sys.exit(0 if ok else 1)
