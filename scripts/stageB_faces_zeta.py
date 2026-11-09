"""Stage B part 1 gate: enumerate face data for zeta's contact set S (stratum 11, p=4) and check that
zeta's own face datum (all 15 contacts in facets) appears among the leaves up to the ambient group.
Usage: python3 scripts/stageB_faces_zeta.py [time_limit_seconds] [--facet-only] [--static] [--determine]
                                            [--prefilter] [--frame-symmetry] [--estimate K] [--resume]
                                            [--heartbeat N] [--out F]
Output: results/stageB_faces_zeta.pkl (full run, gate G6) or results/stageB_faces_zeta_facet.pkl
(--facet-only, gate G6f), or --out F, rewritten atomically at every heartbeat as a checkpoint; --resume
continues from it (the options must match).  --prefilter: certified covector column pre-filter (cp4.covectors).
--frame-symmetry: isomorph rejection under the 24 frame symmetries of zeta's S (static order; frame vectors first;
leaves complete up to G_S x AMB).
--estimate K prints a Knuth tree-size estimate instead (results/stageB_faces_zeta_estimate.json)."""
import sys, os, json, time, pickle, random, argparse
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.faces import enumerate_face_data, estimate_face_data, AMB, act_sign, frame_stabilizer, contact_action
from cp4.lattice import B_of_z, ZETA
from cp4 import lp
from collections import Counter
root = os.path.join(os.path.dirname(__file__), '..')
ap = argparse.ArgumentParser()
ap.add_argument('seconds', type=float, nargs='?')
ap.add_argument('--facet-only', action='store_true')
ap.add_argument('--static', action='store_true')
ap.add_argument('--determine', action='store_true')
ap.add_argument('--prefilter', action='store_true')
ap.add_argument('--frame-symmetry', action='store_true')
ap.add_argument('--out')
ap.add_argument('--estimate', type=int)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--resume', action='store_true')
ap.add_argument('--heartbeat', type=int, default=1000)
a = ap.parse_args()
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
# zeta's X in canonical frame coordinates (convention: row k = sg[k] * X0[pm[k]])
X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
zsig = {i: tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / 4) @ X) for i, c in enumerate(S)}
opts = dict(dynamic=not a.static, signvecs='facet' if a.facet_only else 'all', determine=a.determine,
            prefilter=a.prefilter, frame_symmetry=a.frame_symmetry)
G = frame_stabilizer(S, json.load(open(os.path.join(root, 'data', 'strata.json')))[11]['residues'], 4) \
    if a.frame_symmetry else None
run_opts = dict({k: v for k, v in opts.items() if k != 'frame_symmetry'}, frame_group=G)
print('zeta face datum:', zsig, '| options', opts, '| LP backend', lp.BACKEND, flush=True)

if a.estimate:
    e = estimate_face_data(S, 4, a.estimate, random.Random(a.seed), max_seconds=a.seconds, **run_opts)
    print(f"Knuth estimate: {e}", flush=True)
    f = os.path.join(root, 'results', 'stageB_faces_zeta_estimate.json')
    prev = json.load(open(f)) if os.path.exists(f) else []
    json.dump(prev + [dict(e, seed=a.seed, backend=lp.BACKEND)], open(f, 'w'), indent=1)
    sys.exit(0)

out = a.out or os.path.join(root, 'results', 'stageB_faces_zeta' + ('_facet' if a.facet_only else '') + '.pkl')
FRAME = [contact_action(S, g) for g in (G or [((0, 1, 2, 3), (1, 1, 1, 1))])]
def hits(leaves):
    """Leaves equal to zeta's datum up to AMB (columns), and up to the frame symmetries when they are used."""
    out = []
    for k, L in enumerate(leaves):
        for ca in FRAME:          # frame image: contact i takes eps * (sign vector of contact j)
            t = {i: tuple(e * x for x in L['signs'][j]) for i, (j, e) in enumerate(ca)}
            if any(all(act_sign(q, t[i]) == zsig[i] for i in zsig) for q in AMB):
                out.append(k); break
    return out
def dump(leaves, st):
    tmp = out + '.tmp'
    pickle.dump({'leaves': leaves, 'stats': st, 'zeta_hits': hits(leaves), 'options': opts}, open(tmp, 'wb'))
    os.replace(tmp, out)
resume = None
if a.resume and os.path.exists(out):
    R = pickle.load(open(out, 'rb'))
    old = R.get('options')
    old = dict({'prefilter': False, 'frame_symmetry': False}, **old) if old else old   # older checkpoints
    # A run without the pre-filter may be continued with it ('mixed'): finished subtrees without it contain
    # every leaf the pre-filtered search would find (the filter only drops candidates no NONSINGULAR X completes).
    mixed = bool(old) and a.prefilter and old['prefilter'] in (False, 'mixed') and \
        dict(old, prefilter=True) == dict(opts, prefilter=True) and opts['prefilter'] is True
    assert old == opts or mixed, f"checkpoint options {old} differ from {opts}"
    if mixed:
        opts['prefilter'] = 'mixed'
        print("continuing a checkpoint written without --prefilter with it (options recorded as prefilter='mixed')",
              flush=True)
    if R['stats'].get('complete'):
        sys.exit(f"checkpoint {out} is already complete: zeta hits {R['zeta_hits']} - nothing to do")
    resume = {'done_prefixes': R['stats']['done_prefixes'], 'leaves': R['leaves']}
    print(f"resuming: {len(resume['done_prefixes'])} finished subtrees", flush=True)
run_opts['prefilter'] = opts['prefilter']
leaves, st = enumerate_face_data(S, 4, log_every=a.heartbeat, time_limit=a.seconds, resume=resume,
                                 on_heartbeat=dump, **run_opts)
hit = hits(leaves)
print(f"zeta S: {len(leaves)} leaves, dims {dict(Counter(l['dim'] for l in leaves))}, "
      f"stats { {k: v for k, v in st.items() if k != 'done_prefixes'} }")
print('zeta datum found among leaves:', hit, flush=True)
dump(leaves, st)
