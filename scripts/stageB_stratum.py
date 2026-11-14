"""Stage B1 over all Stage-A contact sets of one stratum: cluster job array, resumable per set.

Usage (from the campaign root):
  python3 scripts/stageB_stratum.py <sid> [--task i --ntasks N] [--seconds T] [--only k,k,...] [--heartbeat N]
                                         [--order filter|static] [--sets F]
  python3 scripts/stageB_stratum.py <sid> --summary
Each set k of results/stageA_<sid>.json (or --sets F) with k = i (mod N), or each k in --only, is searched with the
complete supporting-facet B1 (STATUS §3 Thm B.1', Lemmas B.3-B.5): static frame-first order, frame symmetry,
covector pre-filter, determine-then-check, half-open facet cells, ball rows and the eutaxy leaf prune (Lemma B.6;
--no-eutaxy keeps every leaf).  --order filter (LP-free contact order) is available but measured no gain
(DEAD_ENDS B20).  --slack-tol t (default 1e-9, the campaign's strict-slack threshold): a cell is kept iff its max
common slack exceeds t; 1e-6 halves the tree on zeta's S with the same leaves (STATUS §4) - the owner's call.
Checkpoint per set: results/stageB_<sid>/set_<k>.pkl (leaves, stats, options), rewritten
atomically at every heartbeat; an interrupted set resumes from its finished depth-3 subtrees; complete sets are
skipped.  The run stops itself before T seconds (default: none), so a resubmission continues where it stopped.
Never run two submissions for the same stratum at the same time (no locking).
--summary: complete / partial / untouched sets, nodes, seconds, leaves by dimension; for stratum 11 also whether
zeta's set (if complete) has zeta's face datum among its leaves.  Writes results/stageB_<sid>_summary.json."""
import sys, os, json, time, pickle, argparse
from collections import Counter
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.faces import enumerate_face_data, frame_stabilizer, datum_hits
from cp4 import lp
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
ap = argparse.ArgumentParser()
ap.add_argument('sid', type=int)
ap.add_argument('--task', type=int, default=0)
ap.add_argument('--ntasks', type=int, default=1)
ap.add_argument('--seconds', type=float)
ap.add_argument('--only')
ap.add_argument('--heartbeat', type=int, default=2000)
ap.add_argument('--order', default='static', choices=('static', 'filter'))
ap.add_argument('--sets')
ap.add_argument('--summary', action='store_true')
ap.add_argument('--no-eutaxy', action='store_true')
ap.add_argument('--slack-tol', type=float, default=1e-9)
a = ap.parse_args()
st = json.load(open(os.path.join(root, 'data', 'strata.json')))[a.sid]
p = st['p']
setfile = a.sets or os.path.join(root, 'results', f'stageA_{a.sid}.json')
SETS = [[tuple(c) for c in S] for S in json.load(open(setfile))['contact_sets']]
outdir = os.path.join(root, 'results', f'stageB_{a.sid}')
os.makedirs(outdir, exist_ok=True)
OPTS = dict(dynamic=False, determine=True, prefilter=True, frame_symmetry=True, supporting=True, ball=True,
            contact_order=a.order, eutaxy=not a.no_eutaxy, slack_tol=a.slack_tol, setfile=os.path.basename(setfile))
import cp4.faces as _faces
_faces.SLACK_TOL = a.slack_tol
path = lambda k: os.path.join(outdir, f'set_{k:05d}.pkl')


def load(k):
    return pickle.load(open(path(k), 'rb')) if os.path.exists(path(k)) else None


def zeta_signs(S):
    """zeta's own face datum in the indexing of S, if S is zeta's contact set (stratum 11), else None."""
    from cp4.lattice import B_of_z, ZETA
    z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
    Sz = {tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']}
    if a.sid != z['stratum_id'] or set(S) != Sz:
        return None
    X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
    pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
    X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
    return {i: tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / p) @ X) for i, c in enumerate(S)}


if a.summary:
    rows, dims, nodes, secs = [], Counter(), 0, 0.0
    zeta = None
    for k, S in enumerate(SETS):
        R = load(k)
        if R is None:
            rows.append('untouched'); continue
        stt = R['stats']; nodes += stt.get('nodes_total', stt['nodes']); secs += stt.get('seconds_total', stt.get('seconds', 0))
        rows.append('complete' if stt.get('complete') else 'partial')
        if stt.get('complete'):
            dims.update(L['dim'] for L in R['leaves'])
        zs = zeta_signs(S)
        if zs is not None:
            zeta = {'set': k, 'complete': bool(stt.get('complete')),
                    'hits': datum_hits(S, R['leaves'], zs, frame_stabilizer(S, st['residues'], p))}
    c = Counter(rows)
    out = {'stratum': a.sid, 'sets': len(SETS), 'complete': c['complete'], 'partial': c['partial'],
           'untouched': c['untouched'], 'nodes': nodes, 'cpu_seconds': round(secs), 'leaves_by_dim': dict(sorted(dims.items())),
           'zeta': zeta, 'options': OPTS}
    json.dump(out, open(os.path.join(root, 'results', f'stageB_{a.sid}_summary.json'), 'w'), indent=1)
    print(json.dumps(out))
    if zeta is not None and zeta['complete']:
        assert zeta['hits'], "zeta's set is complete but zeta's face datum is not among its leaves"
    sys.exit(0)

todo = [int(x) for x in a.only.split(',')] if a.only else [k for k in range(len(SETS)) if k % a.ntasks == a.task]
t_start = time.time()
deadline = t_start + a.seconds if a.seconds else None
print(f'stratum {a.sid} (p={p}): task {a.task}/{a.ntasks}, {len(todo)} of {len(SETS)} sets, options {OPTS}, '
      f'LP backend {lp.BACKEND}', flush=True)
for k in todo:
    R = load(k)
    if R is not None:
        old = dict({'slack_tol': 1e-9}, **R['options'])          # checkpoints written before --slack-tol existed
        assert old == OPTS, f'set {k}: checkpoint options {R["options"]} differ from {OPTS}'
        if R['stats'].get('complete'):
            continue
    left = deadline - time.time() if deadline else None
    if left is not None and left < 120:
        print(f'time budget used up before set {k}; stopping (resubmit to continue)', flush=True)
        break
    S = SETS[k]
    G = frame_stabilizer(S, st['residues'], p)
    resume = {'done_prefixes': R['stats']['done_prefixes'], 'leaves': R['leaves']} if R else None
    prev = R['stats'].get('seconds_total', R['stats'].get('seconds', 0)) if R else 0.0
    prev_n = R['stats'].get('nodes_total', R['stats']['nodes']) if R else 0

    def dump(leaves, stt, k=k, prev=prev, prev_n=prev_n):
        stt = dict(stt, seconds_total=prev + stt.get('seconds', 0), nodes_total=prev_n + stt['nodes'])
        tmp = path(k) + '.tmp'
        pickle.dump({'set': k, 'S': S, 'leaves': leaves, 'stats': stt, 'options': OPTS}, open(tmp, 'wb'))
        os.replace(tmp, path(k))
    t = time.time()
    leaves, stt = enumerate_face_data(S, p, log_every=a.heartbeat, time_limit=left, resume=resume, on_heartbeat=dump,
                                      dynamic=False, determine=True, prefilter=True, frame_group=G, supporting=True,
                                      ball=True, contact_order=a.order, eutaxy=OPTS['eutaxy'])
    dump(leaves, stt)
    print(f'set {k}: |S|={len(S)} |G_S|={len(G)} complete={stt["complete"]} nodes={stt["nodes"]} '
          f'leaves={len(leaves)} dims={dict(Counter(L["dim"] for L in leaves))} [{time.time()-t:.0f}s]', flush=True)
    if not stt['complete']:
        print('time budget used up; stopping (resubmit to continue)', flush=True)
        break
print(f'task {a.task}: done in {time.time()-t_start:.0f}s', flush=True)
