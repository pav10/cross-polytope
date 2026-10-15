"""Measure an LP-free covector pre-filter inside row-wise Stage B1 on zeta's contact set (static order).
Usage: python3 scripts/measure_covector_prefilter.py [facet|full] [seconds] [--off]

Pre-filter: with a static contact order the decided contacts at depth k are order[:k]; every column of the
partial sign matrix must be the restriction to order[:k] of a NONZERO covector of the arrangement
{c^perp : c in S} (column j of X is nonzero and realises such a covector).  Candidate sign vectors failing this
are rejected by a set lookup before the LP.  It never rejects the true datum, and any LP-feasible candidate
passes it apart from all-zero partial columns that only the zero covector extends, so the search tree is
(essentially) unchanged and the gain is the LP calls saved.  The original faces._children is replaced in-process
only.  --off runs the same instrumented search without the filter (baseline for the same traversal).
Output: results/covector_prefilter_<mode>[_off].json."""
import sys, os, json, time
import numpy as np
from fractions import Fraction as F
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import cp4.faces as faces
from cp4.lp import solve, BACKEND
root = os.path.join(os.path.dirname(__file__), '..')
mode = sys.argv[1] if len(sys.argv) > 1 else 'facet'
seconds = float(sys.argv[2]) if len(sys.argv) > 2 else None
use_filter = '--off' not in sys.argv
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
A = np.array(S, float) / 4

def cone_ok(sig):
    ub, eq = [], []
    for i, s in enumerate(sig):
        (eq.append(np.r_[A[i], 0.0]) if s == 0 else ub.append(np.r_[-s * A[i], 1.0]))
    c = np.zeros(5); c[-1] = -1
    ok, x, _ = solve(c, A_ub=np.array(ub) if ub else None, b_ub=np.zeros(len(ub)) if ub else None,
                     A_eq=np.array(eq) if eq else None, b_eq=np.zeros(len(eq)) if eq else None,
                     bounds=[(-1, 1)] * 4 + [(None, 1.0)])
    return ok and (x[-1] > 1e-9 or not ub)

cov = [()]
for m in range(len(S)):
    cov = [s + (v,) for s in cov for v in (1, 0, -1) if cone_ok(s + (v,))]
cov = [s for s in cov if any(s)]
order = list(range(len(S)))                       # static order of enumerate_face_data(dynamic=False)
R = [{tuple(s[i] for i in order[:k]) for s in cov} for k in range(len(S) + 1)]
print(f"{len(cov)} nonzero covectors; restriction sets per depth: {[len(r) for r in R]}", flush=True)

stats = {'lp': 0, 'filtered': 0}
_lp = faces.lp_feasible
def counted_lp(eqs, st):
    stats['lp'] += 1
    return _lp(eqs, st)
faces.lp_feasible = counted_lp

def children(C, group, eqs, stricts, signs, dynamic, vecs, st):
    """faces._children for the static order, plus the covector pre-filter before each LP."""
    undecided = [i for i in range(len(C)) if i not in signs]
    if not undecided:
        return None
    ci = undecided[0]
    k = len(signs)
    kids = []
    for s in faces.orbit_reps(group, vecs):
        if use_filter and any(tuple(signs[i][j] for i in order[:k]) + (s[j],) not in R[k + 1] for j in range(4)):
            stats['filtered'] += 1; continue
        e, t = faces._rows_for(C[ci], s)
        ok, _ = faces.lp_feasible(eqs + e, stricts + t)
        if not ok:
            st['lp_dead'] += 1; continue
        kids.append((s, [g for g in group if faces.act_sign(g, s) == s], eqs + e, stricts + t))
    return ci, kids
faces._children = children

leaves, st = faces.enumerate_face_data(S, 4, dynamic=False, signvecs='facet' if mode == 'facet' else 'all',
                                       determine=True, time_limit=seconds)
out = {'mode': mode, 'filter': use_filter, 'backend': BACKEND, 'complete': st['complete'], 'nodes': st['nodes'],
       'leaves': len(leaves), 'leaf_dims': sorted(L['dim'] for L in leaves), 'seconds': st['seconds'],
       'lp_calls': stats['lp'], 'lp_dead': st['lp_dead'], 'filtered': stats['filtered'],
       'nodes_per_s': st['nodes'] / st['seconds'], 'covectors_nonzero': len(cov)}
print(json.dumps(out), flush=True)
json.dump(out, open(os.path.join(root, 'results', f"covector_prefilter_{mode}{'' if use_filter else '_off'}.json"), 'w'),
          indent=1)
