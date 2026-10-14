"""Check that the highspy and scipy LP backends (src/cp4/lp.py) give identical prune verdicts.
Usage: python3 scripts/check_lp_backend.py [n_lps]   (default 2000)
Records every LP issued by Stage-B1 code (faces.lp_feasible, faces.allowed_signs) along random walks on
zeta's contact set, re-solves each with both backends and compares the decisions the code takes from them:
strict-slack feasibility (t > 1e-9) and the closure-range sign class of each optimum (eps = 1e-9)."""
import sys, os, json, time, random
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from fractions import Fraction as F
import cp4.faces as faces
from cp4 import lp
if lp.highspy is None:
    sys.exit('highspy is not installed: nothing to compare (pip install highspy)')
root = os.path.join(os.path.dirname(__file__), '..')
z = json.load(open(os.path.join(root, 'data', 'zeta.json')))
C = [tuple(F(x, 4) for x in c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
want = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
rec = []
faces.lp_solve = lambda *a, **k: (rec.append((a, k)), lp.solve_scipy(*a, **k))[1]
rng = random.Random(0)
while len(rec) < want:                      # random walks through B1 nodes
    eqs, st, todo = [], [], list(range(len(C)))
    rng.shuffle(todo)
    for ci in todo:
        faces.allowed_signs(C[ci], eqs, st)
        for s in rng.sample(faces.SIGNVECS, 6):
            e, t = faces._rows_for(C[ci], s)
            if faces.lp_feasible(eqs + e, st + t)[0]:
                eqs, st = eqs + e, st + t
                break
        else:
            break
def decision(res, slack):
    ok, x, fun = res
    if not ok:
        return 'infeasible'
    if slack:
        return x[-1] > 1e-9
    return -1 if fun < -1e-9 else (0 if fun <= 1e-9 else 1)
mismatch, times = 0, {'scipy': 0.0, 'highspy': 0.0}
for a, k in rec:
    slack = k['bounds'][-1] == (None, 1.0)
    out = {}
    for name, f in (('scipy', lp.solve_scipy), ('highspy', lp.solve_highs)):
        t = time.time(); out[name] = decision(f(*a, **k), slack); times[name] += time.time() - t
    mismatch += out['scipy'] != out['highspy']
n = len(rec)
print(f"{n} LPs: {mismatch} verdict mismatches; scipy {times['scipy'] / n * 1000:.2f} ms/LP, "
      f"highspy {times['highspy'] / n * 1000:.2f} ms/LP ({times['scipy'] / times['highspy']:.1f}x)")
print('BACKENDS AGREE' if mismatch == 0 else 'BACKENDS DISAGREE')
sys.exit(1 if mismatch else 0)
