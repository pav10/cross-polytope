"""Knuth estimate of the G6 tree (zeta's S; static order, determine, prefilter) per depth-1 branch, i.e. per
orbit of contact 0's sign vector.  Usage: python3 scripts/knuth_G6_by_branch.py [probes] [seed]  (~40 s for 150)."""
import sys, os, json, random, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
from fractions import Fraction as F
import numpy as np
import cp4.faces as faces
z = json.load(open('data/zeta.json'))
S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
C = [tuple(F(x, 4) for x in c) for c in S]
st = {'lp_dead': 0, 'prefilter_dead': 0}
flt = faces._column_filter(C, st)
ci, kids = faces._children(C, faces.AMB, [], [], {}, False, faces.SIGNVECS, st, flt)
rng = random.Random(int(sys.argv[2]) if len(sys.argv) > 2 else 0)
probes = int(sys.argv[1]) if len(sys.argv) > 1 else 100
for img0, sub0, e0, s0 in kids:
    t = time.time(); ests = []
    for _ in range(probes):
        group, eqs, stricts, signs = sub0, e0, s0, {ci: img0}
        w, est = 1, 1
        while True:
            res = faces._children(C, group, eqs, stricts, signs, False, faces.SIGNVECS, st, flt)
            if res is None:
                break
            c, ks = res
            if not ks:
                break
            w *= len(ks); est += w
            img, group, eqs, stricts = rng.choice(ks)
            signs = dict(signs); signs[c] = img
            if faces._read_off(C, eqs, signs, faces.SIGNVECS) != 'branch':
                break
        ests.append(est)
    print(img0, f"mean {np.mean(ests):.3g} median {np.median(ests):.3g} max {np.max(ests):.3g} "
          f"({time.time() - t:.0f}s)", flush=True)
