"""Regenerate data/strata.json (39 strata) and cross-check the classical 6 clean types."""
import sys, json, time, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import enumerate_strata
from collections import Counter
t = time.time()
S = enumerate_strata(24, strict=False)
cnt = Counter(s['p'] for s in S)
print(f"{len(S)} strata in {time.time()-t:.0f}s; by p: {dict(sorted(cnt.items()))}")
assert len(S) == 39 and sorted(cnt) == [1,2,3,4,5,6,7,8,9,10,12,16,18]
assert sum(s['n_superlattices_in_orbit'] for s in S) == 1040, sum(s['n_superlattices_in_orbit'] for s in S)
C = enumerate_strata(24, strict=True)
print(f"clean (strict Lee) orbits: {len(C)}, superlattices {sum(s['n_superlattices_in_orbit'] for s in C)}")
assert len(C) == 6 and sum(s['n_superlattices_in_orbit'] for s in C) == 54
for s in S:
    print(f"  id={s['id']:2d} p={s['p']:2d} orbit={s['n_superlattices_in_orbit']:4d} "
          f"|stab|={s['stabilizer_order']:3d} cand={len(s['candidates']):3d} forced={len(s['forced']):2d}")
json.dump(S, open(os.path.join(os.path.dirname(__file__), '..', 'data', 'strata.json'), 'w'))
