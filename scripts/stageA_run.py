"""Stage A driver.  Usage: python3 scripts/stageA_run.py <stratum_id> [node_limit]
Writes results/stageA_<id>.json  (canonical contact sets + search statistics)."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import stabilizer
from cp4.contacts import enumerate_contact_sets, verify_contact_set
root = os.path.join(os.path.dirname(__file__), '..')
strata = json.load(open(os.path.join(root, 'data', 'strata.json')))
sid = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
s = strata[sid]
stab = stabilizer([tuple(r) for r in s['residues']], s['p'])
t = time.time()
res, st = enumerate_contact_sets(s, stab, node_limit=lim, log=20000)
for S in res:
    verify_contact_set(S, s)
from collections import Counter
sizes = Counter(len(S) for S in res)
print(f"stratum {sid} (p={s['p']}): {len(res)} contact sets (mod stabilizer, |stab|={len(stab)}), "
      f"sizes {dict(sorted(sizes.items()))}, {time.time()-t:.1f}s, stats {st}")
json.dump({'stratum': sid, 'p': s['p'], 'contact_sets': [list(map(list, S)) for S in sorted(res)],
           'stats': st, 'seconds': time.time() - t},
          open(os.path.join(root, 'results', f'stageA_{sid}.json'), 'w'))
