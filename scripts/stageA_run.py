"""Stage A driver.  Usage: python3 scripts/stageA_run.py <stratum_id> [node_limit] [--force] [--no-swap]
Writes results/stageA_<id>.json  (canonical contact sets + search statistics).  Strata excluded by proof
(cp4.strata.EXCLUDED, gate G9) are refused unless --force is given.  The frame-swap rule (A6, docs/FRAME_SWAP.md,
cp4.frameswap) is on by default; --no-swap gives the plain (A1)-(A5) search."""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import stabilizer, EXCLUDED
from cp4.contacts import enumerate_contact_sets, verify_contact_set
from cp4.frameswap import FrameSwapRule
root = os.path.join(os.path.dirname(__file__), '..')
strata = json.load(open(os.path.join(root, 'data', 'strata.json')))
force = '--force' in sys.argv
use_swap = '--no-swap' not in sys.argv
args = [x for x in sys.argv[1:] if x not in ('--force', '--no-swap')]
sid = int(args[0]); lim = int(args[1]) if len(args) > 1 else None
if sid in EXCLUDED and not force:
    sys.exit(f"stratum {sid} is excluded by proof (cp4.strata.EXCLUDED; gate G9) - nothing to search (--force to run)")
s = strata[sid]
stab = stabilizer([tuple(r) for r in s['residues']], s['p'])
t = time.time()
swap = FrameSwapRule(strata, sid) if use_swap else None
res, st = enumerate_contact_sets(s, stab, node_limit=lim, log=20000, swap=swap)
for S in res:
    verify_contact_set(S, s)
from collections import Counter
sizes = Counter(len(S) for S in res)
print(f"stratum {sid} (p={s['p']}, frame swap {'on' if use_swap else 'off'}): {len(res)} contact sets (mod stabilizer, |stab|={len(stab)}), "
      f"sizes {dict(sorted(sizes.items()))}, {time.time()-t:.1f}s, stats {st}")
json.dump({'stratum': sid, 'p': s['p'], 'contact_sets': [list(map(list, S)) for S in sorted(res)],
           'stats': st, 'seconds': time.time() - t, 'frame_swap': use_swap},
          open(os.path.join(root, 'results', f'stageA_{sid}.json'), 'w'))
