"""Merged-search driver (src/cp4/merged.py, PLAN §4).
Usage: python3 scripts/merged_run.py <stratum_id> [--estimate K] [--seed N] [--time-limit S] [--node-limit N]
                                      [--adm-prune] [--heartbeat N]
--estimate K: Knuth tree-size estimate from K random probes (no search), printed and saved to
              results/merged_<sid>_estimate.json.  Otherwise: depth-first search with heartbeat; partial
              results are rewritten atomically to results/merged_<sid>.partial.json at every heartbeat and the
              final output goes to results/merged_<sid>.json."""
import sys, os, json, time, random, argparse
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.strata import stabilizer
from cp4.merged import Merged
root = os.path.join(os.path.dirname(__file__), '..')
ap = argparse.ArgumentParser()
ap.add_argument('sid', type=int)
ap.add_argument('--estimate', type=int)
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--time-limit', type=float)
ap.add_argument('--node-limit', type=int)
ap.add_argument('--adm-prune', action='store_true')
ap.add_argument('--heartbeat', type=int, default=2000)
a = ap.parse_args()
s = json.load(open(os.path.join(root, 'data', 'strata.json')))[a.sid]
M = Merged(s, stabilizer([tuple(r) for r in s['residues']], s['p']), adm_prune=a.adm_prune)
print(f"stratum {a.sid} (p={s['p']}): {len(M.cand)} candidates, {len(M.taus)} root frame patterns", flush=True)

def dump(path, extra=None):
    out = {'stratum': a.sid, 'p': s['p'], 'adm_prune': a.adm_prune, 'stats': M.stats,
           'leaves': M.leaves, 'leaves0': M.leaves0, **(extra or {})}
    tmp = path + '.tmp'
    json.dump(out, open(tmp, 'w'), default=str)
    os.replace(tmp, path)

if a.estimate:
    e = M.estimate(a.estimate, random.Random(a.seed), max_seconds=a.time_limit)
    print(f"Knuth estimate ({e['probes']} probes, {e['seconds']:.0f}s): nodes mean {e['nodes_mean']:.3g}, "
          f"median {e['nodes_median']:.3g}, max {e['nodes_max']:.3g}; leaves mean {e['leaves_mean']:.3g}; "
          f"LP {M.stats['lp_s'] / max(1, M.stats['lp']) * 1000:.2f} ms/call", flush=True)
    dump(os.path.join(root, 'results', f'merged_{a.sid}_estimate.json'), {'estimate': e, 'seed': a.seed})
    sys.exit(0)

def heartbeat(m, t):
    st = m.stats
    print(f"  [heartbeat] {t:.0f}s nodes={st['nodes']} ({st['nodes'] / t:.0f}/s) roots {st['roots_done']}/{st['roots']} "
          f"leaves={st['leaves']} dim0={st['leaves0']} band0={st['band0']} lp_share={st['lp_s'] / t:.0%} "
          f"dead lp/det/int/size/excl={st['dead_lp']}/{st['dead_det']}/{st['dead_interior']}/{st['dead_size']}/{st['dead_excl']}",
          flush=True)
    dump(os.path.join(root, 'results', f'merged_{a.sid}.partial.json'))

st = M.run(time_limit=a.time_limit, node_limit=a.node_limit, log_every=a.heartbeat, on_heartbeat=heartbeat)
print(f"done: complete={st['complete']} nodes={st['nodes']} {st['seconds']:.0f}s leaves={st['leaves']} "
      f"dim0={st['leaves0']} canonical S={len(M.canonical_sets())}", flush=True)
dump(os.path.join(root, 'results', f'merged_{a.sid}.json'))
