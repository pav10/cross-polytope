"""Probe of PLAN §2 option 4(a) ("A6"): does the frame-facet LP prune Stage-A contact sets?
Usage: python3 scripts/probe_A6_frame_facets.py <stratum_id> [node_limit]

A6(S): for some choice of closed facets u_i in {±1}^4 of the frame rows X_i (mod the ambient group),
    exists X with  u_i . X_i = 2,  ||c X||_1 <= 2 for all c in S (frame rows included).
Necessary for the true S of a critical lattice (its own X is feasible) and hereditary in S.
The probe runs baseline Stage A (A1-A5), captures every accepted raw set, and counts how many A6 kills.
Measured 14 Oct 2026: 0 of 200 (stratum 0), 0 of 931 (stratum 11, 30 000 nodes) => DEAD_ENDS B12."""
import sys, os, json, time, random
from itertools import product
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import cp4.contacts as contacts
from cp4.strata import stabilizer
from cp4.faces import AMB, act_sign

FACETS = list(product((1, -1), repeat=4))
EPS = 1e-7          # relaxation => an LP is declared infeasible only with margin (conservative)

def facet_pattern_reps():
    """Orbit representatives of (u_1..u_4) under the ambient group (stabilizer chain)."""
    out = []
    def rec(group, prefix):
        if len(prefix) == 4:
            out.append(tuple(prefix)); return
        seen = set()
        for u in FACETS:
            if u in seen:
                continue
            seen.update(act_sign(g, u) for g in group)
            rec([g for g in group if act_sign(g, u) == u], prefix + [u])
    rec(AMB, [])
    return out

def a6_feasible(S_num, p, u):
    rows, rhs = [], []
    for n in S_num:
        c = np.array(n, float) / p
        for v in FACETS:                                   # v . (cX) <= 2
            rows.append(np.outer(c, v).ravel()); rhs.append(2 + EPS)
    for i in range(4):                                     # u_i . X_i >= 2
        r = np.zeros(16); r[4 * i:4 * i + 4] = -np.array(u[i], float)
        rows.append(r); rhs.append(-(2 - EPS))
    res = linprog(np.zeros(16), A_ub=np.array(rows), b_ub=np.array(rhs),
                  bounds=[(-2, 2)] * 16, method='highs')
    return res.status == 0

def alive(S_num, p, reps, first_only=True):
    n = 0
    for u in reps:
        if a6_feasible(S_num, p, u):
            n += 1
            if first_only:
                break
    return n

if __name__ == '__main__':
    root = os.path.join(os.path.dirname(__file__), '..')
    sid = int(sys.argv[1]); lim = int(sys.argv[2]) if len(sys.argv) > 2 else None
    s = json.load(open(os.path.join(root, 'data', 'strata.json')))[sid]; p = s['p']
    reps = facet_pattern_reps()
    got = []
    canon = contacts.canon_set
    contacts.canon_set = lambda S, stab, p_: (got.append(tuple(S)), canon(S, stab, p_))[1]
    t = time.time()
    try:
        contacts.enumerate_contact_sets(s, stabilizer([tuple(r) for r in s['residues']], p), node_limit=lim)
        complete = True
    except RuntimeError:
        complete = False
    contacts.canon_set = canon
    sets = list(dict.fromkeys(got))
    print(f"stratum {sid} (p={p}): {len(sets)} raw sets accepted by Stage A "
          f"({'complete' if complete else f'node limit {lim}'}, {time.time()-t:.0f}s); {len(reps)} facet-pattern reps")
    t = time.time()
    killed = sum(alive(S, p, reps) == 0 for S in sets)
    print(f"A6 kills {killed} of {len(sets)} sets ({time.time()-t:.0f}s)")
    sample = random.Random(1).sample(sets, min(25, len(sets)))
    counts = sorted(alive(S, p, reps, first_only=False) for S in sample)
    print(f"alive patterns on {len(sample)} sampled sets: min {counts[0]}, "
          f"median {counts[len(counts) // 2]}, max {counts[-1]} (of {len(reps)})")
