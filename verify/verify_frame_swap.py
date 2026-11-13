"""Exact check of docs/FRAME_SWAP.md §3-4: single-vector frame swaps for the candidate pools of the remaining
strata.  For a candidate c with a coordinate +-1 (numerator +-p), X' = X with row i replaced by cX is a max-det
frame if c is minimal (Lemma 1); its frame lattice is Gamma T^-1.  Classify: impossible (Lee weight < p, i.e. a
Gamma' point of frame norm < 1), swap into an excluded stratum, swap into another remaining stratum, or only
back into the same stratum.  Run from any directory: python3 verify/verify_frame_swap.py  (prints FRAME SWAP PASS)."""
import json, sys
from fractions import Fraction as F
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from cp4.strata import canon_group, group, EXCLUDED, REMAINING

strata = json.loads((ROOT / 'data/strata.json').read_text())
byid = {s['id']: s for s in strata}
CANON = {}
for s in strata:
    CANON.setdefault((s['p'], canon_group([tuple(r) for r in s['residues']], s['p'])), s['id'])


def inv(M):
    n = len(M); A = [row[:] + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0); A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]; A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fct = A[r][c]; A[r] = [a - fct * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def frame_stratum(sid, Q):
    """Q: 4 numerator rows (over p) with |det| = 1 in frame coordinates.  Returns 'impossible' or a stratum id."""
    st = byid[sid]; p = st['p']
    Ti = inv([[F(x, p) for x in q] for q in Q])
    gens = [tuple(F(x, p) for x in r) for r in st['residues'] if any(r)] + \
           [tuple(F(int(r == k)) for k in range(4)) for r in range(4)]
    imgs = []
    for g in gens:
        w = [sum(g[k] * Ti[k][j] for k in range(4)) for j in range(4)]
        assert all((x * p).denominator == 1 for x in w)
        imgs.append(tuple(int(x * p) % p for x in w))
    G = group(imgs, p)                      # None iff some nonzero element has Lee weight < p
    if G is None:
        return 'impossible'
    return CANON[(p, canon_group(G, p))]


def classify(sid):
    st = byid[sid]; p = st['p']; out = {'impossible': 0, 'excluded': 0, 'cross': 0, 'self': 0, 'none': 0}
    targets = set()
    for c in map(tuple, st['candidates']):
        res = []
        for i in range(4):
            if abs(c[i]) == p:
                Q = [tuple(p * int(r == j) for j in range(4)) for r in range(4)]; Q[i] = c
                res.append(frame_stratum(sid, Q))
        if not res: out['none'] += 1
        elif 'impossible' in res: out['impossible'] += 1
        elif any(r in EXCLUDED for r in res): out['excluded'] += 1
        elif any(r != sid for r in res): out['cross'] += 1; targets |= {r for r in res if r != sid}
        else: out['self'] += 1
    return out, sorted(targets)


EXPECTED = {1: (58, 2, 0, 24), 2: (52, 4, 0, 20), 3: (48, 8, 0, 28), 4: (64, 9, 0, 18), 5: (56, 9, 0, 21),
            8: (84, 12, 32, 24), 9: (82, 10, 18, 32), 10: (74, 10, 22, 16), 11: (68, 14, 15, 9),
            14: (72, 16, 18, 0), 15: (98, 21, 24, 8), 17: (92, 27, 12, 7)}


def check():
    assert all(len(set(map(tuple, byid[s]['candidates']))) == len(byid[s]['candidates']) for s in REMAINING)
    for sid in REMAINING:
        out, targets = classify(sid)
        if sid == 0:
            assert out['impossible'] == out['excluded'] == out['cross'] == 0
            continue
        pool = len(byid[sid]['candidates'])
        assert (pool, out['impossible'], out['excluded'], out['cross']) == EXPECTED[sid], (sid, out)
        assert all(t in REMAINING for t in targets)
        print(f"stratum {sid:2d}: pool {pool}, impossible {out['impossible']}, into excluded {out['excluded']}, "
              f"cross {out['cross']} -> {targets}")
    print("FRAME SWAP PASS")


if __name__ == '__main__':
    check()
