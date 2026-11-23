"""Exact check of docs/COMMON_DESCENT.md §4: the common test set J of all spanning-tree bases.

For each of the 125 labelled spanning trees on five vertices, the roots e_i - e_j of its edges form a Z-basis of
L = {b in Z^5 : sum b_i = 0}.  Transport the single-basis test set of THEORY 1.7 for Gamma = Z^4 (the 1 156 non-basis
pairs of verify/verify_audit.py's finite_tests(), plus the 4 basis vectors) to L through that basis, and intersect
over all trees.  The Wolff set is invariant under signed permutations, so edge orientation and order do not matter.
Expected: 45 pairs up to sign - the 10 roots (-1,0,0,0,1), 15 of type (-1,-1,0,1,1) (disjoint-edge midpoints) and
20 of type (-2,-1,1,1,1) - so 35 inequalities beyond the ten prescribed distances decide full admissibility of an
equilateral five-point configuration.  Run from any directory: python3 verify/verify_common_descent.py
(prints COMMON DESCENT PASS)."""
import sys
from itertools import combinations, product
from collections import Counter
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'verify'))
from verify_audit import finite_tests          # read-only external checker; only its generator is used

T1 = [tuple(q) for q in finite_tests()[0]]
assert len(T1) == 1156
BASE = [tuple(int(i == j) for j in range(4)) for i in range(4)]
W = set()
for q in T1 + BASE:
    W.add(q); W.add(tuple(-x for x in q))


def trees(n=5):
    """All labelled spanning trees on n vertices (edge lists), via all (n-1)-edge subsets that connect."""
    E = list(combinations(range(n), 2))
    for es in combinations(E, n - 1):
        parent = list(range(n))
        def f(x):
            while parent[x] != x:
                x = parent[x]
            return x
        ok = True
        for i, j in es:
            a, b = f(i), f(j)
            if a == b:
                ok = False; break
            parent[a] = b
        if ok:
            yield es


def canon(b):
    return max(b, tuple(-x for x in b))


J = None
count = 0
for es in trees():
    count += 1
    TR = set()
    for q in W:
        b = [0] * 5
        for c, (i, j) in zip(q, es):
            b[i] += c; b[j] -= c
        TR.add(canon(tuple(b)))
    J = TR if J is None else J & TR
assert count == 125
shape = Counter(tuple(sorted(abs(x) for x in b)) for b in J)
print(f'{count} spanning trees; |J| = {len(J)} pairs up to sign; absolute patterns {dict(shape)}')
assert len(J) == 45
assert shape == Counter({(0, 0, 0, 1, 1): 10, (0, 1, 1, 1, 1): 15, (1, 1, 1, 1, 2): 20})
print('COMMON DESCENT PASS')
