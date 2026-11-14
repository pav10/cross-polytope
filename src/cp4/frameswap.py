"""Frame swaps (docs/FRAME_SWAP.md): exact classification of the alternative max-det frames QX of a contact set.

If Q is four contacts (numerator rows over p, frame coordinates) with |det Q| = p^4, i.e. |det| = 1 in frame
coordinates, then QX is again a max-det frame (Lemma 1) with frame lattice Gamma' = Gamma Q^{-1} of the same
index p.  FrameSwapRule(strata, sid).ok(Q) is False iff Gamma' is
  (a) impossible (a nonzero residue of Lee weight < p, i.e. a Gamma' point of frame norm < 1),
  (b) a stratum excluded by proof (cp4.strata.EXCLUDED: then det >= 81/100), or
  (c) a remaining stratum earlier than sid in the canonical order sorted(REMAINING) (the lattice is assigned to,
      and found by, that stratum's search).
All arithmetic is exact (Fractions).  The rules are hereditary (Corollary 2).  verify/verify_frame_swap.py checks
the single-swap counts independently and that this module agrees with them.
"""
from fractions import Fraction as F
from .strata import canon_group, group, EXCLUDED, REMAINING


def inv(M):
    """Exact inverse of a square matrix of Fractions (Gauss-Jordan)."""
    n = len(M)
    A = [list(row) + [F(int(i == j)) for j in range(n)] for i, row in enumerate(M)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [x / pv for x in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                fct = A[r][c]
                A[r] = [a - fct * b for a, b in zip(A[r], A[c])]
    return [row[n:] for row in A]


def canon_index(strata):
    """(p, canonical residue group) -> stratum id."""
    out = {}
    for s in strata:
        out.setdefault((s['p'], canon_group([tuple(r) for r in s['residues']], s['p'])), s['id'])
    return out


def frame_stratum(st, Q, canon):
    """Q: 4 numerator rows (over p = st['p']) with |det| = p^4.  Returns 'impossible' or the stratum id of
    Gamma' = Gamma Q^{-1} (Gamma = st's frame lattice)."""
    p = st['p']
    Ti = inv([[F(x, p) for x in q] for q in Q])
    gens = [tuple(F(x, p) for x in r) for r in st['residues'] if any(r)] + \
           [tuple(F(int(r == k)) for k in range(4)) for r in range(4)]
    imgs = []
    for g in gens:
        w = [sum(g[k] * Ti[k][j] for k in range(4)) for j in range(4)]
        assert all((x * p).denominator == 1 for x in w), 'Gamma Q^-1 must lie in (1/p) Z^4'
        imgs.append(tuple(int(x * p) % p for x in w))
    G = group(imgs, p)                    # None iff some nonzero element has Lee weight < p
    if G is None:
        return 'impossible'
    return canon[(p, canon_group(G, p))]


class FrameSwapRule:
    """The frame-swap pruning rules (a)-(c) for the Stage-A search of stratum sid (cp4.contacts)."""

    def __init__(self, strata, sid):
        self.st = next(s for s in strata if s['id'] == sid)
        self.sid, self.p = sid, self.st['p']
        self.canon = canon_index(strata)
        self.order = {s: k for k, s in enumerate(sorted(REMAINING))}
        self.cache = {}

    def ok(self, Q):
        """False iff the frame Q X is impossible, in an excluded stratum, or in an earlier remaining stratum."""
        key = frozenset(tuple(q) for q in Q)
        if key not in self.cache:
            s2 = frame_stratum(self.st, [tuple(q) for q in Q], self.canon)
            self.cache[key] = (s2 != 'impossible' and s2 not in EXCLUDED
                               and self.order.get(s2, -1) >= self.order[self.sid])
        return self.cache[key]

    def single_swaps(self, c):
        """The frames obtained from the identity by replacing row i with c, for each coordinate with |c_i| = p."""
        p = self.p
        out = []
        for i in range(4):
            if abs(c[i]) == p:
                Q = [tuple(p * int(r == j) for j in range(4)) for r in range(4)]
                Q[i] = tuple(c)
                out.append(Q)
        return out

    def pool_excluded(self, cand):
        """Indices of candidates that no critical lattice of this stratum's canonical assignment has as a
        minimal vector: some single swap with the (forced) frame fails ok()."""
        return {k for k, c in enumerate(cand) if any(not self.ok(Q) for Q in self.single_swaps(c))}
