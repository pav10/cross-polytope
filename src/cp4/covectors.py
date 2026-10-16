"""Covectors of the contact arrangement {c^perp : c in S} in R^4, and the LP-free column pre-filter of Stage B1.

Column view of a face datum (faces.py): column j of the sign matrix (s_cj)_c is the sign vector of
(c . Y_j)_c, Y_j = column j of the frame matrix X, i.e. a covector of the central arrangement {c^perp}.
At a lattice det X != 0, so Y_j != 0 and every column is a covector realised by a NONZERO vector.  Inside
row-wise B1 this gives a prune before any LP: each column of the partial sign matrix (the decided contacts)
must be the restriction of such a covector (PLAN §3; measured 4.5x facet-only, ~10x full on zeta's S).

enumerate_covectors(C): incremental enumeration with float LPs (max common slack), with a witness point
    per covector.  Floats only search here.
certified_covectors(C): exact certificate that the list is complete and correct: (i) every witness is
    rounded to an exact rational point of its cell (exact sign check); (ii) every zero set is a flat
    (exact rank); (iii) for every flat F the number of covectors with zero set exactly F equals the number
    of regions of the restriction to F, r_F = sum_{G >= F} |mu(F, G)| (Zaslavsky).  (i) + distinctness +
    (iii) imply completeness.  Raises AssertionError otherwise.  Returns the covectors realised by nonzero
    vectors: all except the zero covector, which is kept only if S does not span R^4.
ColumnFilter(cov): allowed(signs, ci) -> for each column j, the boolean mask over {-1, 0, 1} of the values
    sigma[ci] of the covectors sigma that agree with the decided entries signs[i][j].  Works for any set of
    decided contacts (static or dynamic order) and is invariant under the ambient signed-permutation group
    (the covector set is closed under negation), so it commutes with the orbit reduction of faces.py.
"""
import numpy as np
from itertools import combinations
from fractions import Fraction as F
from .exact import rank, nullspace
from .lp import solve


def _cone_point(A, sig):
    """Max common slack point y (|y_k| <= 1) of the relatively open cone {sign(A_i . y) = sig_i}, or None."""
    ub, eq = [], []
    for i, s in enumerate(sig):
        (eq.append(np.r_[A[i], 0.0]) if s == 0 else ub.append(np.r_[-s * A[i], 1.0]))
    c = np.zeros(A.shape[1] + 1); c[-1] = -1
    ok, x, _ = solve(c, A_ub=np.array(ub) if ub else None, b_ub=np.zeros(len(ub)) if ub else None,
                     A_eq=np.array(eq) if eq else None, b_eq=np.zeros(len(eq)) if eq else None,
                     bounds=[(-1, 1)] * A.shape[1] + [(None, 1.0)])
    if not ok or (ub and x[-1] <= 1e-9):
        return None
    return x[:-1]


def enumerate_covectors(C):
    """All covectors (incl. the zero covector) of {c^perp : c in C}, as {sign tuple: float witness}."""
    A = np.array([[float(v) for v in c] for c in C])
    A /= np.abs(A).max(axis=1, keepdims=True).clip(min=1e-300)
    cov = {(): np.zeros(A.shape[1])}
    for m in range(len(C)):
        nxt = {}
        for s in cov:
            for v in (1, 0, -1):
                y = _cone_point(A[:m + 1], s + (v,))
                if y is not None:
                    nxt[s + (v,)] = y
        cov = nxt
    return cov


def _exact_point(C, sig, y, N):
    """An exact rational point of the cell of sig near the float witness y (N = exact basis of the
    subspace {c_i . x = 0, sig_i = 0}), or None."""
    if not N:
        return [F(0)] * len(C[0]) if not any(sig) else None
    Nf = np.array([[float(v) for v in b] for b in N]).T
    coef = np.linalg.lstsq(Nf, y, rcond=None)[0]
    for den in (10 ** 4, 10 ** 8, 10 ** 12):
        q = [F(float(a)).limit_denominator(den) for a in coef]
        x = [sum(q[k] * N[k][t] for k in range(len(N))) for t in range(len(C[0]))]
        if all((lambda d: (d > 0) - (d < 0))(sum(c[t] * x[t] for t in range(len(x)))) == s for c, s in zip(C, sig)):
            return x
    return None


def flats(C):
    """Flats of the matroid of C (closed index sets, exact ranks): {frozenset: rank}."""
    n, d = len(C), len(C[0])
    rows = [[F(v) for v in c] for c in C]
    rk = lambda T: rank([rows[i] for i in T]) if T else 0
    out = {}
    def closure(T):
        r = rk(T)
        return frozenset(i for i in range(n) if i in T or rk(list(T) + [i]) == r), r
    for k in range(d + 1):
        for T in combinations(range(n), k):
            if rk(list(T)) == k:                       # independent sets suffice to reach every flat
                Fl, r = closure(set(T))
                out[Fl] = r
    return out


def zaslavsky_counts(fl):
    """For each flat F: r_F = number of faces of the arrangement with zero set exactly F."""
    out = {}
    for Fl in fl:
        up = sorted((G for G in fl if Fl <= G), key=lambda G: fl[G])
        mu = {Fl: 1}
        for G in up[1:]:
            mu[G] = -sum(m for H, m in mu.items() if H < G)
        out[Fl] = sum(abs(v) for v in mu.values())
    return out


def certified_covectors(C, cov=None):
    """Exactly certified covector list of {c^perp : c in C} (see module docstring), restricted to the
    covectors realised by nonzero vectors.  Returns (list of sign tuples, certificate summary)."""
    C = [tuple(F(v) for v in c) for c in C]
    n, d = len(C), len(C[0])
    if cov is None:
        cov = enumerate_covectors(C)
    fl = flats(C)
    r = zaslavsky_counts(fl)
    count = {}
    for sig, y in cov.items():
        Z = frozenset(i for i in range(n) if sig[i] == 0)
        assert Z in fl, f"zero set {sorted(Z)} of {sig} is not a flat"
        N = nullspace([list(C[i]) for i in sorted(Z)], d)
        assert _exact_point(C, sig, y, N) is not None, f"no exact point found in the cell of {sig}"
        count[Z] = count.get(Z, 0) + 1
    bad = [sorted(Fl) for Fl in fl if count.get(Fl, 0) != r[Fl]]
    assert not bad, f"face counts differ from Zaslavsky's r_F on flats {bad[:5]} ({len(bad)} flats)"
    nonzero = [sig for sig in cov if fl[frozenset(i for i in range(n) if sig[i] == 0)] < d]
    by_rank = [sum(1 for Fl in fl if fl[Fl] == k) for k in range(d + 1)]
    summary = {'covectors': len(cov), 'nonzero_realised': len(nonzero), 'topes': sum(1 for s in cov if all(s)),
               'flats_by_rank': by_rank, 'faces_zaslavsky': sum(r.values())}
    return sorted(nonzero), summary


class ColumnFilter:
    """LP-free column test of a partial face datum against a list of covectors (sign tuples over C)."""

    def __init__(self, cov):
        self.M = np.array(cov, dtype=np.int8).reshape(len(cov), -1)
        self.ncov = len(cov)

    def allowed(self, signs, ci):
        """(4, 3) boolean array: allowed[j, v + 1] iff some covector agrees with signs[i][j] for every
        decided contact i and has value v at contact ci."""
        D = sorted(signs)
        out = np.zeros((4, 3), dtype=bool)
        V = np.array([signs[i] for i in D], dtype=np.int8).reshape(len(D), 4)
        MD = self.M[:, D]
        for j in range(4):
            rows = np.all(MD == V[:, j], axis=1) if D else np.ones(self.ncov, dtype=bool)
            out[j, np.unique(self.M[rows, ci]) + 1] = True
        return out

    @staticmethod
    def passes(allowed, s):
        return all(allowed[j, s[j] + 1] for j in range(4))
