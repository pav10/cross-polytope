"""Lattices in R^4 with the l1 norm. Conventions (LOAD-BEARING, see docs/STATUS.md §0):

  * A lattice is Z^4 Y, the ROWS of Y are a basis; a coefficient row q maps to v = qY.
  * Admissible  <=>  ||v||_1 >= 2 for every nonzero v.   Normalisation: lambda_1 = 2.
  * det = |det Y| is the covolume.  density = (2/3)/det.
  * The target constant is D* = det(Lambda_zeta) = 0.80822000090353537...,
    NOT the density 0.8248579173... (confusing the two has wasted cluster runs before).

Two families of functions:
  exact_*   : Fraction arithmetic, certify.  Coefficient box from the rigorous bound
              |q_j| <= R * max_i |(Y^{-1})_{ij}|   (q = v Y^{-1}).
  fast_*    : numpy, same box, for searches only. NEVER use a fixed coefficient box
              (DEAD_ENDS #3: near-singular lattices hide short vectors with large coefficients).
"""
from fractions import Fraction as F
from itertools import product
from math import floor
import numpy as np
from .exact import inv, det, vecmat

def l1(v):
    return sum(abs(x) for x in v)

def canon(q):
    """Representative of +-q with first nonzero entry positive."""
    for x in q:
        if x > 0:
            return tuple(q)
        if x < 0:
            return tuple(-y for y in q)
    return tuple(q)

def coeff_bounds(Y, R):
    Yi = inv(Y)
    return [floor(R * max(abs(Yi[i][j]) for i in range(len(Y)))) for j in range(len(Y))]

def exact_short_vectors(Y, R=2, cap=12):
    """All coefficient pairs q (canonical sign) with ||qY||_1 <= R, exactly. Returns (list of (q, norm))."""
    b = coeff_bounds(Y, F(R))
    if max(b) > cap:
        raise ValueError(f"coefficient box {b} exceeds cap {cap}: lattice too degenerate for this routine")
    out = []
    for q in product(*[range(-k, k + 1) for k in b]):
        if any(q) and canon(q) == q:
            n = l1(vecmat(q, Y))
            if n <= R:
                out.append((q, n))
    return out

def exact_min_and_minimal(Y, cap=12):
    """(lambda_1, list of minimal coefficient pairs), exact.  Uses radius = l1 norm of the
    shortest basis row (an upper bound for lambda_1)."""
    R = min(l1(row) for row in Y)
    sv = exact_short_vectors(Y, R, cap)
    m = min(n for _, n in sv)
    return m, sorted(q for q, n in sv if n == m)

def is_admissible_exact(Y, cap=12):
    """True iff ||qY||_1 >= 2 for all nonzero integer q (exact)."""
    sv = exact_short_vectors(Y, 2, cap)
    return all(n >= 2 for _, n in sv)

# ---------------------------------------------------------------- fast (search only)
def fast_short_vectors(Y, R=2.0, cap=12):
    Y = np.asarray(Y, float)
    Yi = np.linalg.inv(Y)
    b = [int(np.floor(R * np.max(np.abs(Yi[:, j])) + 1e-12)) for j in range(4)]
    if max(b) > cap:
        return None
    Q = np.array([q for q in product(*[range(-k, k + 1) for k in b]) if any(q)], float)
    N = np.abs(Q @ Y).sum(axis=1)
    keep = N <= R + 1e-12
    return Q[keep], N[keep]

def fast_lambda1(Y, cap=12):
    r = fast_short_vectors(Y, float(min(np.abs(np.asarray(Y, float)).sum(axis=1))), cap)
    if r is None:
        return None
    return float(r[1].min())

# ---------------------------------------------------------------- the zeta family
def M_of_z(z):
    """M(z); B(z) = M(z)/4.  Rows have l1 norm 8 for z in [1/2, 2/3]."""
    return [[2 - 2 * z, 2 + z, -2 - z, -2],
            [2 * z - 4, -2 * z, 2, -2],
            [-2 * z, -2 * z, 2 * z - 4, 4 - 2 * z],
            [2 - 2 * z, 2 + z, 2, 2 + z]]

def B_of_z(z):
    return [[x / 4 for x in row] for row in M_of_z(z)]

def zeta_float(prec_iter=200):
    lo, hi = F(1, 2), F(2, 3)
    for _ in range(prec_iter):
        mid = (lo + hi) / 2
        if mid ** 3 - 6 * mid + 3 > 0:
            lo = mid
        else:
            hi = mid
    return lo, hi    # rational isolating interval of width 2^-200 * (1/6)

ZETA = 0.52397639708186596553
DSTAR = 0.80822000090353537301          # det of Lambda_zeta  (= (16 - 9 z + 6 z^2)/16)
DENSITY_STAR = 0.82485791730145054648   # (2/3)/DSTAR
DET_L4_OLD = F(2514935, 3111696)         # Miller's L4: NOT optimal (z = 11/21)
