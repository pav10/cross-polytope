"""Stage A: enumerate candidate contact sets S inside one stratum Gamma.

Frame coordinates: c = num/p with num an integer 4-vector; Gamma given by its residues.
S = the set of ALL minimal vectors of the critical lattice, in frame coordinates, one per +- pair.

Necessary conditions used (all PROVED, docs/STATUS.md §A):
 (A1) S ⊇ forced = {c in Gamma : ||c||_1 = 1}   (includes the frame e_1..e_4).
 (A2) max-det frame: |det(c1,c2,c3,c4)| <= 1 for all c_i in S          [Thm 6.1(ii)]
 (A3) hull closure: let P = conv(±S).  Gamma ∩ int P = {0}, and every nonzero
      Gamma point of P belongs to S.                                       [Lemma A.1]
 (A4) volume: vol(P) <= 16/p   (Minkowski applied to (A3))                [Cor. A.2]
 (A5) 10 <= |S| <= 20                                                      [Thm 5.1, Miller]
 (A6) frame swaps (optional, swap=cp4.frameswap.FrameSwapRule; docs/FRAME_SWAP.md): every quadruple of S with
      |det| = 1 is another max-det frame; S is discarded if that frame is impossible, lies in an excluded
      stratum, or lies in an earlier remaining stratum (canonical assignment).  Candidates whose single swap with
      the frame fails are excluded at the root.
(A2), (A3)-interior, (A4), (A6) and |S|<=20 are HEREDITARY (fail for S => fail for supersets).
(A3)-boundary is a FORCING rule: boundary Gamma points must be added.

Search: binary include/exclude branching over candidates in a fixed order, with forcing
propagation; an excluded candidate that later becomes forced kills the branch.  Results are
reduced modulo the stabilizer of Gamma in the signed permutation group.

Floats are used only for hull facets (coordinates are k/p with p <= 18, tolerance 1e-9);
every OUTPUT S is re-verified in exact arithmetic by verify_contact_set().
"""
import numpy as np
from itertools import combinations
from fractions import Fraction as F
from scipy.spatial import ConvexHull
from .exact import det as exact_det
from .strata import act

TOL = 1e-9

def _hull(points):
    P = np.array(points, float)
    P = np.vstack([P, -P])
    return ConvexHull(P)

def classify_points(S_num, cand_num, p):
    """Return (interior_idx, boundary_idx, volume) for candidates w.r.t. conv(±S)."""
    h = _hull(np.array(S_num, float) / p)
    A = h.equations[:, :4]; b = h.equations[:, 4]
    C = np.array(cand_num, float) / p
    vals = C @ A.T + b                      # <= 0 inside
    mx = vals.max(axis=1)
    interior = np.where(mx < -TOL)[0]
    boundary = np.where(np.abs(mx) <= TOL)[0]
    return interior, boundary, h.volume

def _quad_dets_ok(S_num, new, p):
    """All 4-subsets containing `new` have |det| <= p^4 (i.e. <= 1 in frame coords)."""
    lim = p ** 4 + 0.5
    if len(S_num) < 3:
        return True
    T = np.array(list(combinations(S_num, 3)), float)       # (m,3,4)
    M = np.concatenate([T, np.broadcast_to(np.array(new, float), (T.shape[0], 1, 4))], axis=1)
    d = np.abs(np.linalg.det(M))
    return bool(np.all(d <= lim))

def _swap_ok(S_num, new, p, swap):
    """(A6): every quadruple (T, new), T a triple of S, with |det| = p^4 (exact) passes swap.ok."""
    if swap is None or len(S_num) < 3:
        return True
    lim = p ** 4
    T = np.array(list(combinations(S_num, 3)), float)
    M = np.concatenate([T, np.broadcast_to(np.array(new, float), (T.shape[0], 1, 4))], axis=1)
    d = np.abs(np.linalg.det(M))
    for k in np.nonzero(np.abs(d - lim) < 0.5)[0]:
        Q = [tuple(int(x) for x in row) for row in M[k]]
        if abs(exact_det([list(q) for q in Q])) == lim and not swap.ok(Q):
            return False
    return True

def canon_set(S_num, stab, p):
    best = None
    for g in stab:
        pm, sg = g
        img = []
        for n in S_num:
            m = tuple(sg[k] * n[pm[k]] for k in range(4))
            img.append(max(m, tuple(-x for x in m)))
        key = tuple(sorted(img))
        if best is None or key < best:
            best = key
    return best

def enumerate_contact_sets(stratum, stab, max_size=20, min_size=10, node_limit=None, log=None, swap=None):
    p = stratum['p']
    cand = [tuple(c) for c in stratum['candidates']]
    idx = {c: i for i, c in enumerate(cand)}
    forced = [tuple(c) for c in stratum['forced']]
    stats = {'nodes': 0, 'dead_det': 0, 'dead_interior': 0, 'dead_volume': 0,
             'dead_size': 0, 'dead_excluded_forced': 0, 'leaves': 0, 'dead_swap': 0, 'pool_swap_excluded': 0}
    results = set()
    vol_lim = 16.0 / p + TOL

    def close(S, excluded):
        """Apply forcing until fixpoint. Return closed S or None (dead)."""
        S = list(S)
        while True:
            interior, boundary, vol = classify_points(S, cand, p)
            if len(interior):
                stats['dead_interior'] += 1; return None, None
            if vol > vol_lim:
                stats['dead_volume'] += 1; return None, None
            Sset = set(S)
            new = [cand[i] for i in boundary if cand[i] not in Sset]
            if not new:
                return S, boundary
            for c in new:
                if idx[c] in excluded:
                    stats['dead_excluded_forced'] += 1; return None, None
                if not _quad_dets_ok(S, c, p):
                    stats['dead_det'] += 1; return None, None
                if not _swap_ok(S, c, p, swap):
                    stats['dead_swap'] += 1; return None, None
                S.append(c)
                if len(S) > max_size:
                    stats['dead_size'] += 1; return None, None

    # root: forced set must itself satisfy (A2)
    S0 = []
    for c in forced:
        if not _quad_dets_ok(S0, c, p):
            return results, dict(stats, root_dead='forced set violates det<=1')
        if not _swap_ok(S0, c, p, swap):
            return results, dict(stats, root_dead='forced set fails the frame-swap rule')
        S0.append(c)
    excl0 = set(swap.pool_excluded(cand)) if swap is not None else set()
    stats['pool_swap_excluded'] = len(excl0)
    if any(idx[c] in excl0 for c in S0):
        return results, dict(stats, root_dead='a forced vector fails the frame-swap rule')
    S0, _ = close(S0, excl0)
    if S0 is None:
        return results, dict(stats, root_dead='forced set violates closure')

    def rec(S, excluded, k):
        stats['nodes'] += 1
        if log and stats['nodes'] % log == 0:
            print(f"  [heartbeat] nodes={stats['nodes']} results={len(results)} |S|={len(S)} k={k} {stats}", flush=True)
        if node_limit and stats['nodes'] > node_limit:
            raise RuntimeError('node limit')
        Sset = set(S)
        # next undecided candidate
        while k < len(cand) and (cand[k] in Sset or k in excluded):
            k += 1
        if k == len(cand):
            stats['leaves'] += 1
            if len(S) >= min_size:
                results.add(canon_set(S, stab, p))
            return
        c = cand[k]
        # branch 1: include c
        if len(S) < max_size and _quad_dets_ok(S, c, p):
            if _swap_ok(S, c, p, swap):
                S1, _ = close(S + [c], excluded)
                if S1 is not None:
                    rec(S1, excluded, k + 1)
            else:
                stats['dead_swap'] += 1
        else:
            stats['dead_det'] += 1
        # branch 2: exclude c
        rec(S, excluded | {k}, k + 1)

    rec(S0, frozenset(excl0), 0)
    return results, stats

def verify_contact_set(S_num, stratum):
    """Exact re-verification of (A1),(A2),(A5) for an output set (hull conditions are
    re-checked with exact rational facet arithmetic by verify/test_regression.py)."""
    p = stratum['p']
    S = [tuple(c) for c in S_num]
    assert set(map(tuple, stratum['forced'])) <= set(S)
    for t in combinations(S, 4):
        assert abs(exact_det([list(x) for x in t])) <= p ** 4
    assert 10 <= len(S) <= 20
    return True
