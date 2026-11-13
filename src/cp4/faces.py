"""Stage B (part 1): face data for a fixed contact set S, and the resulting affine families.

Unknown: the frame matrix X (rows = the 4 frame minimal vectors in ambient coordinates),
variables x[4*i+j] = X[i][j].  For a contact c (frame coordinates) put (cX)_j = sum_i c_i X[i][j].

Face datum = for each contact c a sign vector s_c in {-1,0,1}^4 minus {0} (the face of the l1 sphere
of radius 2 containing cX in its RELATIVE INTERIOR):
    sum_j s_cj (cX)_j = 2,      s_cj (cX)_j > 0  (s_cj != 0),      (cX)_j = 0  (s_cj = 0).
A critical lattice lies in the relative interior of exactly one face datum (for its own S), and
is then a critical point of det on the affine set A cut out by the equalities  [STATUS.md §B].

Branching: contacts in a fixed order; for each, the admissible sign vectors are reduced modulo the
subgroup of the ambient signed-permutation group (acting on COLUMNS of X) that fixes all signs
chosen so far (stabilizer chain => complete and non-redundant).  Pruning: LP feasibility of the
strict system (maximise a common slack t; feasible iff t > 0).  Leaves with dim A > 6 are
discarded (perfection, Thm 5.1).

Options of enumerate_face_data (defaults = the original search): static order (dynamic=False, no closure
LPs); signvecs='facet' (facet interiors only: a fast validation pass, not complete); determine=True
(determine-then-check at equation rank 16, PLAN §3); checkpoint/resume by depth-3 subtrees; prefilter=True
(LP-free covector column pre-filter before each LP, cp4.covectors: every column of the partial sign matrix
must extend to a covector of {c^perp : c in S} realised by a nonzero vector -- sound since det X != 0 --
with the covector list certified exactly first).  estimate_face_data gives Knuth tree-size estimates on the
same code path.  LPs go through cp4.lp.  frame_group (static order only): isomorph rejection under the frame
symmetries G_S (signed permutations P of frame coordinates with Gamma P = Gamma and S P = +-S, X -> PX; see
frame_stabilizer); the four frame vectors are decided first and a depth-4 node is kept iff the AMB-orbit of its
frame-row sign matrix contains the lexicographic minimum of its G_S x AMB orbit.
"""
import numpy as np
from itertools import product, permutations
from fractions import Fraction as F
from .exact import solve_affine, rank
from .lp import solve as lp_solve

AMB = [(pm, sg) for pm in permutations(range(4)) for sg in product((1, -1), repeat=4)]
SIGNVECS = [s for s in product((-1, 0, 1), repeat=4) if any(s)]

def act_sign(g, s):
    """Column action: X -> X Q with Q signed permutation; the sign vector transforms the same way."""
    pm, sg = g
    return tuple(sg[k] * s[pm[k]] for k in range(4))

def orbit_reps(group, choices):
    seen = set(); reps = []
    for s in choices:
        if s in seen:
            continue
        reps.append(s)
        for g in group:
            seen.add(act_sign(g, s))
    return reps

def _rows_for(c, s):
    """Equalities and strict inequalities for contact c (Fractions) with sign vector s."""
    eq = []; strict = []
    # (cX)_j as a 16-vector
    cols = []
    for j in range(4):
        r = [F(0)] * 16
        for i in range(4):
            r[4 * i + j] = F(c[i])
        cols.append(r)
    eq.append(([sum(s[j] * cols[j][k] for j in range(4)) for k in range(16)], F(2)))
    for j in range(4):
        if s[j] == 0:
            eq.append((cols[j], F(0)))
        else:
            strict.append([s[j] * x for x in cols[j]])
    return eq, strict

def lp_feasible(eqs, stricts):
    """Max t s.t. A_eq x = b, S x >= t, |x_k| <= 2 (frame rows have l1 norm 2), t <= 1."""
    n = 16
    A_eq = np.array([[float(v) for v in r] + [0.0] for r, _ in eqs]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub = np.array([[-float(v) for v in r] + [1.0] for r in stricts]) if stricts else None
    b_ub = np.zeros(len(stricts)) if stricts else None
    c = np.zeros(n + 1); c[-1] = -1.0
    ok, x, _ = lp_solve(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * n + [(None, 1.0)])
    if not ok:
        return False, None
    return x[-1] > 1e-9, x[:n]

def _closure_ranges(eqs, stricts, row):
    """min and max of row.x over the closure of the current polyhedron (|x|<=2)."""
    n = 16
    A_eq = np.array([[float(v) for v in r] for r, _ in eqs]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub = np.array([[-float(v) for v in r] for r in stricts]) if stricts else None
    b_ub = np.zeros(len(stricts)) if stricts else None
    out = []
    for sgn in (1.0, -1.0):
        ok, _, fun = lp_solve(sgn * np.array([float(v) for v in row]), A_ub=A_ub, b_ub=b_ub,
                              A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * n)
        if not ok:
            return None
        out.append(sgn * fun)
    return out[0], out[1]          # (min, max)

def allowed_signs(c, eqs, stricts, eps=1e-9):
    """Sign vectors for contact c compatible (coordinatewise) with the current polyhedron.
    Conservative: uses the closure, so never discards a feasible sign vector."""
    per = []
    for j in range(4):
        row = [F(0)] * 16
        for i in range(4):
            row[4 * i + j] = F(c[i])
        r = _closure_ranges(eqs, stricts, row)
        if r is None:
            return []
        lo, hi = r
        a = []
        if lo < -eps: a.append(-1)
        if lo <= eps and hi >= -eps: a.append(0)
        if hi > eps: a.append(1)
        per.append(a)
    return [s for s in product(*per) if any(s)]

FACETVECS = [s for s in SIGNVECS if all(s)]     # interiors of the 16 facets (ambient-group invariant)

def _children(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter=None):
    """One branching step: (ci, [(sign vector, subgroup, eqs, stricts)]), or None if all contacts are decided.
    Dynamic order: the undecided contact with the fewest coordinatewise-allowed sign vectors (closure
    LPs, conservative); static order: the first undecided contact.  colfilter (cp4.covectors.ColumnFilter):
    candidates whose partial columns extend to no admissible covector are dropped before their LP; the
    verdict is the same on the whole orbit of the candidate under the subgroup fixing the decided signs."""
    undecided = [i for i in range(len(C)) if i not in signs]
    if not undecided:
        return None
    if dynamic:
        best = None
        for i in undecided:
            al = [s for s in allowed_signs(C[i], eqs, stricts) if vecs is SIGNVECS or all(s)]
            if best is None or len(al) < len(best[1]):
                best = (i, al)
            if len(al) <= 1:
                break
        ci, choices = best
    else:
        ci = undecided[0]; choices = vecs
    # orbit reduction must be over the FULL set of sign vectors, restricted afterwards:
    chs = set(choices)
    allowed = colfilter.allowed(signs, ci) if colfilter is not None else None
    kids = []
    for s in orbit_reps(group, vecs):
        # some image of s under group must be allowed; use that image (same orbit)
        img = next((act_sign(g, s) for g in group if act_sign(g, s) in chs), None)
        if img is None:
            continue
        if allowed is not None and not all(allowed[j, img[j] + 1] for j in range(4)):
            stats['prefilter_dead'] = stats.get('prefilter_dead', 0) + 1; continue
        e, st = _rows_for(C[ci], img)
        ok, _ = lp_feasible(eqs + e, stricts + st)
        if not ok:
            stats['lp_dead'] += 1; continue
        kids.append((img, [g for g in group if act_sign(g, img) == img], eqs + e, stricts + st))
    return ci, kids

def _read_off(C, eqs, signs, vecs):
    """Determine-then-check (PLAN §3): if the face equations have rank 16, X is unique; every undecided
    contact must then satisfy ||cX||_1 = 2 exactly, and its sign vector is read off.  Returns
    'branch' (rank < 16), None (dead) or the leaf (exact, dim 0)."""
    if np.linalg.matrix_rank(np.array([[float(v) for v in r] for r, _ in eqs])) < 16:
        return 'branch'
    sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
    if sol is None or sol[1]:
        return 'branch' if sol else None       # float rank misjudged: keep branching
    x0 = sol[0]
    full = {}
    for i, c in enumerate(C):
        v = [sum(c[k] * x0[4 * k + j] for k in range(4)) for j in range(4)]
        s = tuple((t > 0) - (t < 0) for t in v)
        if sum(abs(t) for t in v) != 2 or (i in signs and s != signs[i]) or (vecs is not SIGNVECS and not all(s)):
            return None
        full[i] = s
    return {'signs': full, 'x0': x0, 'basis': [], 'dim': 0}

def _column_filter(C, stats):
    from .covectors import certified_covectors, ColumnFilter
    cov, cert = certified_covectors(C)
    stats['covector_certificate'] = cert
    return ColumnFilter(cov)

def _canon_pm(v):
    return v if next(x for x in v if x) > 0 else tuple(-x for x in v)

def act_frame(g, c):
    """c -> cP for the frame signed permutation g = (pm, sg) acting by X -> PX, (PX)_i = sg_i X_{pm(i)}."""
    pm, sg = g
    out = [0] * 4
    for i in range(4):
        out[pm[i]] = sg[i] * c[i]
    return tuple(out)

def frame_stabilizer(S, residues, p):
    """Signed permutations g of frame coordinates with Gamma P = Gamma (residues mod p) and S P = +-S."""
    R = {tuple(x % p for x in r) for r in residues}
    Sset = {_canon_pm(tuple(c)) for c in S}
    out = []
    for g in AMB:                                   # the same 384 signed permutations, now acting on rows
        if {tuple(x % p for x in act_frame(g, r)) for r in R} == R and {_canon_pm(act_frame(g, c)) for c in S} == Sset:
            out.append(g)
    return out

def contact_action(S, g):
    """For each contact i: (j, eps) with S[i] P = eps * S[j]."""
    idx = {tuple(c): j for j, c in enumerate(S)}
    out = []
    for c in S:
        d = act_frame(g, tuple(c))
        out.append((idx[d], 1) if d in idx else (idx[tuple(-x for x in d)], -1))
    return out

def _lexmin_amb(Fm):
    return min(tuple(act_sign(q, r) for r in Fm) for q in AMB)

def _frame_canonical(Fm, group, cache):
    """Is the AMB-orbit of the frame-row sign matrix Fm the one holding the lex-minimum of its group x AMB orbit?"""
    if Fm not in cache:
        own = _lexmin_amb(Fm)
        ok = True
        for pm, sg in group:
            if _lexmin_amb(tuple(tuple(sg[i] * x for x in Fm[pm[i]]) for i in range(4))) < own:
                ok = False; break
        cache[Fm] = ok
    return cache[Fm]

def _frame_order(S, p):
    """Internal contact order with the frame vectors p*e_1..p*e_4 first (in coordinate order)."""
    pos = {}
    for i, c in enumerate(S):
        nz = [k for k in range(4) if c[k]]
        if len(nz) == 1 and abs(c[nz[0]]) == p:
            pos[nz[0]] = i
    assert len(pos) == 4, 'frame symmetry needs the four frame vectors in S'
    first = [pos[k] for k in range(4)]
    return first + [i for i in range(len(S)) if i not in first]

def _frame_filter(kids, ci, signs, group, cache, stats):
    """At the step that decides the 4th frame vector (internal indices 0..3), keep only canonical children."""
    if group is None or len(signs) != 3 or ci > 3:
        return kids
    keep = []
    for kid in kids:
        sg = dict(signs); sg[ci] = kid[0]
        if _frame_canonical(tuple(sg[i] for i in range(4)), group, cache):
            keep.append(kid)
        else:
            stats['frame_sym_dead'] = stats.get('frame_sym_dead', 0) + 1
    return keep

def enumerate_face_data(S, p, max_dim=6, log_every=None, time_limit=None, dynamic=True, signvecs='all',
                        determine=False, resume=None, prefix_depth=3, on_heartbeat=None, prefilter=False,
                        frame_group=None):
    """S: list of contacts as numerator tuples over p.  Returns (leaves, stats).
    leaves: dict(signs, x0, basis (affine family), dim, prefix).
    signvecs='facet' restricts every contact to a facet interior (a validation pass, NOT complete).
    determine=True applies determine-then-check at equation rank 16.
    prefilter=True applies the certified covector column pre-filter before each LP (stats['prefilter_dead']).
    Checkpointing: subtrees below each depth-`prefix_depth` path are recorded in stats['done_prefixes'] when
    finished; resume={'done_prefixes', 'leaves'} (a previous result) skips them and keeps only their leaves.
    on_heartbeat(leaves, stats) is called every log_every nodes (e.g. to dump a checkpoint).
    frame_group (from frame_stabilizer; static order only): isomorph rejection under the frame symmetries; leaves
    are then complete up to frame_group x AMB, with 'signs' in the caller's contact indices."""
    import time
    t0 = time.time()
    order = list(range(len(S)))
    if frame_group is not None:
        assert not dynamic, 'frame_group needs the static order'
        order = _frame_order(S, p)
    C = [tuple(F(x, p) for x in S[i]) for i in order]
    back = lambda sg: {order[i]: v for i, v in sg.items()}
    fcache = {}
    vecs = FACETVECS if signvecs == 'facet' else SIGNVECS
    done = set(map(tuple, resume['done_prefixes'])) if resume else set()
    leaves = [L for L in resume['leaves'] if L.get('prefix') in done] if resume else []
    stats = {'nodes': 0, 'lp_dead': 0, 'dim_dead': 0, 'leaves': len(leaves), 'rank16': 0, 'rank16_dead': 0,
             'skipped_subtrees': 0, 'per_depth': {}, 'signvecs': signvecs, 'dynamic': dynamic,
             'determine': determine, 'prefix_depth': prefix_depth, 'prefilter': prefilter, 'prefilter_dead': 0,
             'frame_symmetry': frame_group is not None, 'frame_group_order': len(frame_group or []), 'frame_sym_dead': 0}
    colfilter = _column_filter(C, stats) if prefilter else None

    def rec(group, eqs, stricts, signs, path):
        stats['nodes'] += 1
        d = len(signs)
        stats['per_depth'][d] = stats['per_depth'].get(d, 0) + 1
        if log_every and stats['nodes'] % log_every == 0:
            stats['done_prefixes'] = sorted(done); stats['seconds'] = time.time() - t0
            print(f"  [heartbeat] {time.time()-t0:.0f}s", {k: v for k, v in stats.items() if k != 'done_prefixes'},
                  f"done_subtrees={len(done)}", flush=True)
            if on_heartbeat:
                on_heartbeat(leaves, stats)
        if time_limit and time.time() - t0 > time_limit:
            raise TimeoutError
        if len(path) == prefix_depth and tuple(path) in done:
            stats['skipped_subtrees'] += 1; return
        res = _children(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter)
        if res is None:
            sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
            if sol is None:
                return
            x0, N = sol
            if len(N) > max_dim:
                stats['dim_dead'] += 1; return
            stats['leaves'] += 1
            leaves.append({'signs': back(signs), 'x0': x0, 'basis': N, 'dim': len(N), 'prefix': tuple(path[:prefix_depth])})
            return
        ci, kids = res
        kids = _frame_filter(kids, ci, signs, frame_group, fcache, stats)
        for img, sub, e1, s1 in kids:
            path.append((ci, img))
            if determine:
                signs[ci] = img
                L = _read_off(C, e1, signs, vecs)
                del signs[ci]
                if L != 'branch':
                    stats['rank16'] += 1
                    if L is None:
                        stats['rank16_dead'] += 1
                    else:
                        stats['leaves'] += 1
                        L['signs'] = back(L['signs'])
                        L['prefix'] = tuple(path[:prefix_depth]); leaves.append(L)
                    path.pop(); continue
            signs[ci] = img
            rec(sub, e1, s1, signs, path)
            del signs[ci]
            path.pop()
        if len(path) == prefix_depth:
            done.add(tuple(path))

    try:
        rec(AMB, [], [], {}, [])
        stats['complete'] = True
    except TimeoutError:
        stats['complete'] = False
    stats['done_prefixes'] = sorted(done)
    stats['seconds'] = time.time() - t0
    return leaves, stats

def estimate_face_data(S, p, probes, rng, max_dim=6, dynamic=True, signvecs='all', determine=False, max_seconds=None,
                       prefilter=False, frame_group=None):
    """Knuth's estimator for the enumerate_face_data tree (same branching code): mean over random
    root-to-leaf probes of the tree size and of the number of leaves."""
    import time
    t0 = time.time()
    order = list(range(len(S)))
    if frame_group is not None:
        assert not dynamic, 'frame_group needs the static order'
        order = _frame_order(S, p)
    C = [tuple(F(x, p) for x in S[i]) for i in order]
    fcache = {}
    vecs = FACETVECS if signvecs == 'facet' else SIGNVECS
    stats = {'lp_dead': 0, 'prefilter_dead': 0}
    colfilter = _column_filter(C, stats) if prefilter else None
    t0 = time.time()                                   # the certificate is not part of the per-node cost
    nodes, lvs, depth = [], [], []
    for _ in range(probes):
        group, eqs, stricts, signs = AMB, [], [], {}
        w, est, leaf = 1, 1, 0
        while True:
            res = _children(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter)
            if res is None:
                sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
                leaf = w if sol and len(sol[1]) <= max_dim else 0
                break
            ci, kids = res
            kids = _frame_filter(kids, ci, signs, frame_group, fcache, stats)
            if not kids:
                break
            w *= len(kids); est += w
            img, group, eqs, stricts = rng.choice(kids)
            signs = dict(signs); signs[ci] = img
            if determine:
                L = _read_off(C, eqs, signs, vecs)
                if L != 'branch':
                    leaf = w if L else 0
                    break
        nodes.append(est); lvs.append(leaf); depth.append(len(signs))
        if max_seconds and time.time() - t0 > max_seconds:
            break
    n = sum(d + 1 for d in depth)
    return {'probes': len(nodes), 'nodes_mean': float(np.mean(nodes)), 'nodes_median': float(np.median(nodes)),
            'nodes_max': float(np.max(nodes)), 'leaves_mean': float(np.mean(lvs)), 'depth_mean': float(np.mean(depth)),
            'seconds': time.time() - t0, 's_per_node': (time.time() - t0) / n,
            'signvecs': signvecs, 'dynamic': dynamic, 'determine': determine, 'prefilter': prefilter,
            'frame_symmetry': frame_group is not None, 'lp_dead': stats['lp_dead'],
            'prefilter_dead': stats['prefilter_dead'], 'frame_sym_dead': stats.get('frame_sym_dead', 0)}
