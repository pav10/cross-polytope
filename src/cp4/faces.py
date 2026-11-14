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

class _Row(list):
    """A row of Fractions with a cached float copy, so that each row is converted to floats once, not per LP."""
    __slots__ = ('_f',)

def _fl(r):
    f = getattr(r, '_f', None)
    if f is None:
        f = np.array([float(v) for v in r])
        if isinstance(r, _Row):
            r._f = f
    return f

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
    eq.append((_Row(sum(s[j] * cols[j][k] for j in range(4)) for k in range(16)), F(2)))
    for j in range(4):
        if s[j] == 0:
            eq.append((_Row(cols[j]), F(0)))
        else:
            strict.append(_Row(s[j] * x for x in cols[j]))
    return eq, strict

SLACK_TOL = 1e-9      # strict-slack test of lp_feasible: a cell is kept iff its max common slack t > SLACK_TOL

def _ub_rows(stricts, nonstricts, extra, slack):
    """Inequality rows for the LPs: strict rows S x >= t (slack column, if `slack`), nonstrict rows N x >= 0, and
    the extra rows extra = (A, b), A x <= b (A over 17 columns, last one zero; e.g. ball_rows)."""
    blocks = []
    if stricts:
        B = -np.vstack([_fl(r) for r in stricts])
        blocks.append(np.hstack([B, np.ones((len(B), 1))]) if slack else B)
    if nonstricts:
        B = -np.vstack([_fl(r) for r in nonstricts])
        blocks.append(np.hstack([B, np.zeros((len(B), 1))]) if slack else B)
    A_ub = np.vstack(blocks) if blocks else None
    b_ub = np.zeros(len(A_ub)) if blocks else None
    if extra is not None:
        E = extra[0] if slack else extra[0][:, :16]
        A_ub = E if A_ub is None else np.vstack([A_ub, E])
        b_ub = extra[1] if b_ub is None else np.concatenate([b_ub, extra[1]])
    return A_ub, b_ub

def ball_rows(C):
    """Lemma B.5 (ball rows): -2 <= tau.(cX) <= 2 for every contact c in C (Fractions) and tau in {1} x {+-1}^3,
    i.e. ||cX||_1 <= 2, as float rows (A, b) over the 17 LP columns (x, t) for lp_feasible(extra=...).
    Rows 16i .. 16i+15 belong to contact i (the search uses only those of undecided contacts)."""
    taus = [t for t in product((1, -1), repeat=4) if t[0] == 1]
    A = []
    for c in C:
        for t in taus:
            r = np.zeros(17)
            for i in range(4):
                for j in range(4):
                    r[4 * i + j] = float(c[i]) * t[j]
            A.append(r); A.append(-r)
    A = np.array(A)
    return A, np.full(len(A), 2.0 + 1e-9)

def lp_feasible(eqs, stricts, nonstricts=(), extra=None):
    """Max t s.t. A_eq x = b, S x >= t, N x >= 0 (nonstrict rows), extra rows, |x_k| <= 2 (frame rows have
    l1 norm 2), t <= 1.  Feasible iff t > 0 (strict-slack test)."""
    n = 16
    A_eq = np.hstack([np.vstack([_fl(r) for r, _ in eqs]), np.zeros((len(eqs), 1))]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub, b_ub = _ub_rows(stricts, nonstricts, extra, True)
    c = np.zeros(n + 1); c[-1] = -1.0
    ok, x, _ = lp_solve(c, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * n + [(None, 1.0)])
    if not ok:
        return False, None
    return x[-1] > SLACK_TOL, x[:n]

def _closure_ranges(eqs, stricts, row, nonstricts=(), extra=None):
    """min and max of row.x over the closure of the current polyhedron (|x|<=2)."""
    n = 16
    A_eq = np.vstack([_fl(r) for r, _ in eqs]) if eqs else None
    b_eq = np.array([float(b) for _, b in eqs]) if eqs else None
    A_ub, b_ub = _ub_rows(stricts, nonstricts, extra, False)
    out = []
    for sgn in (1.0, -1.0):
        ok, _, fun = lp_solve(sgn * np.array([float(v) for v in row]), A_ub=A_ub, b_ub=b_ub,
                              A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * n)
        if not ok:
            return None
        out.append(sgn * fun)
    return out[0], out[1]          # (min, max)

def allowed_signs(c, eqs, stricts, eps=1e-9, nonstricts=(), extra=None):
    """Sign vectors for contact c compatible (coordinatewise) with the current polyhedron.
    Conservative: uses the closure, so never discards a feasible sign vector."""
    per = []
    for j in range(4):
        row = [F(0)] * 16
        for i in range(4):
            row[4 * i + j] = F(c[i])
        r = _closure_ranges(eqs, stricts, row, nonstricts, extra)
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

def _children(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter=None, extra=None):
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
            al = [s for s in allowed_signs(C[i], eqs, stricts, extra=extra) if vecs is SIGNVECS or all(s)]
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
        ok, _ = lp_feasible(eqs + e, stricts + st, extra=extra) if extra is not None else lp_feasible(eqs + e, stricts + st)
        if not ok:
            stats['lp_dead'] += 1; continue
        kids.append((img, [g for g in group if act_sign(g, img) == img], eqs + e, stricts + st))
    return ci, kids

def _rows_half(c, s, rho):
    """Lemma B.4: the half-open facet cell of contact c for the facet s in {+-1}^4 and the tie side rho:
    s.(cX) = 2;  s_j (cX)_j >= 0 where s_j = rho_j (zero coordinates go to rho's side), > 0 elsewhere.
    Returns (equalities, strict rows, nonstrict rows)."""
    eq, strict = _rows_for(c, s)
    non = [r for r, j in zip(strict, range(4)) if s[j] == rho[j]]
    strict = [r for r, j in zip(strict, range(4)) if s[j] != rho[j]]
    return eq, strict, non

def _children_supporting(C, ci, eqs, stricts, nons, signs, rhos, stats, colfilter=None, extra=None, eps=1e-9):
    """Branching on a NON-frame contact in supporting mode (Theorem B.1', Lemma B.4): the 16 facets with
    half-open cells.  The tie side rho is chosen per node (any choice is sound; this one avoids empty strict
    sides): from the covector column filter when it rules out one strict side, otherwise from closure-range
    LPs.  No orbit reduction (the group below the frame is taken trivial).  Returns (rho, kids) with kids
    (s, None, eqs, stricts, nons); rho is None if the node is dead."""
    c = C[ci]
    zok = {i: tuple(sg[j] == rhos[i][j] for j in range(4)) for i, sg in signs.items() if i in rhos}
    allowed = colfilter.allowed(signs, ci, zero_ok=zok) if colfilter is not None else None
    rho, pos_ok, neg_ok = [], [], []
    for j in range(4):
        if allowed is not None and not (allowed[j, 0] and allowed[j, 2]):
            p_ok, n_ok = bool(allowed[j, 2]), bool(allowed[j, 0])
            if not (p_ok or n_ok or allowed[j, 1]):
                stats['prefilter_dead'] = stats.get('prefilter_dead', 0) + 1
                return None, []
        else:
            row = [F(0)] * 16
            for i in range(4):
                row[4 * i + j] = F(c[i])
            r = _closure_ranges(eqs, stricts, row, nons, extra)
            stats['range_lps'] = stats.get('range_lps', 0) + 2
            if r is None:
                stats['lp_dead'] += 1
                return None, []
            p_ok, n_ok = r[1] > eps, r[0] < -eps
        rho.append(1 if p_ok else -1); pos_ok.append(p_ok); neg_ok.append(n_ok)
    rho = tuple(rho)
    kids = []
    for s in FACETVECS:
        if any(s[j] != rho[j] and not (pos_ok[j] if s[j] > 0 else neg_ok[j]) for j in range(4)):
            stats['side_dead'] = stats.get('side_dead', 0) + 1; continue
        e, st, no = _rows_half(c, s, rho)
        ok, _ = lp_feasible(eqs + e, stricts + st, nons + no, extra)
        if not ok:
            stats['lp_dead'] += 1; continue
        kids.append((s, None, eqs + e, stricts + st, nons + no))
    return rho, kids

def _read_off(C, eqs, signs, vecs, rhos=None):
    """Determine-then-check (PLAN §3): if the face equations have rank 16, X is unique; every undecided
    contact must then satisfy ||cX||_1 = 2 exactly, and its sign vector is read off.  Returns
    'branch' (rank < 16), None (dead) or the leaf (exact, dim 0).
    Supporting mode (rhos = tie sides of the decided non-frame contacts, a dict): decided non-frame contacts are
    checked against their half-open cells; an undecided non-frame contact only needs ||cX||_1 = 2 (any
    supporting facet will do, Theorem B.1'), and its facet is read off with zeros sent to +1."""
    if rhos is not None:
        return _read_off_supporting(C, eqs, signs, rhos)
    if len(eqs) < 16 or np.linalg.matrix_rank(np.vstack([_fl(r) for r, _ in eqs])) < 16:
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

def _read_off_supporting(C, eqs, signs, rhos):
    if len(eqs) < 16 or np.linalg.matrix_rank(np.vstack([_fl(r) for r, _ in eqs])) < 16:
        return 'branch'
    sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
    if sol is None or sol[1]:
        return 'branch' if sol else None
    x0 = sol[0]
    full = {}
    for i, c in enumerate(C):
        v = [sum(c[k] * x0[4 * k + j] for k in range(4)) for j in range(4)]
        sv = tuple((t > 0) - (t < 0) for t in v)
        if sum(abs(t) for t in v) != 2:
            return None
        if i in rhos:                       # decided non-frame contact: half-open cell
            sg, rho = signs[i], rhos[i]
            if any(sg[j] * v[j] < 0 or (v[j] == 0 and sg[j] != rho[j]) for j in range(4)):
                return None
            full[i] = sg
        elif i in signs:                    # frame row: exact face
            if sv != signs[i]:
                return None
            full[i] = sv
        else:                               # undecided non-frame contact: any supporting facet
            full[i] = tuple(x if x else 1 for x in sv)
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

def eutaxy_residual(C, signs, eqs, ineqs):
    """Leaf prune by l1-eutaxy (THEORY Thm 5.2; STATUS §3 Lemma B.6), RLT relaxation.  At a critical X with
    contacts v_q = qX and ANY supporting facets sigma_q there are lambda >= 0, sum lambda = 2, with
    sum_q lambda_q sigma_q^T v_q = I.  Put Z_q = lambda_q X: then Z_q lies in lambda_q * (the leaf's cell closure),
    i.e. r.Z_q = b lambda_q for every equation (r, b) of the leaf, g.Z_q >= 0 for every sign row g, |Z_q| <= 2 lambda_q,
    and sum_q sigma_q^T (q Z_q) = I.  Returns the minimal max-norm residual of that identity over the relaxation
    (0 if some lambda, Z satisfy it); a leaf with residual > 1e-7 contains no critical X.  C, signs: internal order,
    every contact decided; ineqs: the leaf's strict and nonstrict rows (closure)."""
    n = len(C)
    nv = 17 * n + 1
    sv = nv - 1
    Cf = [np.array([float(x) for x in c]) for c in C]
    Aeq, beq, Aub, bub = [], [], [], []
    for a in range(4):
        for b in range(4):
            r = np.zeros(nv)
            for q in range(n):
                for i in range(4):
                    r[17 * q + 1 + 4 * i + b] += signs[q][a] * Cf[q][i]
            d = 1.0 if a == b else 0.0
            r1 = r.copy(); r1[sv] = -1.0; Aub.append(r1); bub.append(d)
            r2 = -r; r2[sv] = -1.0; Aub.append(r2); bub.append(-d)
    r = np.zeros(nv)
    r[[17 * q for q in range(n)]] = 1.0
    Aeq.append(r); beq.append(2.0)
    E = [(_fl(row), float(b)) for row, b in eqs]
    G = [_fl(row) for row in ineqs]
    for q in range(n):
        z0 = 17 * q + 1
        for row, b in E:
            r = np.zeros(nv); r[z0:z0 + 16] = row; r[17 * q] = -b; Aeq.append(r); beq.append(0.0)
        for row in G:
            r = np.zeros(nv); r[z0:z0 + 16] = -row; Aub.append(r); bub.append(0.0)
        for k in range(16):
            r = np.zeros(nv); r[z0 + k] = 1.0; r[17 * q] = -2.0; Aub.append(r); bub.append(0.0)
            r = np.zeros(nv); r[z0 + k] = -1.0; r[17 * q] = -2.0; Aub.append(r); bub.append(0.0)
    bounds = []
    for q in range(n):
        bounds += [(0, 2)] + [(None, None)] * 16
    bounds.append((0, None))
    c = np.zeros(nv); c[sv] = 1.0
    ok, _, fun = lp_solve(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                          bounds=bounds)
    return fun if ok else float('inf')

EUTAXY_TOL = 1e-7

def _expand(C, group, eqs, stricts, nons, signs, rhos, ctx, stats):
    """One branching step shared by enumerate_face_data and estimate_face_data: None if every contact is decided,
    else (ci, kids) with kids (sign vector, subgroup, eqs, stricts, nonstricts, rho); rho is the tie side of a
    half-open facet cell (supporting mode, non-frame contacts) or None."""
    extra = ctx['extra']
    if extra is not None:
        # ball rows only for undecided contacts: a decided contact's cell already gives ||cY||_1 = 2
        und = [i for i in range(len(C)) if i not in signs]
        if not und:
            extra = None
        else:
            idx = np.concatenate([np.arange(16 * i, 16 * i + 16) for i in und])
            extra = (extra[0][idx], extra[1][idx])
    if ctx['supporting'] and len(signs) >= 4:
        if len(signs) == len(C):
            return None
        if len(signs) == 4 and any(all(signs[i][j] == 0 for i in range(4)) for j in range(4)):
            stats['zero_column_dead'] = stats.get('zero_column_dead', 0) + 1
            return 4, []                    # a zero column of X: det X = 0
        ci = _next_contact(C, signs, rhos, ctx)
        rho, kids = _children_supporting(C, ci, eqs, stricts, nons, signs, rhos, stats, ctx['colfilter'], extra)
        return ci, [(sv, None, e, st, no, rho) for sv, _, e, st, no in kids]
    kw = {'extra': extra} if extra is not None else {}
    res = _children(C, group, eqs, stricts, signs, ctx['dynamic'], ctx['vecs'], stats, ctx['colfilter'], **kw)
    if res is None:
        return None
    ci, kids = res
    kids = _frame_filter(kids, ci, signs, ctx['frame_group'], ctx['fcache'], stats)
    return ci, [(img, sub, e1, s1, nons, None) for img, sub, e1, s1 in kids]

def _next_contact(C, signs, rhos, ctx):
    """Supporting mode, non-frame contacts: the first undecided contact (contact_order='static'), or the undecided
    contact with the fewest facets still allowed by the covector column filter (contact_order='filter'; LP-free).
    Any order is sound: the children of a node partition its cell (Lemma B.4)."""
    und = [i for i in range(len(C)) if i not in signs]
    if ctx['contact_order'] != 'filter' or ctx['colfilter'] is None:
        return und[0]
    zok = {i: tuple(sg[j] == rhos[i][j] for j in range(4)) for i, sg in signs.items() if i in rhos}
    best, bi = None, und[0]
    for i in und:
        a = ctx['colfilter'].allowed(signs, i, zero_ok=zok)
        sc = 1
        for j in range(4):
            sc *= max(1, int(a[j, 0]) + int(a[j, 2]))
        if best is None or sc < best:
            best, bi = sc, i
            if sc == 1:
                break
    return bi

def _setup(S, p, dynamic, signvecs, prefilter, frame_group, supporting, ball, stats, contact_order='static',
           eutaxy=False):
    order = list(range(len(S)))
    if frame_group is not None or supporting:
        assert not dynamic, 'frame_group and supporting need the static order'
        order = _frame_order(S, p)
    if supporting:
        assert signvecs == 'all', 'supporting mode: frame rows on all faces, other contacts on facets'
    assert contact_order in ('static', 'filter') and (contact_order == 'static' or (supporting and prefilter)), \
        "contact_order='filter' needs supporting=True and prefilter=True"
    C = [tuple(F(x, p) for x in S[i]) for i in order]
    ctx = {'dynamic': dynamic, 'vecs': FACETVECS if signvecs == 'facet' else SIGNVECS, 'frame_group': frame_group,
           'fcache': {}, 'supporting': supporting, 'extra': ball_rows(C) if ball else None,
           'colfilter': _column_filter(C, stats) if prefilter else None, 'contact_order': contact_order,
           'eutaxy': eutaxy}
    return order, C, ctx

def enumerate_face_data(S, p, max_dim=6, log_every=None, time_limit=None, dynamic=True, signvecs='all',
                        determine=False, resume=None, prefix_depth=3, on_heartbeat=None, prefilter=False,
                        frame_group=None, supporting=False, ball=False, contact_order='static', eutaxy=False):
    """S: list of contacts as numerator tuples over p.  Returns (leaves, stats).
    leaves: dict(signs, x0, basis (affine family), dim, prefix[, rho]).
    signvecs='facet' restricts every contact to a facet interior (a validation pass, NOT complete).
    determine=True applies determine-then-check at equation rank 16.
    prefilter=True applies the certified covector column pre-filter before each LP (stats['prefilter_dead']).
    Checkpointing: subtrees below each depth-`prefix_depth` path are recorded in stats['done_prefixes'] when
    finished; resume={'done_prefixes', 'leaves'} (a previous result) skips them and keeps only their leaves.
    on_heartbeat(leaves, stats) is called every log_every nodes (e.g. to dump a checkpoint).
    frame_group (from frame_stabilizer; static order only): isomorph rejection under the frame symmetries; leaves
    are then complete up to frame_group x AMB, with 'signs' in the caller's contact indices.
    supporting=True (static order, frame first; STATUS §3 Theorem B.1', Lemma B.4): the four frame vectors take
    relative-interior faces as above; every other contact takes a FACET with a half-open cell (zero coordinates
    on the tie side rho, chosen per node; leaves record rho).  Complete: a critical X is a local minimum of |det|
    on the affine space of any supporting facets.  Leaf families are cut out by one equation per non-frame
    contact.  ball=True adds ||cX||_1 <= 2 for every contact to every LP (Lemma B.5).
    contact_order='filter' (supporting + prefilter): branch next on the non-frame contact with the fewest facets
    allowed by the covector filter (LP-free).  Sound, but no reliable gain (DEAD_ENDS B20): use 'static'.
    eutaxy=True: a leaf is recorded only if the l1-eutaxy relaxation (eutaxy_residual, Lemma B.6) is feasible
    (residual <= EUTAXY_TOL); stats['eutaxy_dead'] counts the others."""
    import time
    t0 = time.time()
    done = set(map(tuple, resume['done_prefixes'])) if resume else set()
    leaves = [L for L in resume['leaves'] if L.get('prefix') in done] if resume else []
    stats = {'nodes': 0, 'lp_dead': 0, 'dim_dead': 0, 'leaves': len(leaves), 'rank16': 0, 'rank16_dead': 0,
             'skipped_subtrees': 0, 'per_depth': {}, 'signvecs': signvecs, 'dynamic': dynamic,
             'determine': determine, 'prefix_depth': prefix_depth, 'prefilter': prefilter, 'prefilter_dead': 0,
             'frame_symmetry': frame_group is not None, 'frame_group_order': len(frame_group or []), 'frame_sym_dead': 0,
             'supporting': supporting, 'ball': ball, 'contact_order': contact_order, 'eutaxy': eutaxy, 'eutaxy_dead': 0}
    order, C, ctx = _setup(S, p, dynamic, signvecs, prefilter, frame_group, supporting, ball, stats, contact_order,
                           eutaxy)
    back = lambda sg: {order[i]: v for i, v in sg.items()}
    vecs = ctx['vecs']

    def rec(group, eqs, stricts, nons, signs, rhos, path):
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
        res = _expand(C, group, eqs, stricts, nons, signs, rhos, ctx, stats)
        if res is None:
            sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
            if sol is None:
                return
            x0, N = sol
            if len(N) > max_dim:
                stats['dim_dead'] += 1; return
            if eutaxy and eutaxy_residual(C, signs, eqs, list(stricts) + list(nons)) > EUTAXY_TOL:
                stats['eutaxy_dead'] += 1; return
            stats['leaves'] += 1
            L = {'signs': back(signs), 'x0': x0, 'basis': N, 'dim': len(N), 'prefix': tuple(path[:prefix_depth])}
            if supporting:
                L['rho'] = back(rhos)
            leaves.append(L)
            return
        ci, kids = res
        for img, sub, e1, s1, n1, rho in kids:
            path.append((ci, img))
            signs[ci] = img
            if rho is not None:
                rhos[ci] = rho
            if determine:
                L = _read_off(C, e1, signs, vecs, rhos if supporting else None)
                if L != 'branch':
                    stats['rank16'] += 1
                    if L is not None and eutaxy and \
                            eutaxy_residual(C, L['signs'], e1, list(s1) + list(n1)) > EUTAXY_TOL:
                        stats['eutaxy_dead'] += 1; L = None
                    if L is None:
                        stats['rank16_dead'] += 1
                    else:
                        stats['leaves'] += 1
                        L['signs'] = back(L['signs'])
                        if supporting:
                            L['rho'] = back(rhos)
                        L['prefix'] = tuple(path[:prefix_depth]); leaves.append(L)
                    del signs[ci]; rhos.pop(ci, None); path.pop(); continue
            rec(sub, e1, s1, n1, signs, rhos, path)
            del signs[ci]; rhos.pop(ci, None)
            path.pop()
        if len(path) == prefix_depth:
            done.add(tuple(path))

    try:
        rec(AMB, [], [], [], {}, {}, [])
        stats['complete'] = True
    except TimeoutError:
        stats['complete'] = False
    stats['done_prefixes'] = sorted(done)
    stats['seconds'] = time.time() - t0
    return leaves, stats

def estimate_face_data(S, p, probes, rng, max_dim=6, dynamic=True, signvecs='all', determine=False, max_seconds=None,
                       prefilter=False, frame_group=None, supporting=False, ball=False, contact_order='static',
                       eutaxy=False):
    """Knuth's estimator for the enumerate_face_data tree (same branching code): mean over random
    root-to-leaf probes of the tree size and of the number of leaves."""
    import time
    stats = {'lp_dead': 0, 'prefilter_dead': 0}
    order, C, ctx = _setup(S, p, dynamic, signvecs, prefilter, frame_group, supporting, ball, stats, contact_order,
                           eutaxy)
    vecs = ctx['vecs']
    t0 = time.time()                                   # the certificate is not part of the per-node cost
    nodes, lvs, depth = [], [], []
    for _ in range(probes):
        group, eqs, stricts, nons, signs, rhos = AMB, [], [], [], {}, {}
        w, est, leaf = 1, 1, 0
        while True:
            res = _expand(C, group, eqs, stricts, nons, signs, rhos, ctx, stats)
            if res is None:
                sol = solve_affine([r for r, _ in eqs], [b for _, b in eqs])
                leaf = w if sol and len(sol[1]) <= max_dim else 0
                if leaf and eutaxy and eutaxy_residual(C, signs, eqs, list(stricts) + list(nons)) > EUTAXY_TOL:
                    leaf = 0
                break
            ci, kids = res
            if not kids:
                break
            w *= len(kids); est += w
            img, group, eqs, stricts, nons, rho = rng.choice(kids)
            signs = dict(signs); signs[ci] = img
            if rho is not None:
                rhos = dict(rhos); rhos[ci] = rho
            if determine:
                L = _read_off(C, eqs, signs, vecs, rhos if supporting else None)
                if L != 'branch':
                    leaf = w if L and not (eutaxy and eutaxy_residual(C, L['signs'], eqs, list(stricts) + list(nons))
                                           > EUTAXY_TOL) else 0
                    break
        nodes.append(est); lvs.append(leaf); depth.append(len(signs))
        if max_seconds and time.time() - t0 > max_seconds:
            break
    n = sum(d + 1 for d in depth)
    return {'probes': len(nodes), 'nodes_mean': float(np.mean(nodes)), 'nodes_median': float(np.median(nodes)),
            'nodes_max': float(np.max(nodes)), 'leaves_mean': float(np.mean(lvs)), 'depth_mean': float(np.mean(depth)),
            'seconds': time.time() - t0, 's_per_node': (time.time() - t0) / n,
            'signvecs': signvecs, 'dynamic': dynamic, 'determine': determine, 'prefilter': prefilter,
            'frame_symmetry': frame_group is not None, 'supporting': supporting, 'ball': ball,
            'contact_order': contact_order, 'eutaxy': eutaxy,
            'lp_dead': stats['lp_dead'], 'range_lps': stats.get('range_lps', 0),
            'prefilter_dead': stats['prefilter_dead'], 'frame_sym_dead': stats.get('frame_sym_dead', 0)}

def datum_hits(S, leaves, sig, frame_group=None):
    """Indices of leaves whose face datum equals `sig` ({contact index of S: sign vector}) up to the ambient group
    (columns) and, if given, the frame symmetries (frame_stabilizer).  Used to check that a known lattice (zeta)
    is found by a complete B1 run."""
    acts = [contact_action(S, g) for g in (frame_group or [((0, 1, 2, 3), (1, 1, 1, 1))])]
    out = []
    for k, L in enumerate(leaves):
        for ca in acts:          # frame image: contact i takes eps * (sign vector of contact j)
            t = {i: tuple(e * x for x in L['signs'][j]) for i, (j, e) in enumerate(ca)}
            if any(all(act_sign(q, t[i]) == sig[i] for i in sig) for q in AMB):
                out.append(k); break
    return out
