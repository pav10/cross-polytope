"""Merged search (PLAN §4): Stage A (contact sets) and Stage B1 (face data) as one tree.

Unknown: the max-det frame X (rows = frame minimal vectors, ambient coordinates), x[4i+j] = X[i][j].
Every stratum candidate c (frame coordinates, docs/STATUS §3) gets a label, in a fixed order:
    in(s): c is a minimal vector and cX lies in the relative interior of the face s of the l1 sphere
           of radius 2 (faces._rows_for: sum_j s_j (cX)_j = 2, s_j (cX)_j > 0, (cX)_j = 0 where s_j = 0);
    out:   c is not minimal, i.e. F(c) = ||cX||_1 > 2 (lazy: no face is chosen for it).
The frame rows e_1..e_4 are 'in'; their sign vectors form the frame pattern.

Prunes (necessary conditions on the true labels of a critical lattice; docs/STATUS §2-§3):
  A2 quadruple dets <= 1; |S| <= 20, and |S| >= 10 still reachable (A5); hull closure (Lemma A.1: no
  candidate in int conv(±S); boundary candidates must be 'in'); strict LP feasibility of the face rows
  (faces.lp_feasible); frame-sign filter (exact: sign((cX)_j) is forced when all nonzero c_i X_ij agree);
  at full leaves dim A <= 6 (Thm 5.1) and, for every 'out' candidate, max F(c) over the cell closure
  is not below 2 (else F(c) > 2 is impossible).
Determine-then-check (PLAN §3): when the face equations have rank 16, X is unique and rational; the
  remaining labels are read off exactly (F < 2 => dead; F = 2 => in; F > 2 => out); dim A = 0.
Symmetry: the root frame pattern is reduced modulo (stratum stabiliser on rows) x (ambient group on
  columns); below the root, sign vectors are reduced modulo the ambient subgroup fixing every chosen
  sign vector (stabiliser chain, as in faces.py).  Complete; not fully non-redundant.

Floats: LP prunes use the strict-slack convention of faces.py and a conservative exclusion margin;
leaves are solved exactly (solve_affine) and dim-0 leaves are checked exactly.  Proof-grade use needs
exact certificates for the LP prunes (PLAN §4).
"""
import time
import numpy as np
from itertools import product, permutations, combinations, combinations_with_replacement
from fractions import Fraction as F
from scipy.optimize import linprog
from .exact import solve_affine, rref, det as exact_det
from .faces import AMB, act_sign, orbit_reps, _rows_for, lp_feasible
from .contacts import classify_points, _quad_dets_ok, canon_set
from .lattice import is_admissible_exact, DSTAR

FACETS = list(product((1, -1), repeat=4))
EXCL_MARGIN = 1e-9


def _colkey(rows):
    """Canonical form of a 4x4 sign pattern under the ambient group (signed column permutations)."""
    cols = []
    for j in range(4):
        col = tuple(r[j] for r in rows)
        if next((v for v in col if v), 0) < 0:
            col = tuple(-v for v in col)
        cols.append(col)
    return tuple(sorted(cols))


def _rows_of(key):
    return tuple(tuple(col[i] for col in key) for i in range(4))


def _row_act(g, rows):
    """Stratum symmetry g on frame coordinates (c'_k = sg_k c_{pm_k}) moves frame rows the same way."""
    pm, sg = g
    return tuple(tuple(sg[k] * x for x in rows[pm[k]]) for k in range(4))


def frame_patterns(stab):
    """Frame sign patterns modulo stab x ambient group.  Patterns forcing det X = 0 (a zero row or
    column, or no nonzero transversal) are dropped (Cor 1.8: det X != 0)."""
    cols = sorted({tuple(-x for x in v) if next(x for x in v if x) < 0 else v
                   for v in product((-1, 0, 1), repeat=4) if any(v)})
    seen, reps = set(), []
    for key in combinations_with_replacement(cols, 4):
        if key in seen:
            continue
        rows = _rows_of(key)
        orbit = {_colkey(_row_act(g, rows)) for g in stab}
        seen |= orbit
        if all(any(r) for r in rows) and any(all(rows[i][pi[i]] for i in range(4))
                                             for pi in permutations(range(4))):
            reps.append(rows)
    return reps


def frame_filter(c, tau):
    """Sign vectors s for contact c compatible with frame pattern tau (exact, no LP)."""
    per = []
    for j in range(4):
        T = {(1 if c[i] > 0 else -1) * tau[i][j] for i in range(4) if c[i] and tau[i][j]}
        per.append((0,) if not T else (1,) if T == {1} else (-1,) if T == {-1} else (-1, 0, 1))
    return [s for s in product(*per) if any(s)]


def gamma_basis(residues, p):
    """Integer basis H (rows) of the Z-span of p*I and the residues; Gamma = Z^4 H / p."""
    A = [[p * int(i == j) for j in range(4)] for i in range(4)] + [list(r) for r in residues]
    basis = []
    for j in range(4):
        while sum(1 for r in A if r[j]) > 1:
            piv = min((r for r in A if r[j]), key=lambda r: abs(r[j]))
            for r in A:
                if r is not piv and r[j]:
                    q = r[j] // piv[j]
                    for t in range(4):
                        r[t] -= q * piv[t]
            A = [r for r in A if any(r)]
        piv = next(r for r in A if r[j])
        basis.append(piv)
        A = [r for r in A if r is not piv]
    return basis


class Node:
    __slots__ = ('k', 'S', 'signs', 'lab', 'eqs', 'st', 'group', 'must', 'x', 'rank', 'root')

    def child(self, **kw):
        n = Node()
        for f in Node.__slots__:
            setattr(n, f, kw.get(f, getattr(self, f)))
        return n


class Merged:
    def __init__(self, stratum, stab, adm_prune=False):
        self.p = p = stratum['p']
        self.residues = stratum['residues']
        self.cand = [tuple(c) for c in stratum['candidates']]
        self.C = [tuple(F(x, p) for x in c) for c in self.cand]
        self.frame = [self.cand.index(tuple(p * int(i == j) for j in range(4))) for i in range(4)]
        self.forced = [self.cand.index(tuple(c)) for c in stratum['forced']]
        self.stab, self.adm_prune = stab, adm_prune
        self.taus = frame_patterns(stab)
        self.leaves, self.leaves0 = [], []
        self.stats = {'nodes': 0, 'roots': len(self.taus), 'roots_done': 0, 'lp': 0, 'lp_s': 0.0,
                      'dead_lp': 0, 'dead_det': 0, 'dead_size': 0, 'dead_interior': 0,
                      'dead_forced_out': 0, 'dead_dim': 0, 'dead_excl': 0, 'dead_rank16': 0,
                      'rank16': 0, 'leaves': 0, 'leaves0': 0, 'band0': 0, 'per_depth': {}}

    # ---------------------------------------------------------------- roots
    def root(self, t):
        return self.root_pattern(self.taus[t])

    def root_pattern(self, tau):
        """Node with the frame rows labelled by pattern tau (no LP needed: rows independent)."""
        allowed = {i: frame_filter(self.C[i], tau) for i in range(len(self.cand))}
        rest = sorted((i for i in range(len(self.cand)) if i not in self.frame),
                      key=lambda i: (len(allowed[i]), sum(map(abs, self.cand[i])), i))
        n = Node()
        n.k, n.S, n.signs, n.lab, n.eqs, n.st = 4, list(self.frame), {}, {}, [], []
        n.group, n.must, n.x, n.root = AMB, frozenset(self.forced) - set(self.frame), None, (tau, rest, allowed)
        for i, ci in enumerate(self.frame):
            e, s = _rows_for(self.C[ci], tau[i])
            n.eqs, n.st = n.eqs + e, n.st + s
            n.signs[ci], n.lab[ci] = tau[i], 'in'
            n.group = [g for g in n.group if act_sign(g, tau[i]) == tau[i]]
        n.rank = int(np.linalg.matrix_rank(np.array([[float(v) for v in r] for r, _ in n.eqs])))
        return n

    def order(self, n):
        return n.root[1]

    # ---------------------------------------------------------------- one label
    def _lp(self, eqs, st):
        t = time.time()
        ok, x = lp_feasible(eqs, st)
        self.stats['lp'] += 1; self.stats['lp_s'] += time.time() - t
        return ok, x

    def pre_include(self, n, ci):
        """Sign-independent checks for adding ci to S: (S1, must) or None.  A2, |S| <= 20, hull closure."""
        p, S = self.p, n.S
        if len(S) >= 20:
            self.stats['dead_size'] += 1; return None
        if not _quad_dets_ok([self.cand[i] for i in S], self.cand[ci], p):
            self.stats['dead_det'] += 1; return None
        S1 = S + [ci]
        interior, boundary, _ = classify_points([self.cand[i] for i in S1], self.cand, p)
        if len(interior):
            self.stats['dead_interior'] += 1; return None
        must = set(n.must)
        for b in boundary:
            b = int(b)
            if b in S1:
                continue
            if n.lab.get(b) == 'out':
                self.stats['dead_forced_out'] += 1; return None
            must.add(b)
        return S1, frozenset(must)

    def include(self, n, ci, s, pre=None):
        pre = pre or self.pre_include(n, ci)
        if pre is None:
            return None
        S1, must = pre
        e, st = _rows_for(self.C[ci], s)
        eqs, sts = n.eqs + e, n.st + st
        ok, x = self._lp(eqs, sts)
        if not ok:
            self.stats['dead_lp'] += 1; return None
        signs = dict(n.signs); signs[ci] = s
        lab = dict(n.lab); lab[ci] = 'in'
        rank = int(np.linalg.matrix_rank(np.array([[float(v) for v in r] for r, _ in eqs])))
        group = [g for g in n.group if act_sign(g, s) == s]
        return n.child(k=n.k + 1, S=S1, signs=signs, lab=lab, eqs=eqs, st=sts, group=group,
                       must=must, x=x, rank=rank)

    def exclude(self, n, ci):
        if ci in n.must:
            self.stats['dead_forced_out'] += 1; return None
        remaining = len(self.order(n)) - (n.k - 4) - 1
        if len(n.S) + remaining < 10:
            self.stats['dead_size'] += 1; return None
        lab = dict(n.lab); lab[ci] = 'out'
        m = n.child(k=n.k + 1, lab=lab)
        if self.adm_prune and n.x is not None and self.cannot_exceed_2(ci, m):
            self.stats['dead_excl'] += 1; return None
        return m

    # ---------------------------------------------------------------- tests on a cell
    def cannot_exceed_2(self, ci, n):
        """True iff max F(c) over the closure of n's cell is < 2 - margin (so c cannot be 'out')."""
        c = np.array([float(v) for v in self.C[ci]])
        if n.x is not None and np.abs(c @ np.array(n.x).reshape(4, 4)).sum() >= 2 - EXCL_MARGIN:
            return False
        A_eq = np.array([[float(v) for v in r] for r, _ in n.eqs])
        b_eq = np.array([float(b) for _, b in n.eqs])
        A_ub = np.array([[-float(v) for v in r] for r in n.st]) if n.st else None
        b_ub = np.zeros(len(n.st)) if n.st else None
        for v in FACETS:
            obj = -np.outer(c, v).ravel()
            t = time.time()
            res = linprog(obj, A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(-2, 2)] * 16, method='highs')
            self.stats['lp'] += 1; self.stats['lp_s'] += time.time() - t
            if res.status == 0 and -res.fun >= 2 - EXCL_MARGIN:
                return False
        return True

    # ---------------------------------------------------------------- terminals
    def full_leaf(self, n, record=True):
        if len(n.S) < 10:
            self.stats['dead_size'] += 1; return 0
        sol = solve_affine([r for r, _ in n.eqs], [b for _, b in n.eqs])
        if sol is None:
            self.stats['dead_lp'] += 1; return 0
        x0, N = sol
        if len(N) > 6:
            self.stats['dead_dim'] += 1; return 0
        for ci, l in n.lab.items():
            if l == 'out' and self.cannot_exceed_2(ci, n):
                self.stats['dead_excl'] += 1; return 0
        self.stats['leaves'] += 1
        if record:
            R, _ = rref([list(r) + [b] for r, b in n.eqs])
            self.leaves.append({'S': [self.cand[i] for i in n.S], 'signs': [n.signs[i] for i in n.S],
                                'dim': len(N), 'family': hash(str(R)), 'frame_pattern': n.root[0]})
        return 1

    def rank16(self, n, record=True):
        """X is unique: read every label off exactly.  Returns 1 for a dim-0 leaf, 0 if dead."""
        self.stats['rank16'] += 1
        sol = solve_affine([r for r, _ in n.eqs], [b for _, b in n.eqs])
        if sol is None or sol[1]:
            return None                      # float rank was wrong: keep branching
        x0 = sol[0]
        X = [x0[4 * i:4 * i + 4] for i in range(4)]
        S = []
        for ci, c in enumerate(self.C):
            v = [sum(c[i] * X[i][j] for i in range(4)) for j in range(4)]
            Fc = sum(abs(t) for t in v)
            sg = tuple((t > 0) - (t < 0) for t in v)
            lab = n.lab.get(ci)
            if Fc < 2 or (lab == 'in' and sg != n.signs[ci]) or (lab == 'out' and Fc == 2) \
                    or (ci in n.must and Fc != 2):
                self.stats['dead_rank16'] += 1; return 0
            if Fc == 2:
                S.append(ci)
        if not 10 <= len(S) <= 20 or any(abs(exact_det([list(self.cand[i]) for i in q])) > self.p ** 4
                                          for q in combinations(S, 4)):
            self.stats['dead_rank16'] += 1; return 0
        d = exact_det(X)
        band = F(2 * self.p, 3) <= abs(d) and abs(d) < self.p * DSTAR
        leaf = {'S': [self.cand[i] for i in S], 'detX': str(d), 'band': band, 'frame_pattern': n.root[0]}
        if band:
            self.stats['band0'] += 1
            H = gamma_basis(self.residues, self.p)
            Y = [[sum(F(H[r][i], self.p) * X[i][j] for i in range(4)) for j in range(4)] for r in range(4)]
            try:
                leaf['admissible'] = is_admissible_exact(Y)
            except ValueError as e:
                leaf['admissible'] = f'undecided: {e}'
        self.stats['leaves0'] += 1
        if record:
            self.leaves0.append(leaf)
        return 1

    # ---------------------------------------------------------------- expansion
    def expand(self, n, record=True):
        """Children of n (labels for the next candidate) and the number of leaves found on the way."""
        order = self.order(n)
        if n.k - 4 == len(order):
            return [], self.full_leaf(n, record)
        ci = order[n.k - 4]
        tau, _, allowed = n.root
        kids, found = [], 0
        pre = self.pre_include(n, ci)
        for s in (orbit_reps(n.group, allowed[ci]) if pre else []):
            m = self.include(n, ci, s, pre)
            if m is None:
                continue
            if m.rank == 16:
                r = self.rank16(m, record)
                if r is not None:
                    found += r
                    continue
            kids.append(m)
        m = self.exclude(n, ci)
        if m is not None:
            kids.append(m)
        return kids, found

    def run(self, time_limit=None, node_limit=None, log_every=None, on_heartbeat=None):
        t0 = time.time()
        st = self.stats
        try:
            for t in range(len(self.taus)):
                stack = [self.root(t)]
                if stack[0].rank == 16:
                    r = self.rank16(stack[0])
                    if r is not None:
                        stack = []
                while stack:
                    n = stack.pop()
                    st['nodes'] += 1
                    d = len(n.lab)
                    st['per_depth'][d] = st['per_depth'].get(d, 0) + 1
                    if log_every and st['nodes'] % log_every == 0 and on_heartbeat:
                        on_heartbeat(self, time.time() - t0)
                    if (time_limit and time.time() - t0 > time_limit) or (node_limit and st['nodes'] > node_limit):
                        raise TimeoutError
                    kids, _ = self.expand(n)
                    stack.extend(reversed(kids))         # 'in' children first, 'out' last
                st['roots_done'] = t + 1
            st['complete'] = True
        except TimeoutError:
            st['complete'] = False
        st['seconds'] = time.time() - t0
        return st

    def estimate(self, probes, rng, max_seconds=None):
        """Knuth's estimator: mean over random root-to-leaf probes of the tree size (nodes) and leaves."""
        t0 = time.time()
        nodes, leaves, done = [], [], 0
        for _ in range(probes):
            w = len(self.taus)
            n = self.root(rng.randrange(w))
            est_n, est_l = 1 + w, 0
            r = self.rank16(n, record=False) if n.rank == 16 else None
            while r is None:
                kids, found = self.expand(n, record=False)
                est_l += w * found
                if not kids:
                    break
                w *= len(kids)
                est_n += w
                n = rng.choice(kids)
            if r:
                est_l += w * r
            nodes.append(est_n); leaves.append(est_l); done += 1
            if max_seconds and time.time() - t0 > max_seconds:
                break
        return {'probes': done, 'nodes_mean': float(np.mean(nodes)), 'nodes_median': float(np.median(nodes)),
                'nodes_max': float(np.max(nodes)), 'leaves_mean': float(np.mean(leaves)),
                'seconds': time.time() - t0}

    def canonical_sets(self):
        return {canon_set(L['S'], self.stab, self.p) for L in self.leaves} | \
               {canon_set(L['S'], self.stab, self.p) for L in self.leaves0}
