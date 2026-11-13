"""Stage B2 over the B1 leaves of a stratum (Cor B.2 / Thm B.1' of STATUS §3): window check of critical values.

Usage (from the campaign root):
  python3 scripts/stageB2_stratum.py <sid> [--max-dim 3] [--pkl F ...] [--out F]
Reads the complete sets of results/stageB_<sid>/set_*.pkl (or the given B1 output pickles, e.g. G6S's
results/stageB_faces_zeta_supporting_ball.pkl).  For every leaf (affine family x0 + span(basis), dim d):
  d = 0         the critical value is |det x0| (exact);
  1 <= d <= max_dim   eliminant R(D) of {grad f = 0, D = f}, f = det on the family (cp4.critical, sympy);
  d > max_dim   'pending' (needs msolve; PLAN §3).
Every real root of R(D) and of R(-D) (the window concerns |det X|; det may be negative on the family) and
|det x0| is compared EXACTLY with the window [2p/3, p*D*):
p*D* is the CRootOf of the minimal polynomial 64 D^3 - 1920 D^2 + 13584 D - 26011 (p = 4; rescaled for other p),
compared by rigorous rational approximation (sympy eval_rational); a root equal to p*D* (same minimal polynomial and
root) is reported as 'zeta' - the boundary, outside the half-open window.  Roots inside the window are Stage C
candidates (cell membership, frame condition, exact admissibility) and are listed, never dropped.
Results are cached by the affine family (RREF).  Output: results/stageB2_<sid>.json (or --out)."""
import sys, os, json, glob, pickle, time, argparse
from collections import Counter
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.exact import det as exact_det, rref
from cp4.critical import critical_value_polynomial
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
ap = argparse.ArgumentParser()
ap.add_argument('sid', type=int)
ap.add_argument('--max-dim', type=int, default=3)
ap.add_argument('--pkl', nargs='*')
ap.add_argument('--out')
a = ap.parse_args()
p = json.load(open(os.path.join(root, 'data', 'strata.json')))[a.sid]['p']
D = sp.Symbol('D')
P4 = 64 * D**3 - 1920 * D**2 + 13584 * D - 26011          # minimal polynomial of 4 D* (PLAN §1, gate G5)
Pp = sp.Poly(sp.numer(sp.together(P4.subs(D, 4 * D / p))), D).primitive()[1]
HI = min(sp.real_roots(Pp), key=lambda r: abs(float(r) - p * 0.80822000090353537301))   # p*D*, exact CRootOf
LO = sp.Rational(2 * p, 3)
EPS = sp.Rational(1, 10**40)
HI_Q = HI.eval_rational(dx=EPS)


def classify_interval(a_, b_, same_minpoly):
    """A real algebraic number known to lie in the rational interval [a_, b_] (width <= EPS), compared EXACTLY with
    the window [LO, p*D*): 'below', 'window', 'zeta' (it is p*D* itself) or 'above'.  A number whose interval
    straddles LO is reported as 'window' (conservative)."""
    if b_ < LO:
        return 'below'
    if b_ < HI_Q - EPS:
        return 'window'
    if a_ > HI_Q + EPS:
        return 'above'
    if same_minpoly:
        return 'zeta'             # a root of p*D*'s minimal polynomial within 2e-40 of p*D*: p*D* itself
    raise ValueError('cannot separate a root from p*D* at 1e-40 although it is not a conjugate')


def classify_value(v, minpoly=None):
    """v: a Rational."""
    return classify_interval(v, v, False)


def family_key(L):
    N = [list(b) for b in L['basis']]
    Rm = [tuple(r) for r in (rref(N)[0] if N else []) if any(r)] if N else []
    x = list(L['x0'])
    for r in Rm:
        piv = next(i for i, v in enumerate(r) if v != 0)
        f = x[piv] / r[piv]
        x = [u - f * w for u, w in zip(x, r)]
    return (tuple(Rm), tuple(x))


def b2_leaf(L):
    if L['dim'] == 0:
        x0 = L['x0']
        v = abs(exact_det([[x0[4 * i + j] for j in range(4)] for i in range(4)]))
        vs = sp.Rational(v.numerator, v.denominator)
        c = classify_value(vs)
        return {'status': c, 'values': [str(vs)] if c in ('window', 'zeta') else []}
    if L['dim'] > a.max_dim:
        return {'status': 'pending', 'values': []}
    try:
        R, _, _ = critical_value_polynomial(L['x0'], L['basis'])
    except ValueError as e:
        return {'status': 'error', 'error': str(e), 'values': []}
    hits, st = [], 'clear'
    for fac0, _ in R[1]:
      for sgn in (1, -1):            # the window concerns |det X|: roots of R(D) and of R(-D)
        fac = sp.expand(fac0.subs(D, sgn * D))
        P = sp.Poly(fac, D)
        if P.degree() == 0:
            continue
        same = P.monic() == Pp.monic()
        for (lo_, hi_), _m in P.intervals(eps=EPS):    # exact isolating intervals of the real roots
            if hi_ < 0:
                continue                               # a negative root of R(sgn D): handled by the other sign
            c = classify_interval(sp.Rational(lo_), sp.Rational(hi_), same)
            if c == 'window':
                st = 'window'; hits.append(f'{sp.N((lo_ + hi_) / 2, 20)} ({fac})')
            elif c == 'zeta':
                st = 'zeta' if st == 'clear' else st; hits.append(f'p*D* ({fac})')
    return {'status': st, 'values': hits}


files = a.pkl or sorted(glob.glob(os.path.join(root, 'results', f'stageB_{a.sid}', 'set_*.pkl')))
cache, out, tot = {}, [], Counter()
t0 = time.time()
for f in files:
    R = pickle.load(open(f, 'rb'))
    if not R['stats'].get('complete'):
        out.append({'file': os.path.basename(f), 'skipped': 'B1 incomplete'}); continue
    res = []
    for k, L in enumerate(R['leaves']):
        key = family_key(L)
        if key not in cache:
            cache[key] = b2_leaf(L)
        r = dict(cache[key], leaf=k, dim=L['dim'])
        res.append(r); tot[r['status']] += 1
    out.append({'file': os.path.basename(f), 'set': R.get('set'), 'leaves': len(res),
                'status': dict(Counter(r['status'] for r in res)),
                'flagged': [r for r in res if r['status'] in ('window', 'zeta', 'error')]})
    print(f"{os.path.basename(f)}: {len(res)} leaves {dict(Counter(r['status'] for r in res))} "
          f"[{time.time()-t0:.0f}s, {len(cache)} families]", flush=True)
summary = {'stratum': a.sid, 'p': p, 'window': f'[{LO}, p*D*)', 'max_dim': a.max_dim, 'totals': dict(tot),
           'families': len(cache), 'seconds': round(time.time() - t0), 'sets': out}
dst = a.out or os.path.join(root, 'results', f'stageB2_{a.sid}.json')
json.dump(summary, open(dst, 'w'), indent=1)
print('B2 totals:', dict(tot), '->', dst, flush=True)
