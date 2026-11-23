"""Stage B2 over the B1 leaves of a stratum (Cor B.2 / Thm B.1' of STATUS §3): window check of critical values.

Usage (from the campaign root):
  python3 scripts/stageB2_stratum.py <sid> [--max-dim 3] [--task i --ntasks N] [--seconds T]   # job-array task
  python3 scripts/stageB2_stratum.py <sid> --merge                  # combine -> results/stageB2_<sid>.json
  python3 scripts/stageB2_stratum.py <sid> --pkl F [F ...] [--out F] # single B1 outputs (gates)
A task processes the B1 outputs results/stageB_<sid>/set_*.pkl with index = i (mod N) and writes
results/stageB2_<sid>/set_*.json (rewritten every 25 leaves; resumed; complete sets skipped).
For every leaf (affine family x0 + span(basis), dim d):
  d = 0               the critical value is |det x0| (exact);
  1 <= d <= max_dim   eliminant R(D) of the ideal (grad f, D - f), f = det on the family (cp4.critical, sympy);
                      a positive-dimensional critical ideal is handled by a lex elimination (f is constant on
                      each component of the critical locus, so R is still nonzero; STATUS §3);
  d > max_dim         'pending' (needs msolve; PLAN §3).
Every real root of R(D) and of R(-D) (the window concerns |det X|; det may be negative on the family) and
|det x0| is compared EXACTLY with the window [2p/3, p*D*), by exact isolating intervals against the CRootOf p*D*
(64 D^3 - 1920 D^2 + 13584 D - 26011 for p = 4, rescaled for other p).  A root equal to p*D* is reported as
'zeta' (the boundary, outside the half-open window) together with the number of complex critical points over its
factor ('zeta_fibre'; = factor degree means one critical point over p*D*: equality classification).  Roots inside
the window are Stage C candidates and are listed, never dropped.  Results are cached by the affine family (RREF)."""
import sys, os, json, glob, pickle, time, argparse
from collections import Counter
from fractions import Fraction as F
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
from cp4.exact import det as exact_det, rref
from cp4.critical import critical_value_polynomial_any, fibre_point_count
root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
ap = argparse.ArgumentParser()
ap.add_argument('sid', type=int)
ap.add_argument('--max-dim', type=int, default=3)
ap.add_argument('--pkl', nargs='*')
ap.add_argument('--out')
ap.add_argument('--task', type=int, default=0)
ap.add_argument('--ntasks', type=int, default=1)
ap.add_argument('--seconds', type=float)
ap.add_argument('--merge', action='store_true')
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
        R, _, _, zerodim = critical_value_polynomial_any(L['x0'], L['basis'])
    except ValueError as e:
        return {'status': 'error', 'error': str(e), 'values': []}
    hits, st = [], 'clear'
    fibres = []
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
                # equality classification: number of complex critical points over the roots of this factor
                # (= degree <=> exactly one critical point over p*D* in this family)
                n = fibre_point_count(L['x0'], L['basis'], fac0) if zerodim else None   # fac0: in det, not |det|
                fibres.append({'factor_degree': P.degree(), 'points': n})
    out = {'status': st, 'values': hits, 'zero_dimensional': zerodim}
    if fibres:
        out['zeta_fibre'] = fibres
    return out


def run_file(f, cache, prev=None, save=None, deadline=None):
    """B2 on one B1 output file; `prev` = partial result to resume, `save(res)` called every 25 leaves."""
    R = pickle.load(open(f, 'rb'))
    if not R['stats'].get('complete'):
        return {'file': os.path.basename(f), 'skipped': 'B1 incomplete'}
    res = list(prev['results']) if prev else []
    for k in range(len(res), len(R['leaves'])):
        if deadline and time.time() > deadline:
            break
        L = R['leaves'][k]
        key = family_key(L)
        if key not in cache:
            cache[key] = b2_leaf(L)
        res.append(dict(cache[key], leaf=k, dim=L['dim']))
        if save and len(res) % 25 == 0:
            save({'file': os.path.basename(f), 'set': R.get('set'), 'leaves': len(R['leaves']), 'results': res,
                  'complete': False})
    return {'file': os.path.basename(f), 'set': R.get('set'), 'leaves': len(R['leaves']), 'results': res,
            'complete': len(res) == len(R['leaves'])}


def digest(r):
    rs = r.get('results', [])
    return {'file': r['file'], 'set': r.get('set'), 'leaves': r.get('leaves'), 'complete': r.get('complete'),
            'skipped': r.get('skipped'), 'status': dict(Counter(x['status'] for x in rs)),
            'flagged': [x for x in rs if x['status'] in ('window', 'zeta', 'error')]}


t0 = time.time()
outdir = os.path.join(root, 'results', f'stageB2_{a.sid}')
if a.pkl:                                    # single files (gates): one summary
    cache, out = {}, []
    for f in a.pkl:
        out.append(digest(run_file(f, cache)))
    tot = Counter()
    for o in out:
        tot.update(o['status'])
elif a.merge:                                # combine the per-set results of a job array
    out, tot = [], Counter()
    for f in sorted(glob.glob(os.path.join(outdir, 'set_*.json'))):
        o = digest(json.load(open(f)))
        out.append(o); tot.update(o['status'])
    n_b1 = len(glob.glob(os.path.join(root, 'results', f'stageB_{a.sid}', 'set_*.pkl')))
    tot['sets_done'] = sum(1 for o in out if o['complete']); tot['sets_total'] = n_b1
else:                                        # job-array task: sets k = task (mod ntasks), resumable per set
    os.makedirs(outdir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(root, 'results', f'stageB_{a.sid}', 'set_*.pkl')))
    deadline = t0 + a.seconds if a.seconds else None
    cache, out, tot = {}, [], Counter()
    for i, f in enumerate(files):
        if i % a.ntasks != a.task:
            continue
        dst1 = os.path.join(outdir, os.path.basename(f).replace('.pkl', '.json'))
        prev = json.load(open(dst1)) if os.path.exists(dst1) else None
        if prev and prev.get('complete'):
            continue
        if prev and prev.get('max_dim', a.max_dim) != a.max_dim:
            sys.exit(f'{dst1}: written with --max-dim {prev["max_dim"]}, not {a.max_dim}')
        if deadline and time.time() > deadline:
            print('time budget used up; resubmit to continue', flush=True); break

        def save(r, dst1=dst1):
            r = dict(r, max_dim=a.max_dim)
            json.dump(r, open(dst1 + '.tmp', 'w')); os.replace(dst1 + '.tmp', dst1)
        r = run_file(f, cache, prev, save, deadline)
        if 'skipped' not in r:
            save(r)
        o = digest(r); out.append(o); tot.update(o['status'])
        print(f"{o['file']}: {o['leaves']} leaves {o['status']} complete={o['complete']} "
              f"[{time.time()-t0:.0f}s, {len(cache)} families]", flush=True)
summary = {'stratum': a.sid, 'p': p, 'window': f'[{LO}, p*D*)', 'max_dim': a.max_dim, 'totals': dict(tot),
           'seconds': round(time.time() - t0), 'sets': out}
dst = a.out or (os.path.join(root, 'results', f'stageB2_{a.sid}.json') if (a.pkl or a.merge)
                else os.path.join(outdir, f'task_{a.task:04d}.json'))
json.dump(summary, open(dst, 'w'), indent=1)
print('B2 totals:', dict(tot), '->', dst, flush=True)
