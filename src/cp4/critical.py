"""Stage B part 2: critical values of det on an affine family (Theorem B.1 / Corollary B.2).

Given a face-data leaf (x0, basis N_1..N_d) in X-space (x[4i+j] = X[i][j]), put
    f(t) = det( X(t) ),   X(t) = x0 + sum_k t_k N_k     (exact rational coefficients).
Every critical lattice with that face datum is a critical point of f (Thm B.1).
critical_value_polynomial() eliminates t from  {df/dt_k = 0, D - f = 0}  and returns a
univariate polynomial R(D) whose roots contain every critical value (Groebner, lex, D last).
Certificate for the leaf (Cor B.2): R has no real root in [2p/3, p*D*), checked by exact
Sturm/real-root isolation (sympy real_roots on rational polynomials is exact).

Caveat: the elimination polynomial may be identically zero if f has a positive-dimensional
critical set with non-constant f -- impossible -- or the critical set is positive-dimensional
(f constant on a component): then R(D) is still nonzero (D - const).  If the Groebner basis
contains no univariate polynomial in D the routine raises; such leaves need manual treatment.
"""
import sympy as sp
from fractions import Fraction as F


def det_poly(x0, basis):
    d = len(basis)
    t = sp.symbols(f't0:{d}')
    X = sp.Matrix(4, 4, lambda i, j: sp.Rational(x0[4 * i + j].numerator, x0[4 * i + j].denominator)
                  + sum(sp.Rational(b[4 * i + j].numerator, b[4 * i + j].denominator) * t[k]
                        for k, b in enumerate(basis)))
    return sp.expand(X.det(method='berkowitz')), t


def critical_value_polynomial(x0, basis):
    f, t = det_poly(x0, basis)
    D = sp.Symbol('D')
    eqs = [sp.diff(f, v) for v in t] + [D - f]
    G = sp.groebner(eqs, *t, D, order='grevlex')
    if not G.is_zero_dimensional:
        raise ValueError('critical ideal not zero-dimensional')
    G = G.fglm('lex')
    uni = [g for g in G.exprs if g.free_symbols <= {D}]
    if not uni:
        raise ValueError('no univariate eliminant (positive-dimensional fibre?)')
    return sp.factor_list(uni[0], D), f, t


def roots_in_window(R, lo, hi):
    """Exact real roots of the rational polynomial R(D) in the half-open window [lo, hi)."""
    D = sp.Symbol('D')
    out = []
    for fac, _ in R[1]:
        P = sp.Poly(fac, D)
        if P.degree() == 0:
            continue
        for r in sp.real_roots(P):
            if sp.Rational(lo) <= r and r < hi:
                out.append(r)
    return out


def critical_value_polynomial_any(x0, basis):
    """As critical_value_polynomial, but also for a POSITIVE-DIMENSIONAL critical ideal.

    Justification (STATUS §3, remark after Lemma B.6): the complex critical locus of f is a variety with finitely
    many irreducible components, and f is constant on each (algebraic Sard; f restricted to a component has zero
    differential).  So the set of critical values is finite, the elimination ideal (grad f, D - f) ∩ Q[D] is
    nonzero, and its generator R(D) vanishes at every critical value, real or complex.  It is computed with a lex
    Groebner basis (D last) - slower than grevlex + FGLM, so only used when the ideal is not zero-dimensional.
    Returns (factor list of R, f, t, zero_dimensional flag)."""
    try:
        R, f, t = critical_value_polynomial(x0, basis)
        return R, f, t, True
    except ValueError:
        pass
    f, t = det_poly(x0, basis)
    D = sp.Symbol('D')
    G = sp.groebner([sp.diff(f, v) for v in t] + [D - f], *t, D, order='lex')
    uni = [g for g in G.exprs if g.free_symbols <= {D}]
    if not uni:
        raise ValueError('no univariate eliminant (impossible by algebraic Sard: check the input)')
    return sp.factor_list(uni[0], D), f, t, False


def fibre_point_count(x0, basis, factor):
    """Number of complex critical points (with multiplicity of the quotient) of f = det on the family lying over
    the roots of the irreducible polynomial `factor` in D, or None if that fibre is not zero-dimensional.  For an
    irreducible factor of degree k, a count of k means exactly one critical point over each root - used to show
    that the equality value p*D* is attained at a single point of a leaf family (equality classification)."""
    f, t = det_poly(x0, basis)
    D = sp.Symbol('D')
    G = sp.groebner([sp.diff(f, v) for v in t] + [D - f, factor], *t, D, order='grevlex')
    if not G.is_zero_dimensional:
        return None
    # vector-space dimension of the quotient = number of standard monomials
    from sympy.polys.orderings import grevlex
    gens = list(t) + [D]
    lead = [sp.Poly(g, *gens).monoms(order='grevlex')[0] for g in G.exprs]
    n = len(gens); count = 0; frontier = [(0,) * n]; seen = set(frontier)
    while frontier:
        m = frontier.pop()
        if any(all(m[i] >= l[i] for i in range(n)) for l in lead):
            continue
        count += 1
        for i in range(n):
            m2 = tuple(m[j] + (j == i) for j in range(n))
            if m2 not in seen:
                seen.add(m2); frontier.append(m2)
    return count
