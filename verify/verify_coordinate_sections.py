"""Exact checks for docs/COORDINATE_SECTIONS.md: strata 6, 7 and 25 excluded.

Run from any directory: python3 verify/verify_coordinate_sections.py   (about 30 s; prints COORDINATE SECTIONS PASS)
  1. strata 6, 7 and 25: the stated residue groups, unit-shell vectors and pairwise support-disjointness of rows
     2, 3, 4 (from the shell, as in THEORY_REVIEW §2); the pair groups Gamma_ij used in the lemma;
  2. stratum 7, residual pattern: the offsets 2(c2, c3) mod 2Z^2 of the four classes c4 = k/4, and the identity
     sum_k dist(offset_k + w, 2Z^2) = 4 for all w (checked exactly on a rational grid and proved in the note);
     the same for stratum 25 on its coset line c1 = 1/2;
  3. stratum 7, the one-dimensional lemma: g(s) = max(|s|/2, |D + e s/2|).  For every D in [0, 81/100],
     e in [0, 1] and delta in [0, 1), some window {delta + j + i : i = 0..3}, j in -4..1, has sum g < 4.
     Certified on the rational grid of step 1/40 in (D, e, delta): the sum is at most 4 - 8/40 at the grid point,
     and each window sum is 8-Lipschitz in the sup-distance to the nearest grid point (|dg/dD| <= 1,
     |dg/de| <= |s|/2 <= 5/2, |dg/ddelta| <= 1/2; half a step in each of D, e, delta; four terms).
     Known answer (non-vacuity): at D = 1, e = 0 some delta satisfies every window (the axis-row limit).
Exact rational arithmetic only (fractions)."""
import json
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
strata = {s['id']: s for s in json.loads((ROOT / 'data/strata.json').read_text())}


def residues(sid):
    return {tuple(r) for r in strata[sid]['residues']}


def disjoint_pairs(sid):
    sh = [tuple(f) for f in strata[sid]['forced']]
    out = set()
    for i, j in combinations(range(4), 2):
        sg = {(n[i] * n[j] > 0) - (n[i] * n[j] < 0) for n in sh}
        if {-1, 1} <= sg:
            out.add((i, j))
    return out


def sub_group(sid, rows):
    """Residues of Gamma supported on the given rows (0-based), as fractions of p."""
    p = strata[sid]['p']
    return {tuple(F(r[i], p) for i in rows) for r in residues(sid) if all(r[k] == 0 for k in range(4) if k not in rows)}


def check_structure():
    # Stratum 6: Gamma = Z^4 + <(0,0,1/2,1/2), (0,1/2,0,1/2)>, rows 2,3,4 (0-based 1,2,3) pairwise disjoint.
    assert strata[6]['p'] == 4 and residues(6) == {(0, 0, 0, 0), (0, 0, 2, 2), (0, 2, 0, 2), (0, 2, 2, 0)}
    assert {(1, 2), (1, 3), (2, 3)} <= disjoint_pairs(6)
    for pair in [(1, 2), (1, 3), (2, 3)]:
        assert sub_group(6, pair) == {(0, 0), (F(1, 2), F(1, 2))}, pair      # Z^2 + Z(1/2,1/2)
    # Stratum 7: Gamma = Z^4 + Z(0,1/2,1/4,1/4).
    assert strata[7]['p'] == 4 and residues(7) == {tuple((k * x) % 4 for x in (0, 2, 1, 1)) for k in range(4)}
    assert {(1, 2), (1, 3), (2, 3)} <= disjoint_pairs(7)
    assert sub_group(7, (2, 3)) == {(0, 0), (F(1, 2), F(1, 2))}               # rows 3,4: Z^2 + Z(1/2,1/2)
    assert sub_group(7, (1, 2)) == {(0, 0)} and sub_group(7, (1, 3)) == {(0, 0)}   # rows 2,3 and 2,4: Z^2
    assert all(r[0] == 0 for r in residues(7))                                  # Gamma = Z (+) Gamma_234
    # Rows 3 <-> 4 symmetry of stratum 7 (used to assume X3 is the axis row in the residual pattern).
    assert {(r[0], r[1], r[3], r[2]) for r in residues(7)} == residues(7)
    # Stratum 25 (p = 8): same unit shell pattern as stratum 7, Gamma_25 contains Gamma_7 with index 2.
    assert strata[25]['p'] == 8 and {(1, 2), (1, 3), (2, 3)} <= disjoint_pairs(25)
    assert sub_group(25, (2, 3)) == {(0, 0), (F(1, 2), F(1, 2))}
    assert {tuple(2 * x % 8 for x in r) for r in residues(7)} <= residues(25)     # (0,1/2,1/4,1/4) in Gamma_25
    assert {(r[0], r[1], r[3], r[2]) for r in residues(25)} == residues(25)
    print("strata 6, 7, 25: residue groups, disjoint rows 2,3,4, pair groups: checked")


def dist2(u):
    """Distance from the rational u to 2Z."""
    r = u % 2
    return min(r, 2 - r)


def check_offsets():
    # Classes of Gamma_234 (coefficients of rows 2,3,4) by c4 = k/4: offsets 2(c2, c3) mod 2Z^2.
    cls = {}
    for r in residues(7):
        c2, c3, c4 = F(r[1], 4), F(r[2], 4), F(r[3], 4)
        cls[int(c4 * 4) % 4] = ((2 * c2) % 2, (2 * c3) % 2)
    assert cls == {0: (0, 0), 1: (1, F(1, 2)), 2: (0, 1), 3: (1, F(3, 2))}, cls
    for a, b in product([F(i, 24) for i in range(-48, 49)], repeat=2):    # sum over the 4 classes is 4
        assert sum(dist2(a + o[0]) + dist2(b + o[1]) for o in cls.values()) == 4
    # Stratum 25, the coset line c1 = 1/2: classes c4 = 1/8 + k/4, offsets 2(c2, c3) mod 2Z^2, and the identity
    # (with the frame-row-1 contribution (alpha, beta)/2 absorbed into the free shift w).
    cls25 = {}
    for r in residues(25):
        if r[0] == 4:
            c2, c3, c4 = F(r[1], 8), F(r[2], 8), F(r[3], 8)
            cls25[c4] = ((2 * c2) % 2, (2 * c3) % 2)
    assert set(cls25) == {F(1, 8), F(3, 8), F(5, 8), F(7, 8)}, cls25
    assert all(r[0] in (0, 4) for r in residues(25))
    assert sorted(r for r in residues(25) if r[0] == 0) == sorted(tuple(2 * x for x in r) for r in residues(7))  # c1 = 0 slice = Gamma_7
    for a, b in product([F(i, 24) for i in range(-48, 49)], repeat=2):
        assert sum(dist2(a + o[0]) + dist2(b + o[1]) for o in cls25.values()) == 4
    print("stratum 7: class offsets and the identity sum f_k = 4: checked (also stratum 25, line c1 = 1/2)")


def g(s, D, e):
    return max(abs(s) / 2, abs(D + e * s / 2))


def min_window(D, e, d):
    return min(sum(g(d + j + i, D, e) for i in range(4)) for j in range(-4, 2))


def check_window_lemma():
    h = F(1, 40)
    worst = F(0)
    for kD in range(0, 34):                    # D = 0 .. 33/40 covers [0, 81/100] with half-step margin
        D = kD * h
        for ke in range(0, 41):
            e = ke * h
            for kd in range(0, 41):
                m = min_window(D, e, kd * h)
                worst = max(worst, m)
    assert worst <= 4 - 8 * h, worst            # Lipschitz margin 8 * (1/40) = 1/5
    # Non-vacuity: at D = 1, e = 0 some delta passes every window (sums >= 4).
    assert any(all(sum(g(F(kd, 40) + j + i, F(1), F(0)) for i in range(4)) >= 4 for j in range(-6, 6))
               for kd in range(40))
    print(f"one-dimensional lemma: max over the grid of the best window sum = {worst} (<= 4 - 1/5); "
          f"hence D >= 81/100 in the residual pattern")


if __name__ == '__main__':
    check_structure()
    check_offsets()
    check_window_lemma()
    print("COORDINATE SECTIONS PASS")
