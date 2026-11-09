"""Exact small certificate for docs/THEORY_REVIEW.md; no search LPs or floats.

Run from any directory: python3 verify/verify_support_reduction.py
The completeness of the 39-stratum list is the separate G1 prerequisite.
"""
import json
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path


def check():
    strata = json.loads((Path(__file__).resolve().parents[1] / 'data/strata.json').read_text())
    excluded = []
    for s in strata:
        p = s['p']
        residues = {tuple(r) for r in s['residues']}
        # Regenerate the unit shell independently of the stored candidate list.
        shell = set()
        for r in residues:
            for shift in product((-1, 0), repeat=4):
                n = tuple(r[i] + p * shift[i] for i in range(4))
                if sum(map(abs, n)) == p:
                    shell.add(max(n, tuple(-a for a in n)))
        assert shell == set(map(tuple, s['forced'])), s['id']
        disjoint = set()
        for i, j in combinations(range(4), 2):
            signs = {(n[i] * n[j] > 0) - (n[i] * n[j] < 0) for n in shell}
            if {-1, 1} <= signs:
                disjoint.add((i, j))
        unseen = set(range(4))
        components = []
        while unseen:
            component = {unseen.pop()}
            while True:
                added = {j for i in component for j in unseen
                         if tuple(sorted((i, j))) not in disjoint}
                if not added:
                    break
                component |= added
                unseen -= added
            components.append(component)
        if any(len(c) > 2 for c in components):
            continue
        # Minkowski gives det >= 2 for minimum-2 l1 lattices in dimensions 1 and 2.
        bound = Fraction(1, p)
        for c in components:
            index = sum(all(r[j] == 0 for j in range(4) if j not in c) for r in residues)
            bound *= 2 * index
        if bound > Fraction(81, 100):
            excluded.append(s['id'])
            print(f"stratum {s['id']:2d}: p={p:2d}, blocks={sorted(map(len, components))}, det >= {bound}")
    expected = [20, 21, 22, 23, 24, 26, 27, 29, 32, 33, 34, 35, 36, 37, 38]
    assert excluded == expected, excluded
    # D* <= d(1/2) = 207/256 < 81/100: rational separation is sufficient.
    assert Fraction(207, 256) < Fraction(81, 100)
    assert 11 not in excluded
    print(f"SUPPORT REDUCTION PASS: {len(excluded)} excluded, {len(strata)-len(excluded)} remain")


if __name__ == '__main__':
    check()
