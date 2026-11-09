"""Measure the facet-class layer prune inside Stage B1 on zeta's S (DEAD_ENDS B16).
Prune: if the decided contacts whose sign vectors put cX on the closed facet +-sigma generate Gamma over Z,
then sigma.Lambda is in 2Z and det >= 5/6 (FACET_LAYERS §1), so the node can be discarded without an LP.
faces._children is patched in-process only.  Static order, determine, covector pre-filter (the G6 mode).
Usage: python3 scripts/measure_layer_prune.py facet|full on|off [seconds]
Output: one JSON line (nodes, finished depth-3 subtrees, leaves, LP and layer rejections, seconds)."""
import sys, os, json, time
from itertools import product
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
import cp4.faces as faces
def hnf_index(rows, n=4):
    """Index in Z^n of the Z-span of integer rows, or 0 if rank < n (integer row reduction)."""
    M = [list(r) for r in rows]
    piv_rows = []
    for col in range(n):
        # gcd-reduce column col among remaining rows
        while True:
            nz = [r for r in M if r[col] != 0]
            if len(nz) <= 1: break
            nz.sort(key=lambda r: abs(r[col]))
            p = nz[0]
            for r in nz[1:]:
                q = r[col] // p[col]
                for k in range(n): r[k] -= q * p[k]
            M = [r for r in M if any(r)]
        nz = [r for r in M if r[col] != 0]
        if not nz: return 0
        piv_rows.append(abs(nz[0][col])); M = [r for r in M if r is not nz[0]]
    out = 1
    for d in piv_rows: out *= d
    return out
SIGMAS = [s for s in product((1, -1), repeat=4) if s[0] == 1]
def below(s, sig):  # sign vector s on closed facet sig or -sig
    return all(a in (0, b) for a, b in zip(s, sig)) or all(a in (0, -b) for a, b in zip(s, sig))
def layer_hit(signs, S, p):
    """signs: {i: sign vector}; S: numerator vectors over p. True if some facet class generates Gamma."""
    target = p ** 3
    for sig in SIGMAS:
        cls = [S[i] for i, s in signs.items() if below(s, sig)]
        if len(cls) >= 4 and hnf_index(cls) == target:
            return sig
    return None

mode = sys.argv[1]; use = sys.argv[2] == 'on'; tl = float(sys.argv[3]) if len(sys.argv) > 3 else None
z = json.load(open('data/zeta.json')); S = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
orig = faces._children; cnt = {'layer_dead': 0}
def children(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter=None):
    res = orig(C, group, eqs, stricts, signs, dynamic, vecs, stats, colfilter)
    if res is None or not use: return res
    ci, kids = res; keep = []
    for kid in kids:
        sg = dict(signs); sg[ci] = kid[0]
        if layer_hit(sg, S, 4) is not None: cnt['layer_dead'] += 1
        else: keep.append(kid)
    return ci, keep
faces._children = children
leaves, st = faces.enumerate_face_data(S, 4, dynamic=False, signvecs='facet' if mode == 'facet' else 'all',
                                       determine=True, prefilter=True, time_limit=tl)
print(json.dumps({"done": len(st["done_prefixes"]), "mode": mode, 'layer': use, 'complete': st['complete'], 'nodes': st['nodes'], 'leaves': len(leaves),
                  'per_depth': st['per_depth'], 'lp_dead': st['lp_dead'], 'layer_dead': cnt['layer_dead'],
                  'seconds': round(st['seconds'], 1)}))
