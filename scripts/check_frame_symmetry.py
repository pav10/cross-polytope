"""Validation of the frame-symmetry option of Stage B1 (enumerate_face_data(frame_group=...)), on known answers.
Usage: python3 scripts/check_frame_symmetry.py facet [on|off]   exact facet-only run on zeta's S (frame vectors
                                                                first; static, determine, prefilter), ~5 min
       python3 scripts/check_frame_symmetry.py compare           compare the two facet-only runs (needs both)
       python3 scripts/check_frame_symmetry.py knuth [probes]    Knuth estimates (full tree) with and without the
                                                                group: zeta's S and the three A4 sets (stratum 0)
'compare' checks that the leaf sets agree up to G_S x AMB, in both directions, and that zeta's datum is found.
Outputs: results/frame_symmetry_{facet_on,facet_off,knuth}.json (+ .pkl for the facet runs)."""
import sys, os, json, pickle, random, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
os.chdir(os.path.join(os.path.dirname(__file__), '..'))
from cp4.faces import (enumerate_face_data, estimate_face_data, frame_stabilizer, contact_action, act_sign, AMB,
                       _frame_order)
from cp4.lattice import B_of_z, ZETA

z = json.load(open('data/zeta.json'))
S0 = [tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']]
order = _frame_order(S0, 4)
S = [S0[i] for i in order]                            # frame vectors first, for both runs
strata = json.load(open('data/strata.json'))
G = frame_stabilizer(S, strata[11]['residues'], 4)
assert len(G) == 24 and ((0, 1, 2, 3), (1, 1, 1, 1)) in G

X0 = np.array(z['maxdet_frame_basis_coords'], float) @ np.array(B_of_z(ZETA), float)
pm = z['signed_perm_to_canonical']['perm']; sg = z['signed_perm_to_canonical']['signs']
X = np.array([sg[k] * X0[pm[k]] for k in range(4)])
zsig = {i: tuple(int(np.sign(round(v, 9))) for v in (np.array(c) / 4) @ X) for i, c in enumerate(S)}
ACT = [contact_action(S, g) for g in G]


def canon(signs):
    """Lex-min of a full face datum over G_S x AMB."""
    best = None
    for ca in ACT:
        t = [tuple(e * x for x in signs[j]) for (j, e) in ca]
        for q in AMB:
            u = tuple(act_sign(q, r) for r in t)
            if best is None or u < best:
                best = u
    return best


mode = sys.argv[1]
if mode == 'facet':
    on = sys.argv[2] == 'on'
    leaves, st = enumerate_face_data(S, 4, dynamic=False, signvecs='facet', determine=True, prefilter=True,
                                     frame_group=G if on else None)
    tag = 'on' if on else 'off'
    pickle.dump({'leaves': leaves, 'stats': st}, open(f'results/frame_symmetry_facet_{tag}.pkl', 'wb'))
    out = {k: st[k] for k in ('nodes', 'complete', 'lp_dead', 'prefilter_dead', 'frame_sym_dead', 'seconds')}
    out['leaves'] = len(leaves)
    json.dump(out, open(f'results/frame_symmetry_facet_{tag}.json', 'w'), indent=1)
    print(json.dumps(out))
elif mode == 'compare':
    A = pickle.load(open('results/frame_symmetry_facet_off.pkl', 'rb'))
    B = pickle.load(open('results/frame_symmetry_facet_on.pkl', 'rb'))
    assert A['stats']['complete'] and B['stats']['complete']
    ca = {canon(L['signs']) for L in A['leaves']}
    cb = {canon(L['signs']) for L in B['leaves']}
    zc = canon(zsig)
    print(f"off: {len(A['leaves'])} leaves, {A['stats']['nodes']} nodes; on: {len(B['leaves'])} leaves, "
          f"{B['stats']['nodes']} nodes; orbit classes off {len(ca)}, on {len(cb)}; equal: {ca == cb}; "
          f"zeta found off {zc in ca}, on {zc in cb}; node ratio {A['stats']['nodes'] / B['stats']['nodes']:.2f}")
    assert ca == cb and zc in cb
    print('FRAME SYMMETRY CHECK OK')
elif mode == 'knuth':
    probes = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    res = {}
    A4 = json.load(open('results/stageA_0.json'))['contact_sets']
    cases = [('zeta', S, 4, strata[11]['residues'])] + [(f'A4_{k}', [tuple(c) for c in T], 1, strata[0]['residues'])
                                                         for k, T in enumerate(A4)]
    for name, T, p, R in cases:
        g = frame_stabilizer(T, R, p)
        for on in (False, True):
            e = estimate_face_data(T, p, probes, random.Random(7), dynamic=False, determine=True, prefilter=True,
                                   frame_group=g if on else None)
            res[f'{name}_{"on" if on else "off"}'] = dict(e, group_order=len(g))
            print(name, 'on' if on else 'off', f"|G|={len(g)}", {k: round(e[k], 1) for k in ('nodes_mean', 'nodes_median')},
                  flush=True)
    json.dump(res, open('results/frame_symmetry_knuth.json', 'w'), indent=1)
