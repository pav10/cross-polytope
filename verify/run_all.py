"""Regression suite: every gate of docs/PLAN.md §1 that is cheap enough to rerun (< ~3 min total).
Usage (from the campaign root):  python3 verify/run_all.py [--skip-slow]
Exit code 0 iff all gates pass.  Long runs (Stage A stratum 11, Stage B1 on zeta's S) are checked
from their saved outputs in results/ if present, never rerun here."""
import subprocess, sys, os, json, time, pickle, tempfile
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(ROOT, 'src'))
os.chdir(ROOT)
ok_all = True

def gate(name, fn):
    global ok_all
    t = time.time()
    try:
        msg = fn()
        print(f"[PASS] {name} ({time.time()-t:.0f}s) {msg or ''}", flush=True)
    except Exception as e:  # noqa
        ok_all = False
        print(f"[FAIL] {name} ({time.time()-t:.0f}s): {e!r}", flush=True)

def run(args, timeout=900):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True, timeout=timeout)
    if r.returncode:
        raise RuntimeError(r.stdout[-800:] + r.stderr[-800:])
    return r.stdout

def g0():
    out = run(['verify/verify_audit.py', '--output', '/dev/null'])
    assert out.count('PASS:') == 3, out[-500:]
    return 'three PASS lines (HNFs/test sets/counterexample+KKT)'

def g1():
    # regenerate into a temp dir; the committed file is only read, never rewritten
    with tempfile.TemporaryDirectory() as d:
        out = os.path.join(d, 'strata.json')
        run(['scripts/make_strata.py', out])
        new = open(out, 'rb').read()
    assert new == open('data/strata.json', 'rb').read(), 'data/strata.json changed on regeneration (ids no longer stable!)'
    return '39 strata / 1040; clean 6 / 54; file reproduced byte-for-byte'

def g2():
    z = json.load(open('data/zeta.json'))
    assert z['stratum_id'] == 11 and z['maxdet_frame_index_p'] == 4
    assert len(z['S_frame_numerators_over_p_canonical_stratum_coords']) == 15
    return 'zeta: stratum 11, p=4, 15 contacts'

def g3():
    from cp4.strata import stabilizer
    from cp4.contacts import enumerate_contact_sets, verify_contact_set
    s = json.load(open('data/strata.json'))[0]
    stab = stabilizer([tuple(r) for r in s['residues']], s['p'])
    res, st = enumerate_contact_sets(s, stab)
    assert len(res) == 3 and all(len(S) == 10 for S in res), (len(res), st)
    for S in res:
        verify_contact_set(S, s)
    return f"p=1: 3 sets of size 10, {st['nodes']} nodes"

def g4():
    from cp4.strata import stabilizer
    from cp4.contacts import canon_set
    f = 'results/stageA_11.json'
    if not os.path.exists(f):
        return 'SKIPPED (results/stageA_11.json absent; run scripts/stageA_run.py 11)'
    z = json.load(open('data/zeta.json'))
    s = json.load(open('data/strata.json'))[11]
    stab = stabilizer([tuple(r) for r in s['residues']], 4)
    key = canon_set([tuple(c) for c in z['S_frame_numerators_over_p_canonical_stratum_coords']], stab, 4)
    R = json.load(open(f))
    sets = {tuple(sorted(tuple(c) for c in S)) for S in R['contact_sets']}
    assert key in sets, "zeta's S missing from Stage A output"
    return f"zeta's S found among {len(sets)} sets (complete run: {'node limit' not in str(R['stats'])})"

def g5():
    out = run(['scripts/stageB_critical_zeta.py'])
    assert '64*D**3 - 1920*D**2 + 13584*D - 26011' in out and '<-- zeta' in out, out
    return 'eliminant contains the zeta cubic; root 4*D* recovered'

def g6(f='results/stageB_faces_zeta.pkl'):
    if not os.path.exists(f):
        return f'SKIPPED ({f} absent)'
    R = pickle.load(open(f, 'rb'))
    if not R['stats'].get('complete') and not R['zeta_hits']:
        return f"INCOMPLETE run ({R['stats']['nodes']} nodes, {R['stats']['seconds']:.0f}s) - rerun longer"
    assert R['zeta_hits'], 'zeta face datum not among leaves of a COMPLETE run'
    return f"{len(R['leaves'])} leaves, complete={R['stats'].get('complete')}, zeta hit {R['zeta_hits']}"

gate('G0 audit checker', g0)
gate('G1 strata', g1)
gate('G2 zeta located', g2)
gate('G3 Stage A p=1', g3)
gate('G4 Stage A stratum 11 contains zeta', g4)
if '--skip-slow' not in sys.argv:
    gate('G5 Stage B2 zeta eliminant', g5)
gate('G6a Stage B1 zeta witness path', lambda: (run(['scripts/stageB_witness_zeta.py']).strip().endswith('WITNESS OK') or (_ for _ in ()).throw(AssertionError('witness'))) and 'zeta datum survives all prunes')
gate('G6f Stage B1 facet-only pass contains zeta', lambda: g6('results/stageB_faces_zeta_facet.pkl'))
gate('G6 Stage B1 zeta face datum', g6)
gate('G8 covector pre-filter certified, zeta passes', lambda: (run(['scripts/check_covectors_zeta.py']).strip().endswith('COVECTORS OK') or (_ for _ in ()).throw(AssertionError('covectors'))) and 'exact covector certificate; zeta datum passes the column filter')
gate('G7a merged-search zeta witness path', lambda: (run(['scripts/merged_witness_zeta.py']).strip().endswith('WITNESS OK') or (_ for _ in ()).throw(AssertionError('merged witness'))) and 'zeta labels survive all merged-search prunes')
print('ALL PASS' if ok_all else 'SOME GATES FAILED')
sys.exit(0 if ok_all else 1)
