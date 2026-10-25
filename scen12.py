"""v12 simulator scenarios: the v11 method with the custom pools asked for at the end, after every ratio is stored.
    python scen12.py A A_typo B noTris C six bad refused range fail
Each scenario prints errors, the acid and Tris validation errors (true pH of the validation wells), every pool's
true final pH against its target, the readings typed, and each time estimate against the simulator's own count."""
import sys, os, json, numpy as np, sim11 as sim, chem
HERE = os.path.dirname(os.path.abspath(__file__))
ESC = os.path.join(HERE, '..', 'esc', 'EDP_Titration_python_ML_hybrid.esc')
SCRIPTS = [os.path.join(HERE, '..', 'python', f) for f in sim.SCRIPTS]
qc = np.load(os.path.join(HERE, 'eluate_C.npy'))
EL_C, T_C = dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4])

def run(prod, seed, **kw):
    el, T = (EL_C, T_C) if prod == 'C' else sim.product(prod)
    sc = dict(eluate=el, T=T, seed=seed, cv=0.007, tip_bias=0.005, read_sd=0.02, answers={})
    sc.update(kw)
    return sim.Sim(ESC, sc, SCRIPTS).run()

def timing(s):
    p = os.path.join(s.data, 'timing_log.csv')
    rows = [l.split(',') for l in open(p).read().splitlines()[1:]] if os.path.exists(p) else []
    out = []
    for blk, sec in s.block_sim.items():
        est = next((float(r[3]) for r in rows if r[1] == 'start' and r[2] == blk), None)
        ok = est is not None and abs(est - sec / 60) <= 0.02 * sec / 60 + 0.05
        out.append(f"{blk} {sec / 60:.1f}/{est}" + ('' if ok else ' MISMATCH'))
    return out

def summary(name, s):
    r = s.report()
    T = s.sc['T']
    acid = [(t, chem.ph_of(s.cont(28, 0, 4, k).comp, T)) for k, t in enumerate((3.6, 3.4, 3.2, 3.5)) if s.cont(28, 0, 4, k).vol > 1]
    tris = [(t, chem.ph_of(s.cont(28, 1, 3, k).comp, T)) for k, t in enumerate((7.2, 7.4, 7.6, 7.4)) if s.cont(28, 1, 3, k).vol > 1]
    try:
        st = json.load(open(os.path.join(s.data, 'run_state.json')))
    except Exception:
        st = {}
    pools = []
    for pl in st.get('pools', []):
        lane = s.cont(*sim.LPH, pl['lane'])
        true = chem.ph_of(lane.comp, T) if lane.vol > 1 else None
        pools.append(f"L{pl['lane']} {pl['target']:g}: " + (f"{true:.3f} ({true - pl['target']:+.3f}) {lane.vol:.1f} mL, readings {pl['readings']}" if true else pl.get('status', '')))
    typed = sum(1 for q in s.questions if q[1].startswith('rd_'))
    print(f"== {name}: errors {len(r['errors'])}, warnings {len(r['warnings'])}, python ERROR lines {len(r['py_errors'])}, messages {len(s.prompts)}, entries typed {typed}")
    for e in (r['errors'] + r['py_errors'])[:5]: print('   !!', e[:170])
    print('   acid validation (true pH - target):', ' '.join(f"{p - t:+.3f}" for t, p in acid) or '-',
          '| Tris validation:', ' '.join(f"{p - t:+.3f}" for t, p in tris) or '-')
    print('   pools:', ' | '.join(pools) or '-')
    print('   robot minutes (simulated / estimate):', ', '.join(timing(s)), '| total', round(s.clock / 60, 1))
    log = open(os.path.join(s.data, 'python_log.txt')).read().splitlines()
    for l in log:
        if any(k in l for k in ('left out', 'cannot be made', 'not usable', 'Custom pools possible', 'Pool request not used')):
            print('   log:', l[l.find(' ', 20) + 1:][:170])
    for p in s.prompts:
        t = next((x for x in (p if isinstance(p, tuple) else (p,)) if isinstance(x, str) and len(x) > 30), '')
        if t.startswith('RATIOS STORED') or t.startswith('CUSTOM POOLS - next') or t.startswith('CUSTOM POOLS - the request'):
            print('   msg:', t.replace('\f\x07', ' / ')[:640])
    return s

if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'A_typo', 'B', 'noTris', 'C', 'six', 'bad', 'refused', 'range', 'fail']
    old = {'pool_lanes': '1,3,5', 'pool_phs': '3.2,3.4,3.6', 'neut': 1}
    if 'A' in which: summary('product A (7 Oct-like), Tris titration, then pools 1,3,5 at 3.2/3.4/3.6', run('A', 1, answers=old))
    if 'A_typo' in which: summary('product A, one Cycle 2 reading typed 0.25 low (as on 7 Oct)', run('A', 2, answers=old, perturb={'rd_acid_c2': (5, -0.25)}))
    if 'B' in which: summary('product B, Tris titration, then lanes 2,4,6 at 3.3, 3.5 and 7.4 (acid then Tris)', run('B', 3, answers={'pool_lanes': '2,4,6', 'pool_phs': '3.3,3.5,7.4', 'neut': 1}))
    if 'noTris' in which: summary('product B, no Tris titration: 7.4 refused, asked again, 2,4 at 3.3/3.5', run('B', 9, answers={'pool_lanes': ['2,4,6', '2,4'], 'pool_phs': ['3.3,3.5,7.4', '3.3,3.5'], 'neut': 0}))
    if 'C' in which: summary('2 Oct product, Tris titration, no pools', run('C', 4, answers={'pool_lanes': '0', 'neut': 1}))
    if 'six' in which: summary('2 Oct product, Tris titration, six pools 3.2-8.0', run('C', 5, answers={'pool_lanes': '1,2,3,4,5,6', 'pool_phs': '3.2,3.4,3.6,4.0,7.4,8.0', 'neut': 1}))
    if 'bad' in which: summary('product A, pool request typed wrong twice', run('A', 6, answers={'pool_lanes': ['1,7', '1,3', '1,3'], 'pool_phs': ['3.2', '3.2', '3.3,3.6'], 'neut': 0}))
    if 'refused' in which: summary('product B, pool at 6.0 (between eluate pH and 6.8) refused, then lane 1 at 3.4', run('B', 7, answers={'pool_lanes': ['1,2', '1'], 'pool_phs': ['3.4,6.0', '3.4'], 'neut': 0}))
    if 'range' in which: summary('product A, pH 2.0 below the acid range refused three times: no pools', run('A', 10, answers={'pool_lanes': ['1', '1', '1'], 'pool_phs': ['2.0', '2.0', '2.0'], 'neut': 0}))
    if 'fail' in which: summary('product A, acid Cycle 1 typed wrong 3 times', run('A', 8, answers=old, typed={'rd_acid_c1': ['1', '2', '3']}))
