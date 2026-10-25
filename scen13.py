"""v13 simulator scenarios: pools above the eluate's pH made with Tris only, from a base titration (Tris straight into
the eluate); the VI neutralisation data (Tris into the acidified eluate) kept as its own optional experiment.
Pipetting: 0.7% CV, 0.5% bias per tip, plus 0.3 uL random and 0.5 uL fixed per tip (matters for 20-50 uL of Tris).
    python scen13.py A A_base A_typo B noBase C six none bad refused range fail
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
    sc = dict(eluate=el, T=T, seed=seed, cv=0.007, tip_bias=0.005, read_sd=0.02, abs_sd=0.0003, abs_bias=0.0005, answers={})
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
        bst = json.load(open(os.path.join(s.data, 'run_state.json'))).get('base', {})
    except Exception:
        bst = {}
    base = [(t, chem.ph_of(s.cont(28, 1, 6, k).comp, T)) for k, (t, r) in enumerate(zip((7.2, 7.4, 7.6, 7.4), bst.get('ratios', [0] * 4)))
            if r and s.cont(28, 1, 6, k).vol > 1]                      # wells at targets above the eluate's pH only
    try:
        st = json.load(open(os.path.join(s.data, 'run_state.json')))
    except Exception:
        st = {}
    pools = []
    for pl in st.get('pools', []):
        lane = s.cont(*sim.LPH, pl['lane'])
        true = chem.ph_of(lane.comp, T) if lane.vol > 1 else None
        pools.append(f"L{pl['lane']} {pl['target']:g} {pl['route']}: " + (f"{true:.3f} ({true - pl['target']:+.3f}) {lane.vol:.1f} mL, readings {pl['readings']}" if true else pl.get('status', '')))
    typed = sum(1 for q in s.questions if q[1].startswith('rd_'))
    print(f"== {name}: errors {len(r['errors'])}, warnings {len(r['warnings'])}, python ERROR lines {len(r['py_errors'])}, messages {len(s.prompts)}, entries typed {typed}")
    for e in (r['errors'] + r['py_errors'])[:5]: print('   !!', e[:170])
    print('   acid validation (true pH - target):', ' '.join(f"{p - t:+.3f}" for t, p in acid) or '-',
          '| VI Tris validation:', ' '.join(f"{p - t:+.3f}" for t, p in tris) or '-',
          '| base validation:', ' '.join(f"{p - t:+.3f}" for t, p in base) or '-')
    print('   pools:', ' | '.join(pools) or '-')
    print('   robot minutes (simulated / estimate):', ', '.join(timing(s)), '| total', round(s.clock / 60, 1))
    log = open(os.path.join(s.data, 'python_log.txt')).read().splitlines()
    for l in log:
        if any(k in l for k in ('left out', 'cannot be made', 'not usable', 'Custom pools possible', 'Pool request not used',
                                'Base fit', 'Base refit', 'aimed', 'than 8 tips', 'Base titration stopped')):
            print('   log:', l[l.find(' ', 20) + 1:][:170])
    for p in s.prompts:
        t = next((x for x in (p if isinstance(p, tuple) else (p,)) if isinstance(x, str) and len(x) > 30), '')
        if t.startswith('RATIOS STORED') or t.startswith('CUSTOM POOLS - next') or t.startswith('CUSTOM POOLS - the request'):
            print('   msg:', t.replace('\f\x07', ' / ')[:640])
    return s

if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'A_base', 'A_typo', 'B', 'noBase', 'C', 'six', 'none', 'bad', 'refused', 'range', 'fail']
    both = {'neut': 1, 'base_on': 1}
    old = dict(both, pool_lanes='1,3,5', pool_phs='3.2,3.4,3.6')
    if 'A' in which: summary('product A (7 Oct-like), VI + base titrations, then pools 1,3,5 at 3.2/3.4/3.6', run('A', 1, answers=old))
    if 'A_base' in which: summary('product A, base titration only, Tris pools 2,4,6 at 7.0/7.4/8.0', run('A', 11, answers={'neut': 0, 'base_on': 1, 'pool_lanes': '2,4,6', 'pool_phs': '7.0,7.4,8.0'}))
    if 'A_typo' in which: summary('product A, acid Cycle 2 and base Cycle 2 readings typed 0.25 off', run('A', 2, answers=dict(both, pool_lanes='1,2', pool_phs='3.2,7.4'), perturb={'rd_acid_c2': (5, -0.25), 'rd_base_c2': (2, 0.25)}))
    if 'B' in which: summary('product B (30 Sep), VI + base, lanes 2,4,6 at 3.3/3.5/7.4 (7.4 with Tris only)', run('B', 3, answers=dict(both, pool_lanes='2,4,6', pool_phs='3.3,3.5,7.4')))
    if 'noBase' in which: summary('product B, no base titration: 7.4 refused, asked again, 2,4 at 3.3/3.5', run('B', 9, answers={'neut': 1, 'base_on': 0, 'pool_lanes': ['2,4,6', '2,4'], 'pool_phs': ['3.3,3.5,7.4', '3.3,3.5']}))
    if 'C' in which: summary('2 Oct product (eluate 7.43, holds Tris), VI + base, pools 1 at 7.2 (acid), 2 at 8.0 (Tris)', run('C', 4, answers=dict(both, pool_lanes='1,2', pool_phs='7.2,8.0')))
    if 'six' in which: summary('2 Oct product, base only, six pools 3.2-8.2', run('C', 5, answers={'neut': 0, 'base_on': 1, 'pool_lanes': '1,2,3,4,5,6', 'pool_phs': '3.2,3.4,3.6,4.0,7.6,8.2'}))
    if 'none' in which: summary('product B, acid titration only, no pools', run('B', 12, answers={'neut': 0, 'base_on': 0, 'pool_lanes': '0'}))
    if 'bad' in which: summary('product A, pool request typed wrong twice', run('A', 6, answers={'neut': 0, 'base_on': 0, 'pool_lanes': ['1,7', '1,3', '1,3'], 'pool_phs': ['3.2', '3.2', '3.3,3.6']}))
    if 'refused' in which: summary('product B, base on: pH 9.9 above the Tris range refused, then lane 1 at 3.4', run('B', 7, answers={'neut': 0, 'base_on': 1, 'pool_lanes': ['1,2', '1'], 'pool_phs': ['3.4,9.9', '3.4']}))
    if 'range' in which: summary('product A, pH 2.0 below the acid range refused three times: no pools', run('A', 10, answers={'neut': 0, 'base_on': 0, 'pool_lanes': ['1', '1', '1'], 'pool_phs': ['2.0', '2.0', '2.0']}))
    if 'fail' in which: summary('product A, acid Cycle 1 typed wrong 3 times', run('A', 8, answers=old, typed={'rd_acid_c1': ['1', '2', '3']}))
