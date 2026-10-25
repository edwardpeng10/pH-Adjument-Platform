"""v17 simulator scenarios (the dilute Tris for the base wells, Product 2W RIGHT well; `E` = the 22 Oct 2026 eluate,
fitted to all 44 readings of that run; `E16` = the same with Dilute_Tris_M 0, all 1 M Tris as v16; `dil_off` = product B
with Dilute_Tris_M 0; `dil_err` = product E with the 0.1 M Tris made 2% weak). v16: tips 5-8 take a second portion on
every trip of two or more; the Acid titration plate in grid 28 position 3; `four` = Use_8_tips 0. v15: the product-only
well and the low acid wells in Cycle 1, the optional base titrations, the D2 drift check, the ratio map, the file-version
check, pools (acid below the eluate's pH, Tris above it).
Pipetting (v17, from the 21-22 Oct runs): 0.7% CV, 0.5% bias per tip, plus 1.4 uL random per dispense and 0.5 uL fixed
per tip (v16: 0.3 uL random - the small Tris dispenses of the real robot scatter ~1.3-1.6 uL).
    python scen17.py A A_base A_typo B C six none bad refused range fail oldparse oldfiles acid_only vi_only base_only drift
                     oct13 four E E16 dil_off dil_err
Each scenario prints errors, the acid, VI and base validation errors (true pH of the validation wells), every pool's
true final pH against its target, the readings typed, and each time estimate against the simulator's own count."""
import sys, os, json, numpy as np, sim17 as sim, chem
HERE = os.path.dirname(os.path.abspath(__file__))
ESC = os.path.join(HERE, '..', 'esc', 'EDP_Titration_python_ML_hybrid.esc')
SCRIPTS = [os.path.join(HERE, '..', 'python', f) for f in sim.SCRIPTS]
qc = np.load(os.path.join(HERE, 'eluate_C.npy'))
EL_C, T_C = dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4])
qe = np.load(os.path.join(HERE, 'eluate_E.npy'))          # v17: 22 Oct 2026 (AZD0292 DurA eluate, pH 4.55), all 44 readings
EL_E, T_E = dict(Ac=qe[0], Na=max(qe[1], 0.0), Cl=max(-qe[1], 0.0), Pc=qe[2], Ph=qe[3], Tr=qe[4]), float(qe[5])

def run(prod, seed, **kw):
    el, T = (EL_C, T_C) if prod == 'C' else (EL_E, T_E) if prod == 'E' else sim.product(prod)
    sc = dict(eluate=el, T=T, seed=seed, cv=0.007, tip_bias=0.005, read_sd=0.02, abs_sd=0.0014, abs_bias=0.0005, answers={})
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
    acid = [(t, chem.ph_of(s.cont(*sim.PLATE_OF['Acid'], 4, k).comp, T)) for k, t in enumerate((3.6, 3.4, 3.2, 3.5)) if s.cont(*sim.PLATE_OF['Acid'], 4, k).vol > 1]
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
        if any(k in l for k in ('Drift check', 'Eluate pH', 'left out', 'cannot be made', 'not usable', 'Custom pools possible', 'Pool request not used',
                                'Base fit', 'Base refit', 'aimed', 'than 8 tips', 'Base titration stopped')):
            print('   log:', l[l.find(' ', 20) + 1:][:170])
    for p in s.prompts:
        t = next((x for x in (p if isinstance(p, tuple) else (p,)) if isinstance(x, str) and len(x) > 30), '')
        if t.startswith('RATIOS STORED') or t.startswith('CUSTOM POOLS - next') or t.startswith('CUSTOM POOLS - the request'):
            print('   msg:', t.replace('\f\x07', ' / ')[:640])
    return s

def map_summary(s):
    """The ratio map against the truth (chemistry model of the simulated eluate): how far each route is, in pH."""
    import pandas as pd
    try:
        df = pd.read_csv(os.path.join(s.data, 'ratio_map_run.csv'))
    except Exception:
        return print('   map: none')
    el, T = s.sc['eluate'], s.sc['T']
    out = []
    for route in ('acid', 'base', 'neutralise'):
        d = df[df.Route == route]
        if d.empty:
            continue
        errs = []
        for _, r in d.iterrows():
            A = r.Acid_per_mL_eluate if route == 'neutralise' else (r.Per_mL_eluate if route == 'acid' else 0.0)
            B = r.Per_mL_eluate if route != 'acid' else 0.0
            liq = chem.mix([(1.0, el), (A, chem.ACID), (B, chem.TRIS)])
            errs.append(chem.ph_of(liq, T) - r.Target_pH)
        e = np.abs(np.array(errs))
        sel = d.Basis.isin(['validated', 'fit']).to_numpy()
        out.append(f"{route} {len(d)} entries: true pH within {np.percentile(e, 50):.3f} (median) / {e.max():.3f} (worst)"
                   + (f", fitted/validated part worst {e[sel].max():.3f}" if sel.any() else ''))
    print('   map:', ' | '.join(out))
    chk = d = None
    log = open(os.path.join(s.data, 'python_log.txt')).read()
    for l in log.splitlines():
        if l.startswith('Ratio map') or l.startswith('Chemistry model'):
            print('   log:', l[:200])

if __name__ == '__main__':
    which = sys.argv[1:] or ['A', 'A_base', 'A_typo', 'B', 'C', 'six', 'none', 'bad', 'refused', 'range', 'fail', 'oldparse', 'oldfiles',
                             'acid_only', 'vi_only', 'base_only', 'drift', 'oct13', 'four', 'E', 'E16', 'dil_off', 'dil_err']
    old = dict(pool_lanes='1,3,5', pool_phs='3.2,3.4,3.6')
    def go(name, prod, seed, **kw):
        s = run(prod, seed, **kw)
        summary(name, s)
        map_summary(s)
        return s
    if 'A' in which: go('product A (7 Oct-like), pools 1,3,5 at 3.2/3.4/3.6', 'A', 1, answers=old)
    if 'A_base' in which: go('product A, Tris pools 2,4,6 at 7.0/7.4/8.0', 'A', 11, answers={'pool_lanes': '2,4,6', 'pool_phs': '7.0,7.4,8.0'})
    if 'A_typo' in which: go('product A, acid Cycle 2 and base Cycle 2 readings typed 0.25 off, pools 3.2 + 7.4 Tris', 'A', 2, answers=dict(pool_lanes='1,2', pool_phs='3.2,7.4'), perturb={'rd_acid_c2': (5, -0.25), 'rd_base_c2': (2, 0.25)})
    if 'B' in which: go('product B (30 Sep / 12 Oct product), lanes 1,2 at 7.0 Tris and 4.0 acid (the 12 Oct request)', 'B', 3, answers=dict(pool_lanes='1,2', pool_phs='7.0,4.0'))
    if 'C' in which: go('2 Oct product (eluate 7.43, holds Tris), pools 1 at 7.2 (acid), 2 at 8.0 (Tris)', 'C', 4, answers=dict(pool_lanes='1,2', pool_phs='7.2,8.0'))
    if 'six' in which: go('2 Oct product, six pools 3.2-8.2', 'C', 5, answers={'pool_lanes': '1,2,3,4,5,6', 'pool_phs': '3.2,3.4,3.6,4.0,7.6,8.2'})
    if 'none' in which: go('product B, no pools (titrations only)', 'B', 12, answers={'pool_lanes': '0'})
    if 'bad' in which: go('product A, pool request typed wrong twice', 'A', 6, answers={'pool_lanes': ['1,7', '1,3', '1,3'], 'pool_phs': ['3.2', '3.2', '3.3,3.6']})
    if 'refused' in which: go('product B: pH 10.5 above the Tris range refused, then lane 1 at 3.4', 'B', 7, answers={'pool_lanes': ['1,2', '1'], 'pool_phs': ['3.4,10.5', '3.4']})
    if 'range' in which: go('product A, pH 2.0 below the acid range refused three times: no pools', 'A', 10, answers={'pool_lanes': ['1', '1', '1'], 'pool_phs': ['2.0', '2.0', '2.0']})
    if 'fail' in which: go('product A, acid Cycle 1 typed wrong 3 times', 'A', 8, answers=old, typed={'rd_acid_c1': ['1', '2', '3']})
    if 'acid_only' in which: go('product B, no base titrations (acid only), pools 3.4 and 5.0', 'B', 15,
                                answers=dict(neut=0, base_on=0, pool_lanes='1,2', pool_phs='3.4,5.0'))
    if 'vi_only' in which: go('product B, VI neutralisation but no base titration: 7.4 refused, then 3.4', 'B', 16,
                              answers=dict(neut=1, base_on=0, pool_lanes=['1,2', '1'], pool_phs=['3.4,7.4', '3.4']))
    if 'base_only' in which: go('product A, base titration but no VI neutralisation, Tris pool 7.4', 'A', 17,
                                answers=dict(neut=0, base_on=1, pool_lanes='3', pool_phs='7.4'))
    if 'drift' in which: go('product B, the probe reads 0.04 high from plate check 3 on (13 Oct pattern): corrected by the column-2 re-read', 'B', 18,
                            answers=dict(pool_lanes='0'), shift={'rd_acid_val': 0.04, 'rd_tris_c1': 0.04, 'rd_base_val': 0.04, 'rd_tris_val': 0.04})
    if 'oct13' in which: go('13 Oct-like (2 Oct product, eluate 7.43): the 5.29 pool now between measured wells, plus 7.6 Tris', 'C', 19,
                            answers=dict(pool_lanes='1,2', pool_phs='5.29,7.6'))
    if 'four' in which: go('product B, Use_8_tips = 0 (tips 1-4 only, as v15), pools 7.0 Tris and 4.0 acid', 'B', 3,
                           answers=dict(pool_lanes='1,2', pool_phs='7.0,4.0'), override={'Use_8_tips': 0})
    if 'E' in which: go('22 Oct eluate (pH 4.55), 0.1 M Tris for the base wells, Tris pools 7.2 and 7.4', 'E', 21,
                        answers=dict(pool_lanes='2,4', pool_phs='7.2,7.4'))
    if 'E16' in which: go('22 Oct eluate, Dilute_Tris_M = 0 (all 1 M Tris, as v16), Tris pools 7.2 and 7.4', 'E', 21,
                          answers=dict(pool_lanes='2,4', pool_phs='7.2,7.4'), override={'Dilute_Tris_M': 0})
    if 'dil_off' in which: go('product B, Dilute_Tris_M = 0 (all 1 M Tris, as v16), pools 7.0 Tris and 4.0 acid', 'B', 3,
                              answers=dict(pool_lanes='1,2', pool_phs='7.0,4.0'), override={'Dilute_Tris_M': 0})
    if 'dil_err' in which: go('22 Oct eluate, the 0.1 M Tris made 2% weak (0.098 M), Tris pool 7.4', 'E', 22,
                              answers=dict(pool_lanes='2', pool_phs='7.4'), dil_true=0.098)
    if 'oldparse' in which:     # the 12 Oct failure: a pool entry the parser cannot read (simulated by a crashing parser)
        go('product B, pool readings parser crashes (12 Oct): must ask again, then stop - never pass silently', 'B', 13,
           answers=dict(pool_lanes='1,2', pool_phs='7.0,4.0'), crash_parse='pools_measured_1.csv')
    if 'oldfiles' in which:     # one helper file from an older version: the run must stop before pipetting
        go('older EDP_parse_readings.py in tecan_data: the version check stops the run at the start', 'B', 14,
           answers=dict(pool_lanes='0'), old_file='EDP_parse_readings.py')
