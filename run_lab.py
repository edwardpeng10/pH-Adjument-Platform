"""v14/v15: fast bench for the three titrations (no EVOware simulation, a few seconds per run). The analysis file's own
commands run in this process, in the script's order (Acid 0-3, Base 0-3, Tris 1-2, Pools range); every well in its
plan files is "made" with the simulator's pipetting errors (per tip: 0.7% CV, 0.5% bias, 0.3 uL random, 0.5 uL fixed)
and "read" from chem.py's truth with the probe's noise (0.02, two decimals).
    python tools/run_lab.py 40 A B C        (runs per product)
Prints, per product: the acid, base and VI-neutralisation validation errors (true pH of those wells - target); the
stored neutralisation ratios of the RATIOS STORED message (from each acid target to the process pH) as the pH they
give; and the ratio map's accuracy per route (true pH of every entry - its target)."""
import sys, os, json, tempfile, shutil, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem, sim11 as sim

qc = np.load(os.path.join(HERE, 'eluate_C.npy'))
PRODS = {'A': sim.product('A'), 'B': sim.product('B'), 'C': (dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4]))}
CV, TIP_BIAS, ABS_SD, ABS_BIAS, READ_SD = 0.007, 0.005, 0.0003, 0.0005, 0.02


class Lab:
    def __init__(self, el, T, seed):
        self.el, self.T = el, T
        self.rng = np.random.default_rng(seed)
        self.bias = self.rng.normal(0, TIP_BIAS, 4)
        self.abs_bias = self.rng.normal(0, ABS_BIAS, 4)

    def disp(self, t, n, ul):
        out = 0.0
        for _ in range(n):
            v = ul / 1000.0
            out += max(v * (1 + self.bias[t] + self.rng.normal(0, CV)) + self.abs_bias[t] + self.rng.normal(0, ABS_SD), 0.0)
        return out

    def plan(self, name):
        """The wells of a plan file: actual (eluate, acid, Tris) mL per well."""
        x = M.nums(name)
        wells = []
        for c in range(len(x) // 15):
            col = x[15 * c:15 * c + 15]
            for t in range(4):
                wells.append(tuple(self.disp(t, int(col[5 * k]), col[5 * k + 1 + t]) for k in range(3)))
        return wells

    def ph(self, w):
        return chem.ph_of(chem.mix([(w[0], self.el), (w[1], chem.ACID), (w[2], chem.TRIS)]), self.T)

    def read(self, wells, n=None):
        return [round(self.ph(w) + self.rng.normal(0, READ_SD), 2) for w in wells[:n] if sum(w) > 0.5]


def write(name, vals):
    with open(M.P(name), 'w') as f:
        f.write(','.join(f'{v:.2f}' for v in vals) + '\n')


def one(prod, seed, keep=None):
    el, T = PRODS[prod]
    work = keep or tempfile.mkdtemp(prefix='runlab_')
    os.makedirs(work, exist_ok=True)
    M.WORKING_DIR, M.LOG_FILE = work, os.path.join(work, 'python_log.txt')
    for name, text in (('titration_params.csv', '1,4,3.6,3.4,3.2,3.5,0.2'), ('tris_params.csv', '7.2,7.4,7.6,7.4'),
                       ('options.csv', '1,1'), ('batch_info.csv', f'lab,{prod},bench')):
        open(M.P(name), 'w').write(text + '\n')
    lab = Lab(el, T, seed)
    M.acid0(); M.base0()
    c1 = lab.plan('plan_c1.csv'); b1 = lab.plan('plan_base_c1.csv')
    if len(c1) >= 8:                                     # v15: A1 product only, then 7 acid wells
        write('acid_c1_measured.csv', lab.read(c1, 8))
    else:                                                # v14: the eluate read in the trough, then 4 wells
        write('acid_c1_measured.csv', [round(chem.ph_of(el, T) + lab.rng.normal(0, READ_SD), 2)] + lab.read(c1, 4))
    write('base_c1_measured.csv', lab.read(b1, 4))
    M.acid1(); M.base1()
    c2w = lab.plan('plan_c2.csv'); c2r = lab.read(c2w, 8); write('cycle2_measured.csv', c2r)
    b2 = lab.plan('plan_base_c2.csv') if M.nums('plan_base_c2.csv') else []
    write('base_c2_measured.csv', lab.read(b2, 4) if b2 else [0, 0, 0, 0])
    M.acid2(); M.base2()
    av = lab.plan('plan_val.csv')
    shift = float(os.environ.get('DRIFT', 0.0))         # probe drift since plate check 2 (v15 drift-check test)
    vr = [round(x + shift, 2) for x in lab.read(av, 4)]
    if M.VERSION >= 'EDP-v15':                           # v15: D2 read again at the end of the acid validation entry
        vr += [round(x + lab.rng.normal(0, 0.01) + shift, 2) for x in c2r[:4]]    # column 2 read again: repeatability 0.01
    write('validation_results.csv', vr)
    tc1 = lab.plan('plan_tris_c1.csv'); write('tris_c1_measured.csv', lab.read(tc1, 8))
    bv = lab.plan('plan_base_val.csv') if M.nums('plan_base_val.csv') else []
    write('base_val_measured.csv', lab.read(bv, 4) if bv else [0, 0, 0, 0])
    M.acid3(); M.tris1(); M.base3()
    tv = lab.plan('plan_tris_val.csv'); write('tris_val_measured.csv', lab.read(tv, 4))
    x = M.nums('plan_tris_val.csv')                     # the planned (error-free) compositions of the VI validation wells
    tv_nom = [tuple(x[5 * k] * x[5 * k + 1 + t] / 1000 for k in range(3)) for t in range(4)] if len(x) >= 15 else []
    M.tris2()
    M.ranges()
    st = M.load_state()
    out = dict(acid_val=[lab.ph(w) - t for w, t in zip(av, (3.6, 3.4, 3.2, 3.5))],
               vi_val=[lab.ph(w) - t for w, t in zip(tv, (7.2, 7.4, 7.6, 7.4))],
               vi_plan=[lab.ph(w) - t for w, t in zip(tv_nom, (7.2, 7.4, 7.6, 7.4))],
               base_val=[lab.ph(w) - t for w, t, r in zip(bv, (7.2, 7.4, 7.6, 7.4), st.get('base', {}).get('ratios', [])) if r],
               vi_ok=(st.get('tris') or {}).get('ok'))
    # the RATIOS STORED neutralisation numbers: per mL of product acidified with the stored acid ratio to each target
    summ = M.nums('run_summary.csv')
    corr = summ[1:5] if len(summ) >= 20 else []
    neut = summ[15:19] if len(summ) >= 20 else []
    out['stored_neut'] = [chem.ph_of(chem.mix([(1.0, el), (R, chem.ACID), (n * (1 + R), chem.TRIS)]), T) - 7.4
                          for R, n in zip(corr, neut) if R > 0 and n > 0]
    import pandas as pd
    try:
        df = pd.read_csv(M.P('ratio_map_run.csv'))
    except Exception:
        df = None
    out['map'] = {}
    if df is not None:
        for route in ('acid', 'base', 'neutralise'):
            d = df[df.Route == route]
            errs = []
            for _, r in d.iterrows():
                A = r.Acid_per_mL_eluate if route == 'neutralise' else (r.Per_mL_eluate if route == 'acid' else 0.0)
                B = r.Per_mL_eluate if route != 'acid' else 0.0
                errs.append((chem.ph_of(chem.mix([(1.0, el), (A, chem.ACID), (B, chem.TRIS)]), T) - r.Target_pH, r.Basis))
            out['map'][route] = errs
    if not keep:
        shutil.rmtree(work, ignore_errors=True)
    return out


def job(args):
    prod, seed = args
    try:
        return prod, seed, one(prod, seed)
    except Exception as e:
        return prod, seed, {'error': repr(e)}


if __name__ == '__main__':
    from multiprocessing import Pool
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    prods = sys.argv[2:] or list(PRODS)
    with Pool(int(os.environ.get('NPROC', 2))) as p:
        res = p.map(job, [(pr, 1000 + k) for pr in prods for k in range(N)])
    for pr in prods:
        rr = [o for p_, k, o in res if p_ == pr]
        bad = [o['error'] for o in rr if 'error' in o]
        rr = [o for o in rr if 'error' not in o]
        def stat(key):
            e = np.abs(np.concatenate([np.array(o[key], float) for o in rr if o[key]] or [np.zeros(0)]))
            return f"RMS {np.sqrt(np.mean(e ** 2)):.3f}, 95% {np.percentile(e, 95):.3f}, max {e.max():.3f}" if len(e) else '-'
        print(f"== {pr}: {len(rr)} runs" + (f", {len(bad)} failed: {bad[:2]}" if bad else '') +
              f"; VI fit not usable in {sum(1 for o in rr if not o['vi_ok'])}")
        print(f"   validation wells: acid {stat('acid_val')} | base {stat('base_val')} | VI {stat('vi_val')}")
        print(f"   RATIOS STORED neutralisation (pH it gives - 7.4): {stat('stored_neut')}")
        for route in ('acid', 'base', 'neutralise'):
            e = [abs(x) for o in rr for x, b in o['map'].get(route, [])]
            ef = [abs(x) for o in rr for x, b in o['map'].get(route, []) if b in ('validated', 'fit')]
            em = [abs(x) for o in rr for x, b in o['map'].get(route, []) if str(b).startswith('chemistry')]
            if e:
                print(f"   map {route}: {len(e)} entries, median {np.median(e):.3f}, 95% {np.percentile(e, 95):.3f}, max {max(e):.3f}"
                      + (f" | fitted/validated 95% {np.percentile(ef, 95):.3f} max {max(ef):.3f}" if ef else '')
                      + (f" | chemistry model 95% {np.percentile(em, 95):.3f} max {max(em):.3f}" if em else ''))
