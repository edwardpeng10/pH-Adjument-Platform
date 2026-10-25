"""v15: fast bench for the acid titration alone (Acid 0-3, the analysis file's own commands, the simulator's pipetting and
probe errors): true pH that the stored (refit) acid ratios give at the 4 targets minus target, and the validation wells.
    python tools/acid_bench.py tools v15 100 A B C          (DRIFT=0.04: the probe 0.04 high at plate check 3)
The first argument is a tools folder whose ../python holds the analysis file; the second says how its Cycle 1 is typed
(v15: A1 product only + 7 wells; v14: the eluate's pH + 4 wells)."""
import sys, os, json, tempfile, shutil, numpy as np
TOOLS, VER = sys.argv[1], sys.argv[2]
sys.path.insert(0, TOOLS); sys.path.insert(0, os.path.join(TOOLS, '..', 'python'))
import run_lab as R, chem
M = R.M
TG = (3.6, 3.4, 3.2, 3.5)

def one(args):
    prod, seed = args
    el, T = R.PRODS[prod]
    work = tempfile.mkdtemp(prefix='ab_')
    M.WORKING_DIR, M.LOG_FILE = work, os.path.join(work, 'python_log.txt')
    try:
        for name, text in (('titration_params.csv', '1,4,3.6,3.4,3.2,3.5,0.2'), ('tris_params.csv', '7.2,7.4,7.6,7.4'),
                           ('options.csv', '0,0'), ('batch_info.csv', f'lab,{prod},bench')):
            open(M.P(name), 'w').write(text + '\n')
        lab = R.Lab(el, T, seed)
        M.acid0()
        c1 = lab.plan('plan_c1.csv')
        if VER == 'v15':
            R.write('acid_c1_measured.csv', lab.read(c1, 8))           # A1 product only first
        else:
            R.write('acid_c1_measured.csv', [round(chem.ph_of(el, T) + lab.rng.normal(0, R.READ_SD), 2)] + lab.read(c1, 4))
        M.acid1()
        c2r = lab.read(lab.plan('plan_c2.csv'), 8); R.write('cycle2_measured.csv', c2r)
        M.acid2()
        D = float(os.environ.get('DRIFT', 0))
        av = lab.plan('plan_val.csv'); vr = [round(x + D, 2) for x in lab.read(av, 4)]
        if VER == 'v15':
            vr += [round(x + lab.rng.normal(0, 0.01) + D, 2) for x in c2r[:4]]
        R.write('validation_results.csv', vr)
        M.acid3()
        corr = M.nums('corr_ratios.csv')[:4]
        stored = [chem.ph_of(chem.mix([(1.0, el), (r, chem.ACID)]), T) - t for r, t in zip(corr, TG)]
        wells = [lab.ph(w) - t for w, t in zip(av, TG)]
        return prod, stored, wells
    except Exception as e:
        return prod, None, repr(e)
    finally:
        shutil.rmtree(work, ignore_errors=True)

if __name__ == '__main__':
    from multiprocessing import Pool
    N = int(sys.argv[3]); prods = sys.argv[4:]
    with Pool(2) as p:
        res = p.map(one, [(pr, 5000 + k) for pr in prods for k in range(N)])
    for pr in prods:
        rr = [r for r in res if r[0] == pr]
        bad = [r[2] for r in rr if r[1] is None]
        s = np.abs(np.array([r[1] for r in rr if r[1] is not None]))
        w = np.abs(np.array([r[2] for r in rr if r[1] is not None]))
        f = lambda x: f"RMS {np.sqrt(np.mean(x**2)):.4f}, median {np.median(x):.4f}, 95% {np.percentile(x, 95):.4f}, max {x.max():.4f}"
        print(f"{VER} {pr}: {len(s)} runs{' (' + str(len(bad)) + ' failed: ' + bad[0][:80] + ')' if bad else ''} | stored ratios {f(s)} | validation wells {f(w)}", flush=True)
