"""v15: fast bench for the VI neutralisation design (no EVOware simulation): recipe wells (eluate + acid + Tris) made
with the simulator's per-tip pipetting errors and read with the probe's noise; the Tris-equation fit of the analysis
file; the validation wells at the process pool; the refit. Compares the v14 design (Cycle 1 at the most and least acid
targets, Tris at 1.06-1.60 x f, validation from the line between them) with v15 variants.
    python tools/vi_lab.py 200 A B C"""
import sys, os, numpy as np
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem
from trim_lab import PRODS, Robot

ACID_T, TRIS_T = (3.6, 3.4, 3.2, 3.5), (7.2, 7.4, 7.6, 7.4)


def true_acid(el, T, ph):
    return brentq(lambda r: chem.ph_of(chem.mix([(1.0, el), (r, chem.ACID)]), T) - ph, 1e-6, 30)


def true_tris(el, T, R, ph):
    """mL Tris per mL of product acidified with R (mL acid per mL eluate) to reach ph."""
    b = brentq(lambda B: chem.ph_of(chem.mix([(1.0, el), (R, chem.ACID), (B, chem.TRIS)]), T) - ph, 1e-6, 30)
    return b / (1 + R)


def make(rob, el, T, rp, t, tip):
    e, a, b = M.recipe_well(rp, t)
    return chem.ph_of(chem.mix([(rob.tip(tip, e), el), (rob.tip(tip, a), chem.ACID), (rob.tip(tip, b), chem.TRIS)]), T)


def te_pk_groups(sets):
    fit = M.fit_tris_joint(sets)
    if not fit or fit[2] > M.TRIS_FIT_OK or not (7.6 <= fit[0] <= 9.2) or any(x is None for x in fit[1]):
        return None
    return fit


def local_res(ts, ys, te, pk, T):
    """Observed minus Tris-equation pH at the readings of one fraction, interpolated (in pH) at T."""
    ts, ys = np.asarray(ts), np.asarray(ys)
    m = (ys >= M.TRIS_REGION) & (ts > te)
    if m.sum() < 2:
        return 0.0
    yh = pk + np.log10(ts[m] / te - 1)
    r = ys[m] - yh
    o = np.argsort(ys[m])
    return float(np.interp(T, ys[m][o], r[o]))


def run(name, seed, design):
    el, T = PRODS[name]
    rng = np.random.default_rng(seed)
    rob = Robot(rng)
    ratios = [true_acid(el, T, p) * (1 + rng.normal(0, 0.01)) for p in ACID_T]          # the found acid ratios
    rp = ratios[3]
    if design['cols'] == 'v14':
        hi, lo = max(ratios), min(ratios)
        cols = [hi, lo]
    else:
        far = max(ratios, key=lambda r: abs(r / (1 + r) - rp / (1 + rp)))
        cols = [rp, far]
    sets, data = [], []
    for k, r in enumerate(cols):
        f = r / (1 + r)
        te0 = f + (1 - f) * design.get('X', 0.0) * (1 + rng.normal(0, 0.3)) if design.get('use_X') else f
        ts = [te0 * x for x in design['levels']]
        ys = [round(make(rob, el, T, r, t, j) + rng.normal(0, 0.02), 2) for j, t in enumerate(ts)]
        sets.append((np.array(ts), np.array(ys))); data.append((f, ts, ys))
    fit = te_pk_groups(sets)
    if fit is None:
        return None
    pk, tes = fit[0], fit[1]
    fs = [d[0] for d in data]
    a_, b_ = M.te_line(list(zip(fs, tes)))
    fp = rp / (1 + rp)
    if design['cols'] == 'v14':
        te = a_ + b_ * fp
    else:
        te = tes[0] + b_ * (fp - fs[0])
    plan = []
    for Tt in TRIS_T:
        corr = design.get('res', 0.0) * local_res(data[0][1], data[0][2], tes[0], pk, Tt) if design['cols'] != 'v14' else 0.0
        plan.append(te * (1 + 10 ** (Tt - corr - pk)))
    val = [round(make(rob, el, T, rp, t, j) + rng.normal(0, 0.02), 2) for j, t in enumerate(plan)]
    plan_err = [chem.ph_of(chem.mix([(1.0, el), (rp, chem.ACID), (t * (1 + rp), chem.TRIS)]), T) - Tt for t, Tt in zip(plan, TRIS_T)]
    # refit with the validation wells (a group of its own) -> stored ratios at the process pool
    fit2 = te_pk_groups(sets + [(np.array(plan), np.array(val))])
    stored_err = []
    if fit2:
        pk2, tes2 = fit2[0], fit2[1]
        for Tt in (7.2, 7.4, 7.6):
            corr = design.get('res2', 0.0) * local_res(plan, val, tes2[-1], pk2, Tt)
            t = tes2[-1] * (1 + 10 ** (Tt - corr - pk2))
            stored_err.append(chem.ph_of(chem.mix([(1.0, el), (rp, chem.ACID), (t * (1 + rp), chem.TRIS)]), T) - Tt)
    return dict(val=[v - Tt for v, Tt in zip(val, TRIS_T)], plan=plan_err, stored=stored_err)


DESIGNS = {
    'v14': dict(cols='v14', levels=(1.06, 1.16, 1.32, 1.60)),
    'v14-L1.09': dict(cols='v14', levels=(1.09, 1.16, 1.32, 1.60)),
    'v14-L1.12': dict(cols='v14', levels=(1.12, 1.20, 1.35, 1.60)),
    'v14-L1.03': dict(cols='v14', levels=(1.03, 1.10, 1.25, 1.60)),
    'proc': dict(cols='proc', levels=(1.06, 1.16, 1.32, 1.60)),
    'proc+dense': dict(cols='proc', levels=(1.04, 1.09, 1.15, 1.30)),
    'proc+dense+res': dict(cols='proc', levels=(1.04, 1.09, 1.15, 1.30), res=1.0),
    'proc+dense+res/2': dict(cols='proc', levels=(1.04, 1.09, 1.15, 1.30), res=0.5, res2=0.5),
    'proc+dense+res2': dict(cols='proc', levels=(1.04, 1.09, 1.15, 1.30), res=1.0, res2=1.0),
}

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    names = [a for a in sys.argv[2:] if a in PRODS] or list(PRODS)
    designs = [a for a in sys.argv[2:] if a in DESIGNS] or list(DESIGNS)
    for name in names:
        print(f"== {name}")
        for d in designs:
            rows = [run(name, s, DESIGNS[d]) for s in range(N)]
            rows = [r for r in rows if r]
            v = np.array([r['val'] for r in rows]); p = np.array([r['plan'] for r in rows])
            s_ = np.array([r['stored'] for r in rows if len(r['stored']) == 3])
            rms = lambda x: np.sqrt(np.mean(x ** 2, axis=0))
            print(f"   {d:18s} n={len(rows):3d} | wells RMS {np.round(rms(v), 3)} max {np.round(np.abs(v).max(0), 2)} | plan RMS {np.round(rms(p), 3)} "
                  f"| stored 7.2/7.4/7.6 RMS {np.round(rms(s_), 3) if len(s_) else '-'}", flush=True)
