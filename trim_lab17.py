"""v17: fast bench for Tris pools (no EVOware simulation, ~1 s per run): the base titration as the script makes it
(v17: base wells from the 0.1 M Tris where one portion of it is enough) -> the run's refit base curve -> a 25 mL pool
made with 8 tips of 1 M Tris, aimed as the script aims it -> up to two trims by the analysis file's trim_for (v17: small
Tris trims from the 0.1 M Tris, with eluate when even its smallest 8-tip addition is too much) -> the true final pH.
Pipetting: per tip 0.7% CV, 0.5% bias, ABS_SD random per dispense (default 1.4 uL - the 21-22 Oct wells), 0.5 uL fixed.
    python tools/trim_lab17.py 200 A B E           (DIL=0: no dilute Tris anywhere, as v16; PYDIR: python folder)
v14's trim_lab.py had 0.3 uL random per dispense."""
import sys, os, tempfile, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('PYDIR', os.path.join(HERE, '..', 'python'))); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem, sim11 as sim

qc, qe = np.load(os.path.join(HERE, 'eluate_C.npy')), np.load(os.path.join(HERE, 'eluate_E.npy'))
PRODS = {'A': sim.product('A'), 'B': sim.product('B'), 'C': (dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4])),
         'E': (dict(Ac=qe[0], Na=max(qe[1], 0.0), Cl=max(-qe[1], 0.0), Pc=qe[2], Ph=qe[3], Tr=qe[4]), float(qe[5]))}
TARGETS = {'A': (7.0, 7.2, 7.4, 8.0), 'B': (7.0, 7.2, 7.4, 8.0), 'C': (7.6, 8.2), 'E': (7.0, 7.2, 7.4, 8.0)}
CV, TIP_BIAS, ABS_BIAS, READ_SD = 0.007, 0.005, 0.0005, 0.02
ABS_SD = float(os.environ.get('ABS_SD', 0.0014))
DIL = float(os.environ.get('DIL', 0.1))
V17 = hasattr(M, 'dil_m')


class Robot:
    def __init__(self, rng):
        self.rng = rng
        self.bias = rng.normal(0, TIP_BIAS, 8)
        self.abs_bias = rng.normal(0, ABS_BIAS, 8)

    def tip(self, t, v):
        out = 0.0
        while v > 1e-9:
            p = min(v, 0.98)
            out += max(p * (1 + self.bias[t] + self.rng.normal(0, CV)) + self.abs_bias[t] + self.rng.normal(0, ABS_SD), 0.0)
            v -= p
        return out

    def lane(self, v):                       # v mL over 8 tips
        return sum(self.tip(t, v / 8) for t in range(8)) if v > 0 else 0.0

    def read(self, ph):
        return round(ph + self.rng.normal(0, READ_SD), 2)


def base_titration(el, T, rob):
    """The run's refit base curve, as the script builds it (Cycle 1 blind, Cycle 2 placed, validation, refit)."""
    def make(t, k):
        if V17:
            w = M.base_well(t)                      # (eluate, acid, 1 M Tris, dilute Tris) mL
            e, b1, bd, s_ = w[0], w[2], w[3], M.base_stock(w)
        else:
            e, b1, bd, s_ = 3.0 / (1 + t), 3.0 - 3.0 / (1 + t), 0.0, 1.0
        parts = [(rob.tip(k, e), el), (rob.tip(k, b1), chem.TRIS)] + ([(rob.tip(k, bd), dict(Tr=DIL))] if bd else [])
        return rob.read(chem.ph_of(chem.mix(parts), T)), s_
    def wells(ts):
        out = [make(t, k) for k, t in enumerate(ts)]
        return [o[0] for o in out], [o[1] for o in out]
    ph0 = rob.read(chem.ph_of(el, T))
    c1 = [M.BASE_MAX / M.BASE_SPAN * M.BASE_SPAN ** (k / 3) for k in range(4)]
    p1, m1 = wells(c1)
    lo, hi = M.base_window(ph0)
    tg2 = [lo + (hi - lo) * f for f in M.BASE_C2_SPLIT]
    curve = (lambda v, p, m: M.BaseCurve(v, p, ph0, m)) if V17 else (lambda v, p, m: M.BaseCurve(v, p, ph0))
    cc = curve(c1, p1, m1)
    c2 = [max((cc.ratio_for(t)[0] if cc.model else None) or M.locate_rising(np.array(c1), np.array(p1), ph0, t), 0.002) for t in tg2]
    for k in range(1, 4):
        c2[k] = max(c2[k], c2[k - 1] + max(0.0005, 0.02 * c2[k - 1]))
    p2, m2 = wells(c2)
    cur = curve(c1 + c2, p1 + p2, m1 + m2)
    rv = [x for x in (cur.ratio_for(t)[0] if t > ph0 + 0.02 else None for t in (7.2, 7.4, 7.6, 7.4)) if x]
    pv, mv = wells(rv)
    return curve(c1 + c2 + rv, p1 + p2 + pv, m1 + m2 + mv), ph0


def pool(el, T, curve, target, rob):
    pl = dict(target=target, route='base', A=0.0, E=0.0, B=0.0, W=0.0, readings=[], hist=[])
    sd = curve.sd_ph(target)
    aim = min(max(2.0 * sd - M.POOL_TOL, 0.0), M.BASE_AIM_MAX) if sd is not None else 0.0
    x = curve.ratio_for(target + aim)[0]
    if x is None:
        aim, x = 0.0, curve.ratio_for(target)[0]
    if x is None:
        return None
    E, B = M.POOL_ML / (1 + x), M.POOL_ML - M.POOL_ML / (1 + x)
    pl.update(E=E, B=B)
    true_e, true_b, true_d = rob.lane(E), rob.lane(B), 0.0          # true_d: mL of the dilute Tris in the lane
    liq = lambda: chem.mix([(true_e, el), (true_b, chem.TRIS)] + ([(true_d, dict(Tr=DIL))] if true_d else []))
    used_dil = 0
    for k in range(1, 4):
        m = rob.read(chem.ph_of(liq(), T))
        pl['readings'].append(m)
        pl['hist'].append((pl['B'] / pl['E'], m))
        if abs(m - target) <= M.POOL_TOL + 1e-9 or k == 3:
            break
        out = M.trim_for(pl, m, None, curve)
        da, de, db, how = out[:4]
        dd = out[4] if len(out) > 4 else 0.0
        if da + de + (dd if dd else db) < M.MIN_TRIM_ML - 1e-9:
            break
        pl['E'] += de; pl['B'] += db
        if dd:
            pl['W'] += dd - db; used_dil += 1
            true_d += rob.lane(dd)
        else:
            true_b += rob.lane(db)
        true_e += rob.lane(de)
    return chem.ph_of(liq(), T) - target, len(pl['readings']) - 1, used_dil


def trial(args):
    name, seed = args
    el, T = PRODS[name]
    work = tempfile.mkdtemp(prefix='trimlab_')
    M.WORKING_DIR, M.LOG_FILE = work, os.path.join(work, 'python_log.txt')
    open(os.path.join(work, 'options.csv'), 'w').write(f'1,1,1,{DIL:g}\n')
    rng = np.random.default_rng(seed)
    rob = Robot(rng)
    curve, ph0 = base_titration(el, T, rob)
    return {t: pool(el, T, curve, t, rob) for t in TARGETS[name] if t > ph0 + 0.05}


if __name__ == '__main__':
    from multiprocessing import Pool
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    names = [a for a in sys.argv[2:] if a in PRODS] or list(PRODS)
    for name in names:
        with Pool(int(os.environ.get('NPROC', 2))) as p:
            outs = p.map(trial, [(name, k) for k in range(N)])
        for t in TARGETS[name]:
            r = [o[t] for o in outs if o.get(t)]
            if not r:
                continue
            e = np.array([x[0] for x in r]); n = np.array([x[1] for x in r]); d = np.array([x[2] for x in r])
            print(f"{name} pH {t}: RMS {np.sqrt(np.mean(e ** 2)):.3f}, within 0.03 {np.mean(np.abs(e) <= 0.03):.0%}, within 0.05 "
                  f"{np.mean(np.abs(e) <= 0.05):.0%}, worst {np.abs(e).max():.3f}, trims {n.mean():.1f} (from the dilute Tris {d.mean():.1f})", flush=True)
