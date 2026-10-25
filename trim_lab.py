"""v14: fast test bench for Tris pools (no EVOware simulation, about 1 s per run): a base titration with the
simulator's pipetting and probe errors -> the run's refit base curve -> a 25 mL pool made with 8 tips, aimed as the
script aims it -> up to two trims by a chosen rule -> the pool's true final pH.
    python tools/trim_lab.py 200 A B                          the script's rule against v13's
    python tools/trim_lab.py 200 A script shape:0.9:0.9       any rules (see trim_rule); DEBUG=0.07 lists bad cases
v14 result (200 runs): A at pH 7.0 - v13 rule 89% within 0.05 (worst 0.47), v14 rule 97% (worst 0.13)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem, sim11 as sim

qc = np.load(os.path.join(HERE, 'eluate_C.npy'))
PRODS = {'A': sim.product('A'), 'B': sim.product('B'), 'C': (dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4]))}
TARGETS = {'A': (7.0, 7.2, 7.4, 8.0), 'B': (7.0, 7.2, 7.4, 8.0), 'C': (7.6, 8.2)}
CV, TIP_BIAS, ABS_SD, ABS_BIAS, READ_SD = 0.007, 0.005, 0.0003, 0.0005, 0.02


class Robot:
    """Pipetting errors as in sim11: per aspiration v(1 + bias_tip + N(0, cv)) + abs_bias_tip + N(0, abs_sd)."""
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


def base_titration(el, T, rng, rob):
    """The run's refit base curve, as the script builds it (Cycle 1 blind, Cycle 2 placed, validation, refit)."""
    def well(t, k):
        e, b = 3.0 / (1 + t), 3.0 - 3.0 / (1 + t)
        return rob.read(chem.ph_of(chem.mix([(rob.tip(k, e), el), (rob.tip(k, b), chem.TRIS)]), T))
    ph0 = rob.read(chem.ph_of(el, T))
    c1 = [M.BASE_MAX / M.BASE_SPAN * M.BASE_SPAN ** (k / 3) for k in range(4)]
    p1 = [well(t, k) for k, t in enumerate(c1)]
    lo, hi = M.base_window(ph0)
    tg2 = [lo + (hi - lo) * f for f in M.BASE_C2_SPLIT]
    cc = M.BaseCurve(c1, p1, ph0)
    c2 = [max((cc.ratio_for(t)[0] if cc.model else None) or M.locate_rising(np.array(c1), np.array(p1), ph0, t), 0.002) for t in tg2]
    for k in range(1, 4):
        c2[k] = max(c2[k], c2[k - 1] + max(0.0005, 0.02 * c2[k - 1]))
    p2 = [well(t, k) for k, t in enumerate(c2)]
    cur = M.BaseCurve(c1 + c2, p1 + p2, ph0)
    rv = [cur.ratio_for(t)[0] if t > ph0 + 0.02 else None for t in (7.2, 7.4, 7.6, 7.4)]
    pv = [well(x, k) if x else None for k, x in enumerate(rv)]
    return M.BaseCurve(c1 + c2 + [x for x in rv if x], p1 + p2 + [y for x, y in zip(rv, pv) if x], ph0), ph0


def trim_rule(rule, pl, m, curve):
    """Trim volumes (eluate, Tris) for a Tris pool. rule 'script' = the analysis file's trim_for; otherwise
    'kind:g1:g2': kind 'secant' (later trims: a straight line through the pool's last two readings - v13, with g1 0.6
    and g2 0.9), 'shape' (later trims: the curve's shift from the pool's last reading, added to its own ratio - the
    v14 rule, g1 0.75, g2 0.9), 'rel' (later trims as the first: in the curve's ratios, relative), 'offset' (the pool
    = the curve shifted by its mean offset); g1 = first-trim gain coming down a steep curve, g2 = later-trim gain."""
    if rule == 'script':
        da, de, db, how, *_ = M.trim_for(pl, m, None, curve)
        return de, db
    kind, g1, g2 = rule.split(':'); g1, g2 = float(g1), float(g2)
    T = pl['target']
    hist = pl['hist']
    ut = curve.ratio_for(T)[0]
    u_last = curve.ratio_for(m)[0]
    if ut is None or u_last is None:
        da, de, db, how, *_ = M.trim_for(pl, m, None, curve)
        return de, db
    if len(hist) < 2:
        sd = curve.sd_ph(T)
        g = g1 if m > T and sd is not None and sd > M.POOL_TOL else M.TRIM_GAIN
        return _amounts(pl, u_last, u_last + g * (ut - u_last))
    r_last = hist[-1][0]
    if kind == 'secant':
        (r1, m1), (r2, m2) = hist[-2], hist[-1]
        if abs(m2 - m1) >= 0.05 and (r2 - r1) / (m2 - m1) > 0:
            return _amounts(pl, r2, r2 + g2 * (r2 - r1) / (m2 - m1) * (T - m))
        return _amounts(pl, u_last, u_last + M.TRIM_GAIN * (ut - u_last))
    if kind == 'shape':
        return _amounts(pl, r_last, r_last + g2 * (ut - u_last))
    if kind == 'rel':
        return _amounts(pl, u_last, u_last + g2 * (ut - u_last))
    if kind == 'offset':
        u = [curve.ratio_for(y)[0] for _, y in hist]
        d = np.mean([ui - r for ui, (r, _) in zip(u, hist)])
        return _amounts(pl, r_last, r_last + g2 * ((ut - d) - r_last))
    raise ValueError(rule)


def _amounts(pl, r_now, r_aim):
    E = pl['E']
    room = max(M.LANE_MAX_ML - pl['A'] - pl['E'] - pl['B'], 0.0)
    if r_aim <= 0:
        return 0.0, 0.0
    de = db = 0.0
    if r_aim > r_now:
        db = (r_aim - r_now) * E
        if db < 0.6 * M.MIN_TRIM_ML:
            db = M.MIN_TRIM_ML
            de = (r_now * E + db) / r_aim - E
    else:
        de = (r_now / r_aim - 1) * E
    de, db = max(de, 0.0), max(db, 0.0)
    tot = de + db
    if 0.6 * M.MIN_TRIM_ML <= tot < M.MIN_TRIM_ML:
        de, db = de * M.MIN_TRIM_ML / tot, db * M.MIN_TRIM_ML / tot
        tot = M.MIN_TRIM_ML
    if tot > room:
        if db and de:
            return 0.0, 0.0
        de, db = de * room / tot, db * room / tot
    return de, db


def pool(el, T, curve, target, rule, rob):
    pl = dict(target=target, route='base', A=0.0, E=0.0, B=0.0, readings=[], hist=[])
    sd = curve.sd_ph(target)
    aim = min(max(float(os.environ.get("AIM_K", 2.0)) * sd - M.POOL_TOL, 0.0), float(os.environ.get("AIM_MAX", M.BASE_AIM_MAX))) if sd is not None else 0.0
    x = curve.ratio_for(target + aim)[0]
    if x is None:
        aim, x = 0.0, curve.ratio_for(target)[0]
    if x is None:
        return None
    E, B = M.POOL_ML / (1 + x), M.POOL_ML - M.POOL_ML / (1 + x)
    pl.update(E=E, B=B, aim=aim)
    true_e, true_b = rob.lane(E), rob.lane(B)
    for k in range(1, 4):
        m = rob.read(chem.ph_of(chem.mix([(true_e, el), (true_b, chem.TRIS)]), T))
        pl['readings'].append(m)
        pl['hist'].append((pl['B'] / pl['E'], m))
        if abs(m - target) <= M.POOL_TOL + 1e-9 or k == 3:
            break
        de, db = trim_rule(rule, pl, m, curve)
        if de + db < M.MIN_TRIM_ML:
            break
        pl['E'] += de; pl['B'] += db
        true_e += rob.lane(de); true_b += rob.lane(db)
    return (chem.ph_of(chem.mix([(true_e, el), (true_b, chem.TRIS)]), T) - target, len(pl['readings']) - 1,
            pl['readings'], round(aim, 3), round(true_b / true_e / (B / E), 4))


def trial(args):
    name, seed, rules = args
    el, T = PRODS[name]
    out = {}
    for r in rules:
        rng = np.random.default_rng(seed)                    # the same run for every rule
        rob = Robot(rng)
        curve, ph0 = base_titration(el, T, rng, rob)
        for t in TARGETS[name]:
            if t > ph0 + 0.05:
                o = pool(el, T, curve, t, r, rob)
                if o:
                    out[(r, t)] = o
    return out


if __name__ == '__main__':
    from multiprocessing import Pool
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    names = [a for a in sys.argv[2:] if a in PRODS] or list(PRODS)
    rules = [a for a in sys.argv[2:] if a not in PRODS] or ['script', 'secant:0.6:0.9']
    for name in names:
        with Pool(int(os.environ.get('NPROC', 2))) as p:
            outs = p.map(trial, [(name, k, rules) for k in range(N)])
        res = {r: {t: [o[(r, t)] for o in outs if (r, t) in o] for t in TARGETS[name]} for r in rules}
        for t in TARGETS[name]:
            line = []
            for r in rules:
                e = np.array([x[0] for x in res[r][t]]); n = np.array([x[1] for x in res[r][t]])
                line.append(f"{r}: RMS {np.sqrt(np.mean(e ** 2)):.3f} w0.03 {np.mean(np.abs(e) <= 0.03):.0%} w0.05 {np.mean(np.abs(e) <= 0.05):.0%} max {np.abs(e).max():.3f} trims {n.mean():.1f}")
            print(f"{name} pH {t}:\n   " + '\n   '.join(line), flush=True)
            if os.environ.get('DEBUG'):
                for r in rules:
                    for k, o in enumerate(outs):
                        if (r, t) in o and abs(o[(r, t)][0]) > float(os.environ['DEBUG']):
                            e_, n_, rd, aim, lf = o[(r, t)]
                            print(f"      {r} seed {k}: error {e_:+.3f}, readings {rd}, aim {aim}, Tris/eluate made vs planned {lf}")
