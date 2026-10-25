"""Base titration alone (Tris straight into the eluate): C1 blind -> C2 placed -> fit -> validation -> refit, with the
simulator's pipetting (per tip: 0.7% CV, 0.5% bias, 0.3 uL random, 0.5 uL fixed) and probe noise (0.02, 2 decimals).
Designs differ in the Cycle-2 wells. Metric: true pH that the stored (refit) ratios give at 7.0-8.0, minus the target.
    python tools/base_bench.py 80          (TRIS_M=0.1: a 0.1 M Tris for the base wells)"""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem, sim11 as sim
qc = np.load(os.path.join(HERE, 'eluate_C.npy')); qd = np.load(os.path.join(HERE, 'eluate_D.npy'))
PRODS = {'A': sim.product('A') + (1.0,), 'B': sim.product('B') + (1.0,),
         'C': (dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4]), 1.0),
         'D': (dict(Ac=qd[0], Na=qd[1], Pc=qd[2], Ph=qd[3], Tr=qd[4]), float(qd[5]), float(qd[6]))}
WELL = 3.0
def designs():
    return {
        'v15: C2 4 wells (5/30/55/85% of 6.8-8.5)': dict(c2=(0.05, 0.30, 0.55, 0.85)),
        'C2 8 wells, same spread': dict(c2=tuple(np.linspace(0.05, 0.85, 8))),
        'C2 8 wells, dense at the targets': dict(c2=(0.03, 0.15, 0.24, 0.33, 0.42, 0.53, 0.68, 0.88)),
        'C2 4 + validation 8 (targets twice)': dict(c2=(0.05, 0.30, 0.55, 0.85), val2=True),
    }
class Lab:
    def __init__(self, el, T, ct, seed, tris_m=1.0):
        self.el, self.T, self.ct, self.tris_m = el, T, ct, tris_m
        self.rng = np.random.default_rng(seed)
        self.bias = self.rng.normal(0, 0.005, 4); self.abs_bias = self.rng.normal(0, 0.0005, 4)
    def disp(self, tip, ml):
        if ml <= 0: return 0.0
        n = int(np.ceil(ml / 0.98 - 1e-9)); out = 0.0
        for _ in range(n):
            v = ml / n
            out += max(v * (1 + self.bias[tip] + self.rng.normal(0, 0.007)) + self.abs_bias[tip] + self.rng.normal(0, 0.0003), 0)
        return out
    def well(self, t, tip):
        """t: mL 1 M-equivalent Tris per mL eluate; with tris_m < 1 a dilute Tris is dispensed (t / tris_m mL)."""
        vt = t / self.tris_m
        e = WELL / (1 + vt); b = WELL - e
        return self.disp(tip, e), self.disp(tip, b)
    def ph(self, e, b):
        return chem.ph_of(chem.mix([(e, self.el), (b, dict(Tr=self.ct * self.tris_m))]), self.T)
    def read(self, ts):
        return [round(self.ph(*self.well(t, k % 4)) + self.rng.normal(0, 0.02), 2) for k, t in enumerate(ts)]
def run(prod, seed, d, tris_m=1.0):
    el, T, ct = PRODS[prod]
    lab = Lab(el, T, ct, seed, tris_m)
    true_ph = lambda x: chem.ph_of(chem.mix([(1, el), (x, dict(Tr=ct))]), T)
    ph0 = round(chem.ph_of(el, T) + lab.rng.normal(0, 0.02), 2)
    c1 = [M.BASE_MAX / M.BASE_SPAN * M.BASE_SPAN ** (k / 3) for k in range(4)]
    p1 = lab.read(c1)
    lo, hi = M.base_window(ph0)
    tg2 = [lo + (hi - lo) * f for f in d['c2']]
    cc = M.BaseCurve(c1, p1, ph0)
    c2 = [max((cc.ratio_for(t)[0] if cc.model else None) or M.locate_rising(np.array(c1), np.array(p1), ph0, t), 0.002) for t in tg2]
    for k in range(1, len(c2)): c2[k] = max(c2[k], c2[k - 1] + max(0.0005, 0.02 * c2[k - 1]))
    p2 = lab.read(c2)
    cur = M.BaseCurve(c1 + c2, p1 + p2, ph0)
    targets = (7.2, 7.4, 7.6, 7.4) + ((7.2, 7.4, 7.6, 7.3) if d.get('val2') else ())
    rv = [cur.ratio_for(t)[0] if t > ph0 + 0.02 else None for t in targets]
    pv = lab.read([x for x in rv if x])
    val_err = [true_ph(x) - t for x, t in zip(rv, targets) if x]
    cur2 = M.BaseCurve(c1 + c2 + [x for x in rv if x], p1 + p2 + pv, ph0)
    stored = []
    for t in (7.0, 7.2, 7.4, 7.6, 7.8, 8.0):
        if t > ph0 + 0.05:
            x = cur2.ratio_for(t)[0]
            if x: stored.append(true_ph(x) - t)
    return val_err, stored
if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    tris_m = float(os.environ.get('TRIS_M', 1.0))
    import contextlib, io
    for prod in ('A', 'B', 'C', 'D'):
        print(f'== {prod} (eluate pH {chem.ph_of(PRODS[prod][0], PRODS[prod][1]):.2f})' + (f', Tris {tris_m} M in the base wells' if tris_m != 1 else ''))
        for name, d in designs().items():
            V, S = [], []
            for s in range(N):
                with contextlib.redirect_stdout(io.StringIO()):
                    v, st = run(prod, 500 + s, d, tris_m)
                V += v; S += st
            V, S = np.abs(V), np.abs(S)
            print(f'   {name:42s} validation wells RMS {np.sqrt(np.mean(V**2)):.3f} | stored 7.0-8.0: RMS {np.sqrt(np.mean(S**2)):.3f}, 95% {np.percentile(S, 95):.3f}, max {S.max():.3f}')
