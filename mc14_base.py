"""v14 (unchanged since v13): Monte Carlo of the base titration alone (C1 blind -> C2 placed -> fit -> validation -> refit) on the modelled eluates."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem, sim11 as sim
from scipy.optimize import brentq
qc = np.load(os.path.join(HERE, 'eluate_C.npy'))
prods = {'A': sim.product('A'), 'B': sim.product('B'), 'C': (dict(Ac=qc[0], Na=qc[1], Pc=qc[2], Ph=qc[3], Tr=qc[5]), float(qc[4]))}
N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
for name, (el, T) in prods.items():
    def well_ph(t, rng, typo=0.0):
        e = 3.0 / (1 + t); b = 3.0 - e
        b_act = b * (1 + rng.normal(0, 0.007) + rng.normal(0, 0.005)) + rng.normal(0, 0.0003) + rng.normal(0, 0.0005)
        return chem.ph_of(chem.mix([(e * (1 + rng.normal(0, 0.007)), el), (max(b_act, 0), chem.TRIS)]), T) + rng.normal(0, 0.02) + typo
    true_ph = lambda x: chem.ph_of(chem.mix([(1, el), (x, chem.TRIS)]), T)
    e_val, e_fit, missed, wrong = [], [], 0, 0
    for seed in range(N):
        rng = np.random.default_rng(seed)
        ph0 = chem.ph_of(el, T) + rng.normal(0, 0.02)
        c1 = [M.BASE_MAX / M.BASE_SPAN * M.BASE_SPAN ** (k / 3) for k in range(4)]
        p1 = [well_ph(t, rng) for t in c1]
        lo, hi = M.base_window(ph0)
        tg2 = [lo + (hi - lo) * f for f in M.BASE_C2_SPLIT]
        cc = M.BaseCurve(c1, p1, ph0)
        c2 = [max((cc.ratio_for(t)[0] if cc.model else None) or M.locate_rising(np.array(c1), np.array(p1), ph0, t), 0.002) for t in tg2]
        for k in range(1, 4): c2[k] = max(c2[k], c2[k - 1] + max(0.0005, 0.02 * c2[k - 1]))
        typo_k = 1 if seed % 4 == 3 else None
        p2 = [well_ph(t, rng, 0.25 if k == typo_k else 0) for k, t in enumerate(c2)]
        cur = M.BaseCurve(c1 + c2, p1 + p2, ph0)
        if typo_k is not None:
            missed += not any(abs(x - c2[typo_k]) < 1e-9 for x, y, r in cur.dropped)
        wrong += sum(1 for x, y, r in cur.dropped if typo_k is None or abs(x - c2[typo_k]) > 1e-9)
        rv = []
        for t in (7.2, 7.4, 7.6, 7.4):
            x = cur.ratio_for(t)[0] if t > ph0 + 0.02 else None
            rv.append(x)
            if x: e_val.append(true_ph(x) - t)
        pv = [well_ph(x, rng) if x else 0 for x in rv]
        cur2 = M.BaseCurve(c1 + c2 + [x for x in rv if x], p1 + p2 + [y for x, y in zip(rv, pv) if x], ph0)
        for t in (7.0, 7.2, 7.4, 7.6, 8.0):
            if t > ph0 + 0.02:
                x = cur2.ratio_for(t)[0]
                e_fit.append(true_ph(x) - t)
    a, b = np.array(e_val), np.array(e_fit)
    print(f"{name}: validation wells (C1+C2 fit) RMS {np.sqrt(np.mean(a**2)):.3f} max {np.abs(a).max():.3f} | after refit, pH 7.0-8.0: RMS {np.sqrt(np.mean(b**2)):.3f}, "
          f"95% within {np.percentile(np.abs(b), 95):.3f}, max {np.abs(b).max():.3f} | typos missed {missed}/{N // 4}, good readings left out {wrong}")
