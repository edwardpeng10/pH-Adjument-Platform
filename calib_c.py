"""Product of the 2 Oct 2026 v7 run (eluate pH 7.43): fit the chemistry model to the acid curve and the three
rebuilt neutralisation Cycle-1 curves (wells made from eluate + acid + Tris from the troughs)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from scipy.optimize import least_squares
import chem
acid = [(0, 7.43), (0.0625, 5.03), (0.125, 4.67), (0.25, 4.35), (0.5, 4.04), (1, 3.74), (2, 3.46), (4, 3.20),
        (1.2419, 3.67), (1.8883, 3.50), (2.9516, 3.31), (4.2592, 3.20),
        (1.5146, 3.60), (2.2966, 3.43), (4.1296, 3.20), (1.8883, 3.51)]
rebuilt = {4.1296: ([0.0000,0.1638,0.3276,0.4914,0.6552,0.8190,0.9828,1.1466], [3.22,4.04,4.40,4.71,5.20,6.80,7.68,7.92]),
           2.2966: ([0.0000,0.1437,0.2873,0.4310,0.5747,0.7183,0.8620,1.0057], [3.44,4.09,4.44,4.76,5.22,6.89,7.67,7.95]),
           1.5146: ([0.0000,0.1261,0.2523,0.3784,0.5046,0.6307,0.7569,0.8830], [3.60,4.11,4.49,4.81,5.26,7.08,7.67,7.92])}
def el(q): return dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3], Tr=q[4])
def res(q):
    e, T = el(q), q[5]
    r = [chem.acid_well(e, x, T) - y for x, y in acid]
    for Rp, (ts, phs) in rebuilt.items():
        pool = chem.pool_liquid(e, Rp)
        # the near-equivalence wells (pH 5.2-7.1) are steep: 1% pipetting moves them ~0.1+, so weight them less
        r += [(chem.base_well(pool, x, T) - y) / (0.1 if 5.0 < y < 7.2 else 0.02) * 0.02 for x, y in zip(ts, phs)]
    return r
best = None
for x0 in ([0.05, 0.1, 0.05, 0.02, 0.02, 20], [0.02, 0.05, 0.1, 0.05, 0.0, 22], [0.1, 0.15, 0.1, 0.01, 0.05, 18]):
    f = least_squares(res, x0, bounds=([0, 0, 0, 0, 0, 10], [2, 2, 1, 0.5, 0.5, 30]))
    if best is None or f.cost < best.cost: best = f
q = best.x; e = el(q)
acid_rms = np.sqrt(np.mean([(chem.acid_well(e, x, q[5]) - y) ** 2 for x, y in acid]))
print('eluate C:', {k: round(float(v), 4) for k, v in e.items()}, 'T', round(q[5], 1), '| acid RMS', round(float(acid_rms), 3))
print('acid fit :', ' '.join(f"{x}:{chem.acid_well(e, x, q[5]):.2f}/{y}" for x, y in acid[:12]))
for Rp, (ts, phs) in rebuilt.items():
    pool = chem.pool_liquid(e, Rp)
    print(f'Rp {Rp}:', ' '.join(f"{x:.3f}:{chem.base_well(pool, x, q[5]):.2f}/{y}" for x, y in zip(ts, phs)))
np.save(os.path.join(HERE, 'eluate_C.npy'), np.r_[q[:4], q[5], q[4]])   # Ac, Na, Pc, Ph, T, Tr
