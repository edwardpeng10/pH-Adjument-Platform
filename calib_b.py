"""Product of the 30 Sep run (eluate pH 5.45, strongly buffered): fit composition to its acid curve and rebuilt neutralisation curves."""
import numpy as np
from scipy.optimize import least_squares
import chem
acid = [(0, 5.45), (0.25, 4.19), (1.0, 3.64), (4.0, 3.14), (0.9645, 3.69), (1.5288, 3.50), (2.3195, 3.35), (3.5112, 3.21),
        (1.1511, 3.61), (2.0559, 3.38), (3.581, 3.17)]
p34 = [(0, 3.41), (0.143, 4.04), (0.286, 4.39), (0.429, 4.79), (0.571, 5.27), (0.714, 7.09), (0.857, 7.67), (1.0, 7.86),
       (0.7139, 7.12), (0.7456, 7.33), (0.7773, 7.46), (0.8815, 7.73)]
p32 = [(0, 3.15), (0.143, 3.94), (0.286, 4.28), (0.429, 4.52), (0.571, 4.99), (0.714, 5.78), (0.857, 7.34), (1.0, 7.73)]
def el(q): return dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3])
def res(q):
    e, T = el(q), q[4]
    r = [chem.acid_well(e, x, T) - y for x, y in acid]
    for Rp, data in ((2.0559, p34), (3.581, p32)):
        pool = chem.pool_liquid(e, Rp)
        r += [chem.base_well(pool, x, T) - y for x, y in data]
    return r
fit = least_squares(res, [0.15, 0.13, 0.02, 0.01, 20.0], bounds=([0, 0, 0, 0, 10], [2, 2, 1, 0.5, 30]))
q = fit.x
print('eluate B:', {k: round(float(v), 4) for k, v in el(q).items()}, 'T', round(q[4], 1), 'RMS', round(float(np.sqrt(np.mean(np.square(fit.fun)))), 3))
e = el(q)
print('acid fit:', ' '.join(f"{x}:{chem.acid_well(e, x, q[4]):.2f}/{y}" for x, y in acid[:8]))
pool = chem.pool_liquid(e, 2.0559)
print('P34 fit:', ' '.join(f"{x}:{chem.base_well(pool, x, q[4]):.2f}/{y}" for x, y in p34))
import os; np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eluate_B.npy'), q)
