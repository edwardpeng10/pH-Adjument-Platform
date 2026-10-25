"""Product of the 21 Oct 2026 v15 run (eluate pH 6.46): fit the chemistry model to every reading of the three titrations,
with and without a free titrant-strength factor (Tris stock relative to the acid)."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from scipy.optimize import least_squares
import chem
acid = [(0, 6.46), (0.01, 5.74), (0.02, 5.36), (0.04, 4.99), (0.08, 4.64), (0.2947, 4.05), (1.0858, 3.53), (4.0, 3.05),
        (0.5671, 3.78), (0.7591, 3.64), (1.0160, 3.52), (1.3600, 3.42), (1.8203, 3.30), (2.4364, 3.20), (3.2611, 3.10), (4.3649, 3.00),
        (0.8606, 3.60), (1.4226, 3.40), (2.4504, 3.19), (1.1018, 3.50)]
base = [(0.005, 7.09), (0.0196, 7.94), (0.0766, 8.60), (0.3, 9.19), (0.0034, 6.82), (0.0068, 7.28), (0.0135, 7.75), (0.0364, 8.27),
        (0.0059, 7.25), (0.0079, 7.42), (0.0108, 7.60), (0.0079, 7.43)]
vi = [(2.4504, 0.7528, 7.43), (2.4504, 0.8238, 7.70), (2.4504, 0.9374, 7.90), (2.4504, 1.1363, 8.13),
      (0.8606, 0.4903, 7.40), (0.8606, 0.5365, 7.67), (0.8606, 0.6105, 7.89), (0.8606, 0.7400, 8.11),
      (1.1020, 0.5276, 7.26), (1.1020, 0.5521, 7.49), (1.1020, 0.5911, 7.63), (1.1020, 0.5521, 7.43)]
def el(q): return dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3], Tr=q[4])
def ph(e, A, B, T, ct):
    return chem.ph_of(chem.mix([(1.0, e), (A, chem.ACID), (B, dict(Tr=ct))]), T)
def residuals(q, free_ct, parts=False):
    e, T = el(q), q[5]
    ct = q[6] if free_ct else 1.0
    ra = [ph(e, x, 0, T, ct) - y for x, y in acid]
    rb = [ph(e, 0, x, T, ct) - y for x, y in base]
    rv = [ph(e, a, t * (1 + a), T, ct) - y for a, t, y in vi]
    if parts: return np.array(ra), np.array(rb), np.array(rv)
    # weights: steep base wells (pipetting of 15-60 uL) count less
    return np.r_[ra, np.array(rb) * 0.5, np.array(rv) * 0.7]
for free_ct in (False, True):
    best = None
    for x0 in ([0.07, 0.02, 0.0, 0.01, 0.04, 20, 1.0], [0.1, 0.05, 0.02, 0.005, 0.02, 22, 1.05], [0.05, 0.0, 0.05, 0.02, 0.06, 18, 1.1]):
        lo = [0, 0, 0, 0, 0, 10, 0.8]; hi = [2, 2, 1, 0.5, 0.5, 30, 1.3]
        x0 = x0 if free_ct else x0[:6]; lo = lo if free_ct else lo[:6]; hi = hi if free_ct else hi[:6]
        f = least_squares(residuals, x0, bounds=(lo, hi), args=(free_ct,))
        if best is None or f.cost < best.cost: best = f
    q = best.x
    ra, rb, rv = residuals(q, free_ct, parts=True)
    rms = lambda r: float(np.sqrt(np.mean(r ** 2)))
    print(('free Tris strength' if free_ct else 'titrants 1 M exactly') + ':', {k: round(float(v), 4) for k, v in el(q).items()},
          'T', round(q[5], 1), ('Tris x %.3f' % q[6]) if free_ct else '',
          f'| RMS acid {rms(ra):.3f} base {rms(rb):.3f} VI {rms(rv):.3f}')
    print('   VI residuals:', ' '.join(f'{x:+.2f}' for x in rv))
    print('   base residuals:', ' '.join(f'{x:+.2f}' for x in rb))
    if free_ct:
        np.save(os.path.join(HERE, 'eluate_D.npy'), q)          # Ac, Na, Pc, Ph, Tr, T, Tris strength
