"""v14: checks the analysis file's vectorised chemistry model (EDP_tecan_gp_master.chem_ph) against the simulator's
chem.py on random eluate + acid + Tris mixtures, and times it.
    python tools/check_chem.py"""
import sys, os, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'python')); sys.path.insert(0, HERE)
import EDP_tecan_gp_master as M, chem

rng = np.random.default_rng(1)
worst = 0.0
for trial in range(300):
    el = dict(Ac=rng.uniform(0, 0.2), Tr=rng.uniform(0, 0.15) * (rng.random() < 0.3), Na=rng.uniform(0, 0.2),
              Cl=rng.uniform(0, 0.05), Pc=rng.uniform(0, 0.1), Ph=rng.uniform(0, 0.02))
    T = rng.uniform(15, 28)
    A, B = rng.uniform(0, 4) * (rng.random() < 0.7), rng.uniform(0, 2) * (rng.random() < 0.6)
    ref = chem.ph_of(chem.mix([(1.0, el), (A, chem.ACID), (B, chem.TRIS)]), T)
    got = float(M.chem_ph(M.chem_mix(el, A, B), T))
    worst = max(worst, abs(got - ref))
print(f"300 random mixtures: worst difference from chem.py {worst:.1e} pH")
el = dict(Ac=0.1, Na=0.1, Pc=0.02, Ph=0.005)
A = np.r_[np.linspace(0, 5, 16), np.zeros(12), np.linspace(0.5, 3.5, 12)]
B = np.r_[np.zeros(16), np.linspace(0.005, 0.3, 12), np.linspace(0.5, 2.0, 12)]
t0 = time.time()
for _ in range(100):
    M.chem_ph(M.chem_mix(el, A, B), 22.0)
print(f"40-mixture batch: {1000 * (time.time() - t0) / 100:.1f} ms per call")
