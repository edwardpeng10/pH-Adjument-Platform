"""v15 pool Monte Carlo (python mc15_pools.py 30 A B C C13; NPROC=cores): the full v15 script in the simulator, many seeds;
final true pH of each pool, trims, volume. C13 = the 13 Oct request on a 2 Oct-type eluate (5.29 acid, 7.6 Tris)."""
import sys, os, json, numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scen15 as S, chem, sim11 as sim
CASES = {'A': ('A', '2,4,6', '7.0,7.4,8.0'), 'B': ('B', '2,4,6', '7.0,7.4,8.0'), 'C': ('C', '1,3,5', '7.2,7.6,8.2'), 'C13': ('C', '1,2', '5.29,7.6')}
def one(args):
    case, seed = args
    prod, lanes, phs = CASES[case]
    s = S.run(prod, 100 + seed, answers={'pool_lanes': lanes, 'pool_phs': phs})
    st = json.load(open(os.path.join(s.data, 'run_state.json')))
    out = []
    for pl in st.get('pools', []):
        lane = s.cont(*sim.LPH, pl['lane'])
        true = chem.ph_of(lane.comp, s.sc['T']) if lane.vol > 1 else None
        out.append((pl['target'], pl['route'], true, pl['trims'], lane.vol, pl.get('note', ''), len(s.errors)))
    return case, seed, out
if __name__ == '__main__':
    N = int(sys.argv[1])
    jobs = [(c, k) for c in sys.argv[2:] for k in range(N)]
    with Pool(int(os.environ.get('NPROC', 8))) as p:
        res = p.map(one, jobs)
    for case in sys.argv[2:]:
        rows = [r for c, k, out in res if c == case for r in out]
        errs = [r[6] for c, k, out in res if c == case for r in out[:1]]
        print(f"== {case} ({N} runs, simulator errors {sum(errs)})")
        for tgt in sorted({r[0] for r in rows}):
            rr = [r for r in rows if r[0] == tgt]
            e = np.array([r[2] - tgt for r in rr if r[2] is not None])
            tr = np.array([r[3] for r in rr])
            stuck = sum(1 for r in rr if 'less than 8 tips' in r[5])
            print(f"   pH {tgt} ({rr[0][1]}): final error RMS {np.sqrt(np.mean(e**2)):.3f}, within 0.03: {np.mean(np.abs(e) <= 0.03):.0%}, within 0.05: {np.mean(np.abs(e) <= 0.05):.0%}, "
                  f"max {np.abs(e).max():.3f}; trims mean {tr.mean():.1f} (0/1/2: {np.mean(tr == 0):.0%}/{np.mean(tr == 1):.0%}/{np.mean(tr == 2):.0%}); low and stuck {stuck}; volume {np.mean([r[4] for r in rr]):.1f} mL")
