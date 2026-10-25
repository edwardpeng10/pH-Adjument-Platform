"""Replay of the real 21 Oct 2026 readings (v15 run, eluate pH 6.46) through the analysis, and what the stored ratios
give on the chemistry model fitted to all 44 readings with a free Tris strength (tools/eluate_D.npy: Tris x 1.099).
    python tools/replay_21oct.py python <empty working folder>"""
import os, sys, json, subprocess, numpy as np, pandas as pd
PY, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]); os.makedirs(work, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import chem
q = np.load(os.path.join(HERE, 'eluate_D.npy')); EL = dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3], Tr=q[4]); T, CT = q[5], q[6]
env = dict(os.environ, TECAN_DATA_DIR=work, TECAN_ARCHIVE_DIR=os.path.join(work, 'archive'))
def w(name, text): open(os.path.join(work, name), 'w').write(text)
def run(*a):
    subprocess.run([sys.executable, os.path.join(work, 'EDP_tecan_gp_master.py'), *a], env=env, check=True)
for f in os.listdir(PY):
    if f.endswith('.py'):
        open(os.path.join(work, f), 'w').write(open(os.path.join(PY, f)).read())
w('titration_params.csv', '1,4,3.6,3.4,3.2,3.5,0.2\n'); w('tris_params.csv', '7.2,7.4,7.6,7.4\n'); w('options.csv', '1,1\n')
w('batch_info.csv', '21Oct,21 Oct product,ProA eluate\n')
run('Acid', '0'); run('Base', '0')
w('acid_c1_measured.csv', '6.46,5.74,5.36,4.99,4.64,4.05,3.53,3.05\n'); run('Acid', '1')
w('base_c1_measured.csv', '7.09,7.94,8.60,9.19\n'); run('Base', '1')
w('cycle2_measured.csv', '3.78,3.64,3.52,3.42,3.30,3.20,3.10,3.00\n'); run('Acid', '2')
w('base_c2_measured.csv', '6.82,7.28,7.75,8.27\n'); run('Base', '2')
w('validation_results.csv', '3.60,3.40,3.19,3.50,3.79,3.65,3.53,3.43\n'); run('Acid', '3')
w('tris_c1_measured.csv', '7.43,7.70,7.90,8.13,7.40,7.67,7.89,8.11\n'); run('Tris', '1')
w('base_val_measured.csv', '7.25,7.42,7.60,7.43\n'); run('Base', '3')
w('tris_val_measured.csv', '7.26,7.49,7.63,7.43\n'); run('Tris', '2')
run('Pools', 'range')
st = json.load(open(os.path.join(work, 'run_state.json')))
log = open(os.path.join(work, 'python_log.txt')).read()
for line in log.splitlines():
    if any(k in line for k in ('Chemistry model', 'Tris refit', 'Tris titration:', 'Base refit', 'Ratio map', 'Tris pH', 'Base pH')):
        print('  ', line.split('] ')[-1][:200])
def ph(A, B):
    return chem.ph_of(chem.mix([(1.0, EL), (A, chem.ACID), (B, dict(Tr=CT))]), T)
df = pd.read_csv(os.path.join(work, 'ratio_map_run.csv'))
print('Stored ratios on the model fitted to all 44 readings (Tris x %.3f, %.1f C): true pH - target' % (CT, T))
for route in ('acid', 'base', 'neutralise'):
    d = df[df.Route == route]
    errs = []
    for _, r in d.iterrows():
        A = r.Acid_per_mL_eluate if route == 'neutralise' else (r.Per_mL_eluate if route == 'acid' else 0.0)
        B = r.Per_mL_eluate if route != 'acid' else 0.0
        errs.append((ph(A, B) - r.Target_pH, r.Basis, r.Start_pH, r.Target_pH))
    e = np.array([x[0] for x in errs])
    print(f'  {route:10s} {len(e)} entries: median |err| {np.median(np.abs(e)):.3f}, 95% {np.percentile(np.abs(e), 95):.3f}, max {np.abs(e).max():.3f}')
    by = {}
    for x, b, s0, t0 in errs:
        by.setdefault(str(b).split(' (')[0], []).append(x)
    for b, xs in by.items():
        xs = np.array(xs)
        print(f'      {b:34s} n={len(xs):3d} mean {xs.mean():+.3f} RMS {np.sqrt(np.mean(xs**2)):.3f} max {np.abs(xs).max():.3f}')
    if route == 'neutralise':
        sel = d[(d.Target_pH.round(2).isin([7.2, 7.4, 7.6])) & (d.Start_pH.round(2).isin([3.2, 3.4, 3.5, 3.6]))]
        for _, r in sel.iterrows():
            print(f'      from {r.Start_pH:.2f} to {r.Target_pH:.2f}: {r.Per_mL_start:.4f} per mL acidified -> true {ph(r.Acid_per_mL_eluate, r.Per_mL_eluate):.3f} [{r.Basis}]')
    if route == 'base':
        sel = d[d.Target_pH.round(2).isin([7.0, 7.2, 7.4, 7.6, 7.8, 8.0])]
        for _, r in sel.iterrows():
            print(f'      to {r.Target_pH:.2f}: {r.Per_mL_eluate:.5f} per mL eluate -> true {ph(0, r.Per_mL_eluate):.3f} [{r.Basis}]')
