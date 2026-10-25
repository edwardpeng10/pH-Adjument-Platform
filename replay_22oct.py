"""Replay of the real 22 Oct 2026 readings (v16 run, AZD0292 DurA eluate, pH 4.55; all base wells with 1 M Tris) through
the analysis, then every base and VI well against the chemistry model the run fits to all 44 readings: how many uL of
Tris (in the 3 mL well) that model needs to land on the reading - the diagnosis behind v17's 0.1 M Tris.
    python tools/replay_22oct.py <python folder> <empty working folder>
BASE_VAL="7.28,7.39,7.60,7.41" replays with another base validation entry (e.g. C6 as its recipe predicts)."""
import os, sys, json, subprocess
import numpy as np
PY, work = os.path.abspath(sys.argv[1]), os.path.abspath(sys.argv[2]); os.makedirs(work, exist_ok=True)
env = dict(os.environ, TECAN_DATA_DIR=work, TECAN_ARCHIVE_DIR=os.path.join(work, 'archive'), EDP_EXPORT_DIR=os.path.join(work, 'export'))
def w(name, text): open(os.path.join(work, name), 'w').write(text)
def run(*a):
    subprocess.run([sys.executable, os.path.join(work, 'EDP_tecan_gp_master.py'), *a], env=env, check=True)
for f in os.listdir(PY):
    if f.endswith('.py'):
        open(os.path.join(work, f), 'w').write(open(os.path.join(PY, f)).read())
w('titration_params.csv', '1,4,3.6,3.4,3.2,3.5,0.2\n'); w('tris_params.csv', '7.2,7.4,7.6,7.4\n')
w('options.csv', '1,1,1\n')                         # as on 22 Oct: no dilute Tris (v16)
w('batch_info.csv', 'GSUL16,AZD0292,DurA Eluate\n')
run('Acid', '0'); run('Base', '0')
w('acid_c1_measured.csv', '4.55,4.28,4.14,3.94,3.71,3.25,2.87,2.60\n'); run('Acid', '1')
w('base_c1_measured.csv', '4.89,7.53,8.69,9.35\n'); run('Base', '1')
w('cycle2_measured.csv', '3.81,3.69,3.57,3.46,3.34,3.24,3.13,3.04\n'); run('Acid', '2')
w('base_c2_measured.csv', '6.87,7.35,7.80,8.26\n'); run('Base', '2')
w('validation_results.csv', '3.61,3.39,3.19,3.49,3.80,3.67,3.57,3.44\n'); run('Acid', '3')
w('tris_c1_measured.csv', '6.97,7.33,7.72,8.02,5.86,6.81,7.49,7.89\n'); run('Tris', '1')
w('base_val_measured.csv', os.environ.get('BASE_VAL', '7.28,7.39,7.39,7.41') + '\n'); run('Base', '3')
w('tris_val_measured.csv', '7.21,7.42,7.62,7.41\n'); run('Tris', '2')
run('Pools', 'range')
log = open(os.path.join(work, 'python_log.txt')).read()
for line in log.splitlines():
    if any(k in line for k in ('Base', 'Chemistry model', 'Tris titration', 'Tris refit', 'left out')):
        print('  ', line[:220])

# every base and VI well against the run's own chemistry model
sys.path.insert(0, work)
os.environ['TECAN_DATA_DIR'] = work
import EDP_tecan_gp_master as M
from scipy.optimize import brentq
st = M.load_state()
q = st['chem']['q']; el = M.chem_eluate(q); T = q[5]
def ph(a, b): return float(M.chem_ph(M.chem_mix(el, np.array([a]), np.array([b])), T)[0])
print(f"\nEach well against the chemistry model fitted to all {st['chem']['n']} readings (uL of Tris in the 3 mL well that "
      "would put the model on the reading):")
for a, b, y, sd, kind, *_ in M.run_points(st):
    if kind not in ('base', 'neutralisation'):
        continue
    try:
        d = brentq(lambda d: ph(a, b + d) - y, -b * 0.6, b * 0.6)
    except ValueError:
        d = float('nan')
    E = 3.0 / (1 + a + b)
    print(f"  {kind:15s} acid {a:6.4f} Tris {b:7.4f}: read {y:5.2f}, model {ph(a, b):5.2f}, Tris dispensed {b * E * 1000:6.1f} uL, "
          f"reading asks for {d * E * 1000:+6.2f} uL")
