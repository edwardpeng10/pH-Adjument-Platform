"""Replay of the real 12 Oct 2026 readings (v12 run, the 30 Sep product) through the analysis: what it stores (v15: the
v12 run had no product-only well and no B1-D1 wells - the eluate pH typed then goes in A1's place, B1-D1 as not read).
    python tools/replay_12oct.py <empty working folder>"""
import os, sys, json, subprocess
work = os.path.abspath(sys.argv[1]); os.makedirs(work, exist_ok=True)
HERE = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(HERE, '..', 'python')
env = dict(os.environ, TECAN_DATA_DIR=work, TECAN_ARCHIVE_DIR=os.path.join(work, 'archive'))
def w(name, text): open(os.path.join(work, name), 'w').write(text)
def run(*a):
    subprocess.run([sys.executable, os.path.join(PY, 'EDP_tecan_gp_master.py'), *a], env=env, check=True)
for f in os.listdir(PY):
    if f.endswith('.py'):
        open(os.path.join(work, f), 'w').write(open(os.path.join(PY, f)).read())
w('titration_params.csv', '1,4,3.6,3.4,3.2,3.5,0.2\n'); w('tris_params.csv', '7.2,7.4,7.6,7.4\n'); w('options.csv', '1,1\n')
w('batch_info.csv', '12Oct,30Sep product,ProA eluate\n')
run('Acid', '0'); run('Base', '0')
w('acid_c1_measured.csv', '5.48,0,0,0,4.68,4.15,3.64,3.15\n'); run('Acid', '1')
w('base_c1_measured.csv', '0,0,0,0\n'); run('Base', '1')                       # no base titration in the v12 run
w('cycle2_measured.csv', '3.80,3.68,3.57,3.45,3.33,3.22,3.12,3.01\n'); run('Acid', '2')
w('validation_results.csv', '3.61,3.41,3.19,3.49\n'); run('Acid', '3')
w('tris_c1_measured.csv', '6.90,7.42,7.73,8.00,6.94,7.39,7.73,8.01\n'); run('Tris', '1')
w('tris_val_measured.csv', '7.26,7.40,7.58,7.43\n'); run('Tris', '2')
run('Pools', 'range')
st = json.load(open(os.path.join(work, 'run_state.json')))
print('acid Cycle 2 ratios (this plan vs the run):', [round(x, 4) for x in st['acid']['c2']], 'vs 0.7451 ... 5.7116')
print('VI Cycle 1 (this plan vs the run):', [(round(c['rp'], 4), [round(t, 4) for t in c['t']]) for c in st['tris']['c1']])
print(open(os.path.join(work, 'python_log.txt')).read().split('pools range')[-1][:3000])
for a, b in ((5.48, 3.4), (3.4, 7.4), (3.5, 7.4), (3.6, 7.4), (3.2, 7.4), (4.0, 7.4), (4.5, 7.2)):
    run('Ratio', str(a), str(b))
    print(f'Ratio {a} -> {b}:', open(os.path.join(work, 'ratio_lookup.csv')).read().split())
