"""v16: trips to the troughs and robot minutes of the three titrations at 3.0, 2.94 and 2.7 mL wells, with and without
tips 5-8 (Edward, 22 Oct: would 2.7 mL save a trip?).    python tools/well_volume.py"""
import sys, os, shutil, collections, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sim16 as sim, scen16 as scen
cnt = collections.Counter()
orig = sim.Sim.pipette
def counted(self, kind, a):
    if kind == 'Aspirate' and (int(a[14]), int(a[15])) in ((2, 7), (2, 8)):
        cnt['trips'] += 1
    return orig(self, kind, a)
sim.Sim.pipette = counted
base = os.path.join(HERE, '..', 'python')
for vol in ('3.0', '2.94', '2.7'):
    d = tempfile.mkdtemp(prefix='wv_')
    for f in os.listdir(base):
        if f.endswith('.py'):
            t = open(os.path.join(base, f)).read()
            if f == 'EDP_tecan_gp_master.py':
                assert 'WELL_ML, PORTION_ML = 3.0, 0.98' in t
                t = t.replace('WELL_ML, PORTION_ML = 3.0, 0.98', f'WELL_ML, PORTION_ML = {vol}, 0.98')
            open(os.path.join(d, f), 'w').write(t)
    scen.SCRIPTS = [os.path.join(d, f) for f in sim.SCRIPTS]
    for eight in (1, 0):
        row = []
        for prod, seed in (('A', 1), ('B', 12), ('C', 4)):
            cnt.clear()
            s = scen.run(prod, seed, answers=dict(pool_lanes='0'), override={'Use_8_tips': eight})
            row.append(f"{prod}: {cnt['trips']} trips, {s.clock / 60:.1f} min" + (f" ERR {s.errors[:1]}" if s.errors else ''))
        print(f"wells {vol} mL, {'tips 1-8' if eight else 'tips 1-4'}: " + ' | '.join(row))
