"""Check a generated .esc stays inside what scripts that have run on the robot already used: per command and argument
position, no longer text, no new non-alphanumeric characters, no new fixed values, no more complex formulas."""
import os, re, sys, collections, sim11 as sim7
REF = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'reference')
OK = [os.path.join(REF, n) for n in ('EDP_v16_22Oct_2026.esc', 'EDP_v15_21Oct_2026.esc', 'EDP_v14_13Oct_2026.esc', 'EDP_v12_12Oct_2026.esc', 'EDP_v7_2Oct_2026.esc', 'EDP_v6_30Sep_2026.esc', 'EDP_28Sep_2026.esc', 'GCB_Titration_phyton_ML_hybrid_PD.esc', 'GCB_Titration_phyton_ML_hybrid.esc')]
FREE = {('Variable', 0), ('Variable', 1), ('Variable', 3), ('Variable', 5), ('Variable', 6), ('Comment', 0), ('Group', 0), ('UserPrompt', 0),
        ('Execute', 0), ('ImportVariable', 0), ('ImportVariable', 1), ('ImportVariable', 3), ('ImportVariable', 5), ('ExportVariable', 0),
        ('ExportVariable', 1), ('If', 0), ('If', 2), ('If', 3), ('IfThen', 0), ('IfThen', 2), ('BeginLoop', 0), ('BeginLoop', 1)}
TIPS58_RUN = {(28, 1), (28, 2), (22, 1), (2, 7), (2, 8)}   # (grid, site) where tips 5-8 have pipetted on this robot

def load(f):
    L = open(f, 'rb').read().decode('latin-1').split('\r\n'); b = L.index('--{ RPG }--') + 1
    return [(m.group(1), sim7.split_args(m.group(2)) if m.group(2) else []) for m in (re.match(r'(\w+)\((.*)\);$', l) for l in L[b:]) if m]
def decode_wells(x):
    """EVOware well selection: 2 hex digits for columns, 2 for rows, then 7 bits per character (offset 48), wells
    numbered down each column. Returns (nx, ny, selected well indices) or None."""
    x = x.strip('"')
    try:
        nx, ny = int(x[0:2], 16), int(x[2:4], 16)
    except ValueError:
        return None
    sel = []
    for k, ch in enumerate(x[4:]):
        v = ord(ch) - 48
        if not 0 <= v < 128:
            return None
        sel += [7 * k + b for b in range(7) if v >> b & 1]
    if any(i >= nx * ny for i in sel):
        return None
    return nx, ny, sel

def well_selection_ok(x, geometries):
    """v11 rule: a well selection that has never run is accepted when it decodes on a plate geometry that HAS run and
    selects whole rows of one column (for example lane 2 of the 6-well plate, rows A-H, where lanes 1, 3, 5 have run).
    It is the same command form with a different selection, not a new parameter."""
    d = decode_wells(x)
    if not d or (d[0], d[1]) not in geometries or not d[2]:
        return False
    nx, ny, sel = d
    cols = {i // ny for i in sel}
    rows = sorted(i % ny for i in sel)
    return len(cols) == 1 and rows == list(range(rows[0], rows[0] + len(rows)))

def check(path):
    mx, chars, vals = collections.defaultdict(int), collections.defaultdict(set), collections.defaultdict(set)
    for f in OK:
        for n, a in load(f):
            for i, x in enumerate(a):
                mx[n, i] = max(mx[n, i], len(x)); vals[n, i].add(x)
                chars[n].update(c for c in x if not c.isalnum())
    geometries = {d[:2] for (n, i), xs in vals.items() if n in ('Aspirate', 'Dispense', 'Mix') and i == 17
                  for d in [decode_wells(x) for x in xs] if d}
    bad = []
    for j, (n, a) in enumerate(load(path)):
        for i, x in enumerate(a):
            if n in ('Aspirate', 'Dispense', 'Mix') and i == 17 and x not in vals[n, i] and well_selection_ok(x, geometries):
                continue                              # v11: a whole-column selection on a geometry that has run
            if len(x) > mx[n, i]: bad.append(f'cmd {j} {n}[{i}] length {len(x)} > {mx[n, i]}')
            new = {c for c in x if not c.isalnum()} - chars[n]
            if new: bad.append(f'cmd {j} {n}[{i}] new characters {sorted(new)}: {x[:70]}')
            if (n, i) not in FREE and not (n in ('Aspirate', 'Dispense', 'Mix') and 2 <= i <= 13) and x not in vals[n, i]:
                bad.append(f'cmd {j} {n}[{i}] value {x} never used')
    # v16: tips 5-8 only where they have run - grid 28 positions 2-3 (v6, 30 Sep: at the 24-well row spacing they sit
    # 77 mm in front of tips 1-4 and cannot reach position 1), the 6-well plate's lanes and the 2W troughs
    for j, (n, a) in enumerate(load(path)):
        if n in ('Aspirate', 'Dispense', 'Mix') and int(a[0]) & 240 and (int(a[14]), int(a[15])) not in TIPS58_RUN:
            bad.append(f'cmd {j} {n} with tips 5-8 at grid {a[14]} site {a[15]}: never run there (out of reach at site 0 of grid 28)')
    stack = []
    for j, (n, a) in enumerate(load(path)):
        if n == 'IfThen': stack.append([j, False])
        elif n == 'Else': stack[-1][1] = True
        elif n == 'EndIf':
            k, has = stack.pop()
            if not has: bad.append(f'cmd {k} IfThen without Else (never used on the robot; crashes EVOware 2.8)')
    return bad
if __name__ == '__main__':
    bad = check(sys.argv[1])
    print('\n'.join(bad) if bad else 'inside the proven envelope: no new lengths, characters or fixed values')
