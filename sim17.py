"""EVOware .esc simulator (v17: the dilute Tris in the Product 2W RIGHT well - the operator fills it at READY with 25 mL of
sc dil_true M Tris (default 0.1; a different value = a dilution made off); v16: the Acid titration plate in grid 28 position 3, tips 5-8 on the plates - reach and
hanging-liquid checks, the second dispense of an 8-tip trip timed, sc override = fixed variables set before the run;
v15: probe-drift hook sc shift; v13: optional absolute pipetting error, sc abs_sd / abs_bias in mL): runs the script command by command with a virtual deck, 8 fixed tips, the real Python
scripts (TECAN_DATA_DIR = a temp folder) and an operator who answers questions and types pH readings computed
from what is actually in each well (chem.py). Collects errors (dry troughs, overfilled wells, zero-count loops,
unknown variables, bad jumps, Python failures) and the true pH of the validation wells."""
import os, re, sys, shutil, subprocess, tempfile, math
import numpy as np
import chem
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'python'))
import EDP_time_estimate as TE          # v8: the same seconds per operation the time estimates use

HERE = os.path.dirname(os.path.abspath(__file__))
LABWARE = {  # (grid, site): name, nx, ny, container per 'column' or 'well', capacity mL
    (2, 7): ('2W Acid & Base', 2, 8, 'column', 150.0), (2, 8): ('Product 2W', 2, 8, 'column', 150.0),
    (22, 0): ('nLpH', 6, 4, 'well', 10.0), (22, 1): ('LpH', 6, 8, 'column', 31.0), (22, 2): ('Pool 2 titration', 6, 4, 'well', 3.4),
    (28, 0): ('Acid titration', 6, 4, 'well', 3.4), (28, 1): ('Pool 1 titration', 6, 4, 'well', 3.4),
    (28, 2): ('Acid titration', 6, 4, 'well', 3.4),                         # v16: the acid plate moved to position 3
    (38, 0): ('Custom pool', 6, 8, 'column', 31.0)}                         # v9: 6-well plate for the custom pool
TROUGH, PRODUCT, LPH, CUSTOM = (2, 7), (2, 8), (22, 1), (38, 0)
PLATE_OF = {'Acid': (28, 2), 'Tris': (28, 1), 'Base titration': (28, 1), 'Pool 1': (28, 1), 'Pool 2': (22, 2), 'Pool 3': (28, 2)}
# v16: where tips 5-8 can work. At the 24-well row spacing (19.3 mm) tips 5-8 sit 77 mm in front of tips 1-4, so they
# cannot reach grid 28 position 1 (site 0); positions 2-3 are where they ran on 30 Sep (v6). Plates in v16: sites 1, 2.
TIPS58_OK = {(28, 1), (28, 2), (22, 1), (2, 7), (2, 8)}
PLATE_SITES = {(28, 1), (28, 2)}

def decode_wells(ws):
    nx, ny = int(ws[0:2], 16), int(ws[2:4], 16)
    bits = []
    for ch in ws[4:]:
        v = ord(ch) - 48
        bits.extend(((v >> b) & 1) for b in range(7))
    return nx, ny, [(k // ny + 1, k % ny) for k, b in enumerate(bits[:nx * ny]) if b]

def split_args(s):
    args, cur, q = [], '', False
    for ch in s:
        if ch == '"':
            q = not q; cur += ch
        elif ch == ',' and not q:
            args.append(cur.strip()); cur = ''
        else:
            cur += ch
    args.append(cur.strip())
    return args

def unq(a): return a[1:-1] if len(a) >= 2 and a[0] == a[-1] == '"' else a

class Liquid:
    def __init__(self, vol=0.0, comp=None):
        self.vol, self.comp = vol, dict(comp or {})
    def add(self, vol, comp):
        if vol <= 0: return
        tot = self.vol + vol
        self.comp = {k: (self.comp.get(k, 0) * self.vol + comp.get(k, 0) * vol) / tot for k in chem.KEYS}
        self.vol = tot
    def take(self, vol):
        self.vol -= vol
        return dict(self.comp)

class Operator:
    def __init__(self, sc, rng):
        self.sc, self.rng, self.attempts, self.log = sc, rng, {}, []
    def answer(self, sim, name, prompt, default):
        if name.startswith('rd_'):
            k = self.attempts.get(name, 0); self.attempts[name] = k + 1
            scripted = self.sc.get('typed', {}).get(name)
            if scripted and k < len(scripted) and scripted[k] is not None:
                return scripted[k]
            return self.perturb(name, self.read(sim, prompt))
        ans = self.sc.get('answers', {}).get(name, default)
        if isinstance(ans, list):                          # v9: successive answers when a question is asked again
            k = self.attempts.get(name, 0); self.attempts[name] = k + 1
            ans = ans[min(k, len(ans) - 1)]
        return str(ans)
    def read(self, sim, prompt):
        head = prompt.split(':')[0]
        if 'Pools in the 6-well plate' in head or 'Pools in the Bulk Pools plate' in head:   # v11/v14: lanes listed for this round
            order = [int(float(x)) for x in open(os.path.join(sim.data, 'pools_order.csv')).read().split() if float(x) > 0]
            wells = [sim.cont(*LPH, c) for c in order]
        elif 'Custom pool' in head:
            wells = [sim.cont(*CUSTOM, 1)]
        elif 'Bulk pools' in head or 'LpH plate' in head:
            wells = [sim.cont(*LPH, c) for c in (1, 3, 5)]
        else:
            plate = next((PLATE_OF[k] for k in ('Base titration', 'Tris', 'Pool 1', 'Pool 2', 'Pool 3', 'Acid') if k in head), None)
            names = re.findall(r'\b([A-H])(\d)\b', head)
            if plate is None or not names:
                sim.err(f'operator cannot tell which wells to read: {head[:80]}'); return '0'
            wells = [sim.cont(plate[0], plate[1], int(c), 'ABCDEFGH'.index(r)) for r, c in names]
        vals = []
        if 'Eluate' in head:                                                # v11: the eluate's own pH first
            vals.append(f"{chem.ph_of(sim.cont(*PRODUCT, 1).comp, self.sc['T']) + self.rng.normal(0, self.sc.get('read_sd', 0.02)):.2f}")
        for w in wells:
            if w.vol <= 1e-9:
                vals.append('0')                                            # v13: a well left empty on purpose
            elif w.vol < 1.0:
                vals.append('0'); sim.warn(f'well read with only {w.vol:.2f} mL - typed 0')
            else:
                vals.append(f"{chem.ph_of(w.comp, self.sc['T']) + self.rng.normal(0, self.sc.get('read_sd', 0.02)):.2f}")
        return ','.join(vals)
    def perturb(self, name, text):
        """v11 test hook: sc['perturb'] = {entry name: (index, delta)} shifts one typed value (a typo or a bad well).
        v15: sc['shift'] = {entry name: delta} shifts every value of that entry (probe drift since the last plate check)."""
        d = self.sc.get('shift', {}).get(name)
        if d:
            text = ','.join(f"{float(x) + d:.2f}" if float(x) > 0 else x for x in text.split(','))
        p = self.sc.get('perturb', {}).get(name)
        if not p or self.attempts.get(name, 0) > 1:
            return text
        vals = text.split(',')
        vals[p[0]] = f"{float(vals[p[0]]) + p[1]:.2f}"
        return ','.join(vals)
    def on_prompt(self, sim, text):
        v = lambda n: float(sim.vars.get(n, 0) or 0)
        if 'After OK the robot works' in text:                             # v11: top up what the message asks for
            if text.startswith('READY TO START') and '0.1 M Tris' in text:   # v17: the dilute Tris, Product 2W RIGHT well
                sim.top_up(*PRODUCT, 2, 25.0, dict(Tr=self.sc.get('dil_true', 0.1)))
            sim.top_up(*PRODUCT, 1, min(v('up_e') + 5, 99), self.sc['eluate'])
            sim.top_up(*TROUGH, 1, min(v('up_a') + 5, 99), chem.ACID)
            sim.top_up(*TROUGH, 2, min(v('up_b') + 5, 99), chem.TRIS)
            return
        if text.startswith('READY TO START (neutralisation only)'):
            for c, liq in zip((1, 3, 5), self.sc['pools']):
                sim.fill(*LPH, c, self.sc.get('pool_ml', 28.0), liq)
            sim.top_up(*PRODUCT, 1, v('up_e') + 1, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('up_a') + 1, chem.ACID)
            sim.top_up(*TROUGH, 2, v('up_b') + 18, chem.TRIS)
        elif text.startswith('READY TO START'):
            sim.top_up(*PRODUCT, 1, v('up_e') + 12, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('up_a') + 10, chem.ACID)
        elif text.startswith('TOP UP NOW (1 of 2)'):
            sim.top_up(*PRODUCT, 1, v('up_e') + 1, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('up_a') + 1, chem.ACID)
        elif text.startswith('TOP UP NOW (2 of 2)'):
            sim.top_up(*PRODUCT, 1, v('up_e') + 1, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('up_a') + 1, chem.ACID)
            sim.top_up(*TROUGH, 2, v('up_b') + 1, chem.TRIS)
        elif text.startswith('TOP UP NOW - until the next plate check'):
            sim.top_up(*PRODUCT, 1, v('up_e') + 1, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('up_a') + 1, chem.ACID)
            sim.top_up(*TROUGH, 2, v('up_b') + 1, chem.TRIS)
        elif text.startswith('CUSTOM POOL - the robot now makes'):
            sim.top_up(*PRODUCT, 1, v('cp_ea') + 2, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('cp_aa') + 2, chem.ACID)
            sim.top_up(*TROUGH, 2, v('cp_ba') + 2, chem.TRIS)
        elif text.startswith('CUSTOM POOL - measured pH'):
            sim.top_up(*PRODUCT, 1, v('cp_tve') + 2, self.sc['eluate'])
            sim.top_up(*TROUGH, 1, v('cp_tva') + 2, chem.ACID)
            sim.top_up(*TROUGH, 2, v('cp_tvb') + 2, chem.TRIS)
        elif text.startswith('TOP UP NOW if needed'):
            sim.top_up(*TROUGH, 2, v('up_b') + 1, chem.TRIS)
            for c in (1, 3, 5):
                if sim.cont(*LPH, c).vol < v('up_e') / 3 + 1:
                    sim.warn(f'LpH lane {c} holds {sim.cont(*LPH, c).vol:.1f} mL before a step needing about {v("up_e") / 3:.1f} mL')

class Sim:
    def __init__(self, esc, sc, scripts):
        self.sc = sc
        self.clock, self.block_start, self.block_sim = 0.0, {}, {}
        raw = open(esc, 'rb').read().decode('latin-1').split('\r\n')
        body = raw[raw.index('--{ RPG }--') + 1:]
        self.cmds = []
        for i, line in enumerate(body):
            if not line.strip(): continue
            m = re.match(r'(\w+)\((.*)\);$', line)
            if not m: raise ValueError(f'unparsable line {i}: {line[:80]}')
            self.cmds.append((m.group(1), split_args(m.group(2)) if m.group(2) else []))
        self.match, stack = {}, []
        for i, (n, a) in enumerate(self.cmds):
            if n in ('IfThen', 'BeginLoop'): stack.append((n, i))
            elif n == 'Else': self.match[stack[-1][1], 'else'] = i; self.match[i, 'start'] = stack[-1][1]
            elif n in ('EndIf', 'EndLoop'):
                k, j = stack.pop()
                assert (k == 'IfThen') == (n == 'EndIf'), f'block mismatch at {i}'
                self.match[j, 'end'] = i; self.match[i, 'start'] = j
                if (j, 'else') in self.match: self.match[self.match[j, 'else'], 'end'] = i
        assert not stack, 'unclosed blocks'
        self.labels = {}
        for i, (n, a) in enumerate(self.cmds):
            if n == 'Comment':
                self.labels.setdefault(unq(a[0]), i)
        self.data = tempfile.mkdtemp(prefix='sim7_'); self.arch = os.path.join(self.data, 'archive'); os.makedirs(self.arch)
        for f in scripts: shutil.copy(f, self.data)
        if sc.get('old_file'):                          # v14 test: one helper file from an older version
            q = os.path.join(self.data, sc['old_file'])
            open(q, 'w').write(re.sub(r'EDP-v1[3-9]-2026-10', 'EDP-v12-2026-10', open(q).read()))   # v16: any current tag
        for f in sc.get('carry', {}): open(os.path.join(self.data, f), 'w').write(sc['carry'][f])
        self.vars, self.errors, self.warnings, self.prompts, self.questions, self.py_calls = {}, [], [], [], [], []
        self.containers, self.tips = {}, [Liquid() for _ in range(8)]
        self.rng = np.random.default_rng(sc.get('seed', 0))
        self.bias = self.rng.normal(0, sc.get('tip_bias', 0.0), 8)
        self.abs_bias = self.rng.normal(0, sc.get('abs_bias', 0.0), 8)   # v13: mL per aspiration, fixed per tip for the run
        self.op = Operator(sc, self.rng)
        self.peak = {}
        self.first_pipette = None

    def err(self, m):
        self.errors.append(f'[{self.pc}] {m}')
    def warn(self, m):
        self.warnings.append(f'[{self.pc}] {m}')
    def cont(self, g, s, col, row=0):
        name, nx, ny, per, cap = LABWARE[(g, s)]
        key = (g, s, col) if per == 'column' else (g, s, col, row)
        return self.containers.setdefault(key, Liquid())
    def fill(self, g, s, col, vol, comp):
        c = self.cont(g, s, col); c.vol, c.comp = vol, dict(comp)
    def top_up(self, g, s, col, want, comp):
        c = self.cont(g, s, col)
        if c.vol < want:
            c.add(want - c.vol, comp)
        self.peak[(g, s, col)] = max(self.peak.get((g, s, col), 0), c.vol)

    def value(self, expr):
        expr = expr.strip()
        if expr == '': return 0.0
        try:
            return float(expr)
        except ValueError:
            pass
        def sub(m):
            n = m.group(0)
            if n not in self.vars:
                raise NameError(n)
            v = self.vars[n]
            if isinstance(v, str):
                try: v = float(v)
                except ValueError: raise TypeError(f'text variable {n} used in arithmetic')
            return repr(float(v))
        try:
            return float(eval(re.sub(r'[A-Za-z_]\w*', sub, expr), {'__builtins__': {}}))
        except NameError as e:
            self.err(f'unknown variable {e} in "{expr}"'); return 0.0
        except Exception as e:
            self.err(f'cannot evaluate "{expr}": {e!r}'); return 0.0

    def text(self, t):
        def sub(m):
            n = m.group(1)
            if n not in self.vars:
                self.err(f'prompt uses unknown ~{n}~'); return '?'
            v = self.vars[n]
            return v if isinstance(v, str) else (f'{v:.4g}' if abs(v) < 1e5 else str(v))
        if t.count('~') % 2:
            self.err('odd number of ~ in prompt')
        return re.sub(r'~(\w*)~', sub, t)

    def cond(self, a, op, b):
        x, y = self.value(a), self.value(b)
        if op == '0': return abs(x - y) < 1e-9
        if op == '2': return x > y
        self.err(f'unsupported comparison operator {op}'); return False

    def path(self, p):
        p = unq(p)
        if 'tecan_data' not in p:
            self.err(f'file outside tecan_data: {p}')
        return os.path.join(self.data, p.split('\\')[-1])

    def pipette(self, kind, a):
        mask = int(a[0])
        ntips = bin(mask).count('1')                     # v8: simulated robot clock, same constants as the estimates
        if kind == 'Aspirate':
            self.clock += TE.T_PORTION8 if ntips == 8 else TE.T_PORTION
        elif kind == 'Mix':
            self.clock += (TE.T_MIX8 + int(a[18]) * TE.T_MIXCYCLE8) if ntips == 8 else (TE.T_MIX4 + TE.T_COLUMN)
        vols = [self.value(unq(x)) if x != '0' else 0.0 for x in a[2:14]]
        g, s, wells = int(a[14]), int(a[15]), unq(a[17])
        if (g, s) not in LABWARE:
            self.err(f'{kind} at grid {g} site {s}: nothing there in the worktable'); return
        name, nx, ny, per, cap = LABWARE[(g, s)]
        wx, wy, sel = decode_wells(wells)
        if (wx, wy) != (nx, ny):
            self.err(f'{kind}: well string is {wx}x{wy} but {name} is {nx}x{ny}')
        if s == 3:
            self.err(f'{kind} at site index 3 - the row closest to the door')
        if mask & 240 and (g, s) not in TIPS58_OK:
            self.err(f'{kind} with tips 5-8 at grid {g} site {s}: out of their reach (v16 rule)')
        if kind in ('Dispense', 'Mix') and (g, s) in PLATE_SITES and mask in (15, 240):
            other = range(4, 8) if mask == 15 else range(0, 4)        # the idle tip group waits one plate away
            over = (g, s + 1) if mask == 15 else (g, s - 1)
            if over in PLATE_SITES and any(getattr(self.tips[t], 'nominal', 0) > 0.002 for t in other):
                self.err(f'{kind} with tips {"1-4" if mask == 15 else "5-8"} at site {s}: the other tips still hold liquid over the plate at site {over[1]}')
            if kind == 'Mix' and over in PLATE_SITES and any(getattr(self.tips[t], 'wet', False) for t in other):
                self.err(f'Mix with tips {"1-4" if mask == 15 else "5-8"} at site {s}: the other tips, wet, wait over the plate at site {over[1]} for the whole mix')
            if kind == 'Dispense' and mask == 240 and per == 'well':
                self.clock += TE.T_PORTION2 - TE.T_PORTION8           # v16: the second dispense of an 8-tip trip
        tips = [t for t in range(8) if mask >> t & 1]
        if len(tips) != len(sel):
            self.err(f'{kind}: {len(tips)} tips but {len(sel)} wells'); return
        if self.first_pipette is None:
            self.first_pipette = self.pc
        for t in tips:
            self.tips[t].wet = True                       # v16: wet until the next wash
        for t, (col, row) in zip(tips, sel):
            v = vols[t] / 1000.0
            if kind == 'Mix':
                c = self.cont(g, s, col, row)
                if c.vol < v + 0.2:
                    self.err(f'mix {vols[t]:.0f} uL in {name} {"ABCDEFGH"[row]}{col} holding {c.vol:.2f} mL')
                continue
            if v <= 0:
                continue                                  # EVOware skips a tip given 0 uL (v6 did this on the instrument)
            if kind == 'Aspirate':
                if v > 1.0:
                    self.err(f'tip {t + 1} asked to hold {vols[t]:.0f} uL')
                c = self.cont(g, s, col, row)
                actual = v * (1 + self.bias[t] + self.rng.normal(0, self.sc.get('cv', 0.0)))
                actual = max(actual + self.abs_bias[t] + self.rng.normal(0, self.sc.get('abs_sd', 0.0)), 0.0)   # v13: small volumes
                if c.vol - actual < (0.5 if per == 'column' else 0.0):
                    self.err(f'{name} {"col" if per == "column" else "well"} {col}{"" if per == "column" else "ABCD"[row]} runs dry (holds {c.vol:.2f} mL, needs {actual:.2f})')
                comp = c.take(actual)
                self.tips[t].add(actual, comp)
                self.tips[t].src = (g, s)
                self.tips[t].nominal = getattr(self.tips[t], 'nominal', 0.0) + v
            else:
                tip = self.tips[t]
                nom = getattr(tip, 'nominal', 0.0)
                if nom < v - 1e-6:
                    self.err(f'tip {t + 1} dispenses {vols[t]:.0f} uL but holds {nom * 1000:.0f}'); continue
                share = tip.vol * v / nom
                if getattr(tip, 'src', None) == LPH and (g, s) in PLATE_OF.values():
                    share *= 1 - self.sc.get('pool_short', 0.0)
                if (g, s) == CUSTOM and getattr(tip, 'src', None) == TROUGH:     # v9 test: lane acid delivered off
                    share *= 1 + self.sc.get('custom_acid_bias', 0.0)
                comp = tip.take(tip.vol * v / nom); tip.nominal = nom - v
                c = self.cont(g, s, col, row)
                c.add(share, comp)
                if c.vol > cap + 1e-6:
                    self.err(f'{name} {"col " + str(col) if per == "column" else "ABCD"[row] + str(col)} overfilled: {c.vol:.2f} mL > {cap}')

    def run(self, max_steps=400000):
        self.pc, steps, loops = 0, 0, []
        while self.pc < len(self.cmds):
            steps += 1
            if steps > max_steps:
                self.err('step limit reached (endless loop?)'); break
            n, a = self.cmds[self.pc]
            nxt = self.pc + 1
            if n in ('Comment', 'Group', 'GroupEnd', 'ROMA'):
                pass
            elif n == 'Wash':
                self.clock += TE.T_WASH
                for t in self.tips:
                    if getattr(t, 'nominal', 0) > 0.002:
                        self.warn(f'wash discards {t.nominal * 1000:.0f} uL left in a tip')
                    t.vol, t.nominal, t.wet = 0.0, 0.0, False
            elif n == 'Variable':
                name, expr, query = a[0], unq(a[1]), a[2] == '1'
                is_text = a[7] == '1'
                if query:
                    prompt = unq(a[3])
                    if '~' in prompt:
                        self.err(f'~ in question text of {name}')
                    ans = self.op.answer(self, name, prompt, expr)
                    self.questions.append((self.pc, name, prompt[:90], ans))
                    if is_text:
                        self.vars[name] = ans
                    else:
                        try:
                            x = float(ans)
                        except ValueError:
                            self.err(f'answer {ans!r} to {name} is not a number'); x = 0.0
                        if a[4] == '1' and not (float(a[5]) - 1e-9 <= x <= float(a[6]) + 1e-9):
                            self.err(f'answer {x} to {name} outside {a[5]}-{a[6]}')
                        self.vars[name] = x
                else:
                    self.vars[name] = expr if is_text else self.value(expr)
                    if name in self.sc.get('override', {}):             # v16 test hook: e.g. Use_8_tips = 0
                        self.vars[name] = self.sc['override'][name]
            elif n == 'ExportVariable':
                vals = []
                for v in a[0].split('#'):
                    if v not in self.vars:
                        self.err(f'export of unknown variable {v}'); vals.append('0'); continue
                    x = self.vars[v]
                    vals.append(x if isinstance(x, str) else f'{x:.6f}')
                with open(self.path(a[1]), 'w') as f: f.write(','.join(vals) + '\n')
            elif n == 'ImportVariable':
                name, p, default, row, col = a[0], self.path(a[1]), unq(a[3]), int(a[5]), int(a[6])
                val = default
                try:
                    lines = [l for l in open(p).read().splitlines()]
                    if len(lines) >= row:
                        fields = lines[row - 1].split(',')
                        if len(fields) >= col and fields[col - 1].strip() != '':
                            val = fields[col - 1].strip()
                except FileNotFoundError:
                    pass
                try:
                    self.vars[name] = float(val)
                except ValueError:
                    self.vars[name] = val
            elif n == 'Execute':
                cmd = unq(a[0]).split()
                script = cmd[1].split('\\')[-1]
                env = dict(os.environ, TECAN_DATA_DIR=self.data, TECAN_ARCHIVE_DIR=self.arch, EDP_EXPORT_DIR=os.path.join(self.data, 'export'))
                if script == 'EDP_time_estimate.py' and len(cmd) >= 4 and cmd[2] in ('start', 'end'):
                    if cmd[2] == 'start':
                        self.block_start[cmd[3]] = self.clock
                    else:
                        self.clock += TE.T_MARKER
                        self.block_sim[cmd[3]] = self.clock - self.block_start.get(cmd[3], self.clock)
                if script == 'EDP_parse_readings.py' and self.sc.get('crash_parse') and self.sc['crash_parse'] in cmd:
                    self.py_calls.append(' '.join([script] + cmd[2:]) + ' (crashed - test)')     # v14 test: the 12 Oct crash
                    self.warn(f'parser crashed on {self.sc["crash_parse"]} (test)')
                    self.pc = nxt
                    continue
                r = subprocess.run([sys.executable, os.path.join(self.data, script)] + cmd[2:], env=env, capture_output=True, text=True, timeout=600)
                self.py_calls.append(' '.join([script] + cmd[2:]))
                if r.returncode:
                    self.err(f'{script} {" ".join(cmd[2:])} failed: {r.stderr[-400:]}')
            elif n == 'UserPrompt':
                t = self.text(unq(a[0]))
                self.prompts.append((self.pc, t))
                self.op.on_prompt(self, t)
            elif n == 'If':
                if self.cond(a[0].strip('"'), a[1], unq(a[2])):
                    label = unq(a[3])
                    if label not in self.labels:
                        self.err(f'jump to missing label {label}')
                    else:
                        if loops:
                            self.err(f'jump to "{label}" from inside a loop')
                        nxt = self.labels[label]
            elif n == 'IfThen':
                if not self.cond(a[0].strip('"'), a[1], unq(a[2])):
                    e = self.match.get((self.pc, 'else'))
                    nxt = (e + 1) if e is not None else self.match[self.pc, 'end'] + 1
            elif n == 'Else':
                nxt = self.match[self.pc, 'end'] + 1
            elif n == 'EndIf':
                pass
            elif n == 'BeginLoop':
                cnt = self.value(unq(a[0]))
                if cnt < 1 or abs(cnt - round(cnt)) > 1e-6:
                    self.err(f'loop {unq(a[1])} with count {cnt}')
                    nxt = self.match[self.pc, 'end'] + 1
                else:
                    loops.append([self.pc, int(round(cnt))])
            elif n == 'EndLoop':
                loops[-1][1] -= 1
                if loops[-1][1] > 0:
                    nxt = loops[-1][0] + 1
                else:
                    loops.pop()
            elif n in ('Aspirate', 'Dispense', 'Mix'):
                self.pipette(n, a)
            else:
                self.err(f'unknown command {n}')
            self.pc = nxt
        return self

    def report(self):
        log = open(os.path.join(self.data, 'python_log.txt')).read() if os.path.exists(os.path.join(self.data, 'python_log.txt')) else ''
        start_q = [q for q in self.questions if self.first_pipette is None or q[0] < self.first_pipette]
        return dict(errors=self.errors, warnings=self.warnings, py_errors=[l for l in log.splitlines() if 'ERROR' in l],
                    start_questions=[q[1] for q in start_q], n_prompts=len(self.prompts), log=log, data=self.data)

    def validation(self):
        """True pH (no reading noise) of every validation well, against its target."""
        out = {}
        def params(f):
            p = os.path.join(self.data, f)
            if not os.path.exists(p) or not open(p).read().strip(): return None
            x = [float(t) for t in open(p).read().split(',')]
            return [x[2], x[3], x[4], x[5]]
        col = 4
        tg = params('titration_params.csv')
        if tg and self.cont(*PLATE_OF['Acid'], col, 0).vol > 1:
            out['acid'] = [(t, round(chem.ph_of(self.cont(*PLATE_OF['Acid'], col, r).comp, self.sc['T']), 3)) for r, t in enumerate(tg)]
        for tag, (g, s) in (('P32', PLATE_OF['Pool 1']), ('P34', PLATE_OF['Pool 2']), ('P36', PLATE_OF['Pool 3'])):
            tg = params(f'titration_params_{tag}.csv')
            if tg and self.cont(g, s, col, 0).vol > 1:
                out[tag] = [(t, round(chem.ph_of(self.cont(g, s, col, r).comp, self.sc['T']), 3)) for r, t in enumerate(tg)]
        return out

def product(name):
    if name == 'A':
        q = np.load(os.path.join(HERE, 'eluate_fit_T.npy')); return dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3]), float(q[5])
    q = np.load(os.path.join(HERE, 'eluate_B.npy')); return dict(Ac=q[0], Na=q[1], Pc=q[2], Ph=q[3]), float(q[4])

SCRIPTS = ['EDP_tecan_gp_master.py', 'EDP_parse_readings.py', 'EDP_archive_run.py', 'EDP_tecan_reset.py', 'EDP_make_workbook.py',
           'EDP_time_estimate.py']

def run(esc, prod='A', pydir=None, **kw):
    el, T = product(prod)
    sc = dict(eluate=el, T=T, seed=0, cv=0.007, tip_bias=0.005, read_sd=0.02, answers={})
    sc.update(kw)
    pydir = pydir or os.path.join(HERE, '..', 'python')
    scripts = [os.path.join(pydir, f) if os.path.exists(os.path.join(pydir, f)) else os.path.join(HERE, f) for f in SCRIPTS]
    return Sim(esc, sc, scripts).run()
