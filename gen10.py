"""Generates EDP_Titration_phyton_ML_hybrid.esc, v10 (v9 plus custom pools above the eluate's own pH: acidified to the
process pH, then Tris - 8-tip Tris from the 2W RIGHT well, a form the GCB scripts ran).
v9 (v8 plus the custom pool: a 6-well plate, 0006-SDWP-Greta, at grid 38
position 1 is added to the worktable, and at the end of a full run the operator can ask for eluate at any pH within
the run's acid curve, made in its column 1, measured and trimmed).
v8 notes follow. Every command form is one that has already run on the robot (the
envelope check proves it); the worktable header is copied verbatim from the reference script (2 Oct deck), with a new
date line and checksum.
v8 changes from v7: time estimates in every message that hands the robot back (EDP_time_estimate.py; step away or stay
nearby), neutralisation Cycle 1 is 4 wells in column 1, the bulk pools skip the wash between acid portions (the lane
holds only acid then), the validation refit (cycle 3) before archiving, and a results summary in the final message."""
import math, zlib, datetime

PY = r'C:\ProgramData\Miniconda3\python.exe'
DATA = 'C:\\ProgramData\\Tecan\\EVOware\\database\\tecan_data\\'
BR = '\f\x07'                                    # line break inside EVOware prompt text
ASP, DSP = 'Titration_tracking_acid_base', 'Titration_tracking'
WASH = 'Wash(255,1,1,1,0,"3",500,"4",500,10,70,30,1,0,1000,0);'
WASH_FIRST = 'Wash(255,1,1,1,0,"50",500,"20",500,10,70,30,1,0,1000,0);'
TROUGH = (2, 7)                                  # 2W - Acid & Base: left well 1 M acetic acid, right well 1 M Tris
PRODUCT = (2, 8)                                 # Product 2W: left well eluate
LPH = (22, 1)                                    # LpH: pools in lanes 1, 3, 5
CUSTOM = (38, 0)                                 # v9: 6-well plate (0006-SDWP-Greta) "Custom pool", grid 38 position 1
PLATES = {'ACID': (28, 0, 'Acid titration'), 'P32': (28, 1, 'Pool 1 titration'),
          'P34': (22, 2, 'Pool 2 titration'), 'P36': (28, 2, 'Pool 3 titration')}   # deck as saved in EVOware 2 Oct 2026
POOLS = [dict(tag='P32', n=1, lane=1, ratio='RATIO_T3', nominal='Target_pH3'),
         dict(tag='P34', n=2, lane=3, ratio='RATIO_T2', nominal='Target_pH2'),
         dict(tag='P36', n=3, lane=5, ratio='RATIO_T1', nominal='Target_pH1')]
ROWS4, ROWS8 = 'ABCD', 'ABCDEFGH'

def encode_wells(nx, ny, names):
    idx = set()
    for n in names:
        r, c = 'ABCDEFGHIJKLMNOP'.index(n[0]), int(n[1:])
        idx.add((c - 1) * ny + r)
    out = f"{nx:02X}{ny:02X}"
    for i in range(math.ceil(nx * ny / 7)):
        out += chr(48 + sum(1 << b for b in range(7) if i * 7 + b in idx))
    return out

def trough_wells(col, rows):   return encode_wells(2, 8, [f'{r}{col}' for r in rows])
def lane_wells(lane, rows):    return encode_wells(6, 8, [f'{r}{lane}' for r in rows])
def plate_wells(col):          return encode_wells(6, 4, [f'{r}{col}' for r in ROWS4])

class Script:
    def __init__(self):
        self.L, self.labels, self.jumps, self.loops = [], set(), set(), set()
    def add(self, s): self.L.append(s)
    def comment(self, text): self.add(f'Comment("{text}");')
    def label(self, text):
        assert text not in self.labels, text
        self.labels.add(text); self.comment(text)
    def group(self, text): self.add(f'Group("{text}");')
    def group_end(self): self.add('GroupEnd();')
    def var(self, name, expr, note='starting value'):
        self.add(f'Variable({name},"{expr}",0,"{note}",0,0.000000,50.000000,0,2,0,0);')
    def ask(self, name, default, prompt, lo=None, hi=None, text=False):
        assert '~' not in prompt, 'EVOware only substitutes ~variables~ in UserPrompt: ' + prompt
        if text:
            self.add(f'Variable({name},"{default}",1,"{prompt}",0,1.000000,50.000000,1,2,1,0);')
        elif lo is not None:
            self.add(f'Variable({name},"{default}",1,"{prompt}",1,{lo:.6f},{hi:.6f},0,2,1,0);')
        else:
            self.add(f'Variable({name},"{default}",1,"{prompt}",0,1.000000,50.000000,0,2,1,0);')
    def prompt(self, *lines):
        text = BR.join(lines)
        assert text.count('~') % 2 == 0, text
        self.add(f'UserPrompt("{text}",0,-1);')
    def export(self, names, fname):
        self.add(f'ExportVariable({"#".join(names)},"{DATA}{fname}",0,1);')
    def imp(self, name, fname, row=1, default='0'):
        self.add(f'ImportVariable({name},"{DATA}{fname}",0,"{default}",1,{row},1,1,0);')
    def py(self, script, args=''):
        self.add(f'Execute("{PY} {DATA}{script}{(" " + args) if args else ""}",2,"",2);')
    def jump_if(self, var, value, target, op=0):
        self.jumps.add(target)
        self.add(f'If("{var}",{op},"{value}","{target}");')
    def if_then(self, var, value, op=0):
        self.__dict__.setdefault('_ifs', []).append(False); self.add(f'IfThen("{var}",{op},"{value}");')
    def els(self):
        self._ifs[-1] = True; self.add('Else();')
    def end_if(self):
        assert self._ifs.pop(), 'IfThen without Else - EVOware 2.8 crashes at run start on this'
        self.add('EndIf();')
    def loop(self, count, name):
        assert name not in self.loops, name
        self.loops.add(name); self.add(f'BeginLoop("{count}","{name}");')
    def end_loop(self): self.add('EndLoop();')
    def wash(self): self.add(WASH)
    def timer(self, event, block): self.py('EDP_time_estimate.py', f'{event} {block}')
    def time_vars(self):
        for k, v in enumerate(('t_min', 't_clock', 't_away', 't_next')):
            self.imp(v, 'time_est.csv', k + 1)
    def handover(self, block, lines, nxt, after=None):
        """A message that hands the robot back: the estimate is planned first, then the message says how long the
        robot works and whether there is time to step away, then the block's start is logged. `after` replaces the
        default 'next: <nxt>' ending (e.g. the first of two top-ups)."""
        self.timer('plan', block)
        self.time_vars()
        away = after or f'You can step away - next: {nxt}.'
        stay = after or f'Please stay nearby - {nxt} follows.'
        self.if_then('t_away', '1')
        self.prompt(*lines, f'After OK the robot works about ~t_min~ min (to about ~t_clock~ if you click OK now). {away}')
        self.els()
        self.prompt(*lines, f'After OK the robot works only about ~t_min~ min (to about ~t_clock~). {stay}')
        self.end_if()
        self.timer('start', block)
    def pip(self, cmd, cls, mask, vols, grid, site, wells, extra='0,0'):
        v = [f'"{x}"' if x != 0 else '0' for x in vols] + ['0'] * (12 - len(vols))
        self.add(f'{cmd}({mask},"{cls}",{",".join(v)},{grid},{site},1,"{wells}",{extra});')

    # ------------------------------------------------------------------ one plate column, tips 1-4
    def column4(self, plan, group, plate, col, src, what):
        """Four wells of one plate column from a dispense plan (15 rows per column: n, uL x4 for eluate/pool,
        acid and Tris). src: 'PRODUCT' or ('LPH', lane). Tips 1-4 only. Each liquid in portions of at most 980 uL,
        tips washed after every liquid, then mixed 5 x 900 uL."""
        pg, ps, _ = PLATES[plate]
        dest = plate_wells(col)
        base = 15 * group
        for k, v in enumerate(['ne', 'e1', 'e2', 'e3', 'e4', 'na', 'a1', 'a2', 'a3', 'a4', 'nb', 'b1', 'b2', 'b3', 'b4']):
            self.imp(v, plan, base + k + 1)
        if src == 'PRODUCT':
            sg, ss, swells, sname = PRODUCT[0], PRODUCT[1], trough_wells(1, ROWS4), 'eluate'
        else:
            sg, ss, swells, sname = LPH[0], LPH[1], lane_wells(src[1], ROWS4), 'pool'
        liquids = [('e', 'ne', (sg, ss, swells), sname),
                   ('a', 'na', (TROUGH[0], TROUGH[1], trough_wells(1, ROWS4)), 'acid'),
                   ('b', 'nb', (TROUGH[0], TROUGH[1], trough_wells(2, ROWS4)), 'Tris')]
        for key, n, (g, s, w), name in liquids:
            vols = [f'{key}{i}' for i in range(1, 5)]
            self.if_then(n, '0', op=2)
            self.loop(n, f'{what}_{name}')
            self.pip('Aspirate', ASP, 15, vols, g, s, w)
            self.pip('Dispense', DSP, 15, vols, pg, ps, dest)
            if key != 'e':
                self.wash()                          # titrant: clean tips before the next portion
            self.end_loop()
            if key == 'e':
                self.wash()
            self.els(); self.comment(f'no {name} here'); self.end_if()
        self.pip('Mix', DSP, 15, ['900'] * 4, pg, ps, dest, extra='5,0,0')
        self.wash()

    # ------------------------------------------------------------------ typed readings with 3 attempts
    def entry(self, var, title, wells, n, out, stop_msg, fail='Run finished'):
        self.var('eok', '0', 'entry accepted flag')
        self.loop('3', f'Entry_{var}')
        self.if_then('eok', '0')
        if n == 1:
            q = (f'{title}: type the pH of the pool. (If this box appears again, the last entry was not a pH between 1 '
                 f'and 14 - please check and type it again.)')
        else:
            q = (f'{title}, {wells}: type the {n} pH values in that order, separated by commas. Type 0 for a well '
                 f'you could not measure. (If this box appears again, the last entry was not {n} values between pH 1 and 14, '
                 f'or was identical to an earlier entry - please check and type it again.)')
        self.ask(var, 'N/A', q, text=True)
        self.export([var], 'entry.csv')
        self.py('EDP_parse_readings.py', f'entry.csv {out} {n}')
        self.imp('eok', 'entry_ok.csv')
        self.els(); self.comment('entry accepted'); self.end_if()
        self.end_loop()
        self.if_then('eok', '0')
        self.prompt(f'The {stop_msg} readings could not be read after 3 attempts.', 'The run will now stop. What was typed is recorded in python_log.txt.')
        self.els(); self.comment('no action needed'); self.end_if()
        self.jump_if('eok', '0', fail)

    # ------------------------------------------------------------------ bulk low-pH pool in one LpH lane, 8 tips
    def bulk_pool(self, lane, ratio):
        dest = lane_wells(lane, ROWS8)
        self.var(f'e_L{lane}', f'Pool_Volume/(1+{ratio})/8', 'eluate per tip (uL)')
        self.var(f'a_L{lane}', f'e_L{lane}*{ratio}', 'acid per tip (uL)')
        for key, (g, s), name in (('a', TROUGH, 'acid'), ('e', PRODUCT, 'eluate')):    # acid first, eluate on top
            self.var('mod', f'{key}_L{lane}%980', 'remainder after 980 uL portions')
            self.var('nloop', f'({key}_L{lane}-mod)/980', 'number of 980 uL portions')
            for kind in ('full', 'rest'):
                if kind == 'full':
                    self.if_then('nloop', '0', op=2); self.loop('nloop', f'Loop_{name}_L{lane}'); vol = '980'
                else:
                    self.if_then('mod', '0', op=2); vol = 'mod'
                self.pip('Aspirate', ASP, 255, [vol] * 8, g, s, trough_wells(1, ROWS8))
                self.pip('Dispense', DSP, 255, [vol] * 8, LPH[0], LPH[1], dest)
                if key == 'e':
                    self.wash()                      # the lane holds acid: never carry it back into the eluate trough
                if kind == 'full':
                    self.end_loop(); self.els(); self.comment('no passes needed'); self.end_if()
                else:
                    self.els(); self.comment('no remainder'); self.end_if()
            if key == 'a':
                self.wash()                          # v8: acid portions go into a lane holding only acid - one wash
        self.pip('Mix', DSP, 255, ['980'] * 8, LPH[0], LPH[1], dest, extra='5,0,0')   # 30 mL lane: 5 + 3 cycles, all 8 tips
        self.pip('Mix', DSP, 255, ['980'] * 8, LPH[0], LPH[1], dest, extra='3,0,0')   # (cycle counts already used on the robot)
        self.wash()

    def lane_add(self, ul_var, src, col, name, tag, wash_each):
        """ul_var uL per tip from column `col` of a 2W trough into lane 1 of the custom pool plate, 8 tips, in portions of
        at most 980 uL. wash_each: a wash after every portion (whenever the lane already holds the pool, so the troughs
        never see it); otherwise one wash after all portions (the first acid, which goes into an empty lane).
        Nothing happens for 0 uL."""
        g, st = src
        self.var('mod', f'{ul_var}%980', 'remainder after 980 uL portions')
        self.var('nloop', f'({ul_var}-mod)/980', 'number of 980 uL portions')
        for kind in ('full', 'rest'):
            if kind == 'full':
                self.if_then('nloop', '0', op=2); self.loop('nloop', f'Loop_{tag}'); vol = '980'
            else:
                self.if_then('mod', '0', op=2); vol = 'mod'
            self.pip('Aspirate', ASP, 255, [vol] * 8, g, st, trough_wells(col, ROWS8))
            self.pip('Dispense', DSP, 255, [vol] * 8, CUSTOM[0], CUSTOM[1], lane_wells(1, ROWS8))
            if wash_each:
                self.wash()
            if kind == 'full':
                self.end_loop(); self.els(); self.comment('no passes needed'); self.end_if()
            else:
                self.els(); self.comment('no remainder'); self.end_if()
        if not wash_each:
            self.if_then(ul_var, '0', op=2); self.wash(); self.els(); self.comment(f'no {name} added'); self.end_if()

    def lane_mix(self):
        dest = lane_wells(1, ROWS8)
        self.pip('Mix', DSP, 255, ['980'] * 8, CUSTOM[0], CUSTOM[1], dest, extra='5,0,0')
        self.pip('Mix', DSP, 255, ['980'] * 8, CUSTOM[0], CUSTOM[1], dest, extra='3,0,0')
        self.wash()

    def needs(self, files, label):
        """Sum eluate/pool, acid and Tris needs from Python's needs_*.csv files into up_e, up_a, up_b (mL)."""
        self.var('up_e', '0', 'eluate or pool needed (mL)'); self.var('up_a', '0', 'acid needed (mL)'); self.var('up_b', '0', 'Tris needed (mL)')
        for f in files:
            for k, v in enumerate(('x_e', 'x_a', 'x_b')):
                self.imp(v, f, k + 1)
            self.var('up_e', 'up_e+x_e', label); self.var('up_a', 'up_a+x_a', label); self.var('up_b', 'up_b+x_b', label)


def build():
    s = Script()
    s.comment('EDP pH titration method v10')
    s.add(WASH_FIRST)
    s.add('ROMA(2,80,75,0,0,0,150,1,0);')
    s.py('EDP_tecan_reset.py')

    # ================================================================== set-up (all questions)
    s.group('Run set-up')
    s.var('one', '1', 'constant')
    s.var('Well_Volume', '3', 'liquid (mL) in every titration well')
    s.var('Pool_Volume', '30000', 'bulk pool volume (uL) in LpH columns 1, 3, 5')
    s.var('Product_Vol', '1', 'ratios are per 1 mL of eluate or pool')
    s.var('pool_ml', 'Pool_Volume/1000', 'pool volume in mL, for messages')
    s.var('Trough_Max', '100', 'mL a 2W trough well is filled to at most; larger needs are split into two top-ups')
    s.export(['Well_Volume'], 'well_limits.csv')
    for name, val in (('Max_Titrant', '4'), ('Target_pH1', '3.6'), ('Target_pH2', '3.4'), ('Target_pH3', '3.2'),
                      ('Process_Target_pH', '3.5'), ('ACCEPT_CRIT', '0.2'), ('nb_max', '0.5'), ('nb_t1', '7.2'),
                      ('nb_t2', '7.4'), ('nb_t3', '7.6'), ('nb_proc', '7.4'), ('RATIO_T1', '0'), ('RATIO_T2', '0'),
                      ('RATIO_T3', '0'), ('RATIO_PROCESS', '0'), ('st', '0'), ('c2v', '0'), ('vpass', '1'),
                      ('ok_P32', '0'), ('ok_P34', '0'), ('ok_P36', '0'), ('split', '0'), ('fin_P32', '0'), ('fin_P34', '0'), ('fin_P36', '0')):
        s.var(name, val)
    s.ask('BatchID', 'N/A', 'Please enter the Batch ID!', text=True)
    s.ask('Product_Name', 'N/A', 'Please enter the Product name!', text=True)
    s.ask('Product_Type', 'N/A', 'Please enter the Product type (e.g. ProA eluate)!', text=True)
    s.export(['BatchID', 'Product_Name', 'Product_Type'], 'batch_info.csv')
    s.ask('run_mode', '0', 'What would you like to run? 0 = acid titration, then neutralisation of the three pools it makes. '
          '1 = neutralisation only (low-pH pools already in LpH columns 1, 3 and 5).', 0, 1)
    s.ask('change', '0', 'Use the default settings? 0 = yes. 1 = change them. Defaults - acid: up to 4 mL acid per mL eluate, '
          'targets pH 3.6, 3.4, 3.2, process pH 3.5. Neutralisation: targets pH 7.2, 7.4, 7.6, process pH 7.4. '
          'Acceptance: 0.2 pH.', 0, 1)
    s.jump_if('change', '0', 'Set-up - settings done')
    s.jump_if('run_mode', '1', 'Set-up - neutralisation settings')
    s.ask('Max_Titrant', '4', 'Acid titration: largest acid ratio in Cycle 1, mL 1 M acetic acid per mL eluate (default 4). Cycle 1 tests 0, then 1/64 of it, doubling up to all of it', 0.05, 10)
    s.ask('Target_pH1', '3.6', 'Acid titration: highest target pH (default 3.6) - Pool 3', 2, 14)
    s.ask('Target_pH2', '3.4', 'Acid titration: middle target pH (default 3.4) - Pool 2', 2, 14)
    s.ask('Target_pH3', '3.2', 'Acid titration: lowest target pH (default 3.2) - Pool 1', 2, 14)
    s.ask('Process_Target_pH', '3.5', 'Acid titration: process target pH (default 3.5)', 2, 14)
    s.label('Set-up - neutralisation settings')
    s.ask('nb_t1', '7.2', 'Neutralisation: lowest target pH (default 7.2)', 2, 14)
    s.ask('nb_t2', '7.4', 'Neutralisation: middle target pH (default 7.4)', 2, 14)
    s.ask('nb_t3', '7.6', 'Neutralisation: highest target pH (default 7.6)', 2, 14)
    s.ask('nb_proc', '7.4', 'Neutralisation: process target pH (default 7.4)', 2, 14)
    s.ask('ACCEPT_CRIT', '0.2', 'Acceptance limit for the validations, in pH units (default 0.2)', 0.05, 0.5)
    s.label('Set-up - settings done')
    s.jump_if('run_mode', '0', 'Set-up - pools confirmed')
    # neutralisation only: the pools' acid ratios are needed to make Cycle 1 in the wells
    for v, row in (('RATIO_T1', 1), ('RATIO_T2', 2), ('RATIO_T3', 3)):
        s.imp(v, 'last_acid_ratios.csv', row)
    s.prompt('Neutralisation only: Cycle 1 of each pool is made in the wells from eluate + acid at the ratio the pool was made with.',
             'The last acid titration gave (mL 1 M acetic acid per mL eluate): Pool 1 (LpH column 1) ~RATIO_T3~, Pool 2 (column 3) ~RATIO_T2~, Pool 3 (column 5) ~RATIO_T1~.')
    s.ask('ratios_ok', '1', 'Were the pools made with those acid ratios? 1 = yes, 0 = no (type them next)', 0, 1)
    s.jump_if('ratios_ok', '1', 'Set-up - pools confirmed')
    s.ask('RATIO_T3', '0.4', 'Pool 1 (LpH column 1): mL 1 M acetic acid per mL eluate it was made with', 0.001, 20)
    s.ask('RATIO_T2', '0.25', 'Pool 2 (LpH column 3): mL 1 M acetic acid per mL eluate it was made with', 0.001, 20)
    s.ask('RATIO_T1', '0.15', 'Pool 3 (LpH column 5): mL 1 M acetic acid per mL eluate it was made with', 0.001, 20)
    s.label('Set-up - pools confirmed')
    s.export(['Product_Vol', 'Max_Titrant', 'Target_pH1', 'Target_pH2', 'Target_pH3', 'Process_Target_pH', 'ACCEPT_CRIT'], 'titration_params.csv')
    s.export(['one', 'one', 'one'], 'pools_selected.csv')
    s.group_end()
    s.jump_if('run_mode', '1', 'Neutralisation - plan Cycle 1')

    # ================================================================== acid titration
    s.py('EDP_tecan_gp_master.py', 'Acid 0')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration could not be planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('acid Cycle 1 planned'); s.end_if()
    s.jump_if('st', '0', 'Run finished')
    s.needs(['needs_c1.csv'], 'acid Cycle 1')
    s.handover('acid_c1', ('READY TO START - check before clicking OK:',
               '- System liquid filled, waste container emptied.',
               '- Product 2W LEFT well: eluate, at least ~up_e~ mL plus about 12 mL (right well not used).',
               '- Acid & Base 2W: LEFT well 1 M acetic acid, at least ~up_a~ mL plus about 10 mL; RIGHT well 1 M Tris.',
               '- Empty 24-well plates: Acid titration, Pool 1, Pool 3 at grid 28 positions 1-3; Pool 2 at grid 22 position 3.',
               '- Empty LpH plate at grid 22 position 2. Nothing in the door row (position 4).'), 'PLATE CHECK 1')
    s.group('Acid titration - Cycle 1 (acid plate, columns 1-2)')
    s.column4('plan_c1.csv', 0, 'ACID', 1, 'PRODUCT', 'AcidC1a')
    s.column4('plan_c1.csv', 1, 'ACID', 2, 'PRODUCT', 'AcidC1b')
    s.group_end()
    s.timer('end', 'acid_c1')
    s.comment('DOOR ACCESS POINT 1 - acid titration Cycle 1')
    s.prompt('PLATE CHECK 1 of 5 - Acid titration, Cycle 1', 'Open the safety door and take out the Acid titration plate.',
             'Measure columns 1 and 2: A1, B1, C1, D1, then A2, B2, C2, D2.', 'Return the plate, close the safety door, then click OK and type the readings.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.entry('rd_acid_c1', 'Acid titration Cycle 1', 'A1 B1 C1 D1 A2 B2 C2 D2', 8, 'acid_c1_measured.csv', 'acid titration Cycle 1')
    s.py('EDP_tecan_gp_master.py', 'Acid 1')
    s.imp('st', 'status.csv')
    s.imp('c2v', 'cycle2_vols.csv')
    s.if_then('c2v', '0')
    s.prompt('The acid titration analysis failed - no Cycle 2 was planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('Cycle 2 planned'); s.end_if()
    s.jump_if('c2v', '0', 'Run finished')
    s.timer('start', 'acid_c2')
    s.group('Acid titration - Cycle 2 (acid plate, column 3)')
    s.column4('plan_c2.csv', 0, 'ACID', 3, 'PRODUCT', 'AcidC2')
    s.group_end()
    s.timer('end', 'acid_c2')
    s.comment('DOOR ACCESS POINT 2 - acid titration Cycle 2')
    s.prompt('PLATE CHECK 2 of 5 - Acid titration, Cycle 2', 'Open the safety door and take out the Acid titration plate.',
             'Measure column 3: A3, B3, C3, D3.', 'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_c2', 'Acid titration Cycle 2', 'A3 B3 C3 D3', 4, 'cycle2_measured.csv', 'acid titration Cycle 2')
    s.py('EDP_tecan_gp_master.py', 'Acid 2')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration analysis failed, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.py('EDP_archive_run.py', 'Acid')
    s.els(); s.comment('acid ratios found'); s.end_if()
    s.jump_if('st', '0', 'Run finished')
    for k, v in enumerate(('RATIO_T1', 'RATIO_T2', 'RATIO_T3', 'RATIO_PROCESS')):
        s.imp(v, 'ratios_clean.csv', k + 1)
    s.jump_if('st', '2', 'Acid - ratios accepted')
    s.prompt('CHECK - some acid ratios lie beyond the measured curve and had to be estimated (less certain).',
             'mL acid per mL eluate: pH ~Target_pH1~: ~RATIO_T1~ | pH ~Target_pH2~: ~RATIO_T2~ | pH ~Target_pH3~: ~RATIO_T3~ | process pH ~Process_Target_pH~: ~RATIO_PROCESS~.',
             'Usually this means Cycle 1 did not go low enough: next time raise the acid ratio in the settings.')
    s.ask('accept', '0', 'Make the pools with these estimated ratios? 1 = yes, 0 = no (the acid titration is archived and the run stops)', 0, 1)
    s.if_then('accept', '0')
    s.py('EDP_archive_run.py', 'Acid')
    s.els(); s.comment('estimated ratios accepted'); s.end_if()
    s.jump_if('accept', '0', 'Run finished')
    s.label('Acid - ratios accepted')

    # ================================================================== neutralisation Cycle 1 plan (both modes)
    s.label('Neutralisation - plan Cycle 1')
    for P in POOLS:
        s.var(f'nom_{P["tag"]}', P['nominal'], 'nominal pool pH')
        s.var(f'pr_{P["tag"]}', P['ratio'], 'acid ratio of the pool')
        s.export(['one', 'nb_max', 'nb_t1', 'nb_t2', 'nb_t3', 'nb_proc', 'ACCEPT_CRIT', f'nom_{P["tag"]}', f'pr_{P["tag"]}'],
                 f'titration_params_{P["tag"]}.csv')
    s.py('EDP_tecan_gp_master.py', 'Base 0 all')
    for P in POOLS:
        s.imp(f'ok_{P["tag"]}', f'status_{P["tag"]}.csv')
    s.var('n_bad', '3-ok_P32-ok_P34-ok_P36', 'pools that could not be planned')
    s.var('stop_now', '0', '1 = stop the run')
    s.if_then('n_bad', '0', op=2)
    s.prompt('The neutralisation could not be planned for every pool (an acid ratio is missing), so the run will stop.',
             'The reason is in python_log.txt in the tecan_data folder.')
    s.var('stop_now', '1', '1 = stop the run')
    s.els(); s.comment('all three pools planned'); s.end_if()
    s.jump_if('stop_now', '1', 'Run finished')
    s.needs([f'needs_c1_{P["tag"]}.csv' for P in POOLS], 'neutralisation Cycle 1')
    s.jump_if('run_mode', '0', 'Neutralisation - full-run top-up')
    s.handover('neut_c1', ('READY TO START (neutralisation only) - check before OK:',
               '- System liquid filled, waste emptied.',
               '- LpH plate at grid 22 position 2: Pools 1, 2, 3 in columns 1, 3, 5, at least 22 mL each.',
               '- Product 2W LEFT well: eluate, at least ~up_e~ mL (Cycle 1 is made from eluate + acid).',
               '- Acid & Base 2W: LEFT 1 M acetic acid, at least ~up_a~ mL; RIGHT 1 M Tris, at least ~up_b~ mL plus about 15 mL.',
               '- Empty 24-well plates: Pool 1, Pool 3 at grid 28 positions 2-3; Pool 2 at grid 22 position 3. Door row empty.'), 'PLATE CHECK 3')
    s.jump_if('one', '1', 'Neutralisation - dispense Cycle 1')
    s.label('Neutralisation - full-run top-up')
    s.var('pe1', 'Pool_Volume/(1+RATIO_T1)', 'eluate for Pool 3 (uL)')
    s.var('pe2', 'Pool_Volume/(1+RATIO_T2)', 'eluate for Pool 2 (uL)')
    s.var('pe3', 'Pool_Volume/(1+RATIO_T3)', 'eluate for Pool 1 (uL)')
    s.var('pe', '(pe1+pe2+pe3)/1000', 'eluate for the pools (mL)')
    s.var('pa', '3*Pool_Volume/1000-pe', 'acid for the pools (mL)')
    s.imp('x_e', 'needs_val.csv', 1); s.imp('x_a', 'needs_val.csv', 2)
    s.var('c1_e', 'up_e', 'eluate for neutralisation Cycle 1 (mL)'); s.var('c1_a', 'up_a', 'acid for neutralisation Cycle 1 (mL)')
    s.var('up_e', 'c1_e+pe+x_e', 'eluate until the next plate check (mL)')
    s.var('up_a', 'c1_a+pa+x_a', 'acid until the next plate check (mL)')
    s.var('split', '0', 'two top-ups needed')
    s.if_then('up_e-Trough_Max', '0', op=2); s.var('split', '1', 'two top-ups needed'); s.els(); s.comment('eluate fits'); s.end_if()
    s.if_then('up_a-Trough_Max', '0', op=2); s.var('split', '1', 'two top-ups needed'); s.els(); s.comment('acid fits'); s.end_if()
    s.jump_if('split', '1', 'Neutralisation - split top-up')
    s.handover('big', ('TOP UP NOW - until the next plate check the robot makes the acid validation (column 4), the three bulk pools '
               '(~pool_ml~ mL each in LpH columns 1, 3, 5) and Cycle 1 of each pool.',
               'It will use about ~up_e~ mL eluate (Product 2W LEFT well), ~up_a~ mL acid (LEFT well) and ~up_b~ mL Tris (RIGHT well).',
               'Top up now if needed (add a few mL spare), then click OK.'), 'PLATE CHECK 3')
    s.jump_if('one', '1', 'Neutralisation - topped up')
    s.label('Neutralisation - split top-up')
    s.var('up_e', 'pe+x_e', 'eluate for acid validation and pools (mL)')
    s.var('up_a', 'pa+x_a', 'acid for acid validation and pools (mL)')
    s.handover('big1', ('TOP UP NOW (1 of 2) - next the robot makes the acid validation (column 4) and the three bulk pools '
               '(~pool_ml~ mL each in LpH columns 1, 3, 5).',
               'They use about ~up_e~ mL eluate (Product 2W LEFT well) and ~up_a~ mL acid (LEFT well). Click OK when topped up.'),
               'TOP UP 2 of 2', after='Then a second top-up is needed (TOP UP 2 of 2) before Cycle 1 of the pools.')
    s.label('Neutralisation - topped up')
    s.group('Acid titration - validation (acid plate, column 4)')
    s.column4('plan_val.csv', 0, 'ACID', 4, 'PRODUCT', 'AcidVal')
    s.group_end()
    s.group('Bulk low-pH pools - LpH columns 1, 3 and 5 (8 tips)')
    for P in POOLS:
        s.bulk_pool(P['lane'], P['ratio'])
    s.group_end()
    s.jump_if('split', '0', 'Neutralisation - dispense Cycle 1')
    s.timer('end', 'big1')
    s.var('up_e', 'c1_e', 'eluate for neutralisation Cycle 1 (mL)'); s.var('up_a', 'c1_a', 'acid for neutralisation Cycle 1 (mL)')
    s.handover('big2', ('TOP UP NOW (2 of 2) - the bulk pools are made. Cycle 1 of the three pools (made in the wells from eluate + acid, then Tris) '
               'uses about ~up_e~ mL eluate (Product 2W LEFT well), ~up_a~ mL acid (LEFT well) and ~up_b~ mL Tris (RIGHT well).',
               'Top up now if needed, then click OK.'), 'PLATE CHECK 3')
    s.label('Neutralisation - dispense Cycle 1')
    s.group('Neutralisation - Cycle 1 (eluate + acid + Tris, each pool plate column 1)')
    for P in POOLS:
        s.column4(f'plan_c1_{P["tag"]}.csv', 0, P['tag'], 1, 'PRODUCT', f'Pool{P["n"]}C1a')
    s.group_end()
    s.if_then('run_mode', '1')
    s.timer('end', 'neut_c1')
    s.els()
    s.if_then('split', '1')
    s.timer('end', 'big2')
    s.els()
    s.timer('end', 'big')
    s.end_if()
    s.end_if()

    # ------------------------------------------------------------------ plate check 3
    s.comment('DOOR ACCESS POINT 3 - acid validation and neutralisation Cycle 1')
    s.jump_if('run_mode', '1', 'Neutralisation - plate check 3 - neutralisation only')
    s.prompt('PLATE CHECK 3 of 5 - acid validation and neutralisation Cycle 1', 'Open the safety door and take out the four plates.',
             'Acid titration plate: measure column 4 (A4, B4, C4, D4).',
             'LpH plate: measure the three bulk pools (columns 1, 3, 5).',
             'Pool 1, Pool 2 and Pool 3 plates: measure column 1 of each (A1, B1, C1, D1).',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_val', 'Acid validation - Acid titration plate', 'A4 B4 C4 D4', 4, 'validation_results.csv', 'acid validation')
    s.entry('rd_pools', 'Bulk pools - LpH plate', 'columns 1 3 5 (Pool 1, Pool 2, Pool 3)', 3, 'pools_measured.csv', 'bulk pool')
    s.jump_if('one', '1', 'Neutralisation - Cycle 1 readings')
    s.label('Neutralisation - plate check 3 - neutralisation only')
    s.prompt('PLATE CHECK 3 of 5 - neutralisation Cycle 1', 'Open the safety door and take out the three pool plates.',
             'Pool 1, Pool 2 and Pool 3 plates: measure column 1 of each (A1, B1, C1, D1).',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.label('Neutralisation - Cycle 1 readings')
    for P in POOLS:
        s.entry(f'rd_{P["tag"]}_c1', f'Pool {P["n"]} Cycle 1 - Pool {P["n"]} titration plate', 'A1 B1 C1 D1', 4,
                f'neutralisation_c1_measured_{P["tag"]}.csv', f'Pool {P["n"]} Cycle 1')
    s.jump_if('run_mode', '1', 'Acid - archived')
    s.py('EDP_tecan_gp_master.py', 'Acid 3')
    s.py('EDP_archive_run.py', 'Acid')
    for k, v in enumerate(('ve1', 've2', 've3', 've4', 'vworst', 'vlimit', 'vpass', 'vn')):
        s.imp(v, 'validation_check.csv', k + 1)
    s.jump_if('vpass', '1', 'Acid - archived')
    s.prompt('ACID VALIDATION OUTSIDE THE LIMIT - the pools may be off target.',
             'Error at pH ~Target_pH1~: ~ve1~ | pH ~Target_pH2~: ~ve2~ | pH ~Target_pH3~: ~ve3~ | process pH ~Process_Target_pH~: ~ve4~ (limit ~vlimit~).')
    s.ask('go_on', '1', 'Continue with the neutralisation? (1 = continue, 0 = stop the run)', 0, 1)
    s.jump_if('go_on', '0', 'Run finished')
    s.label('Acid - archived')
    s.py('EDP_tecan_gp_master.py', 'Base 1 all')
    for P in POOLS:
        t = P['tag']
        s.imp(f'st_{t}', f'status_{t}.csv')
        s.imp(f'c2_{t}', f'cycle2_vols_{t}.csv')
        s.var(f'ok_{t}', '0', 'Cycle 2 planned')
        s.if_then(f'st_{t}', '0')
        s.if_then(f'c2_{t}', '0', op=2)
        s.var(f'ok_{t}', '1', 'Cycle 2 planned')
        s.els(); s.comment('no Cycle 2 volumes'); s.end_if()
        s.els(); s.comment('analysis failed'); s.end_if()
        s.if_then(f'ok_{t}', '0')
        s.prompt(f'Pool {P["n"]}: its Cycle 1 readings could not be analysed, so this pool is skipped.', 'The reason is in python_log.txt in the tecan_data folder.')
        s.els(); s.comment('pool continues'); s.end_if()
    s.var('n_ok', 'ok_P32+ok_P34+ok_P36', 'pools continuing')
    s.jump_if('n_ok', '0', 'Neutralisation - archive')
    s.needs([f'needs_c2_{P["tag"]}.csv' for P in POOLS], 'neutralisation Cycle 2')
    s.handover('neut_c2', ('TOP UP NOW if needed - Cycle 2 of the pools uses about ~up_e~ mL of the pools (LpH) and ~up_b~ mL Tris (RIGHT well).',
               'Click OK to continue.'), 'PLATE CHECK 4')
    s.group('Neutralisation - Cycle 2 (real pool from LpH + Tris, each pool plate column 3)')
    for P in POOLS:
        t = P['tag']
        s.jump_if(f'ok_{t}', '0', f'Neutralisation - Pool {P["n"]} Cycle 2 dispensed')
        s.column4(f'plan_c2_{t}.csv', 0, t, 3, ('LPH', P['lane']), f'Pool{P["n"]}C2')
        s.label(f'Neutralisation - Pool {P["n"]} Cycle 2 dispensed')
    s.group_end()
    s.timer('end', 'neut_c2')

    # ------------------------------------------------------------------ plate check 4
    s.comment('DOOR ACCESS POINT 4 - neutralisation Cycle 2')
    s.prompt('PLATE CHECK 4 of 5 - neutralisation Cycle 2', 'Open the safety door and take out the pool plates.',
             'Measure column 3 (A3, B3, C3, D3) of each pool plate that has it.',
             'Return the plates, close the safety door, then click OK and type the readings.')
    for P in POOLS:
        t = P['tag']
        s.jump_if(f'ok_{t}', '0', f'Neutralisation - Pool {P["n"]} Cycle 2 read')
        s.entry(f'rd_{t}_c2', f'Pool {P["n"]} Cycle 2 - Pool {P["n"]} titration plate', 'A3 B3 C3 D3', 4, f'cycle2_measured_{t}.csv', f'Pool {P["n"]} Cycle 2')
        s.label(f'Neutralisation - Pool {P["n"]} Cycle 2 read')
    s.py('EDP_tecan_gp_master.py', 'Base 2 all')
    for P in POOLS:
        t = P['tag']
        s.var(f'fin_{t}', '0', 'ratios found')
        s.jump_if(f'ok_{t}', '0', f'Neutralisation - Pool {P["n"]} ratios checked')
        s.imp(f'st_{t}', f'status_{t}.csv')
        s.if_then(f'st_{t}', '0', op=2)
        s.var(f'fin_{t}', '1', 'ratios found')
        s.els()
        s.prompt(f'Pool {P["n"]}: the neutralisation analysis failed, so this pool is not validated.', 'The reason is in python_log.txt in the tecan_data folder.')
        s.end_if()
        s.label(f'Neutralisation - Pool {P["n"]} ratios checked')
    s.var('n_fin', 'fin_P32+fin_P34+fin_P36', 'pools with ratios')
    s.jump_if('n_fin', '0', 'Neutralisation - archive')
    s.needs([f'needs_val_{P["tag"]}.csv' for P in POOLS], 'neutralisation validation')
    s.handover('neut_val', ('TOP UP NOW if needed - the validation uses about ~up_e~ mL of the pools (LpH) and ~up_b~ mL Tris (RIGHT well).',
               'Click OK to continue.'), 'PLATE CHECK 5')
    s.group('Neutralisation - validation (real pool from LpH + Tris, each pool plate column 4)')
    for P in POOLS:
        t = P['tag']
        s.jump_if(f'fin_{t}', '0', f'Neutralisation - Pool {P["n"]} validation dispensed')
        s.column4(f'plan_val_{t}.csv', 0, t, 4, ('LPH', P['lane']), f'Pool{P["n"]}Val')
        s.label(f'Neutralisation - Pool {P["n"]} validation dispensed')
    s.group_end()
    s.timer('end', 'neut_val')

    # ------------------------------------------------------------------ plate check 5
    s.comment('DOOR ACCESS POINT 5 - neutralisation validation')
    s.prompt('PLATE CHECK 5 of 5 - neutralisation validation', 'Open the safety door and take out the pool plates.',
             'Measure column 4 (A4, B4, C4, D4) of each pool plate that has it.',
             'Return the plates, close the safety door, then click OK and type the readings.')
    for P in POOLS:
        t = P['tag']
        s.jump_if(f'fin_{t}', '0', f'Neutralisation - Pool {P["n"]} validation read')
        s.entry(f'rd_{t}_val', f'Pool {P["n"]} validation - Pool {P["n"]} titration plate', 'A4 B4 C4 D4', 4, f'validation_results_{t}.csv', f'Pool {P["n"]} validation')
        s.label(f'Neutralisation - Pool {P["n"]} validation read')
    s.label('Neutralisation - archive')
    s.py('EDP_tecan_gp_master.py', 'Base 3 all')
    s.py('EDP_archive_run.py', 'Base all')
    s.label('Run finished')
    s.py('EDP_archive_run.py', 'Session')
    s.py('EDP_make_workbook.py')
    for k, v in enumerate(('sm_va', 'sm_v1', 'sm_v2', 'sm_v3', 'sm_a1', 'sm_a2', 'sm_a3', 'sm_a4',
                           'sm_b1', 'sm_b2', 'sm_b3', 'sm_f1', 'sm_f2', 'sm_f3')):
        s.imp(v, 'run_summary.csv', k + 1)
    s.prompt('Run finished. Validation, worst error in pH (limit ~ACCEPT_CRIT~; -1 = not validated): acid ~sm_va~ | Pool 1 ~sm_v1~ | Pool 2 ~sm_v2~ | Pool 3 ~sm_v3~.',
             'Ratios after the validation refit - acid, mL per mL eluate: pH ~Target_pH1~: ~sm_a1~ | ~Target_pH2~: ~sm_a2~ | ~Target_pH3~: ~sm_a3~ | process: ~sm_a4~.',
             'Tris at the process pH ~nb_proc~, mL per mL pool: Pool 1 ~sm_b1~ | Pool 2 ~sm_b2~ | Pool 3 ~sm_b3~. Bulk pool vs recipe Tris need: ~sm_f1~ | ~sm_f2~ | ~sm_f3~ (1 = as the recipe).',
             'All results: titration_data.xlsx (ML_data_archive); details: python_log.txt (tecan_data).')

    # ================================================================== custom pool (v9)
    s.jump_if('run_mode', '1', 'Custom pool - done')
    s.py('EDP_tecan_gp_master.py', 'Custom check')
    s.imp('cp_can', 'custom_status.csv', 1)
    s.jump_if('cp_can', '0', 'Custom pool - done')
    s.var('cp_ok', '0', 'custom pool planned')
    s.var('cp_ph', '0', 'custom pool pH')
    s.loop('3', 'Custom_ask')
    s.if_then('cp_ok', '0')
    s.ask('cp_ph', '0', 'CUSTOM POOL - type the pH to make about 25 mL of the product at, in column 1 of the 6-well plate '
          '(grid 38 position 1), from everything this run measured. Above the pH of the eluate itself, the eluate is first '
          'taken to the process pH with acid, then up with Tris. Type 0 for no custom pool.', 0, 14)
    s.export(['cp_ph'], 'custom_params.csv')
    s.py('EDP_tecan_gp_master.py', 'Custom 0')
    for k, v in enumerate(('cp_ok', 'cp_lo', 'cp_hi', 'cp_nlo', 'cp_nhi')):
        s.imp(v, 'custom_status.csv', k + 1)
    s.if_then('cp_ok', '0')
    s.prompt('CUSTOM POOL - pH ~cp_ph~ cannot be made from this run. Possible: about pH ~cp_lo~ to ~cp_hi~ with acid, and '
             '~cp_nlo~ to ~cp_nhi~ with acid then Tris (0 to 0 = no neutralisation data).',
             'Click OK and type another pH, or 0 for no custom pool.')
    s.els(); s.comment('custom pool pH accepted or not wanted'); s.end_if()
    s.els(); s.comment('custom pool already planned'); s.end_if()
    s.end_loop()
    s.jump_if('cp_ok', '0', 'Custom pool - done')
    s.jump_if('cp_ok', '2', 'Custom pool - done')
    for k, v in enumerate(('cp_a', 'cp_e', 'cp_ml')):
        s.imp(v, 'custom_plan.csv', k + 1)
    s.imp('cp_b', 'custom_plan.csv', 8)
    s.var('cp_aa', 'cp_a*8/1000', 'acid for the custom pool (mL)')
    s.var('cp_ea', 'cp_e*8/1000', 'eluate for the custom pool (mL)')
    s.var('cp_ba', 'cp_b*8/1000', 'Tris for the custom pool (mL)')
    s.handover('custom', ('CUSTOM POOL - the robot now makes ~cp_ml~ mL at pH ~cp_ph~ in column 1 of the 6-well plate (grid 38 '
               'position 1) from about ~cp_ea~ mL eluate (Product 2W LEFT well), ~cp_aa~ mL 1 M acetic acid (LEFT well) and '
               '~cp_ba~ mL 1 M Tris (RIGHT well; Tris only above the pH of the eluate itself).',
               'Place the empty 6-well plate, top up if needed, then click OK.'), 'measuring the custom pool')
    s.group('Custom pool - 6-well plate column 1 (8 tips)')
    s.lane_add('cp_a', TROUGH, 1, 'acid', 'CustA', wash_each=False)
    s.lane_add('cp_e', PRODUCT, 1, 'eluate', 'CustE', wash_each=True)
    s.lane_mix()
    s.if_then('cp_b', '0', op=2)
    s.lane_add('cp_b', TROUGH, 2, 'Tris', 'CustB', wash_each=True)
    s.lane_mix()
    s.els(); s.comment('no Tris - pool made with acid only'); s.end_if()
    s.group_end()
    s.timer('end', 'custom')
    for k in (1, 2, 3):
        s.prompt(f'CUSTOM POOL - measurement {k} of up to 3: measure the pool in column 1 of the 6-well plate (grid 38 position 1) with the probe.',
                 'Then click OK and type its pH.')
        s.entry(f'rd_cp{k}', 'Custom pool, 6-well plate column 1', 'column 1', 1, 'custom_measured.csv', 'custom pool',
                fail='Custom pool - finished')
        s.py('EDP_tecan_gp_master.py', f'Custom {k}')
        for j, v in enumerate(('cp_ta', 'cp_te', 'cp_done', 'cp_now')):
            s.imp(v, 'custom_trim.csv', j + 1)
        s.imp('cp_tb', 'custom_trim.csv', 7)
        s.jump_if('cp_done', '1', 'Custom pool - finished')
        if k < 3:
            s.var('cp_tva', 'cp_ta*8/1000', 'acid for the trim (mL)')
            s.var('cp_tve', 'cp_te*8/1000', 'eluate for the trim (mL)')
            s.var('cp_tvb', 'cp_tb*8/1000', 'Tris for the trim (mL)')
            s.handover(f'trim{k}', ('CUSTOM POOL - measured pH ~cp_now~ against target ~cp_ph~: the robot adds ~cp_tva~ mL acid, '
                       '~cp_tve~ mL eluate and ~cp_tvb~ mL Tris to the pool and mixes it.', 'Top up if needed, then click OK.'),
                       'the next custom pool measurement')
            s.group(f'Custom pool - trim {k}')
            s.lane_add('cp_ta', TROUGH, 1, 'acid', f'CustTA{k}', wash_each=True)
            s.lane_add('cp_te', PRODUCT, 1, 'eluate', f'CustTE{k}', wash_each=True)
            s.lane_add('cp_tb', TROUGH, 2, 'Tris', f'CustTB{k}', wash_each=True)
            s.lane_mix()
            s.group_end()
            s.timer('end', f'trim{k}')
    s.label('Custom pool - finished')
    s.py('EDP_archive_run.py', 'Custom')
    s.py('EDP_make_workbook.py')
    for k, v in enumerate(('cp_fph', 'cp_fv', 'cp_fa', 'cp_fe')):
        s.imp(v, 'custom_final.csv', k + 1)
    s.imp('cp_fb', 'custom_final.csv', 8)
    s.prompt('CUSTOM POOL DONE - column 1 of the 6-well plate holds about ~cp_fv~ mL at pH ~cp_fph~ (target ~cp_ph~): '
             '~cp_fe~ mL eluate + ~cp_fa~ mL 1 M acetic acid + ~cp_fb~ mL 1 M Tris.',
             'Recorded in titration_data.xlsx (Custom Pools sheet).')
    s.label('Custom pool - done')
    s.comment('End of script')
    missing = s.jumps - s.labels
    assert not missing, missing
    return s.L


def add_custom_plate(head):
    """v9: put 0006-SDWP-Greta (label "Custom pool") on the empty 4-position carrier at grid 38, site index 0. The
    worktable lists, after the carrier line (14;...), one line per grid without a carrier and two per grid with one
    (labware per site, then labels) - the same two lines EVOware rewrote when the deck changed on 2 Oct."""
    i = next(k for k, l in enumerate(head) if l.startswith('14;'))
    carriers = head[i].split(';')[1:]
    j = i + 1
    for g in range(CUSTOM[0]):
        j += 2 if carriers[g] not in ('-1', '') else 1
    assert carriers[CUSTOM[0]] == '151' and head[j] == '998;4;;;;;' and head[j + 1] == '998;;;;;', (head[j], head[j + 1])
    head[j], head[j + 1] = '998;4;0006-SDWP-Greta;;;;', '998;Custom pool;;;;'
    return head

def assemble(body, header_src, out):
    raw = open(header_src, 'rb').read().decode('latin-1').split('\r\n')
    head = add_custom_plate(raw[1:raw.index('--{ RPG }--') + 1])
    head[0] = datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + head[0][15:]
    rest = '\r\n'.join(head + body + ['']).encode('latin-1')
    crc = format((~zlib.crc32(rest)) & 0xffffffff, '08X')
    open(out, 'wb').write(crc.encode('latin-1') + b'\r\n' + rest)
    return crc

if __name__ == '__main__':
    import sys
    body = build()
    crc = assemble(body, sys.argv[1], sys.argv[2])
    print(len(body), 'lines, checksum', crc)
