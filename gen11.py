"""Generates EDP_Titration_python_ML_hybrid.esc (the old file name's "phyton" typo fixed), v11: the redesigned pH adjustment method (8 October 2026).

Flow: set-up (pools: lanes + target pHs, optional Tris titration) -> acid Cycle 1 (4 locate wells) -> plate check 1
(eluate pH + Cycle 1) -> acid Cycle 2 (8 wells across the located window) -> plate check 2 -> acid validation, pools
and Tris titration -> plate check 3 -> trims, new pools, Tris validation -> plate check 4 -> second trims -> plate
check 5 -> archive and results.

Liquid handling: eluate first, then acid, then Tris into each well or lane, with no wash between them (each liquid's
film stays in the same well instead of being swapped for wash water) and one wash after each column or pool; the
pool lanes are filled eluate-first, so acid never meets the eluate trough. Only command forms that have already run
on the robot are used (envelope.py), and the worktable header is the 2 Oct deck, verbatim: the pools go into the
6-well plate at grid 22 position 2 (the former LpH plate); grid 38 is not used."""
import datetime
import zlib
import gen10 as G
from gen10 import (WASH, WASH_FIRST, TROUGH, PRODUCT, LPH, ROWS4, ROWS8, ASP, DSP,
                   trough_wells, lane_wells, plate_wells)

PLATES = {'ACID': (28, 0, 'Acid titration'), 'TRIS': (28, 1, 'Tris titration')}
LANES = range(1, 7)


class Script(G.Script):
    def column(self, plan, group, plate, col, what):
        """Four wells of one plate column (tips 1-4) from a dispense plan: eluate, then acid, then Tris, each in its
        own portions of at most 980 uL, no wash between them; then 5 x 900 uL mix and one wash."""
        pg, ps, _ = PLATES[plate]
        dest = plate_wells(col)
        base = 15 * group
        for k, v in enumerate(['ne', 'e1', 'e2', 'e3', 'e4', 'na', 'a1', 'a2', 'a3', 'a4', 'nb', 'b1', 'b2', 'b3', 'b4']):
            self.imp(v, plan, base + k + 1)
        for key, n, (g, s, w), name in (('e', 'ne', (PRODUCT[0], PRODUCT[1], trough_wells(1, ROWS4)), 'eluate'),
                                         ('a', 'na', (TROUGH[0], TROUGH[1], trough_wells(1, ROWS4)), 'acid'),
                                         ('b', 'nb', (TROUGH[0], TROUGH[1], trough_wells(2, ROWS4)), 'Tris')):
            vols = [f'{key}{i}' for i in range(1, 5)]
            self.if_then(n, '0', op=2)
            self.loop(n, f'{what}_{name}')
            self.pip('Aspirate', ASP, 15, vols, g, s, w)
            self.pip('Dispense', DSP, 15, vols, pg, ps, dest)
            self.end_loop()
            self.els(); self.comment(f'no {name} here'); self.end_if()
        self.pip('Mix', DSP, 15, ['900'] * 4, pg, ps, dest, extra='5,0,0')
        self.wash()

    def lane_add(self, ul_var, src, col, lane, tag, wash_each=False):
        """ul_var uL per tip from column `col` of a 2W trough into one lane of the 6-well plate, 8 tips, portions of
        at most 980 uL. wash_each: a wash after every portion (eluate going into a pool that already holds acid)."""
        g, st = src
        self.var('mod', f'{ul_var}%980', 'remainder after 980 uL portions')
        self.var('nloop', f'({ul_var}-mod)/980', 'number of 980 uL portions')
        for kind in ('full', 'rest'):
            if kind == 'full':
                self.if_then('nloop', '0', op=2); self.loop('nloop', f'Loop_{tag}'); vol = '980'
            else:
                self.if_then('mod', '0', op=2); vol = 'mod'
            self.pip('Aspirate', ASP, 255, [vol] * 8, g, st, trough_wells(col, ROWS8))
            self.pip('Dispense', DSP, 255, [vol] * 8, LPH[0], LPH[1], lane_wells(lane, ROWS8))
            if wash_each:
                self.wash()
            if kind == 'full':
                self.end_loop(); self.els(); self.comment('no passes needed'); self.end_if()
            else:
                self.els(); self.comment('no remainder'); self.end_if()

    def lane_mix(self, lane, wash=True):
        dest = lane_wells(lane, ROWS8)
        self.pip('Mix', DSP, 255, ['980'] * 8, LPH[0], LPH[1], dest, extra='5,0,0')
        self.pip('Mix', DSP, 255, ['980'] * 8, LPH[0], LPH[1], dest, extra='3,0,0')
        if wash:
            self.wash()

    def lane_block(self, plan, lane, kind, p):
        """One lane of the 6-well plate from a lane plan (per lane: flag, acid, eluate, Tris uL per tip).
        make: eluate, acid, mix - then, for pools above the eluate's pH, Tris and a second mix - and one wash.
        trim: eluate (washed after each portion), acid, Tris, mix, wash."""
        names = [f'{p}f{lane}', f'{p}a{lane}', f'{p}e{lane}', f'{p}b{lane}']
        for k, v in enumerate(names):
            self.imp(v, plan, 4 * (lane - 1) + k + 1)
        f, a, e, b = names
        tag = f"{plan.split('.')[0].replace('pool_', 'P')}L{lane}"     # loop names must be unique in the script
        self.if_then(f, '0', op=2)
        if kind == 'make':
            self.lane_add(e, PRODUCT, 1, lane, f'{tag}E')
            self.lane_add(a, TROUGH, 1, lane, f'{tag}A')
            self.if_then(b, '0', op=2)
            self.lane_mix(lane, wash=False)
            self.lane_add(b, TROUGH, 2, lane, f'{tag}B')
            self.els(); self.comment('no Tris in this pool'); self.end_if()
            self.lane_mix(lane)
        else:
            self.lane_add(e, PRODUCT, 1, lane, f'{tag}E', wash_each=True)
            self.lane_add(a, TROUGH, 1, lane, f'{tag}A')
            self.lane_add(b, TROUGH, 2, lane, f'{tag}B')
            self.lane_mix(lane)
        self.els(); self.comment(f'lane {lane} - nothing to add'); self.end_if()

    def pool_order(self):
        for k in LANES:
            self.imp(f'po{k}', 'pools_order.csv', k)

    def entry(self, var, title, wells, n, out, stop_msg, fail='Run finished'):
        self.var('eok', '0', 'entry accepted flag')
        self.loop('3', f'Entry_{var}')
        self.if_then('eok', '0')
        if n == 'pools':
            q = (f'{title}, lanes as listed in the message: type one pH per pool in that order, separated by commas. Type 0 '
                 f'for a pool you could not measure. (If this box appears again, the number of values or a value was not '
                 f'right - please check and type it again.)')
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


def build():
    s = Script()
    s.comment('EDP pH adjustment method v11')
    s.add(WASH_FIRST)
    s.add('ROMA(2,80,75,0,0,0,150,1,0);')
    s.py('EDP_tecan_reset.py')

    # ================================================================== set-up
    s.group('Run set-up')
    s.var('one', '1', 'constant')
    s.var('Well_Volume', '3', 'liquid (mL) in every titration well')
    s.var('Product_Vol', '1', 'ratios are per 1 mL of eluate or pool')
    s.export(['Well_Volume'], 'well_limits.csv')
    for name, val in (('Max_Titrant', '4'), ('Target_pH1', '3.6'), ('Target_pH2', '3.4'), ('Target_pH3', '3.2'),
                      ('Process_Target_pH', '3.5'), ('ACCEPT_CRIT', '0.2'), ('nb_t1', '7.2'), ('nb_t2', '7.4'),
                      ('nb_t3', '7.6'), ('nb_proc', '7.4'), ('st', '0'), ('c2v', '0'), ('pool_ok', '0'), ('n_pools', '0'),
                      ('neut', '0'), ('n_ref', '0')):
        s.var(name, val)
    s.ask('BatchID', 'N/A', 'Please enter the Batch ID!', text=True)
    s.ask('Product_Name', 'N/A', 'Please enter the Product name!', text=True)
    s.ask('Product_Type', 'N/A', 'Please enter the Product type (e.g. ProA eluate)!', text=True)
    s.export(['BatchID', 'Product_Name', 'Product_Type'], 'batch_info.csv')
    s.loop('3', 'Pools_ask')
    s.if_then('pool_ok', '0')
    s.ask('pool_lanes', '0', 'POOLS in the 6-well plate (grid 22 position 2): type the lanes to fill, 1 to 6, separated by '
          'commas - for example 1,3,5 for the old bulk pools. Type 0 for no pools.', text=True)
    s.ask('pool_phs', '0', 'POOLS: type one target pH per lane, in the same order, separated by commas - for example '
          '3.2,3.4,3.6. Up to the pH of the eluate itself a pool is made with acid; from 6.8 to 8.5 with acid to the '
          'process pH, then Tris. Every pool is measured and trimmed.', text=True)
    s.export(['pool_lanes'], 'pool_lanes.csv')
    s.export(['pool_phs'], 'pool_phs.csv')
    s.py('EDP_tecan_gp_master.py', 'Setup')
    s.imp('pool_ok', 'setup_status.csv', 1)
    s.imp('n_pools', 'setup_status.csv', 2)
    s.if_then('pool_ok', '0')
    s.prompt('POOLS could not be read: type the lanes (1 to 6, each once), then the same number of target pH values '
             '(1 to 14), each list separated by commas.', 'Click OK to type them again. After 3 tries the run goes on without pools.')
    s.els(); s.comment('pool request read'); s.end_if()
    s.els(); s.comment('pool request already read'); s.end_if()
    s.end_loop()
    s.ask('neut', '0', 'TRIS TITRATION: 1 = also find the Tris needed to neutralise the process pool and the low-pH pools '
          '(VI method, 12 wells on the Tris plate). 0 = no. It runs anyway when a pool is above pH 6.8.', 0, 1)
    s.ask('change', '0', 'Use the default settings? 0 = yes. 1 = change them. Defaults - acid: up to 4 mL acid per mL eluate, '
          'targets pH 3.6, 3.4, 3.2, process pH 3.5. Tris: targets pH 7.2, 7.4, 7.6, process pH 7.4. '
          'Acceptance: 0.2 pH.', 0, 1)
    s.jump_if('change', '0', 'Set-up - settings done')
    s.ask('Max_Titrant', '4', 'Acid titration: largest acid ratio in Cycle 1, mL 1 M acetic acid per mL eluate (default 4). Cycle 1 tests 1/50 of it up to all of it in 4 steps', 0.05, 10)
    s.ask('Target_pH1', '3.6', 'Acid titration: highest target pH (default 3.6)', 2, 14)
    s.ask('Target_pH2', '3.4', 'Acid titration: middle target pH (default 3.4)', 2, 14)
    s.ask('Target_pH3', '3.2', 'Acid titration: lowest target pH (default 3.2)', 2, 14)
    s.ask('Process_Target_pH', '3.5', 'Acid titration: process target pH (default 3.5)', 2, 14)
    s.ask('nb_t1', '7.2', 'Tris titration: lowest target pH (default 7.2)', 2, 14)
    s.ask('nb_t2', '7.4', 'Tris titration: middle target pH (default 7.4)', 2, 14)
    s.ask('nb_t3', '7.6', 'Tris titration: highest target pH (default 7.6)', 2, 14)
    s.ask('nb_proc', '7.4', 'Tris titration: process target pH (default 7.4)', 2, 14)
    s.ask('ACCEPT_CRIT', '0.2', 'Acceptance limit for the validations, in pH units (default 0.2)', 0.05, 0.5)
    s.label('Set-up - settings done')
    s.export(['Product_Vol', 'Max_Titrant', 'Target_pH1', 'Target_pH2', 'Target_pH3', 'Process_Target_pH', 'ACCEPT_CRIT'], 'titration_params.csv')
    s.export(['nb_t1', 'nb_t2', 'nb_t3', 'nb_proc'], 'tris_params.csv')
    s.export(['neut'], 'options.csv')
    s.group_end()

    # ================================================================== acid titration: locate, then dense
    s.py('EDP_tecan_gp_master.py', 'Acid 0')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration could not be planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('acid Cycle 1 planned'); s.end_if()
    s.jump_if('st', '0', 'Run finished')
    s.needs(['needs_c1.csv'], 'acid Cycle 1')
    s.handover('acid_c1', ('READY TO START - check before clicking OK:',
               '- System liquid filled, waste container emptied.',
               '- Product 2W LEFT well: eluate, at least ~up_e~ mL (Cycles 1 and 2; pools need more later).',
               '- Acid & Base 2W: LEFT 1 M acetic acid, at least ~up_a~ mL; RIGHT 1 M Tris if needed.',
               '- Empty 24-well plates: Acid titration at grid 28 position 1; Tris titration at grid 28 position 2 if used.',
               '- Empty 6-well plate at grid 22 position 2 if pools are made. Nothing in the door row (position 4).'), 'PLATE CHECK 1')
    s.group('Acid titration - Cycle 1, locate (acid plate, column 1)')
    s.column('plan_c1.csv', 0, 'ACID', 1, 'AcidC1')
    s.group_end()
    s.timer('end', 'acid_c1')
    s.comment('DOOR ACCESS POINT 1 - eluate pH and acid titration Cycle 1')
    s.prompt('PLATE CHECK 1 - eluate pH and acid titration Cycle 1',
             'Measure the eluate itself with the probe (a few mL in a tube, or the Product 2W LEFT well).',
             'Open the safety door, take out the Acid titration plate and measure column 1 (A1, B1, C1, D1).',
             'Return the plate, close the safety door, then click OK and type the eluate pH followed by the four readings.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.entry('rd_acid_c1', 'Eluate and Acid titration Cycle 1', 'eluate then A1 B1 C1 D1', 5, 'acid_c1_measured.csv', 'acid titration Cycle 1')
    s.py('EDP_tecan_gp_master.py', 'Acid 1')
    s.imp('c2v', 'cycle2_vols.csv')
    s.imp('n_ref', 'pools_refused.csv')
    s.if_then('c2v', '0')
    s.prompt('The acid titration analysis failed - no Cycle 2 was planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('Cycle 2 planned'); s.end_if()
    s.jump_if('c2v', '0', 'Run finished')
    s.if_then('n_ref', '0', op=2)
    s.prompt('POOLS - ~n_ref~ requested pool(s) cannot be made from this eluate and will be skipped: their pH lies between '
             'the pH of the eluate itself and 6.8, or above 8.5. Details in python_log.txt.')
    s.els(); s.comment('every requested pool possible so far'); s.end_if()
    s.timer('start', 'acid_c2')
    s.group('Acid titration - Cycle 2, dense (acid plate, columns 2-3)')
    s.column('plan_c2.csv', 0, 'ACID', 2, 'AcidC2a')
    s.column('plan_c2.csv', 1, 'ACID', 3, 'AcidC2b')
    s.group_end()
    s.timer('end', 'acid_c2')
    s.comment('DOOR ACCESS POINT 2 - acid titration Cycle 2')
    s.prompt('PLATE CHECK 2 - acid titration Cycle 2', 'Open the safety door and take out the Acid titration plate.',
             'Measure columns 2 and 3: A2, B2, C2, D2, then A3, B3, C3, D3.',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_c2', 'Acid titration Cycle 2', 'A2 B2 C2 D2 A3 B3 C3 D3', 8, 'cycle2_measured.csv', 'acid titration Cycle 2')
    s.py('EDP_tecan_gp_master.py', 'Acid 2')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration analysis failed, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('acid ratios found'); s.end_if()
    s.jump_if('st', '0', 'Run finished')

    # ================================================================== acid validation, pools, Tris titration
    for k, v in enumerate(('n_mk1', 'n_basic', 'n_ref', 'tris_on', 'n_m1')):
        s.imp(v, 'pools_status.csv', k + 1)
    s.needs(['needs_main.csv'], 'acid validation, pools and Tris titration')
    s.pool_order()
    s.handover('main', ('NEXT - the robot makes the acid validation (Acid titration plate, column 4), ~n_mk1~ pool(s) in the '
               '6-well plate (grid 22 position 2) and, if requested, the Tris titration (Tris plate, columns 1 and 2).',
               'It uses about ~up_e~ mL eluate (Product 2W LEFT well), ~up_a~ mL acid (LEFT well) and ~up_b~ mL Tris '
               '(RIGHT well). Top up now if needed, then click OK.'), 'PLATE CHECK 3')
    s.group('Acid titration - validation (acid plate, column 4)')
    s.column('plan_val.csv', 0, 'ACID', 4, 'AcidVal')
    s.group_end()
    s.group('Pools - 6-well plate, eluate then acid (8 tips)')
    for lane in LANES:
        s.lane_block('pool_make1.csv', lane, 'make', 'mk')
    s.group_end()
    s.jump_if('tris_on', '0', 'Main - Tris titration done')
    s.group('Tris titration - Cycle 1, two acid fractions (Tris plate, columns 1-2)')
    s.column('plan_tris_c1.csv', 0, 'TRIS', 1, 'TrisC1a')
    s.column('plan_tris_c1.csv', 1, 'TRIS', 2, 'TrisC1b')
    s.group_end()
    s.label('Main - Tris titration done')
    s.timer('end', 'main')
    s.comment('DOOR ACCESS POINT 3 - acid validation, pools, Tris titration')
    s.prompt('PLATE CHECK 3 - acid validation, pools and Tris titration', 'Open the safety door and take out the plates.',
             'Acid titration plate: measure column 4 (A4, B4, C4, D4).',
             '6-well plate: measure the pools in lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Tris plate, if used: measure columns 1 and 2 (A1-D1, then A2-D2).',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_val', 'Acid validation - Acid titration plate', 'A4 B4 C4 D4', 4, 'validation_results.csv', 'acid validation')
    s.jump_if('n_m1', '0', 'Round 1 - pools read')
    s.entry('rd_pools1', 'Pools in the 6-well plate', '', 'pools', 'pools_measured_1.csv', 'pool')
    s.label('Round 1 - pools read')
    s.jump_if('tris_on', '0', 'Round 1 - Tris read')
    s.entry('rd_tris_c1', 'Tris titration Cycle 1 - Tris plate', 'A1 B1 C1 D1 A2 B2 C2 D2', 8, 'tris_c1_measured.csv', 'Tris titration Cycle 1')
    s.label('Round 1 - Tris read')
    s.py('EDP_tecan_gp_master.py', 'Acid 3')
    s.py('EDP_archive_run.py', 'Acid')
    s.jump_if('tris_on', '0', 'Round 1 - Tris fitted')
    s.py('EDP_tecan_gp_master.py', 'Tris 1')
    s.label('Round 1 - Tris fitted')
    s.py('EDP_tecan_gp_master.py', 'Pools 1')
    for k, v in enumerate(('r_work', 'n_m2', 'n_mk2')):
        s.imp(v, 'round2.csv', k + 1)
    s.jump_if('r_work', '0', 'Run finished')

    # ================================================================== round 2: trims, new pools, Tris validation
    s.needs(['needs_round2.csv'], 'pool trims, new pools and Tris validation')
    s.pool_order()
    s.handover('round2', ('NEXT - the robot trims the pools that missed their target, makes ~n_mk2~ new pool(s) and, if '
               'requested, the Tris validation (Tris plate, column 3).',
               'It uses about ~up_e~ mL eluate, ~up_a~ mL acid and ~up_b~ mL Tris. Top up if needed, then click OK.'), 'PLATE CHECK 4')
    s.group('Pools - trims and new pools (8 tips)')
    for lane in LANES:
        s.lane_block('pool_trim1.csv', lane, 'trim', 'tr')
    for lane in LANES:
        s.lane_block('pool_make2.csv', lane, 'make', 'mk')
    s.group_end()
    s.jump_if('tris_on', '0', 'Round 2 - Tris validation done')
    s.group('Tris titration - validation at the process pool (Tris plate, column 3)')
    s.column('plan_tris_val.csv', 0, 'TRIS', 3, 'TrisVal')
    s.group_end()
    s.label('Round 2 - Tris validation done')
    s.timer('end', 'round2')
    s.comment('DOOR ACCESS POINT 4 - pools and Tris validation')
    s.prompt('PLATE CHECK 4 - pools and Tris validation', 'Open the safety door and take out the plates.',
             '6-well plate: measure the pools in lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Tris plate, if used: measure column 3 (A3, B3, C3, D3).',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.jump_if('n_m2', '0', 'Round 2 - pools read')
    s.entry('rd_pools2', 'Pools in the 6-well plate', '', 'pools', 'pools_measured_2.csv', 'pool')
    s.label('Round 2 - pools read')
    s.jump_if('tris_on', '0', 'Round 2 - Tris read')
    s.entry('rd_tris_val', 'Tris validation - Tris plate', 'A3 B3 C3 D3', 4, 'tris_val_measured.csv', 'Tris validation')
    s.py('EDP_tecan_gp_master.py', 'Tris 2')
    s.py('EDP_archive_run.py', 'Base TR')
    s.label('Round 2 - Tris read')
    s.py('EDP_tecan_gp_master.py', 'Pools 2')
    for k, v in enumerate(('r_work', 'n_m3')):
        s.imp(v, 'round3.csv', k + 1)
    s.jump_if('r_work', '0', 'Run finished')

    # ================================================================== round 3: second trims
    s.needs(['needs_round3.csv'], 'second pool trims')
    s.pool_order()
    s.handover('round3', ('NEXT - the robot trims the pools that still missed their target.',
               'It uses about ~up_e~ mL eluate, ~up_a~ mL acid and ~up_b~ mL Tris. Top up if needed, then click OK.'), 'PLATE CHECK 5')
    s.group('Pools - second trims (8 tips)')
    for lane in LANES:
        s.lane_block('pool_trim2.csv', lane, 'trim', 'tr')
    s.group_end()
    s.timer('end', 'round3')
    s.comment('DOOR ACCESS POINT 5 - pools')
    s.prompt('PLATE CHECK 5 - pools', 'Open the safety door and take out the 6-well plate.',
             'Measure the pools in lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_pools3', 'Pools in the 6-well plate', '', 'pools', 'pools_measured_3.csv', 'pool')
    s.py('EDP_tecan_gp_master.py', 'Pools 3')

    # ================================================================== archive and results
    s.label('Run finished')
    s.py('EDP_archive_run.py', 'Pools')
    s.py('EDP_archive_run.py', 'Session')
    s.py('EDP_make_workbook.py')
    s.py('EDP_tecan_gp_master.py', 'Summary')
    for k, v in enumerate(('sm_va', 'sm_a1', 'sm_a2', 'sm_a3', 'sm_a4', 'sm_vb', 'sm_b1', 'sm_b2', 'sm_b3', 'sm_b4')):
        s.imp(v, 'run_summary.csv', k + 1)
    s.prompt('Run finished. Acid titration, mL 1 M acetic acid per mL eluate after the validation refit: pH ~Target_pH1~: '
             '~sm_a1~ | ~Target_pH2~: ~sm_a2~ | ~Target_pH3~: ~sm_a3~ | process ~Process_Target_pH~: ~sm_a4~. Worst validation '
             'error ~sm_va~ (limit ~ACCEPT_CRIT~; -1 = not validated).',
             'Tris titration, mL 1 M Tris per mL of the process pool: pH ~nb_t1~: ~sm_b1~ | ~nb_t2~: ~sm_b2~ | ~nb_t3~: '
             '~sm_b3~ | process ~nb_proc~: ~sm_b4~. Worst validation error ~sm_vb~ (-1 = not run).',
             'All results: titration_data.xlsx (ML_data_archive); details: python_log.txt (tecan_data).')
    s.imp('n_any', 'pools_any.csv')
    s.jump_if('n_any', '0', 'Run done')
    for lane in LANES:
        s.imp(f'pf{lane}', 'pool_final.csv', 2 * lane - 1)
        s.imp(f'pt{lane}', 'pool_final.csv', 2 * lane)
    s.prompt('POOLS in the 6-well plate - last measured pH (target); 0 = no pool, -1 = not made:',
             'lane 1: ~pf1~ (~pt1~) | lane 2: ~pf2~ (~pt2~) | lane 3: ~pf3~ (~pt3~) | lane 4: ~pf4~ (~pt4~) | lane 5: ~pf5~ '
             '(~pt5~) | lane 6: ~pf6~ (~pt6~).',
             "Volumes, ratios, every reading and each low-pH pool's Tris ratios: titration_data.xlsx, Pools sheet.")
    s.label('Run done')
    s.comment('End of script')
    missing = s.jumps - s.labels
    assert not missing, missing
    return s.L


def assemble(body, header_src, out):
    """Worktable header copied verbatim from an EVOware-saved reference (the 2 Oct deck), new date line, checksum."""
    raw = open(header_src, 'rb').read().decode('latin-1').split('\r\n')
    head = raw[1:raw.index('--{ RPG }--') + 1]
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
