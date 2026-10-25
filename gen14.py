"""Generates EDP_Titration_python_ML_hybrid.esc, v14 (12 October 2026).

The script's first job is to measure and store the product's acid, base and neutralisation ratios so its pH can be
adjusted later from any pH it is at (Edward, 12 Oct); all three titrations run in every run and the ratio map is
archived (ratio_map.csv). Custom pools at the end stay optional.

v14 against v13:
  - worktable header from Edward's EVOware-saved script of the 12 Oct run (labels "Bulk Pools" at grid 22 position 2,
    "Base titration" at grid 28 position 2); the messages use those names
  - a file-version check at the start: on 12 Oct the robot had v12 analysis files but older helpers, the old
    EDP_parse_readings.py crashed on the pool entry, and the previous entry's 'accepted' flag let the run go on
  - every reading entry resets that flag first, so a failed parse asks again instead of passing
  - the VI-neutralisation and base-titration questions are gone: both always run
  - pools holding Tris get 5 more mix cycles; the pool checks ask for each pool to be stirred before reading
  - two RATIOS STORED messages: acid and base; neutralisation from each acid pH
  - 'EDP_archive_run.py Map' archives the ratio map and the chemistry fingerprint

Flow: set-up -> version check -> acid Cycle 1 + base Cycle 1 -> plate check 1 -> acid Cycle 2 + base Cycle 2 -> plate
check 2 -> acid validation + VI Cycle 1 + base validation -> plate check 3 -> VI validation -> plate check 4 -> ratio
map -> ratios stored (2 messages) -> custom pool ranges -> pool request -> pools -> pool checks -> archive, results."""
import gen11
from gen11 import Script as Script11, LANES, assemble, WASH_FIRST, PRODUCT, TROUGH, LPH, ROWS8, DSP
from gen10 import lane_wells


class Script(Script11):
    def lane_block(self, plan, lane, kind, p):
        """One lane of the 6-well plate from a lane plan (per lane: flag, acid, eluate, Tris uL per tip).
        make (v13): eluate, then acid or Tris (no pool holds both any more), mix, one wash.
        trim: eluate (washed after each portion), acid, Tris, mix, wash."""
        names = [f'{p}f{lane}', f'{p}a{lane}', f'{p}e{lane}', f'{p}b{lane}']
        for k, v in enumerate(names):
            self.imp(v, plan, 4 * (lane - 1) + k + 1)
        f, a, e, b = names
        tag = f"{plan.split('.')[0].replace('pool_', 'P')}L{lane}"     # loop names must be unique in the script
        self.if_then(f, '0', op=2)
        self.lane_add(e, PRODUCT, 1, lane, f'{tag}E', wash_each=(kind == 'trim'))
        self.lane_add(a, TROUGH, 1, lane, f'{tag}A')
        self.lane_add(b, TROUGH, 2, lane, f'{tag}B')
        self.if_then(b, '0', op=2)            # v14: Tris (denser than eluate) gets 5 more mix cycles
        self.pip('Mix', DSP, 255, ['980'] * 8, LPH[0], LPH[1], lane_wells(lane, ROWS8), extra='5,0,0')
        self.els(); self.comment('no Tris in this lane'); self.end_if()
        self.lane_mix(lane)
        self.els(); self.comment(f'lane {lane} - nothing to add'); self.end_if()

    def entry(self, var, title, wells, n, out, stop_msg, fail='Run finished'):
        """v14: as v11, but the 'accepted' flag file is reset to 0 before every parse, so a parse that fails (on 12 Oct
        an older EDP_parse_readings.py crashed on the pool entry) asks again instead of passing on the last entry's 1."""
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
        self.export(['eok'], 'entry_ok.csv')
        self.py('EDP_parse_readings.py', f'entry.csv {out} {n}')
        self.imp('eok', 'entry_ok.csv')
        self.els(); self.comment('entry accepted'); self.end_if()
        self.end_loop()
        self.if_then('eok', '0')
        self.prompt(f'The {stop_msg} readings could not be read after 3 attempts.', 'The run will now stop. What was typed is recorded in python_log.txt.')
        self.els(); self.comment('no action needed'); self.end_if()
        self.jump_if('eok', '0', fail)

    def skip_unless(self, var, label):
        """Jump forward to `label` when var is 0 (EVOware conditions: expression left, literal right)."""
        self.jump_if(var, '0', label)


def build():
    s = Script()
    s.comment('EDP pH adjustment method v14')
    s.add(WASH_FIRST)
    s.add('ROMA(2,80,75,0,0,0,150,1,0);')
    s.py('EDP_tecan_reset.py')
    s.var('vchk', '0', 'file version check')
    s.var('vbad', '0', 'first file with the wrong version')
    s.export(['vchk'], 'check_status.csv')
    s.py('EDP_tecan_gp_master.py', 'Check')
    s.imp('vchk', 'check_status.csv', 1)
    s.imp('vbad', 'check_status.csv', 2)
    s.if_then('vchk', '0')
    s.prompt('The Python files in tecan_data are not all v14, so the run will stop before anything is pipetted.',
             'First file found older: ~vbad~ (1 EDP_tecan_gp_master.py, 2 EDP_parse_readings.py, 3 EDP_archive_run.py, '
             '4 EDP_tecan_reset.py, 5 EDP_make_workbook.py, 6 EDP_time_estimate.py; 0 = the analysis file itself is older).',
             'Copy all six EDP Python files from the python folder of the zip into tecan_data, then start the run again.')
    s.els(); s.comment('all six EDP files are v14'); s.end_if()
    s.jump_if('vchk', '0', 'Run done')

    # ================================================================== set-up
    s.group('Run set-up')
    s.var('one', '1', 'constant')
    s.var('Well_Volume', '3', 'liquid (mL) in every titration well')
    s.var('Product_Vol', '1', 'ratios are per 1 mL of eluate or pool')
    s.export(['Well_Volume'], 'well_limits.csv')
    for name, val in (('Max_Titrant', '4'), ('Target_pH1', '3.6'), ('Target_pH2', '3.4'), ('Target_pH3', '3.2'),
                      ('Process_Target_pH', '3.5'), ('ACCEPT_CRIT', '0.2'), ('nb_t1', '7.2'), ('nb_t2', '7.4'),
                      ('nb_t3', '7.6'), ('nb_proc', '7.4'), ('st', '0'), ('c2v', '0'), ('pool_ok', '0'), ('n_pools', '0'),
                      ('neut', '1'), ('base_on', '1'), ('tris_on', '0')):
        s.var(name, val)
    s.ask('BatchID', 'N/A', 'Please enter the Batch ID!', text=True)
    s.ask('Product_Name', 'N/A', 'Please enter the Product name!', text=True)
    s.ask('Product_Type', 'N/A', 'Please enter the Product type (e.g. ProA eluate)!', text=True)
    s.export(['BatchID', 'Product_Name', 'Product_Type'], 'batch_info.csv')
    s.ask('change', '0', 'Use the default settings? 0 = yes. 1 = change them. Defaults - acid: up to 4 mL acid per mL eluate, '
          'targets pH 3.6, 3.4, 3.2, process pH 3.5. Tris (base and neutralisation): targets pH 7.2, 7.4, 7.6, process pH 7.4. '
          'Acceptance: 0.2 pH.', 0, 1)
    s.jump_if('change', '0', 'Set-up - settings done')
    s.ask('Max_Titrant', '4', 'Acid titration: largest acid ratio in Cycle 1, mL 1 M acetic acid per mL eluate (default 4). Cycle 1 tests 1/50 of it up to all of it in 4 steps', 0.05, 10)
    s.ask('Target_pH1', '3.6', 'Acid titration: highest target pH (default 3.6)', 2, 14)
    s.ask('Target_pH2', '3.4', 'Acid titration: middle target pH (default 3.4)', 2, 14)
    s.ask('Target_pH3', '3.2', 'Acid titration: lowest target pH (default 3.2)', 2, 14)
    s.ask('Process_Target_pH', '3.5', 'Acid titration: process target pH (default 3.5)', 2, 14)
    s.ask('nb_t1', '7.2', 'Tris titrations: lowest target pH (default 7.2)', 2, 14)
    s.ask('nb_t2', '7.4', 'Tris titrations: middle target pH (default 7.4)', 2, 14)
    s.ask('nb_t3', '7.6', 'Tris titrations: highest target pH (default 7.6)', 2, 14)
    s.ask('nb_proc', '7.4', 'Tris titrations: process target pH (default 7.4)', 2, 14)
    s.ask('ACCEPT_CRIT', '0.2', 'Acceptance limit for the validations, in pH units (default 0.2)', 0.05, 0.5)
    s.label('Set-up - settings done')
    s.export(['Product_Vol', 'Max_Titrant', 'Target_pH1', 'Target_pH2', 'Target_pH3', 'Process_Target_pH', 'ACCEPT_CRIT'], 'titration_params.csv')
    s.export(['nb_t1', 'nb_t2', 'nb_t3', 'nb_proc'], 'tris_params.csv')
    s.export(['neut', 'base_on'], 'options.csv')
    s.group_end()

    # ================================================================== Cycle 1: acid (locate) + base (locate)
    s.py('EDP_tecan_gp_master.py', 'Acid 0')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration could not be planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('acid Cycle 1 planned'); s.end_if()
    s.jump_if('st', '0', 'Run finished')
    s.py('EDP_tecan_gp_master.py', 'Base 0')
    s.needs(['needs_c1.csv', 'needs_base_c1.csv'], 'acid and base Cycle 1')
    s.handover('acid_c1', ('READY TO START - check before clicking OK:',
               '- System liquid filled, waste container emptied.',
               '- Product 2W LEFT well: eluate, at least ~up_e~ mL (Cycles 1 and 2; pools need more later).',
               '- 2W - Acid & Base: LEFT 1 M acetic acid, at least ~up_a~ mL; RIGHT 1 M Tris, at least ~up_b~ mL.',
               '- Empty 24-well plates: Acid titration (grid 28 position 1) and Base titration (grid 28 position 2).',
               '- Empty Bulk Pools plate (6-well, grid 22 position 2) for custom pools. Nothing in the door row.'), 'PLATE CHECK 1')
    s.group('Cycle 1 - acid (Acid titration plate, column 1) and base (Base titration plate, column 4)')
    s.column('plan_c1.csv', 0, 'ACID', 1, 'AcidC1')
    s.skip_unless('base_on', 'Cycle 1 - base done')
    s.column('plan_base_c1.csv', 0, 'TRIS', 4, 'BaseC1')
    s.label('Cycle 1 - base done')
    s.group_end()
    s.timer('end', 'acid_c1')
    s.comment('DOOR ACCESS POINT 1 - eluate pH, acid and base Cycle 1')
    s.if_then('base_on', '1')
    s.prompt('PLATE CHECK 1 - eluate pH, acid and base titration Cycle 1',
             'Measure the eluate itself with the probe (a few mL in a tube, or the Product 2W LEFT well).',
             'Open the safety door. Acid titration plate: measure column 1 (A1-D1). Base titration plate: measure column 4 (A4-D4).',
             'Return the plates, close the safety door, click OK, then type the eluate pH and A1-D1, then A4-D4.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.els()
    s.prompt('PLATE CHECK 1 - eluate pH and acid titration Cycle 1',
             'Measure the eluate itself with the probe (a few mL in a tube, or the Product 2W LEFT well).',
             'Open the safety door, take out the Acid titration plate and measure column 1 (A1, B1, C1, D1).',
             'Return the plate, close the safety door, then click OK and type the eluate pH followed by the four readings.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.end_if()
    s.entry('rd_acid_c1', 'Eluate and Acid titration Cycle 1', 'eluate then A1 B1 C1 D1', 5, 'acid_c1_measured.csv', 'acid titration Cycle 1')
    s.py('EDP_tecan_gp_master.py', 'Acid 1')
    s.imp('c2v', 'cycle2_vols.csv')
    s.if_then('c2v', '0')
    s.prompt('The acid titration analysis failed - no Cycle 2 was planned, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('Cycle 2 planned'); s.end_if()
    s.jump_if('c2v', '0', 'Run finished')
    s.skip_unless('base_on', 'Base Cycle 1 entered')
    s.entry('rd_base_c1', 'Base titration Cycle 1 - Base titration plate', 'A4 B4 C4 D4', 4, 'base_c1_measured.csv', 'base titration Cycle 1')
    s.py('EDP_tecan_gp_master.py', 'Base 1')
    s.imp('base_on', 'base_status.csv')
    s.label('Base Cycle 1 entered')

    # ================================================================== Cycle 2: acid (dense) + base (across the window)
    s.timer('start', 'acid_c2')
    s.group('Cycle 2 - acid (Acid titration plate, cols 2-3) and base (Base titration plate, col 5)')
    s.column('plan_c2.csv', 0, 'ACID', 2, 'AcidC2a')
    s.column('plan_c2.csv', 1, 'ACID', 3, 'AcidC2b')
    s.skip_unless('base_on', 'Cycle 2 - base done')
    s.column('plan_base_c2.csv', 0, 'TRIS', 5, 'BaseC2')
    s.label('Cycle 2 - base done')
    s.group_end()
    s.timer('end', 'acid_c2')
    s.comment('DOOR ACCESS POINT 2 - acid and base Cycle 2')
    s.if_then('base_on', '1')
    s.prompt('PLATE CHECK 2 - acid and base titration Cycle 2',
             'Open the safety door. Acid titration plate: measure columns 2 and 3 (A2-D2, then A3-D3). Base titration plate: measure column 5 (A5-D5).',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.els()
    s.prompt('PLATE CHECK 2 - acid titration Cycle 2', 'Open the safety door and take out the Acid titration plate.',
             'Measure columns 2 and 3: A2, B2, C2, D2, then A3, B3, C3, D3.',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.end_if()
    s.entry('rd_acid_c2', 'Acid titration Cycle 2', 'A2 B2 C2 D2 A3 B3 C3 D3', 8, 'cycle2_measured.csv', 'acid titration Cycle 2')
    s.py('EDP_tecan_gp_master.py', 'Acid 2')
    s.imp('st', 'status.csv')
    s.if_then('st', '0')
    s.prompt('The acid titration analysis failed, so the run will stop.', 'The reason is in python_log.txt in the tecan_data folder.')
    s.els(); s.comment('acid ratios found'); s.end_if()
    s.jump_if('st', '0', 'Run finished')
    s.imp('tris_on', 'val_status.csv')
    s.skip_unless('base_on', 'Base Cycle 2 entered')
    s.entry('rd_base_c2', 'Base titration Cycle 2 - Base titration plate', 'A5 B5 C5 D5', 4, 'base_c2_measured.csv', 'base titration Cycle 2')
    s.py('EDP_tecan_gp_master.py', 'Base 2')
    s.imp('base_on', 'base_status.csv')
    s.label('Base Cycle 2 entered')

    # ================================================================== validations (+ VI Tris Cycle 1)
    s.needs(['needs_val.csv', 'needs_tris_c1.csv', 'needs_base_val.csv'], 'validations')
    s.handover('val', ('NEXT - the robot makes the acid validation (Acid titration plate, column 4), the VI neutralisation '
               'Cycle 1 (Base titration plate, columns 1 and 2) and the base validation (Base titration plate, column 6).',
               'It uses about ~up_e~ mL eluate (Product 2W LEFT well), ~up_a~ mL acid (LEFT well) and ~up_b~ mL Tris '
               '(RIGHT well). Top up now if needed, then click OK.'), 'PLATE CHECK 3')
    s.group('Validations - acid (col 4), VI Tris Cycle 1 (Base titration plate cols 1-2), base (col 6)')
    s.column('plan_val.csv', 0, 'ACID', 4, 'AcidVal')
    s.skip_unless('tris_on', 'Validation - VI Tris Cycle 1 done')
    s.column('plan_tris_c1.csv', 0, 'TRIS', 1, 'TrisC1a')
    s.column('plan_tris_c1.csv', 1, 'TRIS', 2, 'TrisC1b')
    s.label('Validation - VI Tris Cycle 1 done')
    s.skip_unless('base_on', 'Validation - base done')
    s.column('plan_base_val.csv', 0, 'TRIS', 6, 'BaseVal')
    s.label('Validation - base done')
    s.group_end()
    s.timer('end', 'val')
    s.comment('DOOR ACCESS POINT 3 - acid validation, VI neutralisation Cycle 1, base validation')
    s.prompt('PLATE CHECK 3 - validations', 'Open the safety door and take out the plates.',
             'Acid titration plate: measure column 4 (A4, B4, C4, D4).',
             'Base titration plate: measure columns 1 and 2 (VI neutralisation), then column 6 (base validation) - each entry box names its wells.',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_val', 'Acid validation - Acid titration plate', 'A4 B4 C4 D4', 4, 'validation_results.csv', 'acid validation')
    s.py('EDP_tecan_gp_master.py', 'Acid 3')
    s.py('EDP_archive_run.py', 'Acid')
    s.skip_unless('tris_on', 'VI Tris Cycle 1 entered')
    s.entry('rd_tris_c1', 'VI Cycle 1 - Base titration plate', 'A1 B1 C1 D1 A2 B2 C2 D2', 8, 'tris_c1_measured.csv', 'VI neutralisation Cycle 1')
    s.py('EDP_tecan_gp_master.py', 'Tris 1')
    s.label('VI Tris Cycle 1 entered')
    s.skip_unless('base_on', 'Base validation entered')
    s.entry('rd_base_val', 'Base validation - Base titration plate', 'A6 B6 C6 D6', 4, 'base_val_measured.csv', 'base validation')
    s.py('EDP_tecan_gp_master.py', 'Base 3')
    s.py('EDP_archive_run.py', 'Base EL')
    s.label('Base validation entered')
    s.jump_if('tris_on', '0', 'Ratios stored')
    s.needs(['needs_tris_val.csv'], 'VI neutralisation validation')
    s.handover('tris_val', ('NEXT - the robot makes the VI neutralisation validation (Base titration plate, column 3) at the process pool.',
               'It uses about ~up_e~ mL eluate, ~up_a~ mL acid and ~up_b~ mL Tris. Top up if needed, then click OK.'), 'PLATE CHECK 4')
    s.group('VI neutralisation - validation at the process pool (Base titration plate, column 3)')
    s.column('plan_tris_val.csv', 0, 'TRIS', 3, 'TrisVal')
    s.group_end()
    s.timer('end', 'tris_val')
    s.comment('DOOR ACCESS POINT 4 - VI neutralisation validation')
    s.prompt('PLATE CHECK 4 - VI neutralisation validation', 'Open the safety door and take out the Base titration plate.',
             'Measure column 3 (A3, B3, C3, D3).', 'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_tris_val', 'VI neutralisation validation - Base titration plate', 'A3 B3 C3 D3', 4, 'tris_val_measured.csv', 'VI neutralisation validation')
    s.py('EDP_tecan_gp_master.py', 'Tris 2')
    s.py('EDP_archive_run.py', 'Base TR')

    # ================================================================== ratios stored -> ranges -> custom pool request
    s.label('Ratios stored')
    s.py('EDP_tecan_gp_master.py', 'Pools range')
    s.py('EDP_archive_run.py', 'Map')
    for k, v in enumerate(('sm_va', 'sm_a1', 'sm_a2', 'sm_a3', 'sm_a4', 'sm_vb', 'sm_b1', 'sm_b2', 'sm_b3', 'sm_b4',
                           'sm_vc', 'sm_c1', 'sm_c2', 'sm_c3', 'sm_c4', 'sm_n1', 'sm_n2', 'sm_n3', 'sm_n4', 'sm_ph0')):
        s.imp(v, 'run_summary.csv', k + 1)
    for k, v in enumerate(('rg_lo', 'rg_hi', 'rg_dlo', 'rg_dhi', 'rg_blo', 'rg_bhi', 'rg_bdlo', 'rg_bdhi')):
        s.imp(v, 'pool_range.csv', k + 1)
    s.prompt('RATIOS STORED (1 of 2) - archived with the full ratio map: mL to add from any start pH (titration_data.xlsx, RatioMap).',
             'Acid, mL 1 M acetic acid per mL eluate (eluate pH ~sm_ph0~): pH ~Target_pH1~: ~sm_a1~ | ~Target_pH2~: ~sm_a2~ | '
             '~Target_pH3~: ~sm_a3~ | process ~Process_Target_pH~: ~sm_a4~. Worst validation error ~sm_va~.',
             'Base, mL 1 M Tris per mL eluate: pH ~nb_t1~: ~sm_c1~ | ~nb_t2~: ~sm_c2~ | ~nb_t3~: ~sm_c3~ (0 = at or below the eluate pH). Worst validation error ~sm_vc~.')
    s.prompt('RATIOS STORED (2 of 2) - neutralisation (VI method), mL 1 M Tris per mL of acidified product to reach pH ~nb_proc~:',
             'from pH ~Target_pH1~: ~sm_n1~ | from ~Target_pH2~: ~sm_n2~ | from ~Target_pH3~: ~sm_n3~ | from process ~Process_Target_pH~: ~sm_n4~.',
             'From the process pool to pH ~nb_t1~ / ~nb_t3~: ~sm_b1~ / ~sm_b3~. Worst validation error ~sm_vb~ (-1 = not measured).')
    s.prompt('CUSTOM POOLS (optional) - next you can fill lanes of the Bulk Pools plate (grid 22 position 2), about 25 mL each, at any pH in:',
             '- with acid: pH ~rg_lo~ to ~rg_hi~ (the eluate pH) - most accurate from ~rg_dlo~ to ~rg_dhi~;',
             '- with Tris: pH ~rg_blo~ to ~rg_bhi~ - most accurate from ~rg_bdlo~ to ~rg_bdhi~ (0 to 0 = the base titration did not finish).',
             'Each pool is measured and trimmed until within 0.03. Click OK, then type the lanes and the pH values (0 for none).')
    s.var('bad_ln', '0', 'lane of a refused pool pH')
    s.var('bad_ph', '0', 'refused pool pH')
    s.loop('3', 'Pools_ask')
    s.if_then('pool_ok', '0')
    s.ask('pool_lanes', '0', 'CUSTOM POOLS in the Bulk Pools plate: type the lanes to fill, 1 to 6, separated by commas - for '
          'example 1,3,5 for the old bulk pools. Type 0 for no pools.', text=True)
    s.ask('pool_phs', '0', 'CUSTOM POOLS: type one target pH per lane, in the same order, separated by commas - for example '
          '3.2,3.4,3.6 - within the ranges shown before. Every pool is measured and trimmed.', text=True)
    s.export(['pool_lanes'], 'pool_lanes.csv')
    s.export(['pool_phs'], 'pool_phs.csv')
    s.py('EDP_tecan_gp_master.py', 'Pools 0')
    for k, v in enumerate(('pool_ok', 'n_pools', 'bad_ln', 'bad_ph', 'n_mk1', 'n_m1')):
        s.imp(v, 'pool_status.csv', k + 1)
    s.if_then('pool_ok', '0')
    s.prompt('CUSTOM POOLS - the request could not be used. Lanes are 1 to 6, each once, with one pH each; every pH must lie '
             'in the ranges: pH ~rg_lo~ to ~rg_hi~ with acid, ~rg_blo~ to ~rg_bhi~ with Tris.',
             'Outside the ranges: lane ~bad_ln~, pH ~bad_ph~ (0 = the lists did not match). Click OK and type them again; after 3 tries the run ends without pools.')
    s.els(); s.comment('pool request accepted'); s.end_if()
    s.els(); s.comment('pool request already accepted'); s.end_if()
    s.end_loop()
    s.jump_if('pool_ok', '0', 'Run finished')
    s.jump_if('n_pools', '0', 'Run finished')

    # ================================================================== custom pools: make, measure, trim
    s.needs(['needs_pools.csv'], 'custom pools')
    s.pool_order()
    s.handover('pools', ('NEXT - the robot makes ~n_mk1~ custom pool(s) in the Bulk Pools plate (grid 22 position 2): eluate, '
               'then acid - or, above the pH of the eluate itself, Tris.',
               'It uses about ~up_e~ mL eluate (Product 2W LEFT well), ~up_a~ mL acid (LEFT well) and ~up_b~ mL Tris '
               '(RIGHT well). Place the empty Bulk Pools plate, top up if needed, then click OK.'), 'POOL CHECK 1')
    s.group('Custom pools - Bulk Pools plate (8 tips)')
    for lane in LANES:
        s.lane_block('pool_make1.csv', lane, 'make', 'mk')
    s.group_end()
    s.timer('end', 'pools')
    s.comment('DOOR ACCESS POINT 5 - custom pools, check 1')
    s.prompt('POOL CHECK 1 of up to 3 - custom pools', 'Open the safety door and take out the Bulk Pools plate.',
             'Stir each pool gently before reading (Tris is denser than the eluate). Measure lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_pools1', 'Pools in the Bulk Pools plate', '', 'pools', 'pools_measured_1.csv', 'pool')
    s.py('EDP_tecan_gp_master.py', 'Pools 1')
    for k, v in enumerate(('r_work', 'n_m2', 'n_mk2')):
        s.imp(v, 'round2.csv', k + 1)
    s.jump_if('r_work', '0', 'Run finished')
    s.needs(['needs_round2.csv'], 'pool trims and deferred pools')
    s.pool_order()
    s.handover('round2', ('NEXT - the robot trims the pools that missed their target and makes ~n_mk2~ pool(s) that did not '
               'fit in the troughs before.',
               'It uses about ~up_e~ mL eluate, ~up_a~ mL acid and ~up_b~ mL Tris. Top up if needed, then click OK.'), 'POOL CHECK 2')
    s.group('Custom pools - trims and deferred pools (8 tips)')
    for lane in LANES:
        s.lane_block('pool_trim1.csv', lane, 'trim', 'tr')
    for lane in LANES:
        s.lane_block('pool_make2.csv', lane, 'make', 'mk')
    s.group_end()
    s.timer('end', 'round2')
    s.comment('DOOR ACCESS POINT 6 - custom pools, check 2')
    s.prompt('POOL CHECK 2 of up to 3 - custom pools', 'Open the safety door and take out the Bulk Pools plate.',
             'Stir each pool gently before reading (Tris is denser than the eluate). Measure lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_pools2', 'Pools in the Bulk Pools plate', '', 'pools', 'pools_measured_2.csv', 'pool')
    s.py('EDP_tecan_gp_master.py', 'Pools 2')
    for k, v in enumerate(('r_work', 'n_m3')):
        s.imp(v, 'round3.csv', k + 1)
    s.jump_if('r_work', '0', 'Run finished')
    s.needs(['needs_round3.csv'], 'second pool trims')
    s.pool_order()
    s.handover('round3', ('NEXT - the robot trims the pools that still missed their target.',
               'It uses about ~up_e~ mL eluate, ~up_a~ mL acid and ~up_b~ mL Tris. Top up if needed, then click OK.'), 'POOL CHECK 3')
    s.group('Custom pools - second trims (8 tips)')
    for lane in LANES:
        s.lane_block('pool_trim2.csv', lane, 'trim', 'tr')
    s.group_end()
    s.timer('end', 'round3')
    s.comment('DOOR ACCESS POINT 7 - custom pools, check 3')
    s.prompt('POOL CHECK 3 of 3 - custom pools', 'Open the safety door and take out the Bulk Pools plate.',
             'Stir each pool gently before reading (Tris is denser than the eluate). Measure lanes ~po1~ ~po2~ ~po3~ ~po4~ ~po5~ ~po6~, in that order (0 = no more).',
             'Return the plate, close the safety door, then click OK and type the readings.')
    s.entry('rd_pools3', 'Pools in the Bulk Pools plate', '', 'pools', 'pools_measured_3.csv', 'pool')
    s.py('EDP_tecan_gp_master.py', 'Pools 3')

    # ================================================================== archive and results
    s.label('Run finished')
    s.py('EDP_archive_run.py', 'Pools')
    s.py('EDP_archive_run.py', 'Session')
    s.py('EDP_make_workbook.py')
    s.py('EDP_tecan_gp_master.py', 'Summary')
    s.imp('n_any', 'pools_any.csv')
    s.if_then('n_any', '0')
    s.prompt('Run finished. All results are in titration_data.xlsx (ML_data_archive); details in python_log.txt (tecan_data).')
    s.els(); s.comment('custom pools made'); s.end_if()
    s.jump_if('n_any', '0', 'Run done')
    for lane in LANES:
        s.imp(f'pf{lane}', 'pool_final.csv', 2 * lane - 1)
        s.imp(f'pt{lane}', 'pool_final.csv', 2 * lane)
    s.prompt('Run finished. CUSTOM POOLS in the Bulk Pools plate - last measured pH (target); 0 = no pool, -1 = not made:',
             'lane 1: ~pf1~ (~pt1~) | lane 2: ~pf2~ (~pt2~) | lane 3: ~pf3~ (~pt3~) | lane 4: ~pf4~ (~pt4~) | lane 5: ~pf5~ '
             '(~pt5~) | lane 6: ~pf6~ (~pt6~).',
             "Volumes, ratios, every reading and each low-pH pool's Tris ratios: titration_data.xlsx, Pools sheet.")
    s.label('Run done')
    s.comment('End of script')
    missing = s.jumps - s.labels
    assert not missing, missing
    return s.L


if __name__ == '__main__':
    import sys
    body = build()
    crc = assemble(body, sys.argv[1], sys.argv[2])
    print(len(body), 'lines, checksum', crc)
