"""Generates EDP_Titration_python_ML_hybrid.esc, v16 (22 October 2026).

v16 against v15 (Edward, 22 Oct, after the 21 Oct v15 run):
  - tips 5-8 work too. They were washed with tips 1-4 at every wash but never pipetted on the plates since v7. Now,
    whenever a liquid needs two or more portions per well, all 8 tips aspirate together and tips 5-8 put a second
    portion into the same four wells, so each trip to a trough carries twice as much (v6 did this on 30 Sep). An odd
    last portion goes with tips 1-4 alone; one-portion liquids (the small acid and Tris volumes) are unchanged, so no
    dispense is split into smaller ones. Use_8_tips (Run set-up) = 0 puts everything back on tips 1-4.
  - the Acid titration plate moves from grid 28 position 1 to position 3. At the 24-well row spacing tips 5-8 sit
    77 mm in front of tips 1-4, so they cannot reach the rearmost position; positions 2 and 3 are where tips 5-8 ran on
    30 Sep (the tip-5 crash that day was at grid 38, a carrier no longer used). Position 1 stays empty. The worktable
    labels follow (position 1 no label, position 3 "Acid titration"); nothing else in the header changes.
  - order of the two dispenses of a trip: at the back plate of the two (Base titration, position 2) tips 5-8 go first,
    so tips 1-4 still holding liquid wait over the empty position 1; at the front plate (Acid titration, position 3)
    tips 1-4 go first and tips 5-8 wait over the door row. Full tips never wait over a plate.
  - the 5 x 900 uL mix of a column: tips 5-8 at the Base titration plate, tips 1-4 at the Acid titration plate - the
    idle tips, wet from the column, wait over position 1 or the door row, never over the other plate (v6 mixed with
    tips 5-8 at grid 28 on 30 Sep, Edward's test script of 29 Sep at position 2)
  - wells stay at 3 mL: with two portions per trip any 3 mL well takes at most two trips per liquid (2.7 mL would too)
  - the run's results workbook (EDP_make_workbook.py), named after the batch and the date, is written for export.

Everything else as v15: flow, plate checks, the product of a block first, washes at every change of liquid."""
import datetime
import zlib
import gen15
from gen11 import LANES, WASH_FIRST, PRODUCT, TROUGH, LPH, ROWS8, DSP, ROWS4, ASP, trough_wells, plate_wells
from gen10 import lane_wells

PLATES = {'ACID': (28, 2, 'Acid titration'), 'TRIS': (28, 1, 'Base titration')}   # v16: acid plate in position 3
BACK_SITE = 1                                  # the rear plate of the two: tips 5-8 dispense first there
LABELS_V15 = '998;Acid titration;Base titration;Pool 3 titration;;'
LABELS_V16 = '998;;Base titration;Acid titration;;'


class Script(gen15.Script):
    def trips(self, n, vols, src, src_col, pg, ps, dest, tag):
        """n portions of `vols` (uL per tip, tips 1-4 -> rows A-D of `dest`). v16: with Use_8_tips = 1 two portions go
        in one trip - tips 1-8 aspirate from rows A-H of the trough column, then tips 5-8 and tips 1-4 each put their
        portion into the same four wells (back plate: tips 5-8 first; front plate: tips 1-4 first, so tips still holding
        liquid never wait over a plate). An odd last portion goes with tips 1-4 alone; with Use_8_tips = 0, all of them."""
        self.var('d8', f'({n}-({n}%2))/2*Use_8_tips', 'trips with tips 1-8: two portions per well')
        self.var('d1', f'{n}-2*d8', 'trips with tips 1-4 alone')
        hi = [0, 0, 0, 0] + list(vols)
        order = ((240, hi), (15, vols)) if ps == BACK_SITE else ((15, vols), (240, hi))
        self.if_then('d8', '0', op=2)
        self.loop('d8', f'{tag}8')
        self.pip('Aspirate', ASP, 255, list(vols) * 2, src[0], src[1], trough_wells(src_col, ROWS8))
        for mask, v in order:
            self.pip('Dispense', DSP, mask, v, pg, ps, dest)
        self.end_loop()
        self.els(); self.comment('no trips with tips 1-8'); self.end_if()
        self.if_then('d1', '0', op=2)
        self.loop('d1', f'{tag}4')
        self.pip('Aspirate', ASP, 15, vols, src[0], src[1], trough_wells(src_col, ROWS4))
        self.pip('Dispense', DSP, 15, vols, pg, ps, dest)
        self.end_loop()
        self.els(); self.comment('no trip with tips 1-4 alone'); self.end_if()

    def product_part(self, plan, group, plate, col, what):
        """The eluate portions of one plate column, from its 15 plan numbers (portions, uL per tip x 4 for eluate, then
        acid, then Tris); two portions per trip with tips 1-8 (v16)."""
        pg, ps, _ = PLATES[plate]
        vols = ['e1', 'e2', 'e3', 'e4']
        for k, v in enumerate(['ne'] + vols):
            self.imp(v, plan, 15 * group + k + 1)
        self.if_then('ne', '0', op=2)
        self.trips('ne', vols, PRODUCT, 1, pg, ps, plate_wells(col), f'{what}_eluate')
        self.els(); self.comment('no eluate here'); self.end_if()

    def reagent_part(self, plan, group, plate, col, what):
        """The acid and Tris portions of one plate column (a wash between them if both; two portions per trip with tips
        1-8), its 5 x 900 uL mix, one wash."""
        pg, ps, _ = PLATES[plate]
        dest = plate_wells(col)
        acid, tris = ['a1', 'a2', 'a3', 'a4'], ['b1', 'b2', 'b3', 'b4']
        for k, v in enumerate(['na'] + acid + ['nb'] + tris):
            self.imp(v, plan, 15 * group + k + 6)
        self.if_then('na', '0', op=2)
        self.trips('na', acid, TROUGH, 1, pg, ps, dest, f'{what}_acid')
        self.els(); self.comment('no acid here'); self.end_if()
        self.if_then('nb', '0', op=2)
        self.if_then('na', '0', op=2); self.wash(); self.els(); self.comment('no acid before the Tris'); self.end_if()
        self.trips('nb', tris, TROUGH, 2, pg, ps, dest, f'{what}_Tris')
        self.els(); self.comment('no Tris here'); self.end_if()
        # v16: the mix by the tips whose idle partners wait over an empty position - tips 5-8 at the back plate (tips
        # 1-4 over position 1), tips 1-4 at the front plate (tips 5-8 over the door row): wet tips never wait over a plate
        if ps == BACK_SITE:
            self.pip('Mix', DSP, 240, ['0'] * 4 + ['900'] * 4, pg, ps, dest, extra='5,0,0')
        else:
            self.pip('Mix', DSP, 15, ['900'] * 4, pg, ps, dest, extra='5,0,0')
        self.wash()


def build():
    s = Script()
    s.comment('EDP pH adjustment method v16')
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
    s.prompt('The Python files in tecan_data are not all v16, so the run will stop before anything is pipetted.',
             'First file found older: ~vbad~ (1 EDP_tecan_gp_master.py, 2 EDP_parse_readings.py, 3 EDP_archive_run.py, '
             '4 EDP_tecan_reset.py, 5 EDP_make_workbook.py, 6 EDP_time_estimate.py; 0 = the analysis file itself is older).',
             'Copy all six EDP Python files from the python folder of the zip into tecan_data, then start the run again.')
    s.els(); s.comment('all six EDP files are v16'); s.end_if()
    s.jump_if('vchk', '0', 'Run done')

    # ================================================================== set-up
    s.group('Run set-up')
    s.var('one', '1', 'constant')
    s.var('Well_Volume', '3', 'liquid (mL) in every titration well')
    s.var('Product_Vol', '1', 'ratios are per 1 mL of eluate or pool')
    s.var('Use_8_tips', '1', '1 = tips 5-8 take a second portion into the same wells on every trip of two or more '
          '(plates at grid 28 positions 2 and 3); 0 = tips 1-4 only, as v15')
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
    s.ask('neut', '1', 'VI NEUTRALISATION TITRATION - allows the VI method: Tris added to the acidified eluate, storing the Tris '
          'that brings the product back to neutral from each acid pH (12 wells, Base titration plate columns 1-3; about 10 min '
          'more robot time and a 4th plate check). 1 = yes, 0 = no.', 0, 1)
    s.ask('base_on', '1', 'BASE TITRATION OF THE ELUATE - allows custom pools above the eluate pH, e.g. pH 7.4 straight from '
          'the eluate: Tris added straight to the eluate (12 wells, Base titration plate columns 4-6; about 8 min more robot time). 1 = yes, '
          '0 = no.', 0, 1)
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
    s.export(['neut', 'base_on', 'Use_8_tips'], 'options.csv')
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
               '- Empty 24-well plates, grid 28: Acid titration in position 3, Base titration in position 2 if used; position 1 empty.',
               '- Empty Bulk Pools plate (6-well, grid 22 position 2) for custom pools. Nothing in the door row.'), 'PLATE CHECK 1')
    s.group('Cycle 1 - acid (Acid titration plate, cols 1 and 5), base (Base titration plate, col 4)')
    s.plate_block('Cycle 1', [('plan_c1.csv', 0, 'ACID', 1, 'AcidC1a', None), ('plan_c1.csv', 1, 'ACID', 5, 'AcidC1b', None),
                              ('plan_base_c1.csv', 0, 'TRIS', 4, 'BaseC1', 'base_on')])
    s.group_end()
    s.timer('end', 'acid_c1')
    s.comment('DOOR ACCESS POINT 1 - acid Cycle 1 (A1 product only) and base Cycle 1')
    s.if_then('base_on', '1')
    s.prompt('PLATE CHECK 1 - acid and base titration Cycle 1',
             'Open the safety door. Acid titration plate: measure column 1 (A1 is product only - the eluate pH), then column 5 (A5-D5).',
             'Base titration plate: measure column 4 (A4-D4).',
             'Return the plates, close the safety door, click OK, then type A1-D1 and A5-D5, then A4-D4.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.els()
    s.prompt('PLATE CHECK 1 - acid titration Cycle 1',
             'Open the safety door, take out the Acid titration plate and measure column 1 (A1 is product only - the eluate pH), then column 5.',
             'Return the plate, close the safety door, then click OK and type A1, B1, C1, D1, then A5, B5, C5, D5.',
             'After the readings the robot makes Cycle 2 in about ~t_next~ min. Please stay nearby - PLATE CHECK 2 follows.')
    s.end_if()
    s.entry('rd_acid_c1', 'Acid titration Cycle 1', 'A1 B1 C1 D1 A5 B5 C5 D5', 8, 'acid_c1_measured.csv', 'acid titration Cycle 1')
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
    s.plate_block('Cycle 2', [('plan_c2.csv', 0, 'ACID', 2, 'AcidC2a', None), ('plan_c2.csv', 1, 'ACID', 3, 'AcidC2b', None),
                              ('plan_base_c2.csv', 0, 'TRIS', 5, 'BaseC2', 'base_on')])
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
    s.plate_block('Validations', [('plan_val.csv', 0, 'ACID', 4, 'AcidVal', None),
                                  ('plan_tris_c1.csv', 0, 'TRIS', 1, 'TrisC1a', 'tris_on'), ('plan_tris_c1.csv', 1, 'TRIS', 2, 'TrisC1b', 'tris_on'),
                                  ('plan_base_val.csv', 0, 'TRIS', 6, 'BaseVal', 'base_on')])
    s.group_end()
    s.timer('end', 'val')
    s.comment('DOOR ACCESS POINT 3 - acid validation, VI neutralisation Cycle 1, base validation')
    s.prompt('PLATE CHECK 3 - validations', 'Open the safety door and take out the plates.',
             'Acid titration plate: measure column 4 (A4-D4), then column 2 again (A2-D2) - a drift check against plate check 2.',
             'Base titration plate, if used: columns 1 and 2 (VI neutralisation), then column 6 (base validation) - each entry box names its wells.',
             'Return the plates, close the safety door, then click OK and type the readings.')
    s.entry('rd_acid_val', 'Acid validation', 'A4 B4 C4 D4 then A2 B2 C2 D2 again', 8, 'validation_results.csv', 'acid validation')
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
    s.plate_block('VI validation', [('plan_tris_val.csv', 0, 'TRIS', 3, 'TrisVal', None)])
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
                           'sm_vc', 'sm_c1', 'sm_c2', 'sm_c3', 'sm_c4', 'sm_n1', 'sm_n2', 'sm_n3', 'sm_n4', 'sm_ph0', 'sm_dr')):
        s.imp(v, 'run_summary.csv', k + 1)
    for k, v in enumerate(('rg_lo', 'rg_hi', 'rg_dlo', 'rg_dhi', 'rg_blo', 'rg_bhi', 'rg_bdlo', 'rg_bdhi')):
        s.imp(v, 'pool_range.csv', k + 1)
    s.prompt('RATIOS STORED (1 of 2) - archived with the full ratio map: mL to add from any start pH (results workbook, RatioMap sheet).',
             'Acid, mL 1 M acetic acid per mL eluate (eluate pH ~sm_ph0~): pH ~Target_pH1~: ~sm_a1~ | ~Target_pH2~: ~sm_a2~ | '
             '~Target_pH3~: ~sm_a3~ | process ~Process_Target_pH~: ~sm_a4~. Worst validation error ~sm_va~ (probe drift taken off: ~sm_dr~).',
             'Base, mL 1 M Tris per mL eluate: pH ~nb_t1~: ~sm_c1~ | ~nb_t2~: ~sm_c2~ | ~nb_t3~: ~sm_c3~ (0 = at or below the eluate pH, -1 = no base titration). Worst validation error ~sm_vc~.')
    s.if_then('tris_on', '1')
    s.prompt('RATIOS STORED (2 of 2) - neutralisation (VI method), mL 1 M Tris per mL of acidified product to reach pH ~nb_proc~:',
             'from pH ~Target_pH1~: ~sm_n1~ | from ~Target_pH2~: ~sm_n2~ | from ~Target_pH3~: ~sm_n3~ | from process ~Process_Target_pH~: ~sm_n4~.',
             'From the process pool to pH ~nb_t1~ / ~nb_t3~: ~sm_b1~ / ~sm_b3~. Worst validation error ~sm_vb~ (-1 = not measured).')
    s.els(); s.comment('no VI neutralisation titration this run'); s.end_if()
    s.prompt('CUSTOM POOLS (optional) - next you can fill lanes of the Bulk Pools plate (grid 22 position 2), about 25 mL each, at any pH in:',
             '- with acid: pH ~rg_lo~ to ~rg_hi~ (the eluate pH) - most accurate from ~rg_dlo~ to ~rg_dhi~;',
             '- with Tris: pH ~rg_blo~ to ~rg_bhi~ - most accurate from ~rg_bdlo~ to ~rg_bdhi~ (0 to 0 = no base titration this run).',
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
    s.make_block('pool_make1.csv', 'mk', 'n_mk1')
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
    s.make_block('pool_make2.csv', 'mk', 'n_mk2')
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
    s.py('EDP_tecan_gp_master.py', 'Summary')
    s.py('EDP_make_workbook.py')                  # v16: after Summary, so the run's workbook has the pools' last readings
    s.imp('n_any', 'pools_any.csv')
    s.if_then('n_any', '0')
    s.prompt('Run finished. This run is in its own workbook, EDP_results folder in Public Documents (file named after the batch '
             'and the date); the whole archive is in titration_data.xlsx (ML_data_archive); details in python_log.txt (tecan_data).')
    s.els(); s.comment('custom pools made'); s.end_if()
    s.jump_if('n_any', '0', 'Run done')
    for lane in LANES:
        s.imp(f'pf{lane}', 'pool_final.csv', 2 * lane - 1)
        s.imp(f'pt{lane}', 'pool_final.csv', 2 * lane)
    s.prompt('Run finished. CUSTOM POOLS in the Bulk Pools plate - last measured pH (target); 0 = no pool, -1 = not made:',
             'lane 1: ~pf1~ (~pt1~) | lane 2: ~pf2~ (~pt2~) | lane 3: ~pf3~ (~pt3~) | lane 4: ~pf4~ (~pt4~) | lane 5: ~pf5~ '
             '(~pt5~) | lane 6: ~pf6~ (~pt6~).',
             "Volumes, readings and trims: this run's workbook (EDP_results folder in Public Documents, file named after the batch and the date), Pools sheet.")
    s.label('Run done')
    s.comment('End of script')
    missing = s.jumps - s.labels
    assert not missing, missing
    return s.L


def assemble(body, header_src, out):
    """Worktable header copied verbatim from Edward's EVOware-saved script of 12 Oct, a new date line, and the grid 28
    labels moved with the Acid titration plate (position 1 no label, position 3 "Acid titration"); checksum."""
    raw = open(header_src, 'rb').read().decode('latin-1').split('\r\n')
    head = raw[1:raw.index('--{ RPG }--') + 1]
    head[0] = datetime.datetime.now().strftime('%Y%m%d_%H%M%S') + head[0][15:]
    assert head.count(LABELS_V15) == 1, 'grid 28 label line not found in the header source'
    head[head.index(LABELS_V15)] = LABELS_V16
    rest = '\r\n'.join(head + body + ['']).encode('latin-1')
    crc = format((~zlib.crc32(rest)) & 0xffffffff, '08X')
    open(out, 'wb').write(crc.encode('latin-1') + b'\r\n' + rest)
    return crc


if __name__ == '__main__':
    import sys
    body = build()
    crc = assemble(body, sys.argv[1], sys.argv[2])
    print(len(body), 'lines, checksum', crc)
