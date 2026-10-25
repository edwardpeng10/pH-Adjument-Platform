import os
import re
import sys
import json
import math
from datetime import datetime

if sys.stdout is None: sys.stdout = open(os.devnull, "w")      # started with pythonw: there is no console
if sys.stderr is None: sys.stderr = open(os.devnull, "w")

import numpy as np
import pandas as pd

# =====================================================================================
# EDP pH adjustment platform - analysis, v17 (22 October 2026; v14 notes below).
#
# v17 (Edward, 22 Oct, after the v16 run that day: "base validation very inaccurate ... the whole point of this stage is
# to collect accurate data for machine learning"): the base wells are made with a 0.1 M Tris wherever one portion of it
# is enough. That run's eluate (pH 4.55) needed only 47-65 uL of 1 M Tris per 3 mL base well near pH 7, where 1 uL moved
# the pH 0.05-0.09; the chemistry model fitted to all 44 readings puts those wells 0.04 to 5.5 uL of Tris off their
# recipes (about 2 uL typical; 21 Oct about 1.3 uL), while the acid and VI wells (300-1500 uL dispenses) were within
# 0.02 pH. The same few uL in 10x the volume of a 0.1 M Tris are 10x less Tris: simulated base validation wells on that
# eluate 0.08 -> 0.035 pH RMS (weakest eluates 0.26-0.29 -> 0.05-0.06). Ratios stay in mL 1 M Tris (mmol) per mL
# eluate; every well records its Tris stock (Tris_Stock_M) and the chemistry model includes the water that comes with
# the dilute Tris. ABS_UL (the uL scatter of a dispense that does not shrink with its volume) is now 1.5 uL, from those
# two runs (0.5 before). Dilute_Tris_M = 0 in the script's Run set-up turns it off (all 1 M Tris, as v16).
#
# EDP pH adjustment platform - analysis, v14 (12 October 2026).
#
# The script's first job is to measure and STORE the product's acid, base and neutralisation ratios, so that its pH
# can be adjusted later from any pH it is at (Edward, 12 Oct): acid from the eluate, Tris straight into the eluate, and
# Tris into the acidified product (the VI method, its own script). All three titrations run in every run; custom
# pools at the end are optional. v14 stores a RATIO MAP per run (ratio_map.csv in the archive, RatioMap sheet):
#   acid path    mL acid per mL eluate to reach each pH below the eluate's own
#   base path    mL Tris per mL eluate to reach each pH above it
#   neutralise   mL Tris per mL of acidified product, from every pH on the acid path to every pH 6.8-8.5
# From pH A to pH B on the same path: (ratio(B) - ratio(A)) / (1 + ratio(A)) per mL of the product at A; the command
# 'Ratio <from pH> <to pH> [product]' does this for another script and writes ratio_lookup.csv.
#
#   EDP_tecan_gp_master.py Check          all six EDP files the same version? (check_status.csv; the run stops if not)
#                          Acid 0|1|2|3   acid titration: Cycle 1 plan | Cycle 2 plan | ratios, validation and VI Tris
#                                         titration plans | refit with the validation wells
#                          Base 0|1|2|3   base titration, Tris straight into the eluate: Cycle 1 plan | Cycle 2 plan |
#                                         fit and validation plan | refit with the validation wells
#                          Tris 1|2       VI neutralisation data, Tris into the acidified eluate (recipe wells): fit and
#                                         validation plan | refit
#                          Pools range    after every ratio is found: the ratio map, the chemistry cross-check, the
#                                         ratios and the pH ranges custom pools can be made in, for the messages
#                          Pools 0        read the pool request (lanes, target pHs), check it against the ranges, plan
#                          Pools 1|2|3    pool readings of round 1, 2, 3 -> trims (and deferred pools) for the next round
#                          Ratio a b [p]  mL titrant per mL of product at pH a to reach pH b (ratio_lookup.csv)
#                          Summary        numbers for the final message
#   v14 (12 Oct run, v12 files on the robot except the helpers): the ratio map; the neutralisation line anchored at
#   zero acid by the base titration (so neutralising from any acid-path pH is interpolation, not extrapolation); a
#   chemistry-model fit to all three titrations as an independent cross-check and a product fingerprint for the ML;
#   a file-version check; a lane check on the first pool reading (12 Oct: a lane pool read 7.51 where its recipe
#   gives 7.0 by every model - the validation well with MORE Tris read 7.26 - so the lane, not the ratio, was off).
#   The RATIOS STORED messages show the map's own entries. Later trims of a Tris pool follow the run's curve from the
#   pool's last reading (v13's secant overshot near the equivalence point). Without a usable Tris-equation fit, the
#   VI validation wells are planned from the chemistry model when it reproduces the run's readings.
#   v13: pools above the eluate's own pH are eluate + Tris only, planned from a base titration of the eluate itself.
#   v10-v12 made them as eluate taken to the process pH with acid, then Tris - the only base data those versions
#   measured was the VI neutralisation (Tris into the acidified eluate), so that was the only composition they could
#   predict. For a pool that just has to be at pH 7.4 that is the wrong product: on the 30 Sep eluate it held 24% eluate
#   and ~0.4 M acetate + ~0.4 M Tris, against 99% eluate and 8 mM Tris for Tris alone, and the protein sat at pH 3.5
#   on the way. Acidify-hold-neutralise is the VI method's job (its own script); this script measures the data it needs
#   (acid ratios, and the Tris that neutralises the acidified product) and keeps the two Tris experiments apart.
#   v12: the pool request comes after all calculations - the operator sees the stored ratios and the ranges first.
#   Acid Cycle 2 is therefore placed around the titration targets only; a pool outside that window starts from
#   fewer readings and its trims do more of the work (the message shows where the titration measured densely).
#
# What changed from v10 and why (7 Oct run: the 3.2 target missed by +0.16 and Pool 1 came out at 3.41):
#   * The eluate's own pH is typed with the Cycle-1 readings, so no well is spent on ratio 0.
#   * Acid Cycle 1 is 4 wells log-spaced from Max/50 to Max (0.08-4): only to LOCATE the window around the targets.
#     v10's 8 doubling wells put 4 wells at pH 2.6-2.9 for that product, far beyond every target, and left
#     0.125-0.25 (where three targets were) without an intermediate point.
#   * Acid Cycle 2 is 8 wells (two columns) log-spaced across the located window, widened 0.8x/1.25x: dense data
#     exactly where the ratios are read.
#   * Ratios come from a robust fit (Huber-weighted polynomial in log ratio over the readings around the targets),
#     not from a curve drawn through every reading. A reading more than 0.10 pH (and 4 noise SD) off the fit is left
#     out and logged (at most 2). On 7 Oct one Cycle-1 reading (3.08 at ratio 0.25, ~0.25 below its neighbours)
#     decided the 3.2 ratio; a fit averages the reading noise instead of copying it.
#   * Pools: up to 6, in any lanes of the 6-well plate (grid 22 position 2), each at its own target pH. Up to the
#     eluate's own pH a pool is eluate + acid; from 6.8 to 8.5 it is eluate taken to the process pH with acid, then
#     Tris. Every pool is measured and trimmed (at most twice) until within POOL_TOL. The old bulk pools are lanes
#     1, 3, 5 at 3.2, 3.4, 3.6.
#   * Tris titration (optional; automatic for pools above 6.8) on RECIPE wells - eluate + acid + Tris dispensed from
#     the troughs, never drawn back out of a lane - at two acid fractions. The Tris needed per mL of low-pH pool at
#     equivalence is linear in the pool's acid fraction f = Rp/(1+Rp) (te = X + (1 - X) f, X = the eluate's own net
#     acid demand), so two fractions give it for every pool. 12 wells instead of v10's 36, the bulk-pool gap (real
#     pools needing ~25% less Tris on 30 Sep, 2 Oct and 7 Oct) cannot enter, and the pools stay untouched.
# =====================================================================================
WORKING_DIR = os.environ.get("TECAN_DATA_DIR", r"C:\ProgramData\Tecan\EVOware\database\tecan_data")
LOG_FILE = os.path.join(WORKING_DIR, "python_log.txt")
DECIMAL_SEPARATOR = '.'
VERSION = "EDP-v17-2026-10"

PH_VALID = (1.0, 14.0)
READ_SD, VOL_CV = 0.02, 0.01           # probe reading noise (pH); pipetting scatter (fraction)
TRIS_REGION, TRIS_PKA, TRIS_FIT_OK, PK_PRIOR_SD = 6.8, 8.40, 0.06, 0.03
MIN_STEP, MAX_RATIO = 0.003, 20.0
WELL_ML, PORTION_ML = 3.0, 0.98
C1_SPAN = 50.0                         # acid Cycle 1: 4 wells log-spaced from Max/50 to Max (column 5 since v15)
C1_LOW = (400.0, 200.0, 100.0)         # v15: acid Cycle 1 also at Max/400, Max/200, Max/100 (B1-D1, under the product-only
                                       # well A1): the first locate well (Max/50) left no reading between a pH-7.56 eluate
                                       # and 4.94 on 13 Oct, where the 5.29 pool was asked for
GAP_PH = 0.35                          # v15: map entries between readings this close are labelled "between readings"
DRIFT_MIN = 0.03                       # v15: plate check 3 re-reads column 2 (pH ~3.45-3.8): a mean change this large is drift
C2_N, C2_PAD = 8, (0.80, 1.25)         # acid Cycle 2: 8 wells across the located window, widened by these factors
WINDOW_MARGIN = 0.10                   # pH beyond the lowest and highest acid target the window must reach
FIT_MARGIN = 0.25                      # readings within this pH of the targets enter the robust fit
CLUSTER_PH = 0.6                       # targets this close share one fit
OUTLIER_PH = 0.10                      # a reading this far (5 x the probe noise) off the others' fit is left out, at most 2
EXTEND_PH = 0.15                       # pH a ratio may be extended beyond the readings
POOL_ML, LANE_MAX_ML = 25.0, 29.0      # pools are made at 25 mL; trims stop at 29 mL planned (30 mL lanes)
POOL_TOL, MIN_TRIM_ML = 0.03, 0.08     # trim a pool missing by more than 0.03; smallest addition with 8 tips
BASE_MAX, BASE_SPAN = 0.3, 60.0        # base Cycle 1: 4 wells log-spaced from 0.005 to 0.3 mL 1 M Tris per mL eluate
BASE_WINDOW = (6.8, 8.5)               # base Cycle 2 spreads 4 wells across this pH window (from the eluate's pH if higher)
BASE_REGION = 6.0                      # from this pH up the base curve follows its 3-term equation (acetate, carboxyls spent)
ABS_UL = 1.3                           # uL of pipetting scatter per dispense that does not shrink with the volume (v17: 1.3 uL
                                       # RMS on the 10-65 uL Tris dispenses of 21 and 22 Oct against the chemistry model fitted
                                       # to all readings, 22 Oct's C6 aside; 0.5 until v16). A 3 mL well takes only 20-60 uL
                                       # of 1 M Tris for a lightly buffered eluate, where 1 uL is 0.05-0.2 pH
DIL_TRIS_M = 0.1                       # v17: the dilute Tris (Product 2W RIGHT well) for base wells, unless options.csv says
                                       # otherwise (4th value, Dilute_Tris_M in the script; 0 = not used)
DIL_MAX_ML = 0.98                      # v17: a base well takes the dilute Tris when one portion of it is enough (<= 980 uL,
                                       # i.e. <= 0.485 x its molarity in mL 1 M Tris per mL eluate) - so its dispense is
                                       # always one 4-tip trip, a form that has run with a tip given 0 uL
DIL_TRIM_MAX_ML = 1.0                  # v17: a Tris pool's trim that needs at most this much dilute Tris takes it (8 tips, up to
                                       # 125 uL each; at most ~4% water in a 25 mL pool) - else 1 M Tris, as v16
BASE_FIT_OK = 2.5                      # Tris-equation fit kept if its residuals average within 2.5 x their expected scatter
BASE_OUT_SD = 3.0                      # v17: a base reading is left out if off by more than OUTLIER_PH and this many expected
                                       # scatters (4 until v16, with ABS_UL 0.5: 22 Oct's C6, 7.39 for 7.6, was 3.6 of the
                                       # larger v17 scatter - a well 5 uL of 1 M Tris off its recipe, or a reading slip)
BASE_K_PRIOR = (0.002, 0.004)          # weak-acid term k: typical value and spread (acetate 20-100 mM gives 0.001-0.006);
                                       # without it a few noisy readings low on the curve can push k 10x too high
BASE_C2_SPLIT = (0.05, 0.30, 0.55, 0.85)   # base Cycle 2 wells at these fractions of the window (most where pools usually are)
MAP_STEP, NEUT_STEP = 0.05, 0.1        # ratio map grids: acid and base paths; neutralisation start and target pHs
NEUT_TARGETS = (6.8, 8.5)              # neutralisation targets covered by the map (the Tris buffer region)
NEUT_START_MAX = 5.0                   # neutralisation map from every acid-path pH up to this (or the eluate's pH)
ANCHOR_SD_EXTRA = 0.002                # mL/mL added to the base titration's equivalence SD when it anchors the VI line
LANE_FLAG = 0.15                       # a pool reading this far from what its recipe gives on the run's curve is flagged
TITRANT_SD, TITRANT_RANGE, TITRANT_FLAG = 0.10, (0.75, 1.33), 0.05   # v16: fitted Tris/acid stock strength: prior SD, bounds,
                                       # and the departure from 1 that the log flags (check the stocks)
TITRANT_CHI2 = 9.0                     # v16: the free strength is kept only if it lowers the chi-square this much (3 SD)
CHEM_ACID_OK, CHEM_FIT_OK = 0.04, 2.5  # chemistry model trusted in the gaps of a route if it fits that route's readings to
                                       # 0.04 pH (acid) / 2.5 x their expected scatter (base, neutralisation)
FILES_CHECKED = ("EDP_tecan_gp_master.py", "EDP_parse_readings.py", "EDP_archive_run.py", "EDP_tecan_reset.py",
                 "EDP_make_workbook.py", "EDP_time_estimate.py")
TROUGH_ML = 95.0                       # what one block may take from a 2W trough well (100 mL, 5 mL spare)
TRIM_GAIN = 0.8                        # each trim applies 80% of the computed correction: no overshoot on steep curves
BASE_TRIM_FIRST, BASE_TRIM_LATER = 0.75, 0.9   # v14: Tris pool coming down a steep curve (first trim); later Tris trims
BASE_AIM_MAX = 0.30                    # Tris pools: the first make is aimed above target by up to this (see pool_plan)
TRIS_POOL_MAX = 5.0                    # pools up to this pH count as low-pH pools for the Tris titration
TRIS_C1_FACTORS = (1.06, 1.16, 1.32, 1.60)   # Tris titration: x the Tris the pool's acid needs (stoichiometry)
LANES = 6


def P(name):
    return os.path.join(WORKING_DIR, name)

def read_text(name):
    try:
        with open(P(name)) as f:
            return f.read().strip()
    except Exception:
        return ""

def load_state():
    try:
        with open(P("run_state.json")) as f:
            return json.load(f)
    except Exception:
        return {}

def save_state(st):
    with open(P("run_state.json"), "w") as f:
        json.dump(st, f, indent=1)

def acid_params():
    row = nums(("titration_params.csv"))   # Product_Vol, Max_Titrant, T1, T2, T3, Process, ACCEPT
    return dict(max=min(row[1], MAX_RATIO), targets=[row[2], row[3], row[4], row[5]], accept=row[6]) if len(row) >= 7 else \
        dict(max=4.0, targets=[3.6, 3.4, 3.2, 3.5], accept=0.2)

def tris_params():
    row = nums(("tris_params.csv"))        # nb_t1, nb_t2, nb_t3, nb_proc
    return dict(targets=row[:4])

def tris_wanted():
    o = nums(("options.csv"))                # neut (VI neutralisation data), base_on (base titration)
    return bool(o and o[0] == 1)

def base_wanted():
    o = nums(("options.csv"))
    return bool(len(o) > 1 and o[1] == 1)

def dil_m():
    """v17: molarity of the dilute Tris for the base wells (options.csv: neut, base_on, Use_8_tips, Dilute_Tris_M);
    0 = not used. An options.csv without the 4th value (older script) means not used."""
    o = nums(("options.csv"))
    if len(o) > 3:
        return o[3] if 0.01 <= o[3] <= 0.5 else 0.0
    return 0.0

# ----------------------------------------------------------------------------- proven helpers (unchanged from v10)
_NUM = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")

def nlog(msg):
    try:
        with open(LOG_FILE, "a") as f: f.write(msg + "\n")
    except Exception: pass

def read_numbers(path, expected=None):
    """Numeric values from an EVOware export, whatever the layout (header row or not, one row or one per line).
    If the export appends rather than overwrites, keep only the most recent `expected` values."""
    with open(path, "r") as f:
        tokens = re.split(r"[,;=\s]+", f.read())
    vals = [float(t.strip('"')) for t in tokens if _NUM.fullmatch(t.strip('"'))]
    if expected and len(vals) > expected and len(vals) % expected == 0:
        nlog(f"WARNING: {os.path.basename(path)} holds {len(vals)} values (export appends?) - using the last {expected}")
        vals = vals[-expected:]
    return vals

def write_lines(path, values, fmt="{:.4f}"):
    with open(path, "w") as f:
        for v in values: f.write(fmt.format(v).replace(".", DECIMAL_SEPARATOR) + "\n")

def sorted_arrays(v, p):
    v = np.asarray(v, float); p = np.asarray(p, float)
    o = np.argsort(v, kind="stable")
    return v[o], p[o]

def isotonic(v, y, increasing=True):
    """Least-squares fit that only moves one way (pool-adjacent-violators); v sorted, repeated v values averaged."""
    y = np.asarray(y, float) if increasing else -np.asarray(y, float)
    _, inv, cnt = np.unique(v, return_inverse=True, return_counts=True)
    means = np.bincount(inv, weights=y) / cnt
    blocks = []                                        # [mean, weight, number of distinct v values]
    for m_, w_ in zip(means, cnt):
        blocks.append([m_, float(w_), 1])
        while len(blocks) > 1 and blocks[-2][0] > blocks[-1][0]:
            b2, b1 = blocks.pop(), blocks.pop()
            w = b1[1] + b2[1]
            blocks.append([(b1[0] * b1[1] + b2[0] * b2[1]) / w, w, b1[2] + b2[2]])
    fit = np.concatenate([[b[0]] * b[2] for b in blocks])[inv]
    return fit if increasing else -fit

def monotone(v, p, rising):
    """Sorted ratios and readings plus the best-fitting curve that only moves one way with added titrant."""
    v, p = sorted_arrays(v, p)
    return v, p, isotonic(v, p, bool(rising))

def pooled(v, p, rising):
    """Monotone curve as distinct points: readings the monotone fit had to average become one point at their mean ratio."""
    v, p, fit = monotone(v, p, rising)
    bv, bp, start = [], [], 0
    for k in range(1, len(v) + 1):
        if k == len(v) or not np.isclose(fit[k], fit[start]):
            bv.append(v[start:k].mean()); bp.append(fit[start]); start = k
    return np.array(bv), np.array(bp)

def pchip_ratio(v, p, target, rising):
    """Ratio at which the monotone curve through the readings reaches the target (monotone cubic interpolation),
    or None if the target lies outside the readings."""
    bv, bp = pooled(v, p, rising)
    if len(bv) < 2:
        return None
    s, t = (bp, target) if rising else (-bp, -target)
    if not (s[0] <= t <= s[-1]):
        return None
    if np.any(np.diff(s) <= 0):                        # flat steps: straight line between the bracketing points
        i = int(np.flatnonzero(s >= t)[0])
        if i == 0: return float(bv[0])
        return float(bv[i - 1] + (t - s[i - 1]) * (bv[i] - bv[i - 1]) / (s[i] - s[i - 1]))
    from scipy.interpolate import PchipInterpolator
    return float(PchipInterpolator(s, bv)(t))

def well_sd(t, te):
    """Expected pH scatter of a well at Tris ratio t: the reading, plus the pipetting error magnified by the steepness
    of the curve (just past the equivalence point 1% more Tris moves pH by 0.1 or more), so those wells count less."""
    t = np.asarray(t, float)
    steep = t / np.maximum(t - te, 1e-6) / np.log(10)
    return np.sqrt(READ_SD ** 2 + (VOL_CV * steep) ** 2)

TRIS_LEFT_OUT = []                     # v17: (ratio, pH) of the readings the last fit_tris_joint left out (wells table)

def fit_tris_joint(sets, _retry=True, pk_prior=None):
    """One Tris-equation fit over several groups of wells (readings above TRIS_REGION): a shared pKa' and one
    equivalence point per group. sets: list of (ratios, pH). Returns (pk, [te per group], rms, n) or None.
    v8: a reading at or just below the equivalence point can read above pH 6.8 (the curve is nearly vertical there)
    but cannot follow the equation; it shows up as the group's equivalence point pinned at that reading's ratio, the
    upper bound of the fit. Such a reading is dropped and the fit repeated once.
    v17: each reading's weight comes from the fitted equivalence point (its expected scatter grows without bound as te
    nears its ratio), so a fit can also pin te on a reading that DOES follow the equation. On the 22 Oct eluate in the
    simulator, with both columns' first readings below 6.8, the retry then left out both columns' next readings and too
    few were left: the fit was not usable. When that happens (or the retry still misses), the fit is repeated with each
    reading weighted at its own measured pH ((1 + 10^(pKa' - pH)) / ln 10 times the pipetting error - no vanishing
    weights), leaving out the lowest reading of a group only if needed: the fewest that bring it within TRIS_FIT_OK.
    On the real readings of 12, 21 and 22 Oct the first fit (or its retry) is usable, so nothing changes there."""
    from scipy.optimize import least_squares
    from itertools import combinations
    groups = []
    for v, p in sets:
        v, p = np.asarray(v, float), np.asarray(p, float)
        m = p >= TRIS_REGION
        groups.append((v[m], p[m]))

    def solve(gs, at_reading=False):
        n = sum(len(g[0]) for g in gs)
        used = [i for i, g in enumerate(gs) if len(g[0])]
        if n < len(used) + 2 or not used:
            return None
        def res(q):
            out = []
            for j, i in enumerate(used):
                tv, tp = gs[i]
                pred = q[0] + np.log10(np.maximum(tv - q[1 + j], 1e-9) / q[1 + j])
                sd = (np.sqrt(READ_SD ** 2 + (VOL_CV * (1 + 10 ** (q[0] - tp)) / np.log(10)) ** 2) if at_reading
                      else well_sd(tv, q[1 + j]))
                out.append((pred - tp) / sd)
            if pk_prior is not None:
                out.append(np.array([(q[0] - pk_prior) / PK_PRIOR_SD]))
            return np.concatenate(out)
        his = [gs[i][0].min() * (1 - 1e-4) for i in used]
        best = None
        for frac in (0.9, 0.97, 0.995):
            try:
                r = least_squares(res, [TRIS_PKA if pk_prior is None else pk_prior] + [frac * h for h in his],
                                  bounds=([7.3] + [1e-4] * len(used), [9.5] + his))
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
        if best is None:
            return None
        te = [None] * len(gs)
        raw = []
        for j, i in enumerate(used):
            te[i] = float(best.x[1 + j])
            tv, tp = gs[i]
            raw.append(best.x[0] + np.log10(np.maximum(tv - te[i], 1e-9) / te[i]) - tp)
        return float(best.x[0]), te, float(np.sqrt(np.mean(np.concatenate(raw) ** 2))), n, his, used

    lowest = lambda g: int(np.argmin(g[0]))
    drop_low = lambda gs, idx: [(np.delete(g[0], lowest(g)), np.delete(g[1], lowest(g))) if i in idx else g for i, g in enumerate(gs)]
    def done(fit, dropped):
        del TRIS_LEFT_OUT[:]
        for i in dropped:
            g = groups[i]
            nlog(f"Tris fit: reading pH {g[1][lowest(g)]:.2f} at ratio {g[0][lowest(g)]:.4f} sits at the equivalence point - left out")
            TRIS_LEFT_OUT.append([float(g[0][lowest(g)]), float(g[1][lowest(g)])])
        return fit[:4] if fit else None

    full = solve(groups)
    if full is None or full[2] <= TRIS_FIT_OK or not _retry:
        return done(full, ())
    best, best_drop = full, ()
    pinned = [i for j, i in enumerate(full[5]) if full[1][i] >= full[4][j] * (1 - 5e-3) and len(groups[i][0]) >= 2]
    if pinned:                                                   # v8: the reading at the equivalence point
        again = solve(drop_low(groups, pinned))
        if again and again[2] < best[2]:
            best, best_drop = again, tuple(pinned)
            if best[2] <= TRIS_FIT_OK:
                return done(best, best_drop)
    # v17: weights at each reading's own pH; the fewest lowest readings left out that bring the fit within TRIS_FIT_OK
    cands = [i for i, g in enumerate(groups) if len(g[0]) >= 2]
    for k in range(0, len(cands) + 1):
        found = None
        for drop in combinations(cands, k):
            r = solve(drop_low(groups, drop), at_reading=True)
            if r is not None and (found is None or r[2] < found[0][2]):
                found = (r, drop)
        if found and found[0][2] < best[2]:
            best, best_drop = found
            if best[2] <= TRIS_FIT_OK:
                break
    return done(best, best_drop)


# ----------------------------------------------------------------------------- plans and recipes
def acid_well(r):
    """3 mL well: eluate + 1 M acetic acid at ratio r (mL acid per mL eluate)."""
    e = WELL_ML / (1 + r)
    return (e, WELL_ML - e, 0.0)

def recipe_well(rp, t):
    """3 mL well: a low-pH pool made from eluate + acid at ratio rp, then t mL 1 M Tris per mL of that pool."""
    pool = WELL_ML / (1 + t)
    e = pool / (1 + rp)
    return (e, pool - e, WELL_ML - pool)

def write_plan(name, wells, needs=None):
    """Dispense plan, 15 numbers per 4-well plate column: for eluate (or pool), acid and Tris in turn the number of
    portions and the uL per portion for tips 1-4. v17: base plans (wells of 4 liquids) have 20 numbers per column, the
    4th liquid the dilute Tris (Product 2W RIGHT well). Returns the mL of each liquid."""
    wells = list(wells)
    nliq = max((len(w) for w in wells), default=3)
    while len(wells) % 4:
        wells.append((0.0,) * nliq)
    rows = []
    for i in range(0, len(wells), 4):
        col = wells[i:i + 4]
        for k in range(nliq):
            vols = [w[k] if k < len(w) else 0.0 for w in col]
            top = max(vols)
            n = int(math.ceil(top / PORTION_ML - 1e-9)) if top > 1e-6 else 0
            rows += [n] + [round(1000.0 * v / n, 1) if n else 0.0 for v in vols]
    write_lines(P(name), rows, fmt="{:g}")
    tot = [sum(w[k] if k < len(w) else 0.0 for w in wells) for k in range(nliq)]
    if needs:
        write_lines(P(needs), tot, fmt="{:.1f}")
    return tot

def lane_rows(per_lane):
    """6 lanes x [flag, acid uL/tip, eluate uL/tip, Tris uL/tip] for 8-tip additions into the 6-well plate. v17: flag 2 =
    the Tris comes from the dilute Tris (Product 2W RIGHT well), its uL per tip in the Tris place (trims only)."""
    rows = []
    for lane in range(1, LANES + 1):
        a, e, b, *kind = per_lane.get(lane, (0.0, 0.0, 0.0))
        flag = (kind[0] if kind else 1) if a + e + b > 0 else 0
        rows += [flag, round(a * 1000 / 8, 1), round(e * 1000 / 8, 1), round(b * 1000 / 8, 1)]
    return rows

# ----------------------------------------------------------------------------- the acid curve
class AcidCurve:
    """pH against acid ratio (mL 1 M acetic acid per mL eluate) for this run's eluate.
    Targets within CLUSTER_PH of each other share one fit (the usual 3.2-3.6 set is one cluster; a pool at 4.2 or 7.4
    gets its own): the readings within FIT_MARGIN of the cluster (at least the 5 nearest), first checked one by one
    against a fit to all the others (leave-one-out) - a reading off by more than OUTLIER_PH (1.5 x at the window's
    ends) is left out, worst first, at most 2 per cluster - then a Huber-weighted polynomial in log10(ratio) (degree 3
    with 8+ readings, 2 with 5+, else 1) that must fall steadily and fit within twice the probe noise. Elsewhere:
    monotone interpolation through the kept readings and the eluate's own pH at ratio 0, then a short log-linear
    extension. A threshold set by the probe's noise, not by how well the curve fits, keeps a poor fit from hiding a
    typo (7 Oct: 3.08 at ratio 0.25, 0.26 below its neighbours, decided v9's 3.2 ratio)."""

    def __init__(self, v, p, ph0=None, targets=None):
        v, p = np.asarray(v, float), np.asarray(p, float)
        ok = (v > 0) & (p >= PH_VALID[0]) & (p <= PH_VALID[1])
        order = np.argsort(v[ok], kind="stable")
        self.v, self.p = v[ok][order], p[ok][order]
        self.ph0 = ph0 if ph0 is not None and PH_VALID[0] <= ph0 <= PH_VALID[1] else None
        self.keep = np.ones(len(self.v), bool)
        self.dropped, self.models = [], []
        clusters = []
        for t in sorted(set(targets or [])):
            if clusters and t - clusters[-1][-1] <= CLUSTER_PH:
                clusters[-1].append(t)
            else:
                clusters.append([t])
        for cl in clusters:
            if len(self.v) >= 3:
                self._fit(min(cl), max(cl))

    def _window(self, lo, hi):
        kept = [i for i in range(len(self.v)) if self.keep[i]]
        sel = [i for i in kept if lo - FIT_MARGIN <= self.p[i] <= hi + FIT_MARGIN]
        if len(sel) < 5:
            rest = sorted((i for i in kept if i not in sel), key=lambda i: min(abs(self.p[i] - lo), abs(self.p[i] - hi)))
            sel += rest[:5 - len(sel)]
        return np.array(sorted(sel, key=lambda i: self.v[i]), int)

    @staticmethod
    def _robust_poly(u, y):
        n = len(u)
        for deg in ([3, 2, 1] if n >= 8 else [2, 1] if n >= 5 else [1]):
            w = np.ones(n)
            for _ in range(10):
                coef = np.polyfit(u, y, deg, w=np.sqrt(w))
                res = y - np.polyval(coef, u)
                s = max(READ_SD, 1.4826 * float(np.median(np.abs(res))))
                w = np.minimum(1.0, 1.5 * s / np.maximum(np.abs(res), 1e-12))
            grid = np.linspace(u.min() - 0.1, u.max() + 0.1, 80)
            if np.all(np.polyval(np.polyder(coef), grid) < 0):
                return coef, res, deg
        return None, None, 0

    def _fit(self, lo, hi):
        for _ in range(2):                            # leave-one-out check, at most 2 readings left out
            idx = self._window(lo, hi)
            if len(idx) < 6:
                break
            u, y = np.log10(self.v[idx]), self.p[idx]
            worst, worst_x, worst_r = None, 1.0, 0.0
            for j in range(len(idx)):
                m = np.ones(len(idx), bool); m[j] = False
                coef, _, _ = self._robust_poly(u[m], y[m])
                if coef is None:
                    continue
                r = y[j] - float(np.polyval(coef, u[j]))
                limit = OUTLIER_PH if 0 < j < len(idx) - 1 else 1.5 * OUTLIER_PH
                if abs(r) / limit > worst_x:
                    worst, worst_x, worst_r = j, abs(r) / limit, r
            if worst is None:
                break
            k = idx[worst]
            self.keep[k] = False
            self.dropped.append((float(self.v[k]), float(self.p[k]), float(worst_r)))
        idx = self._window(lo, hi)
        if len(idx) < 3:
            return
        u, y = np.log10(self.v[idx]), self.p[idx]
        coef, res, deg = self._robust_poly(u, y)
        if coef is None:
            return
        rms = float(np.sqrt(np.mean(res ** 2)))
        if rms <= 2 * READ_SD:
            self.models.append(dict(lo=lo, hi=hi, coef=coef, urange=(float(u.min()), float(u.max())), deg=deg, rms=rms, n=len(idx)))

    def points(self):
        v, p = self.v[self.keep], self.p[self.keep]
        if self.ph0 is not None:
            v, p = np.r_[0.0, v], np.r_[self.ph0, p]
        return v, p

    def gap_at(self, ph):
        """v15: pH distance between the two readings (the eluate's own pH counts) that bracket ph; inf outside them."""
        v, p = self.points()
        if len(p) < 2:
            return float("inf")
        bv, bp = pooled(v, p, False)
        for k in range(len(bp) - 1):
            if bp[k + 1] - 1e-9 <= ph <= bp[k] + 1e-9:
                return float(bp[k] - bp[k + 1])
        return float("inf")

    def ratio_for(self, ph):
        """(ratio, how, sure) - or (None, reason, False)."""
        for m in sorted(self.models, key=lambda m: 0 if m["lo"] <= ph <= m["hi"] else min(abs(ph - m["lo"]), abs(ph - m["hi"]))):
            a, b = m["urange"][0] - 0.1, m["urange"][1] + 0.1
            if np.polyval(m["coef"], b) <= ph <= np.polyval(m["coef"], a):
                lo, hi = a, b
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    if np.polyval(m["coef"], mid) > ph:
                        lo = mid
                    else:
                        hi = mid
                return float(10 ** (0.5 * (lo + hi))), f"robust fit ({m['n']} readings, degree {m['deg']}, RMS {m['rms']:.3f})", True
        v, p = self.points()
        bv, bp = pooled(v, p, False)
        if self.ph0 is not None and len(bv) >= 2 and bv[0] == 0 and bp[1] <= ph <= bp[0] and bp[0] > bp[1]:
            # v13: between the eluate and the first reading there is no data to shape a curve; near the eluate's own
            # pH its buffering is roughly constant, so a straight line in ratio (2 Oct eluate, pool at 7.2: the
            # monotone cubic gave 0.0006 mL acid per mL where ~0.009 was needed)
            return float(bv[1] * (bp[0] - ph) / (bp[0] - bp[1])), "straight line from the eluate's own pH to the first reading", False
        x = pchip_ratio(v, p, ph, False)
        if x is not None:
            return float(x), "monotone interpolation", True
        if ph > bp[0]:
            return None, "above the eluate's own pH", False
        if ph < bp[-1] - EXTEND_PH:
            return None, f"more than {EXTEND_PH} pH below the lowest reading", False
        if len(bv) >= 2 and bv[-2] > 0 and bp[-1] != bp[-2]:
            slope = (bp[-1] - bp[-2]) / (math.log10(bv[-1]) - math.log10(bv[-2]))
            if slope < 0:
                return float(min(10 ** (math.log10(bv[-1]) + (ph - bp[-1]) / slope), MAX_RATIO)), "extended beyond the readings", False
        return None, "beyond the readings", False

def acid_readings(st, with_val):
    """(ratios, pH, cycle) of the acid titration so far; the eluate's own pH is kept separately (state)."""
    a = st.get("acid", {})
    v, p, c = [], [], []
    for key, fname, cyc in (("c1", "acid_c1_measured.csv", 1), ("c2", "cycle2_measured.csv", 2), ("val", "validation_results.csv", 3)):
        if key == "val" and not with_val:
            continue
        ratios = (a.get("ratios") if key == "val" else a.get(key)) or []      # the validation wells sit at the found ratios
        vals = nums((fname))
        if key == "c1" and vals:
            vals = vals[1:]                            # the first value typed is the eluate's own pH (A1, product only)
        if key == "val" and a.get("drift"):            # v15: probe drift between plate checks 2 and 3 (column 2 re-read)
            vals = [x - a["drift"] if PH_VALID[0] <= x <= PH_VALID[1] else x for x in vals]
        if ratios and len(vals) >= len(ratios):
            v += list(ratios); p += list(vals[:len(ratios)]); c += [cyc] * len(ratios)
    return np.array(v, float), np.array(p, float), np.array(c, int)

def acid_targets(st):
    """Every acid pH the curve must be good at: the four titration targets, the acid pools, and the process pH."""
    ap = acid_params()
    t = list(ap["targets"])
    for pl in st.get("pools", []):
        if pl["route"] == "acid":
            t.append(pl["target"])
    return t

def make_curve(st, with_val):
    v, p, c = acid_readings(st, with_val)
    return AcidCurve(v, p, st.get("acid", {}).get("ph0"), acid_targets(st)), v, p, c

def nums(name):
    """Numbers in a working file; [] if it is missing or empty."""
    try:
        return read_numbers(P(name)) if os.path.exists(P(name)) else []
    except Exception:
        return []

def write_curve_table(curve, v, p, c, ph0):
    """current_ph_data.csv for the archive: the eluate (ratio 0), then every reading with its cycle and whether the
    fit kept it."""
    kept = {(round(x, 6), round(y, 4)) for x, y in zip(curve.v[curve.keep], curve.p[curve.keep])}
    rows = [dict(Volume_Added=0.0, pH=ph0 if ph0 is not None else 0.0, Cycle=1, Valid=1 if ph0 is not None else 0)]
    for x, y, cyc in zip(v, p, c):
        ok = PH_VALID[0] <= y <= PH_VALID[1] and ((round(x, 6), round(y, 4)) in kept)
        rows.append(dict(Volume_Added=x, pH=y, Cycle=int(cyc), Valid=int(ok)))
    df = pd.DataFrame(rows)
    df["GPR_Predicted"], df["GPR_Variance"] = np.nan, np.nan
    df[["Volume_Added", "pH", "GPR_Predicted", "GPR_Variance", "Cycle", "Valid"]].to_csv(P("current_ph_data.csv"), index=False)

def log_dropped(curve):
    for x, y, r in curve.dropped:
        nlog(f"Reading pH {y:.2f} at ratio {x:.4f} is {r:+.2f} off the curve through its neighbours - left out")

# ----------------------------------------------------------------------------- set-up: the pool request
def parse_numbers(text):
    out = []
    for tok in re.split(r"[,;\s/]+", text.replace('"', " ").strip()):
        if tok:
            out.append(float(tok))
    return out

# ----------------------------------------------------------------------------- start-up: are all six files this version?
def check():
    """Every EDP Python file in tecan_data must carry this VERSION. On 12 Oct the robot ran v12 analysis files with
    older helpers: the old EDP_parse_readings.py could not read the pool entry, crashed, and the previous entry's
    'accepted' flag let the run go on - the pool readings were never stored. Now a mismatch stops the run at the start."""
    bad = 0
    for k, name in enumerate(FILES_CHECKED, start=1):
        try:
            with open(P(name), encoding="latin-1") as f:
                ok = VERSION in f.read()
        except Exception:
            ok = False
        if not ok and not bad:
            bad = k
            nlog(f"ERROR: {name} in {WORKING_DIR} is not {VERSION} - copy all six EDP_*.py files from the zip's python folder")
    write_lines(P("check_status.csv"), [0 if bad else 1, bad], fmt="{}")
    if not bad:
        nlog(f"All six EDP files are {VERSION}")

# ----------------------------------------------------------------------------- acid titration
def acid0():
    """Acid Cycle 1 (v15), 8 wells: column 1 = A1 product only (the eluate's own pH - no probe in the product trough)
    and acid at Max/400, Max/200, Max/100; column 5 = the 4 locate wells, Max/50 to Max."""
    st = {"version": VERSION, "started": datetime.now().strftime("%Y-%m-%d %H:%M:%S"), "pools": []}
    ap = acid_params()
    R = ap["max"]
    low = [R / d for d in C1_LOW]
    ratios = [R / C1_SPAN * C1_SPAN ** (k / 3) for k in range(4)]
    need = write_plan("plan_c1.csv", [acid_well(0.0)] + [acid_well(r) for r in low] + [acid_well(r) for r in ratios])
    # nothing asks for a top-up between plate check 1 and Cycle 2, so READY covers Cycle 2's worst case too
    write_lines(P("needs_c1.csv"), [need[0] + 8 * WELL_ML, need[1] + 8 * WELL_ML * R / (1 + R), 0.0], fmt="{:.1f}")
    st.setdefault("acid", {})["c1"] = low + ratios
    save_state(st)
    nlog("Acid Cycle 1: A1 product only; acid ratios " + ", ".join(f"{r:.4f}" for r in low) + " (B1-D1), "
         + ", ".join(f"{r:.4f}" for r in ratios) + " (A5-D5, locate)")
    write_lines(P("status.csv"), [1], fmt="{}")

def locate(v, p, ph0, ph):
    """Ratio where the Cycle-1 readings (and the eluate's own pH at 0) reach `ph`: straight lines in log10(ratio)
    between readings, linear in ratio between the eluate and the first reading, extended in log ratio beyond them."""
    bv, bp = pooled(v, p, False)
    if ph0 is not None and ph >= bp[0]:
        return bv[0] * max(ph0 - ph, 0.0) / (ph0 - bp[0]) if ph0 > bp[0] else bv[0] * 0.5
    u = np.log10(bv)
    if ph >= bp[0]:
        k = 0
    elif ph <= bp[-1]:
        k = len(bv) - 2
    else:
        return float(10 ** np.interp(ph, bp[::-1], u[::-1]))
    if len(bv) < 2 or bp[k] == bp[k + 1]:
        return float(bv[k])
    slope = (bp[k + 1] - bp[k]) / (u[k + 1] - u[k])
    return float(min(10 ** (u[k] + (ph - bp[k]) / slope), MAX_RATIO)) if slope < 0 else float(bv[k])

def acid1():
    st = load_state()
    a = st["acid"]
    vals = nums(("acid_c1_measured.csv"))
    ph0 = vals[0] if vals and PH_VALID[0] <= vals[0] <= PH_VALID[1] else None
    a["ph0"] = ph0
    v, p = np.array(a["c1"], float), np.array((vals[1:1 + len(a["c1"])] + [0.0] * len(a["c1"]))[:len(a["c1"])], float)
    good = (p >= PH_VALID[0]) & (p <= PH_VALID[1])
    if good.sum() < 3:
        write_lines(P("status.csv"), [0], fmt="{}"); write_lines(P("cycle2_vols.csv"), [0], fmt="{}")
        save_state(st)
        return nlog(f"ERROR: only {int(good.sum())} usable acid Cycle-1 readings")
    v, p = v[good], p[good]
    if np.any(np.diff(p) > 0.05):
        nlog("Acid Cycle 1: the pH does not fall steadily with acid - check the readings (pooled for placement)")
    if ph0 is None:
        nlog("Acid Cycle 1: no reading for the product-only well A1 - the eluate's own pH is unknown this run")
    targets = acid_targets(st)
    thi, tlo = max(targets) + WINDOW_MARGIN, min(targets) - WINDOW_MARGIN
    r_hi, r_lo = locate(v, p, ph0, thi), locate(v, p, ph0, tlo)
    lo, hi = max(r_hi * C2_PAD[0], 0.003), min(max(r_lo * C2_PAD[1], r_hi * 1.5), MAX_RATIO)
    pts = [lo * (hi / lo) ** (k / (C2_N - 1)) for k in range(C2_N)]
    for k in range(1, C2_N):
        pts[k] = max(pts[k], pts[k - 1] + MIN_STEP)
    a["c2"] = pts
    write_plan("plan_c2.csv", [acid_well(r) for r in pts], needs="needs_c2.csv")
    write_lines(P("cycle2_vols.csv"), [pts[0]])
    write_lines(P("status.csv"), [1], fmt="{}")
    save_state(st)
    nlog(f"Eluate pH {ph0} (A1, product only); acid Cycle 1 readings " + ", ".join(f"{x:.4f}:{y:.2f}" for x, y in zip(v, p)))
    nlog(f"Acid window pH {tlo:.2f}-{thi:.2f} located at ratios {r_lo:.4f}-{r_hi:.4f}; Cycle 2 (8 wells): "
         + ", ".join(f"{x:.4f}" for x in pts))

def acid2():
    st = load_state()
    a = st["acid"]
    ap = acid_params()
    curve, v, p, c = make_curve(st, with_val=False)
    log_dropped(curve)
    ratios, methods, sure = [], [], True
    for t in ap["targets"]:
        x, how, ok = curve.ratio_for(t)
        if x is None:
            x, how, ok = 0.0, how, False
        ratios.append(x); methods.append(how); sure &= ok
        nlog(f"pH {t:g}: {how} = ratio {x:.4f}")
    if not all(r > 0 for r in ratios[:4]):
        nlog("WARNING: not every titration target could be reached - check the acid settings")
    status = 2 if sure else 3
    write_lines(P("ratios_clean.csv"), ratios)
    write_lines(P("corr_ratios.csv"), ratios)
    write_lines(P("last_acid_ratios.csv"), ratios)
    with open(P("method_used.txt"), "w") as f:
        f.write("\n".join(methods))
    write_lines(P("status.csv"), [status], fmt="{}")
    write_curve_table(curve, v, p, c, a.get("ph0"))
    a["ratios"] = ratios
    write_plan("plan_val.csv", [acid_well(r) for r in ratios], needs="needs_val.csv")
    tris_on = tris_wanted()
    if tris_on:
        tris_plan_c1(st, curve, ratios)
    else:
        for name in ("plan_tris_c1.csv", "needs_tris_c1.csv"):
            open(P(name), "w").close()
    st["tris_on"] = tris_on
    write_lines(P("val_status.csv"), [int(tris_on)], fmt="{}")
    save_state(st)
    nlog(f"Status {status}")

def drift_check(st):
    """v15: plate check 3 ends the acid validation entry with column 2 read again (A2-D2: Cycle 2, pH ~3.45-3.8, read at
    plate check 2). An acetate solution at pH 3.5 does not change in a few hours, so a different reading is the probe (or
    its temperature). On 13 Oct all four validation wells read 0.02-0.04 high after 4.5 h, while Cycle 2 had fitted to
    0.006 and the run's own chemistry model put Cycle 1 and 2 within 0.005 of it - the validation wells 0.02 above it.
    The mean change of the four wells (one re-read alone is as noisy as the drift it looks for) is taken off the
    validation readings if it is DRIFT_MIN or more, so the refit compares them with Cycle 2 in the same frame. The raw
    readings stay in validation_results.csv. v16: the drift also as % of the reading (Edward, 22 Oct), per well and as a
    mean, in run_state.json (acid.drift_check), drift_check.csv (one row per well, then the mean) and the run's
    workbook (Probe drift sheet); the archive's sessions.csv gets the mean in pH and %."""
    a = st["acid"]
    vals = nums("validation_results.csv")
    a["drift"] = 0.0
    if len(vals) <= 4:
        return vals
    c2 = nums("cycle2_measured.csv")
    with open(P("validation_results.csv"), "w") as f:              # the four validation wells only, for the archive
        f.write(",".join(f"{x:.2f}" for x in vals[:4]) + "\n")
    ok = lambda x: PH_VALID[0] <= x <= PH_VALID[1]
    rows = [(w, t, n) for w, n, t in zip(("A2", "B2", "C2", "D2"), vals[4:8], c2[:4]) if ok(n) and ok(t)]
    if len(rows) < 2:
        nlog("Drift check: column 2 not read again (or not at plate check 2) - no correction")
        return vals[:4]
    d = round(float(np.mean([n - t for _, t, n in rows])), 3)
    pct = round(float(np.mean([100.0 * (n - t) / t for _, t, n in rows])), 2)
    applied = d if abs(d) >= DRIFT_MIN else 0.0
    a["drift_check"] = dict(wells=[w for w, _, _ in rows], pc2=[t for _, t, _ in rows], pc3=[n for _, _, n in rows],
                            mean_ph=d, mean_pct=pct, applied_ph=applied, limit=DRIFT_MIN)
    with open(P("drift_check.csv"), "w") as f:
        f.write("Well,Plate_check_2_pH,Plate_check_3_pH,Change_pH,Change_pct\n")
        for w, t, n in rows:
            f.write(f"{w},{t:.2f},{n:.2f},{n - t:+.3f},{100.0 * (n - t) / t:+.2f}\n")
        f.write(f"mean,,,{d:+.3f},{pct:+.2f}\n")
    txt = ", ".join(f"{t:.2f}->{n:.2f}" for _, t, n in rows)
    if applied:
        a["drift"] = d
        nlog(f"Drift check: column 2 read again {txt} (mean {d:+.3f} pH, {pct:+.2f}%) - probe drift: {d:+.2f} taken off the "
             f"validation readings ({', '.join(f'{x:.2f}' for x in vals[:4])} -> {', '.join(f'{x - d:.2f}' for x in vals[:4])})")
    else:
        nlog(f"Drift check: column 2 read again {txt} (mean {d:+.3f} pH, {pct:+.2f}%) - within the probe's repeatability, no correction")
    return vals[:4]

def acid3():
    st = load_state()
    a = st["acid"]
    ap = acid_params()
    drift_check(st)
    curve, v, p, c = make_curve(st, with_val=True)
    log_dropped(curve)
    vals = [x - a.get("drift", 0.0) if PH_VALID[0] <= x <= PH_VALID[1] else x for x in nums(("validation_results.csv"))]
    corr = []
    for t, r0, m in zip(ap["targets"], a.get("ratios", [0] * 4), vals + [0] * 4):
        x, how, ok = curve.ratio_for(t)
        x = x if x is not None else r0
        corr.append(x)
        nlog(f"pH {t:g}: validated at {m:.2f} with ratio {r0:.4f} -> refit ratio {x:.4f} ({100 * (x / r0 - 1) if r0 else 0:+.1f}%) - {how}")
    write_lines(P("corr_ratios.csv"), corr)
    write_curve_table(curve, v, p, c, a.get("ph0"))
    a["corr"] = corr
    a["dropped"] = [[x, y] for x, y, r in curve.dropped]          # v17: readings the refit left out (wells table)
    save_state(st)

# ----------------------------------------------------------------------------- base titration: Tris straight into the eluate
def base_well(t, dil=None):
    """3 mL well: eluate + Tris at ratio t (mL 1 M Tris - i.e. mmol Tris - per mL eluate), no acid. v17: from the dilute
    Tris (molarity `dil`, default dil_m()) when one portion of it is enough (DIL_MAX_ML), else from the 1 M Tris.
    Returns (eluate, acid, 1 M Tris, dilute Tris) in mL."""
    c = dil_m() if dil is None else dil
    if c > 0 and t > 0:
        e = WELL_ML / (1 + t / c)
        if WELL_ML - e <= DIL_MAX_ML + 1e-9:
            return (e, 0.0, 0.0, WELL_ML - e)
    e = WELL_ML / (1 + t)
    return (e, 0.0, WELL_ML - e, 0.0)

def base_stock(w, dil=None):
    """The Tris stock a base well was made with (M): 1, or the dilute Tris's molarity."""
    return (dil_m() if dil is None else dil) if len(w) > 3 and w[3] > 0 else 1.0

def base_ratio(w, dil=None):
    """mL 1 M Tris (mmol) per mL eluate in a base well (eluate, acid, 1 M Tris, dilute Tris)."""
    c = dil_m() if dil is None else dil
    return (w[2] + (c * w[3] if len(w) > 3 else 0.0)) / w[0]

class BaseCurve:
    """pH against Tris ratio t (mL 1 M Tris per mL eluate) for Tris added straight to the eluate.
    From pH BASE_REGION up, the Tris needed to reach a pH is, exactly in form:
        t(pH) = e + a 10^(pH - 8) - k 10^(6 - pH)
    e: the Tris the eluate's acids take up once fully spent (negative if the eluate already holds Tris base);
    a 10^(pH - 8): the Tris left as free base (the Tris buffer); k 10^(6 - pH): the protons the eluate's weaker acids
    (acetate, protein carboxyls) still hold - for any mix of them, while they are mostly spent (pH > pKa + 1).
    Leaving out that last term (v13 drafts) put the curve 0.1-0.5 pH off at 6.5-7.4 for lightly buffered eluates.
    Linear in (e, a, k): a weighted least-squares fit of the Tris ratios, each reading weighted by its expected scatter
    in Tris - the pipetting error (VOL_CV, and ABS_UL that does not shrink with the volume: a 3 mL well takes only
    20-50 uL of Tris for a lightly buffered eluate) plus the probe noise times the slope dt/dpH. k >= 0, with a weak
    prior (BASE_K_PRIOR) so a few noisy readings low on the curve cannot push it far from what real eluates hold. Before the fit, each reading is checked against a fit to all the others
    (leave-one-out): one off by more than OUTLIER_PH and BASE_OUT_SD x its expected scatter is left out, worst first, at most 2,
    and only where at least two other readings lie within 0.5 pH (a lone reading is what pins the curve there).
    Below BASE_REGION, or if the fit is poor (a buffer near the pool pHs, e.g. citrate or histidine): monotone
    interpolation through the readings and the eluate's own pH at ratio 0, then a short log-linear extension."""

    def __init__(self, v, p, ph0=None, m=None):
        v, p = np.asarray(v, float), np.asarray(p, float)
        m = np.ones(len(v)) if m is None else np.asarray(m, float)       # v17: Tris stock of each reading (M)
        ok = (v > 0) & (p >= PH_VALID[0]) & (p <= PH_VALID[1])
        order = np.argsort(v[ok], kind="stable")
        self.v, self.p, self.m = v[ok][order], p[ok][order], m[ok][order]
        self.ph0 = ph0 if ph0 is not None and PH_VALID[0] <= ph0 <= PH_VALID[1] else None
        self.keep = np.ones(len(self.v), bool)
        self.dropped, self.model = [], None
        self._fit()

    @staticmethod
    def _cols(ph, with_k=True):
        ph = np.asarray(ph, float)
        cols = [np.ones_like(ph), 10 ** (ph - 8.0)] + ([-(10.0 ** (6.0 - ph))] if with_k else [])
        return np.stack(cols, axis=1)

    @staticmethod
    def t_of(q, ph):
        e, a, k = q
        return e + a * 10 ** (np.asarray(ph, float) - 8.0) - k * 10.0 ** (6.0 - np.asarray(ph, float))

    @staticmethod
    def dt_dph(q, ph):
        e, a, k = q
        return np.log(10) * (a * 10 ** (np.asarray(ph, float) - 8.0) + k * 10.0 ** (6.0 - np.asarray(ph, float)))

    @staticmethod
    def _sd_t(q, t, ph, m=None):
        """Expected scatter of a reading, in Tris ratio: pipetting (relative + absolute) and the probe noise. v17: the
        absolute part is ABS_UL of the stock the well was made with (m, M) - 10x less Tris for a 0.1 M Tris - over
        the eluate in the well (3 mL / (1 + t/m))."""
        t = np.asarray(t, float)
        m = np.ones_like(t) if m is None else np.asarray(m, float)
        pip = np.sqrt((VOL_CV * t) ** 2 + (ABS_UL / 1000.0 * m * (1 + t / m) / WELL_ML) ** 2) * (t > 0)
        return np.sqrt(pip ** 2 + (READ_SD * BaseCurve.dt_dph(q, ph)) ** 2)

    @staticmethod
    def _system(t, y, w, with_k):
        """Weighted design matrix and targets, with the k prior as one more row when k is fitted."""
        X = BaseCurve._cols(y, with_k) * w[:, None]
        z = t * w
        if with_k:
            k0, sk = BASE_K_PRIOR
            X, z = np.vstack([X, [0.0, 0.0, 1.0 / sk]]), np.r_[z, k0 / sk]
        return X, z

    @staticmethod
    def _solve(t, y, m=None):
        """(e, a, k) or None: iteratively reweighted linear least squares, k >= 0 with a weak prior, a > 0."""
        with_k = len(t) >= 3
        q = None
        w = np.ones(len(t))
        for _ in range(6):
            coef, *_ = np.linalg.lstsq(*BaseCurve._system(t, y, w, with_k), rcond=None)
            if with_k and coef[2] < 0:                           # no weak-acid term needed (or noise): fit without it
                with_k = False
                coef, *_ = np.linalg.lstsq(*BaseCurve._system(t, y, w, False), rcond=None)
            q = (float(coef[0]), float(coef[1]), float(coef[2]) if with_k else 0.0)
            if q[1] <= 0:
                return None
            w = 1.0 / BaseCurve._sd_t(q, t, y, m)
        return q

    def _region(self):
        idx = [i for i in range(len(self.v)) if self.keep[i] and self.p[i] >= BASE_REGION]
        t, y, m = list(self.v[idx]), list(self.p[idx]), list(self.m[idx])
        if self.ph0 is not None and self.ph0 >= BASE_REGION:            # the eluate itself is in range: t = 0
            t, y, m, idx = [0.0] + t, [self.ph0] + y, [1.0] + m, [-1] + idx
        return np.array(idx, int), np.array(t, float), np.array(y, float), np.array(m, float)

    def _fit(self):
        for _ in range(2):                                                # leave-one-out check, at most 2 left out
            idx, t, y, ms = self._region()
            if len(t) < 6:
                break
            worst, worst_x, worst_r = None, 1.0, 0.0
            for j in range(len(t)):
                if idx[j] < 0 or np.sum(np.abs(y - y[j]) <= 0.5) < 3:
                    continue              # the eluate's own pH, or fewer than 2 other readings within 0.5 pH
                m = np.ones(len(t), bool); m[j] = False
                q = self._solve(t[m], y[m], ms[m])
                if q is None:
                    continue
                rt = t[j] - float(self.t_of(q, y[j]))
                rp = rt / max(float(self.dt_dph(q, y[j])), 1e-12)        # the same miss, in pH
                x = min(abs(rp) / OUTLIER_PH, abs(rt) / (BASE_OUT_SD * float(self._sd_t(q, [t[j]], [y[j]], [ms[j]])[0])))
                if x > worst_x:
                    worst, worst_x, worst_r = j, x, -rp
            if worst is None:
                break
            k = idx[worst]
            self.keep[k] = False
            self.dropped.append((float(self.v[k]), float(self.p[k]), float(worst_r)))
        idx, t, y, ms = self._region()
        if len(t) < 3:
            return
        q = self._solve(t, y, ms)
        if q is None:
            return
        rt = t - self.t_of(q, y)
        rp = rt / np.maximum(self.dt_dph(q, y), 1e-12)
        sd = self._sd_t(q, t, y, ms)
        rmsw = float(np.sqrt(np.mean((rt / sd) ** 2)))
        if rmsw <= BASE_FIT_OK:
            with_k = q[2] > 0
            X = self._system(t, y, 1.0 / sd, with_k)[0]
            try:                                  # parameter covariance, scaled up if the readings scatter more than expected
                cov = np.linalg.inv(X.T @ X) * max(rmsw ** 2, 1.0)
            except np.linalg.LinAlgError:
                cov = None
            self.model = dict(q=q, n=len(t), rms=float(np.sqrt(np.mean(rp ** 2))), rmsw=rmsw, with_k=with_k, cov=cov,
                              pmin=float(y.min()), pmax=float(y.max()))

    def points(self):
        v, p = self.v[self.keep], self.p[self.keep]
        if self.ph0 is not None:
            v, p = np.r_[0.0, v], np.r_[self.ph0, p]
        return v, p

    def top(self):
        """Highest pH the readings support (with the short extension)."""
        kept = self.p[self.keep]
        return float(kept.max()) + EXTEND_PH if len(kept) else 0.0

    def sd_ph(self, ph):
        """How uncertain the fitted curve is at this pH (1 SD, in pH), from the fit's parameter covariance."""
        m = self.model
        if not m or m["cov"] is None:
            return None
        x = self._cols([ph], m["with_k"])[0]
        sd_t = math.sqrt(max(float(x @ m["cov"] @ x), 0.0))
        return sd_t / max(float(self.dt_dph(m["q"], ph)), 1e-12)

    def describe(self):
        m = self.model
        return (f"t = {m['q'][0]:.5f} + {m['q'][1]:.5f} x 10^(pH-8) - {m['q'][2]:.5f} x 10^(6-pH), RMS {m['rms']:.3f} pH "
                f"({m['rmsw']:.1f} x the expected scatter) over {m['n']} readings") if m else "no Tris-equation fit"

    def ratio_for(self, ph):
        """(ratio, how, sure) - or (None, reason, False)."""
        m = self.model
        if self.ph0 is not None and ph <= self.ph0:
            return None, "at or below the eluate's own pH", False
        if m and BASE_REGION - 0.05 <= ph <= m["pmax"] + EXTEND_PH:
            t = float(self.t_of(m["q"], ph))
            if t > 0:
                return t, f"Tris-equation fit ({m['n']} readings, RMS {m['rms']:.3f})", bool(m["pmin"] - 0.1 <= ph <= m["pmax"] + 0.05)
        v, p = self.points()
        bv, bp = pooled(v, p, True)
        if self.ph0 is not None and len(bv) >= 2 and bv[0] == 0 and bp[0] <= ph <= bp[1] and bp[1] > bp[0]:
            return float(bv[1] * (ph - bp[0]) / (bp[1] - bp[0])), "straight line from the eluate's own pH to the first reading", False
        x = pchip_ratio(v, p, ph, True)
        if x is not None:
            return float(x), "monotone interpolation", True
        if len(bv) and ph > bp[-1] + EXTEND_PH:
            return None, f"more than {EXTEND_PH} pH above the highest reading", False
        if len(bv) >= 2 and bv[-2] > 0 and bp[-1] != bp[-2]:
            slope = (bp[-1] - bp[-2]) / (math.log10(bv[-1]) - math.log10(bv[-2]))
            if slope > 0:
                return float(min(10 ** (math.log10(bv[-1]) + (ph - bp[-1]) / slope), MAX_RATIO)), "extended beyond the readings", False
        return None, "beyond the readings", False

def base_readings(st, with_val, stocks=False):
    """(ratios, pH, cycle) of the base titration so far (the eluate's own pH is the acid titration's first value);
    v17: with stocks=True also the Tris stock each well was made with (M; 1 for runs before v17)."""
    b = st.get("base", {})
    v, p, c, m = [], [], [], []
    for key, fname, cyc in (("c1", "base_c1_measured.csv", 1), ("c2", "base_c2_measured.csv", 2), ("val", "base_val_measured.csv", 3)):
        if key == "val" and not with_val:
            continue
        ratios = (b.get("val_wells") if key == "val" else b.get(key)) or []
        ms = b.get(f"{key}_m") or [1.0] * len(ratios)
        vals = nums(fname)
        if ratios and len(vals) >= len(ratios):
            for x, y, s_ in zip(ratios, vals[:len(ratios)], ms):
                if x > 0:
                    v.append(x); p.append(y); c.append(cyc); m.append(s_)
    out = (np.array(v, float), np.array(p, float), np.array(c, int))
    return out + (np.array(m, float),) if stocks else out

def make_base_curve(st, with_val=True):
    v, p, c, m = base_readings(st, with_val, stocks=True)
    return BaseCurve(v, p, st.get("acid", {}).get("ph0"), m), v, p, c

def base_window(ph0):
    lo, hi = BASE_WINDOW
    if ph0 is not None:
        lo = max(lo, ph0 + 0.2)
    return lo, hi

def base0():
    """Base Cycle 1 (locate): 4 wells log-spaced from BASE_MAX/BASE_SPAN to BASE_MAX - blind, made with acid Cycle 1
    before the eluate's pH is known. Measured at plate check 1 with no extra stop."""
    st = load_state()
    if not base_wanted():
        for name in ("plan_base_c1.csv", "needs_base_c1.csv"):
            open(P(name), "w").close()
        st["base"] = {"on": 0}
        save_state(st)
        return
    ratios = [BASE_MAX / BASE_SPAN * BASE_SPAN ** (k / 3) for k in range(4)]
    wells = [base_well(t) for t in ratios]
    need = write_plan("plan_base_c1.csv", wells)
    # nothing asks for a top-up before base Cycle 2: READY covers its worst case too (the dilute Tris: READY asks for 25 mL; the base wells take at most ~8.5 mL of it, the Tris pool trims up to 1 mL each)
    write_lines(P("needs_base_c1.csv"), [need[0] + 4 * WELL_ML, 0.0, need[2] + 4 * WELL_ML * BASE_MAX / (1 + BASE_MAX), need[3]],
                fmt="{:.1f}")
    st["base"] = {"on": 1, "c1": ratios, "c1_m": [base_stock(w) for w in wells], "dil": dil_m()}
    save_state(st)
    nlog("Base Cycle 1 (locate, Tris straight into the eluate): ratios " + ", ".join(f"{t:.4f}" for t in ratios) + stock_note(wells))

def stock_note(wells):
    """v17, for the log: which Tris each base well gets and how much (uL)."""
    c = dil_m()
    if not c:
        return " (1 M Tris)"
    return " - Tris: " + ", ".join(f"{1000 * w[3]:.0f} uL {c:g} M" if w[3] > 0 else f"{1000 * w[2]:.0f} uL 1 M" for w in wells)

def locate_rising(v, p, ph0, ph):
    """Ratio where the base Cycle-1 readings (and the eluate's own pH at 0) reach `ph`: straight lines in log10(ratio)
    between readings, linear in ratio from the eluate to the first reading, extended in log ratio beyond them."""
    bv, bp = pooled(v, p, True)
    if ph0 is not None and ph <= bp[0]:
        return bv[0] * max(ph - ph0, 0.0) / (bp[0] - ph0) if bp[0] > ph0 else bv[0] * 0.5
    u = np.log10(bv)
    if bp[0] < ph < bp[-1]:
        return float(10 ** np.interp(ph, bp, u))
    k = 0 if ph <= bp[0] else len(bv) - 2
    if len(bv) < 2 or bp[k] == bp[k + 1]:
        return float(bv[k])
    slope = (bp[k + 1] - bp[k]) / (u[k + 1] - u[k])
    return float(min(10 ** (u[k] + (ph - bp[k]) / slope), 1.0)) if slope > 0 else float(bv[k])

def base1():
    """After plate check 1: base Cycle 2, 4 wells across the window, most of them low in it (where the curve is steep and
    pools are usually wanted), placed with the Tris equation through the Cycle-1 readings when 3 are in its range."""
    st = load_state()
    b = st.get("base", {})
    ph0 = st.get("acid", {}).get("ph0")
    vals = nums("base_c1_measured.csv")
    v, p = np.array(b.get("c1", []), float), np.array(vals[:4], float)
    good = (p >= PH_VALID[0]) & (p <= PH_VALID[1]) if len(p) == len(v) else np.zeros(len(v), bool)
    lo, hi = base_window(ph0)
    if good.sum() < 2 or lo > hi - 0.3:
        b["on"] = 0
        write_lines(P("base_status.csv"), [0], fmt="{}")
        open(P("plan_base_c2.csv"), "w").close(); open(P("needs_base_c2.csv"), "w").close()
        save_state(st)
        return nlog(f"Base titration stopped: {int(good.sum())} usable Cycle-1 readings" if good.sum() < 2 else
                    f"Base titration stopped: the eluate (pH {ph0}) is already at the top of the window")
    ms = np.array(b.get("c1_m") or [1.0] * len(v), float)[good]
    v, p = v[good], p[good]
    nlog("Base Cycle 1 readings " + ", ".join(f"{x:.4f}:{y:.2f}" for x, y in zip(v, p)))
    if np.any(np.diff(p) < -0.05):
        nlog("Base Cycle 1: the pH does not rise steadily with Tris - check the readings (pooled for placement)")
    targets = [lo + (hi - lo) * f for f in BASE_C2_SPLIT]
    c1 = BaseCurve(v, p, ph0, ms)              # the Tris equation through Cycle 1, if 3 readings are in its range
    how = "Tris-equation fit to Cycle 1" if c1.model else "straight lines between the Cycle-1 readings"
    pts = []
    for t in targets:
        x = c1.ratio_for(t)[0] if c1.model else None
        pts.append(max(x if x else locate_rising(v, p, ph0, t), 0.002))
    for k in range(1, 4):                      # at least 1.5 uL more Tris than the previous well (3 mL wells)
        pts[k] = max(pts[k], pts[k - 1] + max(0.0005, 0.02 * pts[k - 1]))
    wells = [base_well(t) for t in pts]
    b["c2"], b["c2_m"] = pts, [base_stock(w) for w in wells]
    write_plan("plan_base_c2.csv", wells, needs="needs_base_c2.csv")
    write_lines(P("base_status.csv"), [1], fmt="{}")
    save_state(st)
    nlog(f"Base window pH {lo:.2f}-{hi:.2f}; Cycle 2 (4 wells, {how}) at pH " + ", ".join(f"{t:.2f}" for t in targets)
         + ": ratios " + ", ".join(f"{x:.4f}" for x in pts) + stock_note(wells))

def base_targets():
    return (tris_params()["targets"] or [7.2, 7.4, 7.6, 7.4])[:4]

def base2():
    """After plate check 2: the base curve from Cycles 1 and 2, and the validation wells at the Tris targets (a target
    at or below the eluate's own pH gets no well)."""
    st = load_state()
    b = st.get("base", {})
    curve = make_base_curve(st, with_val=False)[0]
    for x, y, r in curve.dropped:
        nlog(f"Base reading pH {y:.2f} at ratio {x:.4f} is {r:+.2f} off the curve through the others - left out")
    ratios, wells, extra = [], [], []
    targets = base_targets()
    below = [t for t in targets if curve.ph0 is not None and t <= curve.ph0 + 0.02]
    lo, hi = base_window(curve.ph0)
    spare = [lo + (hi - lo) * (j + 0.5) / max(len(below), 1) for j in range(len(below))]   # extra curve points instead
    for t in targets:
        x, how, sure = (None, "at or below the eluate's own pH", False) if t in below else curve.ratio_for(t)
        ratios.append(x or 0.0)
        if x is None:                      # a well is never left empty: mixing an empty well can stop the robot
            t2 = spare.pop(0) if spare else hi
            x2 = curve.ratio_for(t2)[0] or locate_rising(*base_readings(st, False)[:2], curve.ph0, t2)
            extra.append(x2)
            wells.append(base_well(max(x2, 0.002)))
            nlog(f"Base pH {t:g}: {how} - its well becomes a curve point at pH {t2:.2f} (ratio {x2:.4f})")
        else:
            wells.append(base_well(x))
            nlog(f"Base pH {t:g}: {how} = {x:.4f} mL Tris per mL eluate")
    b["ratios"], b["val_wells"], b["val_m"] = ratios, [base_ratio(w) for w in wells], [base_stock(w) for w in wells]
    nlog("Base fit (Cycles 1 and 2): " + curve.describe())
    nlog("Base validation wells" + stock_note(wells))
    write_plan("plan_base_val.csv", wells, needs="needs_base_val.csv")
    write_lines(P("base_status.csv"), [1], fmt="{}")
    write_lines(P("ratios_clean_EL.csv"), ratios)
    save_state(st)

def base3():
    """After plate check 3: refit with the validation wells -> the base ratios at the Tris targets, archived as EL."""
    st = load_state()
    b = st.get("base", {})
    if not b.get("ratios"):
        return
    curve, v, p, c = make_base_curve(st, with_val=True)
    for x, y, r in curve.dropped:
        nlog(f"Base reading pH {y:.2f} at ratio {x:.4f} is {r:+.2f} off the curve through the others - left out")
    val = nums("base_val_measured.csv") + [0.0] * 4
    corr, methods = [], []
    for t, r0, m in zip(base_targets(), b["ratios"], val):
        if not r0:
            corr.append(0.0); methods.append("at or below the eluate's own pH")
            continue
        x, how, ok = curve.ratio_for(t)
        x = x if x is not None else r0
        corr.append(x); methods.append(how)
        nlog(f"Base pH {t:g}: validated at {m:.2f} with {r0:.4f} mL Tris per mL eluate -> refit {x:.4f} ({100 * (x / r0 - 1):+.1f}%) - {how}")
    b["corr"] = corr
    b["dropped"] = [[x, y] for x, y, r in curve.dropped]          # v17: readings the refit left out (wells table)
    if curve.model:                                    # the eluate's own Tris equivalence: anchors the VI line at zero acid
        b["e"] = float(curve.model["q"][0])
        cov = curve.model.get("cov")
        b["sd_e"] = float(math.sqrt(max(cov[0][0], 0.0))) if cov is not None else None
    nlog("Base refit with the validation wells: " + curve.describe())
    write_lines(P("base_corr.csv"), corr)
    # archive files, tagged EL (EDP_archive_run.py Base EL)
    ap = acid_params()
    tt = base_targets()
    with open(P("titration_params_EL.csv"), "w") as f:      # a target at or below the eluate's pH: 0 (not archived)
        f.write(",".join(f"{x:.6f}" for x in [1, BASE_MAX, *[t if r else 0.0 for t, r in zip(tt, b["ratios"])], ap["accept"]]) + "\n")
    kept = {(round(x, 6), round(y, 4)) for x, y in zip(curve.v[curve.keep], curve.p[curve.keep])}
    stock = base_readings(st, True, stocks=True)[3]
    rows = [dict(Volume_Added=0.0, pH=curve.ph0 if curve.ph0 is not None else 0.0, Cycle=1, Valid=1 if curve.ph0 is not None else 0,
                 Stock_M=1.0)]
    for x, y, cyc, s_ in zip(v, p, c, stock):
        rows.append(dict(Volume_Added=x, pH=y, Cycle=int(cyc), Valid=int((round(x, 6), round(y, 4)) in kept), Stock_M=s_))
    df = pd.DataFrame(rows)
    df["GPR_Predicted"], df["GPR_Variance"] = np.nan, np.nan
    # v17: Stock_M - the Tris stock the well was made with (1 M, or the dilute Tris); Volume_Added stays mL 1 M Tris
    # (mmol) per mL eluate, so the curves of every version compare directly
    df[["Volume_Added", "pH", "GPR_Predicted", "GPR_Variance", "Cycle", "Valid", "Stock_M"]].to_csv(P("current_ph_data_EL.csv"), index=False)
    with open(P("validation_results_EL.csv"), "w") as f:
        f.write(",".join(f"{(m if r0 else 0.0):.6f}" for m, r0 in zip(val[:4], b["ratios"])) + "\n")
    write_lines(P("corr_ratios_EL.csv"), corr)
    with open(P("method_used_EL.txt"), "w") as f:
        f.write("\n".join(methods))
    write_lines(P("status_EL.csv"), [2 if curve.model else 3], fmt="{}")
    save_state(st)

# ----------------------------------------------------------------------------- VI neutralisation data: Tris titration (recipe wells)
def tris_plan_c1(st, curve, ratios):
    """Two plate columns of recipe wells (eluate + acid at one acid ratio, then Tris), one per acid fraction: the
    most and the least acid of the run's targets and pools."""
    rs = [r for r in ratios if r > 0]
    for pl in st.get("pools", []):                    # low-pH pools only: what the Tris titration is for
        if pl["route"] == "acid" and pl["target"] <= TRIS_POOL_MAX:
            x, _, _ = curve.ratio_for(pl["target"])
            if x:
                rs.append(x)
    hi_r, lo_r = max(rs), min(rs)
    if hi_r / (1 + hi_r) - lo_r / (1 + lo_r) < 0.03:
        lo_r = hi_r * 0.5
    wells, cols = [], []
    for rp in (hi_r, lo_r):
        f = rp / (1 + rp)
        ts = [f * k for k in TRIS_C1_FACTORS]
        wells += [recipe_well(rp, t) for t in ts]
        cols.append(dict(rp=rp, t=ts))
    need = write_plan("plan_tris_c1.csv", wells, needs="needs_tris_c1.csv")
    st["tris"] = dict(c1=cols, need_c1=need)
    nlog("Tris titration Cycle 1: acid ratios " + ", ".join(f"{c_['rp']:.4f}" for c_ in cols) + "; Tris ratios "
         + " | ".join(", ".join(f"{t:.4f}" for t in c_["t"]) for c_ in cols))

def te_line(points):
    """Least-squares line te = a + b f through (f, te) points."""
    f = np.array([x for x, _ in points]); te = np.array([y for _, y in points])
    if len(points) == 1:
        return 0.0, float(te[0] / f[0])
    b, a = np.polyfit(f, te, 1)
    return float(a), float(b)

def process_ratio(st):
    a = st.get("acid", {})
    r = (a.get("corr") or a.get("ratios") or [0, 0, 0, 0])[3]
    return r

def tris1():
    st = load_state()
    tr = st.get("tris")
    if not tr:
        return
    vals = nums(("tris_c1_measured.csv"))
    sets = []
    for k, col in enumerate(tr["c1"]):
        ph = vals[4 * k:4 * k + 4]
        ok = [(t, y) for t, y in zip(col["t"], ph) if PH_VALID[0] <= y <= PH_VALID[1]]
        sets.append(sorted_arrays([t for t, _ in ok], [y for _, y in ok]))
    fit = fit_tris_joint(sets)
    tr["dropped"] = list(TRIS_LEFT_OUT)                # v17: per mL of pool (ratio, pH) - for the wells table
    fs = [c_["rp"] / (1 + c_["rp"]) for c_ in tr["c1"]]
    if fit and fit[2] <= TRIS_FIT_OK and 7.6 <= fit[0] <= 9.2 and all(x is not None for x in fit[1]):
        pk, tes, rms, n = fit
        a_, b_ = te_line(list(zip(fs, tes)))
        tr.update(pk=pk, a=a_, b=b_, ok=1, rms=rms)
        nlog(f"Tris titration: pKa' {pk:.3f}, equivalence {tes[0]:.4f} at acid fraction {fs[0]:.3f} and {tes[1]:.4f} at "
             f"{fs[1]:.3f} -> te = {a_:.4f} + {b_:.4f} f (RMS {rms:.3f} pH over {n} readings)")
    else:
        tr.update(pk=8.30, a=0.0, b=1.0, ok=0, rms=None)
        nlog("Tris titration fit not usable" + (f" (RMS {fit[2]:.3f}, pKa' {fit[0]:.2f})" if fit else "")
             + " - stoichiometry (te = f, pKa' 8.30) used; pools are trimmed on their readings")
    rp = process_ratio(st)
    fp = rp / (1 + rp)
    te = tr["a"] + tr["b"] * fp
    targets = tris_params()["targets"]
    ts = [te * (1 + 10 ** (t - tr["pk"])) for t in targets]
    how = "Tris equation"
    near = [k for k, t in enumerate(targets) if t < tr["pk"] - 1.0]
    if tr["ok"] and near:
        # v16: within 1 pH below Tris's pKa' (7.2 on most products) the Tris equation leaves out the protons acetate and
        # the eluate's weak groups still hold. On 21 Oct its 7.2 well read 7.26, where the run's chemistry model (with the
        # Tris strength fitted) put that same recipe at 7.26; the map takes those entries from the model already ("near
        # equivalence"), so the validation wells there are planned from it too when it reproduces the VI readings.
        try:
            ch = fit_chemistry(st)
        except Exception:
            ch = None
        if ch and ch.get("good", {}).get("neutralise"):
            el, T = chem_eluate(ch["q"]), ch["q"][5]
            x = chem_titrant_for(el, T, np.full(len(near), rp), np.zeros(len(near)), np.array([targets[k] for k in near], float), "Tris")
            if np.all(np.isfinite(x)):
                for k, xi in zip(near, x):
                    ts[k] = float(xi) / (1 + rp)
                how = f"Tris equation; below pH {tr['pk'] - 1.0:.2f} the chemistry model"
    if not tr["ok"]:
        # v14: without a usable Tris-equation fit, the validation wells are planned from the chemistry model if it
        # reproduces the run's acid and VI readings (stoichiometry alone put them up to 0.5 pH off in the simulator)
        st["tris"] = tr
        try:
            ch = fit_chemistry(st)
        except Exception:
            ch = None
        if ch and ch.get("good", {}).get("neutralise"):
            el, T = chem_eluate(ch["q"]), ch["q"][5]
            x = chem_titrant_for(el, T, np.full(len(targets), rp), np.zeros(len(targets)), np.array(targets, float), "Tris")
            if np.all(np.isfinite(x)):
                ts, how = [float(v) / (1 + rp) for v in x], "chemistry model"
    write_plan("plan_tris_val.csv", [recipe_well(rp, t) for t in ts], needs="needs_tris_val.csv")
    tr.update(val_rp=rp, val_t=ts)
    write_lines(P("tris_ratios.csv"), ts)
    write_lines(P("tris_status.csv"), [tr["ok"]], fmt="{}")
    save_state(st)
    nlog(f"Tris validation at the process pool (acid ratio {rp:.4f}, {how}): " + ", ".join(f"pH {t:g}: {x:.4f}" for t, x in zip(targets, ts)))

def tris2():
    st = load_state()
    tr = st.get("tris")
    if not tr or "val_t" not in tr:
        return
    c1 = nums(("tris_c1_measured.csv"))
    val = nums(("tris_val_measured.csv"))
    groups, fs = [], []
    for k, col in enumerate(tr["c1"]):
        ph = c1[4 * k:4 * k + 4]
        ok = [(t, y) for t, y in zip(col["t"], ph) if PH_VALID[0] <= y <= PH_VALID[1]]
        groups.append(sorted_arrays([t for t, _ in ok], [y for _, y in ok])); fs.append(col["rp"] / (1 + col["rp"]))
    ok = [(t, y) for t, y in zip(tr["val_t"], val) if PH_VALID[0] <= y <= PH_VALID[1]]
    if ok:
        groups.append(sorted_arrays([t for t, _ in ok], [y for _, y in ok])); fs.append(tr["val_rp"] / (1 + tr["val_rp"]))
    fit = fit_tris_joint(groups)
    tr["dropped"] = list(TRIS_LEFT_OUT)                # v17: the refit's own left-out readings (wells table)
    targets = tris_params()["targets"]
    if fit and fit[2] <= TRIS_FIT_OK and 7.6 <= fit[0] <= 9.2 and all(x is not None for x in fit[1]):
        pk, tes, rms, n = fit
        a_, b_ = te_line(list(zip(fs, tes)))
        if not tr.get("ok"):                           # v14: usable with the validation wells, though Cycle 1 alone was not
            tr.update(ok=1, ok_from="refit")
        tr.update(pk2=pk, a2=a_, b2=b_, rms2=rms)
        nlog(f"Tris refit with the validation wells: pKa' {pk:.3f}, te = {a_:.4f} + {b_:.4f} f, RMS {rms:.3f} pH over {n} readings")
    else:
        pk, a_, b_ = tr["pk"], tr["a"], tr["b"]
        tr.update(pk2=pk, a2=a_, b2=b_, rms2=None)
        nlog("Tris refit not usable - Cycle-1 fit kept")
    fp = tr["val_rp"] / (1 + tr["val_rp"])
    corr = [(a_ + b_ * fp) * (1 + 10 ** (t - pk)) for t in targets]
    write_lines(P("tris_corr.csv"), corr)
    for t, r0, r1, m in zip(targets, tr["val_t"], corr, val + [0] * 4):
        nlog(f"Tris pH {t:g}: validated at {m:.2f} with {r0:.4f} mL Tris per mL pool -> refit {r1:.4f} ({100 * (r1 / r0 - 1) if r0 else 0:+.1f}%)")
    write_tris_archive(st, tr, targets, val, corr, pk, a_ + b_ * fp, fit)
    save_state(st)

def write_tris_archive(st, tr, targets, val, corr, pk, te_p, fit):
    """Files the archive reads for a neutralisation experiment, tagged TR (one per run, at the process pool)."""
    ap = acid_params()
    with open(P("titration_params_TR.csv"), "w") as f:              # one comma-separated row, as EVOware exports it
        f.write(",".join(f"{x:.6f}" for x in [1, 0, *targets, ap["accept"], ap["targets"][3], tr["val_rp"]]) + "\n")
    rows = []
    c1 = nums(("tris_c1_measured.csv"))
    for k, col in enumerate(tr["c1"]):
        for t, y in zip(col["t"], c1[4 * k:4 * k + 4]):
            rows.append(dict(Volume_Added=t, pH=y, Cycle=1 if k == 0 else 2, Valid=int(PH_VALID[0] <= y <= PH_VALID[1])))
    for t, y in zip(tr["val_t"], val):
        rows.append(dict(Volume_Added=t, pH=y, Cycle=3, Valid=int(PH_VALID[0] <= y <= PH_VALID[1])))
    df = pd.DataFrame(rows)
    df["GPR_Predicted"], df["GPR_Variance"] = np.nan, np.nan
    df[["Volume_Added", "pH", "GPR_Predicted", "GPR_Variance", "Cycle", "Valid"]].to_csv(P("current_ph_data_TR.csv"), index=False)
    write_lines(P("ratios_clean_TR.csv"), tr["val_t"])
    with open(P("validation_results_TR.csv"), "w") as f:
        f.write(",".join(f"{x:.6f}" for x in (val[:4] if len(val) >= 4 else val + [0] * (4 - len(val)))) + "\n")
    write_lines(P("corr_ratios_TR.csv"), corr)
    with open(P("method_used_TR.txt"), "w") as f:
        f.write("\n".join(["Recipe Tris fit (two acid fractions)"] * 4))
    write_lines(P("status_TR.csv"), [2 if tr.get("ok") else 3], fmt="{}")
    with open(P("tris_fit_TR.csv"), "w") as f:
        f.write("tris_pka,te_recipe,te_real,real_recipe_factor,fit_rms_pH,n_points\n")
        f.write(f"{pk:.4f},{te_p:.5f},{te_p:.5f},1.0000,{(fit[2] if fit else float('nan')):.4f},{fit[3] if fit else 0}\n")

# ----------------------------------------------------------------------------- chemistry model: one physical model of the product
# The same physics as the simulator (charge balance, Davies activities): acetate (pKa 4.756), protein carboxyls (4.0),
# histidine-like groups (6.3), strong ions, Tris (8.07 at 25 C, -0.028 per C). Fitted to every reading of the run's
# three titrations, it is an independent cross-check of the ratio map and a compact fingerprint of the product for the
# ML (acetate, strong ions, carboxyls, histidines, Tris, temperature). It does not drive any decision in the run.
# On the 12 Oct readings (16 acid wells, 12 VI wells) it fitted the acid wells to 0.009 pH and predicted the VI
# validation wells from Cycle 1 alone within 0.06, as the Tris equation did.
CHEM_KEYS = ("Ac", "Tr", "Na", "Cl", "Pc", "Ph")

def chem_ph(c, T, it=200):
    """pH of many mixtures at once (c: dict of arrays, M; T in C): Newton on the charge balance (steps limited to 1 pH),
    Davies coefficients recomputed from the ionic strength as it converges. Matches the simulator's chem.py to 1e-7."""
    c = {k: np.asarray(c.get(k, 0.0), float) for k in CHEM_KEYS}
    n = np.broadcast(*c.values()).shape
    c = {k: np.broadcast_to(v, n).astype(float) for k, v in c.items()}
    pka_ac, pka_tr, pka_c, pka_h = 4.756, 8.07 - 0.028 * (T - 25.0), 4.0, 6.3
    kw = 10 ** -(14.00 + 0.034 * (25.0 - T))
    L = np.log(10)
    ph, I = np.full(n, 7.0), np.full(n, 0.1)
    for _ in range(it):
        sq = np.sqrt(I)
        g = 10 ** (-0.51 * (sq / (1 + sq) - 0.3 * I))
        aH = 10 ** -ph
        H, OH = aH / g, kw / aH / g
        x_ac, x_c = 1 / (1 + aH * g / 10 ** -pka_ac), 1 / (1 + aH * g / 10 ** -pka_c)
        x_tr, x_h = 1 / (1 + 10 ** -pka_tr * g / aH), 1 / (1 + 10 ** -pka_h * g / aH)
        f = H + c["Na"] + c["Tr"] * x_tr + c["Ph"] * x_h - OH - c["Cl"] - c["Ac"] * x_ac - c["Pc"] * x_c
        df = -L * (H + OH + c["Ac"] * x_ac * (1 - x_ac) + c["Pc"] * x_c * (1 - x_c) + c["Tr"] * x_tr * (1 - x_tr)
                   + c["Ph"] * x_h * (1 - x_h))
        step = np.clip(-f / df, -1.0, 1.0)
        I_new = 0.5 * (H + OH + c["Na"] + c["Cl"] + c["Ac"] * x_ac + c["Pc"] * x_c + c["Tr"] * x_tr + c["Ph"] * x_h)
        ph = np.clip(ph + step, -1.0, 15.0)
        dI = I_new - I
        I = I + 0.5 * dI
        if np.all(np.abs(step) < 1e-9) and np.all(np.abs(dI) < 1e-9):
            break
    return ph

def chem_eluate(q):
    """Parameters -> eluate composition: Ac, strong-ion difference (+ Na, - Cl), Pc, Ph, Tr; q[5] is the temperature;
    v16: q[6], if there, the strength of the Tris stock relative to the acid stock (1 = both exactly 1 M), kept with the
    eluate as "_s" so every mixture made from it (chem_mix) uses it."""
    return dict(Ac=q[0], Na=max(q[1], 0.0), Cl=max(-q[1], 0.0), Pc=q[2], Ph=q[3], Tr=q[4], _s=q[6] if len(q) > 6 else 1.0)

def chem_mix(el, A, B, W=0.0):
    """1 mL eluate + A mL 1 M acetic acid + B mL 1 M Tris (arrays) -> concentrations (Tris x the fitted strength, v16).
    v17: + W mL of water - what comes with B mmol of Tris given as a dilute Tris (B / molarity - B mL)."""
    A, B, W = np.asarray(A, float), np.asarray(B, float), np.asarray(W, float)
    V = 1.0 + A + B + W
    s_ = el.get("_s", 1.0)
    return {k: (el.get(k, 0.0) + (A if k == "Ac" else 0.0) + (B * s_ if k == "Tr" else 0.0)) / V for k in CHEM_KEYS}

def run_points(st):
    """Every reading of the run as (acid mL, Tris mL per mL eluate, pH, expected scatter, experiment, water mL per mL
    eluate - v17: what came with a dilute Tris)."""
    pts = []
    ph0 = st.get("acid", {}).get("ph0")
    if ph0 is not None:
        pts.append((0.0, 0.0, ph0, 0.02, "eluate", 0.0))
    v, p, c = acid_readings(st, with_val=True)
    pts += [(x, 0.0, y, 0.025, "acid", 0.0) for x, y in zip(v, p) if PH_VALID[0] <= y <= PH_VALID[1]]
    bc = make_base_curve(st)[0] if st.get("base", {}).get("c1") else None
    if bc is not None:
        q = bc.model["q"] if bc.model else None
        for x, y, s_ in zip(bc.v, bc.p, bc.m):
            sd = float(np.sqrt(READ_SD ** 2 + (BaseCurve._sd_t(q, [x], [y], [s_])[0] / max(float(BaseCurve.dt_dph(q, y)), 1e-9)) ** 2)) \
                if q else 0.1
            pts.append((0.0, x, y, min(max(sd, 0.02), 0.3), "base", x / s_ - x))
    tr = st.get("tris") or {}
    groups = [(col["rp"], col["t"], nums("tris_c1_measured.csv")[4 * k:4 * k + 4]) for k, col in enumerate(tr.get("c1", []))]
    if tr.get("val_t"):
        groups.append((tr["val_rp"], tr["val_t"], nums("tris_val_measured.csv")[:4]))
    for rp, ts, ys in groups:
        f = rp / (1 + rp)
        te = tr.get("a2", tr.get("a", 0.0)) + tr.get("b2", tr.get("b", 1.0)) * f
        for t, y in zip(ts, ys):
            if PH_VALID[0] <= y <= PH_VALID[1] and t > 0:
                sd = math.sqrt(READ_SD ** 2 + (VOL_CV * t / (math.log(10) * max(t - te, 1e-3))) ** 2)
                pts.append((rp, t * (1 + rp), y, min(sd, 0.3), "neutralisation", 0.0))
    return pts

def fit_chemistry(st):
    """Least-squares fit of the chemistry model to every reading (robust loss, several starts). Returns a dict or None.
    v16: with VI neutralisation wells (the only wells holding both titrants in bulk) the model is also fitted with the
    strength of the Tris stock relative to the acid stock free (a weak prior at 1, SD TITRANT_SD), and that fit is kept
    only if the readings clearly ask for it - its chi-square (residuals in units of their expected scatter) at least
    TITRANT_CHI2 lower than with 1:1 stocks; otherwise the stocks stay 1:1, as v15 (in the simulator, with exactly 1:1
    stocks, a free strength lands within ~0.5% of 1, which alone moves the neutralisation entries ~0.02 pH). On 21 Oct
    the VI wells said the Tris behaved ~10% stronger than the acid (te = 0.0136 + 0.900 f; 12 Oct: 0.0004 + 1.012 f):
    with 1:1 stocks the model could not reproduce them (RMS 0.229 pH, at an implausible 15 C) and was not used in the
    map's gaps; with the strength free it fits all 44 readings (VI 0.025, acid 0.022, base 0.034 pH; 23 C, Tris x1.10)."""
    from scipy.optimize import least_squares
    pts = run_points(st)
    if len(pts) < 8:
        return None
    A = np.array([x[0] for x in pts]); B = np.array([x[1] for x in pts]); y = np.array([x[2] for x in pts])
    sd = np.array([x[3] for x in pts]); kind = [x[4] for x in pts]; W = np.array([x[5] for x in pts])
    def model(q):
        return chem_ph(chem_mix(chem_eluate(q), A, B, W), q[5])
    def fit(starts, lo, hi, prior):
        def res(q):
            r = (model(q) - y) / sd
            return np.r_[r, (q[6] - 1.0) / TITRANT_SD] if prior else r
        best = None
        for x0 in starts:
            try:
                r = least_squares(res, x0, bounds=(lo, hi), loss="soft_l1", f_scale=2.0, max_nfev=300)
            except Exception:
                continue
            if best is None or r.cost < best.cost:
                best = r
        return best
    lo, hi = [0, -0.5, 0, 0, 0, 10], [1.5, 1.5, 0.6, 0.2, 0.6, 35]
    starts = ([0.05, 0.05, 0.05, 0.005, 0.0, 21.0], [0.15, 0.12, 0.02, 0.01, 0.0, 22.0], [0.05, 0.0, 0.05, 0.005, 0.1, 21.0])
    best = fit(starts, lo, hi, False)
    strength, chi2 = None, None
    if sum(1 for k in kind if k == "neutralisation") >= 4:
        s0 = [list(best.x) + [1.0]] if best is not None else []
        free = fit(s0 + [x + [1.0] for x in starts[:2]] + [starts[0] + [v] for v in (1.08, 0.93)],
                   lo + [TITRANT_RANGE[0]], hi + [TITRANT_RANGE[1]], True)
        chi = lambda r: float(np.sum(((model(r.x) - y) / sd) ** 2)) if r is not None else float("inf")
        chi2 = (chi(best), chi(free))
        if free is not None and chi2[0] - chi2[1] >= TITRANT_CHI2:
            best, strength = free, float(free.x[6])
    if best is None:
        return None
    q = best.x
    fitted = chem_ph(chem_mix(chem_eluate(q), A, B, W), q[5])
    out = dict(q=[float(x) for x in q], n=len(pts), strength=strength, chi2=chi2)
    for k in ("acid", "base", "neutralisation"):
        m = [i for i, kk in enumerate(kind) if kk == k]
        out[f"rms_{k}"] = float(np.sqrt(np.mean((fitted[m] - y[m]) ** 2))) if m else None
        out[f"rmsw_{k}"] = float(np.sqrt(np.mean(((fitted[m] - y[m]) / sd[m]) ** 2))) if m else None
    out["rms"] = float(np.sqrt(np.mean((fitted - y) ** 2)))
    # trusted to fill the gaps between a route's readings only where it reproduces that route's own readings
    out["good"] = {"acid": out["rms_acid"] is not None and out["rms_acid"] <= CHEM_ACID_OK,
                   "base": out["rmsw_base"] is not None and out["rmsw_base"] <= CHEM_FIT_OK,
                   "neutralise": out["rmsw_neutralisation"] is not None and out["rmsw_neutralisation"] <= CHEM_FIT_OK
                   and out["rms_acid"] is not None and out["rms_acid"] <= CHEM_ACID_OK}
    return out

def chem_titrant_for(el, T, A0, B0, target, titrant, top=30.0):
    """Vectorised: mL of titrant (acid or Tris) per mL eluate to add to (A0, B0) to reach `target` (NaN if out of reach)."""
    A0, B0, target = (np.asarray(x, float) for x in (A0, B0, target))
    lo, hi = np.zeros_like(target), np.full_like(target, top)
    sign = -1.0 if titrant == "acid" else 1.0           # acid lowers the pH, Tris raises it
    ph_lo = chem_ph(chem_mix(el, A0, B0), T)
    ph_hi = chem_ph(chem_mix(el, A0 + (hi if titrant == "acid" else 0), B0 + (hi if titrant == "Tris" else 0)), T)
    reach = (sign * (target - ph_lo) >= 0) & (sign * (ph_hi - target) >= 0)
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        ph = chem_ph(chem_mix(el, A0 + (mid if titrant == "acid" else 0), B0 + (mid if titrant == "Tris" else 0)), T)
        below = sign * (ph - target) < 0                  # not there yet
        lo, hi = np.where(below, mid, lo), np.where(below, hi, mid)
    return np.where(reach, 0.5 * (lo + hi), np.nan)

# ----------------------------------------------------------------------------- the ratio map
def vi_line(st):
    """Neutralisation (VI) data: te(f) = a + b f per mL of acidified product, pKa', the acid fractions measured, and the
    base titration's equivalence (f = 0 anchor) if it ran."""
    tr = st.get("tris") or {}
    if "pk" not in tr:
        return None
    fs = [c_["rp"] / (1 + c_["rp"]) for c_ in tr.get("c1", [])]
    if tr.get("val_rp"):
        fs.append(tr["val_rp"] / (1 + tr["val_rp"]))
    line = dict(a=tr.get("a2", tr["a"]), b=tr.get("b2", tr["b"]), pk=tr.get("pk2", tr["pk"]), fmin=min(fs), fmax=max(fs), ok=tr.get("ok", 0))
    b = st.get("base", {})
    if b.get("e") is not None:
        line["e0"], line["sd_e0"] = b["e"], math.hypot(b.get("sd_e") or 0.0, ANCHOR_SD_EXTRA)
    return line

def vi_te(line, f):
    """Tris equivalence per mL of acidified product at acid fraction f, and the basis: measured range, anchored at zero
    acid by the base titration, or extended."""
    if f >= line["fmin"] - 1e-9 or "e0" not in line:
        basis = "fit" if line["fmin"] - 0.02 <= f <= line["fmax"] + 0.02 else "extended"
        return line["a"] + line["b"] * f, basis
    t_min = line["a"] + line["b"] * line["fmin"]
    return line["e0"] + (t_min - line["e0"]) * f / line["fmin"], "anchored by the base titration"

def grid(lo, hi, step):
    """Grid points from lo to hi on multiples of step (both ends included if on the grid)."""
    k0, k1 = math.ceil(lo / step - 1e-9), math.floor(hi / step + 1e-9)
    return [round(k * step, 4) for k in range(k0, k1 + 1)]

def build_map(st):
    """The run's ratio map (rows of dicts) from the validated curves; the chemistry model's value as a check column."""
    rows = []
    a = st.get("acid", {})
    ph0 = a.get("ph0")
    curve, v, p, c = make_curve(st, with_val=True)
    kept = curve.p[curve.keep]
    if ph0 is None or not len(kept):
        return rows
    val_ph = acid_params()["targets"]
    def basis_of(how, ph, targets, sure=True, gap=float("inf")):
        if any(abs(ph - t) <= 0.1 for t in targets) and sure:
            return "validated"
        if how.startswith(("robust", "Tris-equation")):
            return "fit" if sure else "extended"      # a fitted curve beyond its readings is an extrapolation
        if how == "monotone interpolation" and gap < GAP_PH:
            return "between readings"                 # v15: readings under GAP_PH apart (e.g. the Cycle-1 wells under pH 7)
        return "extended" if "extend" in how else "interpolated"
    lo = float(kept.min()) - EXTEND_PH
    for ph in grid(lo, ph0, MAP_STEP)[::-1]:
        x, how, sure = curve.ratio_for(ph)
        if x is None:
            continue
        rows.append(dict(Route="acid", Start="eluate", Start_pH=ph0, Target_pH=ph, Titrant="1 M acetic acid",
                         Per_mL_start=x, Per_mL_eluate=x, Basis=basis_of(how, ph, val_ph, sure, curve.gap_at(ph))))
    acid_r = {r["Target_pH"]: r["Per_mL_eluate"] for r in rows}
    b = st.get("base", {})
    if b.get("c2"):
        bc = make_base_curve(st)[0]
        tt = [t for t, r in zip(base_targets(), b.get("ratios", [])) if r]
        for ph in grid(ph0 + MAP_STEP / 2, bc.top(), MAP_STEP):
            x, how, sure = bc.ratio_for(ph)
            if x is None:
                continue
            rows.append(dict(Route="base", Start="eluate", Start_pH=ph0, Target_pH=ph, Titrant="1 M Tris",
                             Per_mL_start=x, Per_mL_eluate=x, Basis=basis_of(how, ph, tt, sure)))
    line = vi_line(st)
    if line:
        starts = sorted(set([s_ for s_ in grid(lo, min(ph0 - NEUT_STEP / 2, NEUT_START_MAX), NEUT_STEP)]
                            + [round(t, 4) for t in val_ph if t < ph0]))
        targets = sorted(set(grid(*NEUT_TARGETS, NEUT_STEP) + [round(t, 4) for t in base_targets()]))
        for s_ in starts:
            R_, how, sure = curve.ratio_for(s_)
            if R_ is None:
                continue
            f = R_ / (1 + R_)
            te, basis = vi_te(line, f)
            for t in targets:
                per_pool = te * (1 + 10 ** (t - line["pk"]))
                # within 1 pH of Tris's pKa' below it, the acetate still holds protons the Tris equation leaves out
                # (0.05-0.16 pH at 6.8-7.1 in the simulator): marked, and filled from the chemistry model if trusted
                b_ = basis if t >= line["pk"] - 1.0 or basis != "fit" else "near equivalence"
                rows.append(dict(Route="neutralise", Start="acidified product", Start_pH=s_, Target_pH=t, Titrant="1 M Tris",
                                 Per_mL_start=per_pool, Per_mL_eluate=per_pool * (1 + R_), Acid_per_mL_eluate=R_,
                                 Basis=b_ if line.get("ok") else "uncertain (VI fit not usable)"))
    # how far each entry is from the nearest reading it rests on (pH): acid and base paths - that path's readings and the
    # eluate itself; neutralisation - the further of (start to the nearest acid reading, target to the nearest VI reading)
    acid_ph = list(kept) + [ph0]
    base_ph = (list(make_base_curve(st)[0].p) + [ph0]) if b.get("c1") else [ph0]
    vi_ph = [y for y in nums("tris_c1_measured.csv") + nums("tris_val_measured.csv") if PH_VALID[0] <= y <= PH_VALID[1]]
    near = lambda x, pts: float(min(abs(x - y) for y in pts)) if pts else None
    for r in rows:
        if r["Route"] == "acid":
            r["Nearest_reading_pH"] = near(r["Target_pH"], acid_ph)
        elif r["Route"] == "base":
            r["Nearest_reading_pH"] = near(r["Target_pH"], base_ph)
        else:
            d1, d2 = near(r["Start_pH"], acid_ph), near(r["Target_pH"], vi_ph)
            r["Nearest_reading_pH"] = max(d1, d2) if d2 is not None else d1
    # the chemistry model's value for each row: an independent check where the readings are dense, and the map's value
    # where the curves only bridge a gap between readings (straight line, monotone interpolation, extension, the zero-acid
    # anchor) - if the model reproduces that route's own readings (CHEM_ACID_OK / CHEM_FIT_OK). In the simulator the
    # curve-only map was up to 0.3-0.9 pH off in such gaps (eluate to first acid well; eluate to the Tris region; a
    # Tris-holding eluate neutralised from mild acid), the model within 0.04.
    ch = st.get("chem")
    if ch and rows:
        el, T = chem_eluate(ch["q"]), ch["q"][5]
        for route, titrant in (("acid", "acid"), ("base", "Tris"), ("neutralise", "Tris")):
            idx = [i for i, r in enumerate(rows) if r["Route"] == route]
            if not idx:
                continue
            A0 = np.array([rows[i].get("Acid_per_mL_eluate", 0.0) for i in idx])
            tg = np.array([rows[i]["Target_pH"] for i in idx])
            x = chem_titrant_for(el, T, A0, np.zeros(len(idx)), tg, titrant)
            for i, xi in zip(idx, x):
                r = rows[i]
                r["Model_Per_mL_eluate"] = None if not np.isfinite(xi) else float(xi)
                r["Curve_Per_mL_eluate"] = r["Per_mL_eluate"]
                if np.isfinite(xi) and r["Per_mL_eluate"] > 0:
                    r["Model_vs_curve_pct"] = 100.0 * (float(xi) / r["Per_mL_eluate"] - 1)
                if r["Basis"] not in ("validated", "fit") and ch.get("good", {}).get(route) and np.isfinite(xi):
                    r["Per_mL_eluate"] = float(xi)
                    r["Per_mL_start"] = float(xi) / (1 + r.get("Acid_per_mL_eluate", 0.0)) if route == "neutralise" else float(xi)
                    r["Basis"] = ("chemistry model (between readings)" if r["Basis"] == "between readings"
                                  else f"chemistry model (gap: {r['Basis']})")
    return rows

def write_map(st, rows):
    cols = ["Route", "Start", "Start_pH", "Target_pH", "Titrant", "Per_mL_start", "Per_mL_eluate", "Basis",
            "Nearest_reading_pH", "Acid_per_mL_eluate", "Curve_Per_mL_eluate", "Model_Per_mL_eluate", "Model_vs_curve_pct"]
    df = pd.DataFrame(rows).reindex(columns=cols)
    for k in ("Per_mL_start", "Per_mL_eluate", "Acid_per_mL_eluate", "Curve_Per_mL_eluate", "Model_Per_mL_eluate"):
        df[k] = df[k].astype(float).round(5)
    for k in ("Nearest_reading_pH", "Model_vs_curve_pct"):
        df[k] = df[k].astype(float).round(2)
    df.to_csv(P("ratio_map_run.csv"), index=False)
    return df

def map_check(df, st):
    """How far the chemistry model is from the map at the validated points, as pH (log)."""
    ch = st.get("chem")
    if not ch or df.empty or "Model_Per_mL_eluate" not in df:
        return
    for route in ("acid", "base", "neutralise"):
        d = df[(df.Route == route) & (df.Basis.isin(["validated", "fit", "between readings"])) & df.Model_Per_mL_eluate.notna()]
        if len(d):
            dev = (d.Model_Per_mL_eluate / d.Per_mL_eluate - 1).abs()
            nlog(f"Ratio map {route}: chemistry model within {100 * dev.median():.1f}% (median) / {100 * dev.max():.1f}% (worst) of the fitted ratios over {len(d)} entries")
        n_gap = int((df.Route.eq(route) & df.Basis.str.startswith("chemistry model")).sum())
        if n_gap:
            nlog(f"Ratio map {route}: {n_gap} entries between readings taken from the chemistry model")

def lookup_ratio(df, ph_from, ph_to):
    """mL titrant per mL of product at ph_from to reach ph_to, from a run's map rows (DataFrame). Returns
    (ratio, titrant code 1 acid / 2 Tris, how, distance in pH to the nearest reading the entries rest on) or
    (None, 0, why, None). Same path: (ratio(to) - ratio(from)) / (1 + ratio(from)); acidified -> 6.8-8.5: the
    neutralisation table (per mL of acidified product)."""
    if df is None or df.empty or not (df.Route == "acid").any():
        return None, 0, "no ratio map", None
    ph0 = float(df.Start_pH[df.Route == "acid"].iloc[0])
    def table(route):
        d = df[df.Route == route].sort_values("Target_pH")
        return d.Target_pH.to_numpy(float), d.Per_mL_eluate.to_numpy(float), d.Nearest_reading_pH.to_numpy(float) if "Nearest_reading_pH" in d else None
    def at(route, ph):
        """(ratio per mL eluate, distance) at ph on a path; the eluate itself is 0."""
        if abs(ph - ph0) < 1e-6:
            return 0.0, 0.0
        x, y, dist = table(route)
        if not len(x) or not (x[0] - 1e-9 <= ph <= x[-1] + 1e-9):
            return None, None
        k = int(np.argmin(np.abs(x - ph)))
        return float(np.interp(ph, x, y)), (float(dist[k]) if dist is not None else None)
    def worst(*d):
        d = [x for x in d if x is not None]
        return max(d) if d else None
    if ph_to <= ph_from and ph_from <= ph0 + 1e-6:                   # down the acid path
        (r0, d0), (r1, d1) = at("acid", ph_from), at("acid", ph_to)
        if r0 is None or r1 is None:
            return None, 0, "outside the acid path of the map", None
        return (r1 - r0) / (1 + r0), 1, f"acid path: ({r1:.4f} - {r0:.4f}) / (1 + {r0:.4f})", worst(d0, d1)
    if ph_to >= ph_from and ph_from >= ph0 - 1e-6:                   # up the base path
        (r0, d0), (r1, d1) = at("base", ph_from), at("base", ph_to)
        if r0 is None or r1 is None:
            return None, 0, "outside the base path of the map", None
        return (r1 - r0) / (1 + r0), 2, f"base path: ({r1:.4f} - {r0:.4f}) / (1 + {r0:.4f})", worst(d0, d1)
    if ph_from < ph0 < ph_to:                                        # neutralise an acidified product
        d = df[df.Route == "neutralise"]
        if d.empty:
            return None, 0, "no neutralisation data in the map", None
        starts = np.unique(d.Start_pH.to_numpy(float))
        if not (starts[0] - 1e-9 <= ph_from <= starts[-1] + 1e-9):
            return None, 0, f"start pH outside the neutralisation map ({starts[0]:g}-{starts[-1]:g})", None
        k = int(np.searchsorted(starts, ph_from))
        s0, s1 = starts[max(k - 1, 0)], starts[min(k, len(starts) - 1)]
        vals, dists = [], []
        for s_ in (s0, s1):
            e = d[np.isclose(d.Start_pH, s_)].sort_values("Target_pH")
            x, y = e.Target_pH.to_numpy(float), e.Per_mL_start.to_numpy(float)
            if not (x[0] - 1e-9 <= ph_to <= x[-1] + 1e-9):
                return None, 0, f"target outside the neutralisation map ({x[0]:g}-{x[-1]:g})", None
            vals.append(float(np.exp(np.interp(ph_to, x, np.log(y)))))
            if "Nearest_reading_pH" in e:
                dists.append(float(e.Nearest_reading_pH.iloc[int(np.argmin(np.abs(x - ph_to)))]))
        w = 0.0 if s1 == s0 else (ph_from - s0) / (s1 - s0)
        return (1 - w) * vals[0] + w * vals[1], 2, f"neutralisation from pH {ph_from:g}, per mL of acidified product", worst(*dists)
    return None, 0, "this direction is not in the map (e.g. acid into a product already raised with Tris)", None

def lookup(args):
    """Ratio <from pH> <to pH> [product name]: from the current run's map, or the latest run of that product in the
    archive. Writes ratio_lookup.csv: mL titrant per mL of product at the start pH, titrant (1 acid, 2 Tris, 0 none),
    status (1 found, 0 not in the map), and how far (pH) the entry is from the nearest reading it rests on."""
    try:
        ph_from, ph_to = float(args[0]), float(args[1])
    except Exception:
        write_lines(P("ratio_lookup.csv"), [0, 0, 0, 0], fmt="{:g}")
        return nlog("ERROR: Ratio needs a start pH and a target pH")
    product = " ".join(args[2:]).strip()
    df = None
    if product:
        arch = os.environ.get("TECAN_ARCHIVE_DIR", os.path.join(os.path.dirname(WORKING_DIR), "ML_data_archive"))
        try:
            m = pd.read_csv(os.path.join(arch, "ratio_map.csv"))
            m = m[m.Product_Name.astype(str).str.strip().str.lower() == product.lower()]
            if len(m):
                df = m[m.Run_ID == m.Run_ID.iloc[-1]]
        except Exception:
            df = None
    else:
        try:
            df = pd.read_csv(P("ratio_map_run.csv"))
        except Exception:
            df = None
    x, code, how, dist = lookup_ratio(df, ph_from, ph_to)
    write_lines(P("ratio_lookup.csv"), [x if x is not None else 0.0, code, 1 if x is not None else 0, dist or 0.0], fmt="{:.5f}")
    nlog(f"Ratio pH {ph_from:g} -> {ph_to:g}" + (f" ({product})" if product else "") + ": "
         + (f"{x:.5f} mL {'acid' if code == 1 else 'Tris'} per mL - {how}; nearest reading {dist if dist is not None else float('nan'):.2f} pH away"
            if x is not None else f"not found - {how}"))

# ----------------------------------------------------------------------------- pools
def pool_ranges(st):
    """pH ranges custom pools can be made in, from this run's stored results:
    acid - from the lowest kept acid reading minus EXTEND_PH up to the eluate's own pH; most accurate across the pH the
    dense readings (Cycle 2 and validation) span; Tris - from the eluate's own pH up to the highest kept base reading
    plus EXTEND_PH, only if the base titration ran; most accurate across its Cycle 2 and validation readings."""
    curve, v, p, c = make_curve(st, with_val=True)
    kept = curve.p[curve.keep]
    ph0 = st.get("acid", {}).get("ph0")
    hi = ph0 if ph0 is not None else (float(kept.max()) if len(kept) else 0.0)
    lo = float(kept.min()) - EXTEND_PH if len(kept) else 0.0
    dense = [y for y, cyc in zip(p, c) if cyc >= 2 and PH_VALID[0] <= y <= PH_VALID[1]]
    dlo, dhi = (min(dense), max(dense)) if dense else (lo, hi)
    blo = bhi = bdlo = bdhi = 0.0
    if st.get("base", {}).get("c2"):
        bc, bv, bp, bcyc = make_base_curve(st, with_val=True)
        if bc.top() > hi + POOL_TOL:
            blo, bhi = hi, bc.top()
            bd = [y for y, cyc in zip(bp, bcyc) if cyc >= 2 and y > hi and PH_VALID[0] <= y <= PH_VALID[1]]
            bdlo, bdhi = (min(bd), max(bd)) if bd else (blo, bhi)
    return dict(lo=lo, hi=hi, dlo=dlo, dhi=dhi, blo=blo, bhi=bhi, bdlo=bdlo, bdhi=bdhi)

def ranges():
    """After every ratio is found: the chemistry cross-check, the ratio map, and the numbers for the messages before
    the pool request."""
    st = load_state()
    try:
        ch = fit_chemistry(st)
    except Exception as e:
        ch = None
        nlog(f"Chemistry model fit failed ({e!r}) - the ratio map is written without its check column")
    if ch:
        st["chem"] = ch
        q = ch["q"]
        nlog(f"Chemistry model: acetate {q[0]:.4f} M, strong-ion difference {q[1]:+.4f} M, carboxyls {q[2]:.4f} M, histidine-like "
             f"{q[3]:.4f} M, Tris {q[4]:.4f} M, {q[5]:.1f} C"
             + (f", Tris stock x{ch['strength']:.3f} the strength of the acid stock (the VI wells reject 1:1 stocks: chi-square "
                f"{ch['chi2'][0]:.0f} -> {ch['chi2'][1]:.0f})" if ch.get("strength") is not None else
                (f", 1:1 stocks (a free Tris strength would gain chi-square {ch['chi2'][0] - ch['chi2'][1]:.1f})" if ch.get("chi2") else ""))
             + f"; RMS {ch['rms']:.3f} pH over {ch['n']} readings ("
             + ", ".join(f"{k} {ch['rms_' + k]:.3f}" for k in ("acid", "base", "neutralisation") if ch.get('rms_' + k) is not None) + ")")
        if ch.get("strength") is not None and abs(ch["strength"] - 1.0) >= TITRANT_FLAG:
            nlog(f"WARNING (titrants): in the VI wells the 1 M Tris behaves as {ch['strength']:.2f} x the strength of the 1 M acetic "
                 f"acid ({100 * (ch['strength'] - 1):+.0f}%): one of the stocks is likely not 1 M. This run's ratios are measured with "
                 f"these stocks and hold for them; ratios stored from runs with other stocks differ by about as much")
    try:
        rows = build_map(st)
        df = write_map(st, rows)
        map_check(df, st)
        nlog(f"Ratio map: {len(df)} entries (" + ", ".join(f"{r} {int((df.Route == r).sum())}" for r in ("acid", "base", "neutralise")) + ")")
    except Exception as e:
        nlog(f"ERROR: ratio map not written ({e!r})")
    save_state(st)
    summary()
    g = pool_ranges(st)
    st["ranges"] = g
    save_state(st)
    write_lines(P("pool_range.csv"), [g["lo"], g["hi"], g["dlo"], g["dhi"], g["blo"], g["bhi"], g["bdlo"], g["bdhi"]], fmt="{:.2f}")
    nlog(f"Custom pools possible: pH {g['lo']:.2f}-{g['hi']:.2f} with acid (dense readings {g['dlo']:.2f}-{g['dhi']:.2f})"
         + (f", {g['blo']:.2f}-{g['bhi']:.2f} with Tris (dense readings {g['bdlo']:.2f}-{g['bdhi']:.2f})" if g["bhi"]
            else ", none with Tris (no base titration)"))

def model_ratio(st, route, ph):
    """The chemistry model's ratio from the eluate (acid or Tris per mL eluate), if trusted for that route; else None."""
    ch = st.get("chem")
    if not ch or not ch.get("good", {}).get(route):
        return None
    x = chem_titrant_for(chem_eluate(ch["q"]), ch["q"][5], np.zeros(1), np.zeros(1), np.array([ph]), "acid" if route == "acid" else "Tris")[0]
    return float(x) if np.isfinite(x) else None

def in_gap(how, sure=True):
    """Is a curve value only bridging readings (not a fit through them)? Then the chemistry model plans the pool if it
    is trusted for the route - also between close readings (v15 bench: the model 2-3x closer there than the monotone
    curve through noisy readings)."""
    return not sure or not how.startswith(("robust", "Tris-equation"))

def pool_plan(pl, acid_curve, base_curve, ph0, st=None):
    """First make of one pool: (acid, eluate, Tris) mL and the ratio - or None with the reason. A target between the
    readings (where the curve only bridges a gap) is planned from the chemistry model if it is trusted for the route."""
    T = pl["target"]
    if pl["route"] == "acid":
        x, how, sure = acid_curve.ratio_for(T)
        if x is not None and in_gap(how, sure) and st is not None:
            xm = model_ratio(st, "acid", T)
            if xm is not None:
                x, how = xm, f"chemistry model (between readings; the curve gave {x:.4f})"
        if x is None and ph0 is not None and T >= ph0 - 1e-9:
            x, how = 0.0, "at the eluate's own pH: eluate only"
        if x is None:
            return None, how
        e = POOL_ML / (1 + x)
        A = POOL_ML - e if POOL_ML - e >= MIN_TRIM_ML else 0.0
        return ((A, POOL_ML - A, 0.0), x, how), ""
    if not base_curve:
        return None, "no base titration this run"
    # A Tris pool that comes out high is brought down with eluate (mL volumes, precise); one that comes out low needs a
    # few uL more Tris for a lightly buffered eluate - less than 8 tips can add. So the first make is aimed above the
    # target by about 2 SD of the curve's uncertainty there, less the tolerance: 0 for a well-buffered eluate, up to
    # BASE_AIM_MAX for one where pH 7.2 is only a few uL of Tris past the eluate's equivalence point.
    sd = base_curve.sd_ph(T)
    aim = min(max(2 * sd - POOL_TOL, 0.0), BASE_AIM_MAX) if sd is not None else 0.0
    x, how, sure = base_curve.ratio_for(T + aim)
    if x is None and aim:
        aim = 0.0
        x, how, sure = base_curve.ratio_for(T)
    if x is not None and in_gap(how, sure) and st is not None:      # below the Tris region, between the readings
        xm = model_ratio(st, "base", T)
        if xm is not None:
            x, how, aim = xm, f"chemistry model (between readings; the curve gave {x:.4f})", 0.0
    if x is None:
        return None, how
    pl["aim"] = aim
    if aim:
        how += f", aimed {aim:+.2f} above target (curve +-{sd:.2f} there)"
    e = POOL_ML / (1 + x)
    return ((0.0, e, POOL_ML - e), x, how), ""

def pools0():
    """The pool request (lanes and target pHs), checked against the ranges, and the first make of every pool that
    fits in the troughs (the rest are made with the round-2 block). Up to the eluate's own pH a pool is eluate + acid;
    above it, eluate + Tris."""
    st = load_state()
    g = st.get("ranges") or pool_ranges(st)
    st["pools"], st["measure"] = [], {}
    bad_lane, bad_ph, ok = 0, 0.0, True
    try:
        lanes = [int(round(x)) for x in parse_numbers(read_text("pool_lanes.csv"))]
        phs = parse_numbers(read_text("pool_phs.csv"))
    except ValueError:
        lanes, phs, ok = [], [], False
    if ok and lanes and lanes != [0]:
        if any(l < 1 or l > LANES for l in lanes) or len(set(lanes)) != len(lanes) or len(phs) != len(lanes):
            ok = False
        else:
            for l, t in zip(lanes, phs):
                acid = g["lo"] <= t <= g["hi"] + POOL_TOL
                base = g["bhi"] > 0 and g["blo"] < t <= g["bhi"]
                if not (acid or base):
                    ok, bad_lane, bad_ph = False, l, t
                    break
    elif not ok or (lanes != [0] and lanes):
        ok = False
    if not ok:
        write_lines(P("pool_status.csv"), [0, 0, bad_lane, bad_ph, 0, 0], fmt="{:g}")
        nlog(f"Pool request not used (lanes {read_text('pool_lanes.csv')!r}, pH {read_text('pool_phs.csv')!r})"
             + (f": lane {bad_lane} pH {bad_ph:g} is outside the ranges" if bad_lane else "") + " - asked again")
        save_state(st)
        return
    if not lanes or lanes == [0]:
        write_lines(P("pool_status.csv"), [1, 0, 0, 0, 0, 0], fmt="{:g}")
        nlog("No custom pools requested")
        save_state(st)
        return
    ph0 = st.get("acid", {}).get("ph0")
    for l, t in zip(lanes, phs):
        base = g["bhi"] > 0 and t > g["hi"] + POOL_TOL
        st["pools"].append(dict(lane=l, target=float(t), route="base" if base else "acid", A=0.0, E=0.0, B=0.0,
                                readings=[], trims=0, status="requested", note=""))
    base_curve = make_base_curve(st)[0] if st.get("base", {}).get("c2") else None
    curve = make_curve(st, with_val=True)[0]
    curve = AcidCurve(curve.v, curve.p, curve.ph0, acid_targets(st))   # acid pool pHs get their own fit clusters
    make1, order, used = {}, [], [0.0, 0.0, 0.0]
    for pl in st["pools"]:
        plan, why = pool_plan(pl, curve, base_curve, ph0, st)
        if plan is None:
            pl["route"], pl["status"], pl["note"] = "refused", "refused", why
            nlog(f"Pool lane {pl['lane']} pH {pl['target']:g}: not made - {why}")
            continue
        (A, E, B), pl["planned_ratio"], pl["how"] = plan
        pl["plan"] = [A, E, B]
        if used[0] + E > TROUGH_ML or used[1] + A > TROUGH_ML or used[2] + B > TROUGH_ML:
            pl["status"] = "deferred"                                       # made with the round-2 block
            nlog(f"Pool lane {pl['lane']} pH {pl['target']:g}: deferred to the next block (trough volume)")
            continue
        used = [used[0] + E, used[1] + A, used[2] + B]
        pl.update(A=A, E=E, B=B, status="made")
        make1[pl["lane"]] = (A, E, B); order.append(pl["lane"])
        nlog(f"Pool lane {pl['lane']} pH {pl['target']:g} ({'eluate + Tris' if pl['route'] == 'base' else 'eluate + acid'}): "
             f"ratio {pl['planned_ratio']:.4f} ({pl['how']}) -> {A:.2f} mL acid + {E:.2f} mL eluate + {B:.3f} mL Tris")
    write_lines(P("pool_make1.csv"), lane_rows(make1), fmt="{:g}")
    write_lines(P("needs_pools.csv"), used, fmt="{:.1f}")
    st["measure"] = {"1": order}
    write_measure(order)
    write_lines(P("pool_status.csv"), [1, len(st["pools"]), 0, 0, len(order), len(order)], fmt="{:g}")
    pool_tris(st)
    save_state(st)

def pool_tris(st):
    """Each low-pH pool's own Tris need (mL 1 M Tris per mL pool) at the Tris targets, from its composition."""
    tr = st.get("tris") or {}
    if not tr.get("ok"):
        return
    pk, a_, b_ = tr.get("pk2", tr["pk"]), tr.get("a2", tr["a"]), tr.get("b2", tr["b"])
    targets = tris_params()["targets"]
    for pl in st.get("pools", []):
        if pl["route"] == "acid" and pl["E"] > 0 and pl["target"] <= TRIS_POOL_MAX:
            f = pl["A"] / (pl["A"] + pl["E"])
            pl["tris_ratios"] = [(a_ + b_ * f) * (1 + 10 ** (t - pk)) for t in targets]

def write_measure(order):
    write_lines(P("pools_expected.csv"), [len(order)], fmt="{}")
    write_lines(P("pools_order.csv"), list(order) + [0] * (LANES - len(order)), fmt="{}")

def trim_for(pl, m, curve, base_curve):
    """Acid, eluate and Tris (mL) that move a pool from its measured pH m towards its target, and how it was worked out.
    First trim: along the run's curve - its ratios at the measured and the target pH give the shift, so an error in the
    pool's own make does not enter - applying TRIM_GAIN of it, or BASE_TRIM_FIRST for a Tris pool coming down a steep
    curve (where the curve's slope is uncertain, so the pool stays above target and eluate can finish it).
    Later trims of a Tris pool (v14): the run's curve again, its shift from the pool's last reading to the target added
    to the pool's own ratio, applying BASE_TRIM_LATER. v13 used a straight line through the pool's last two readings
    (a secant), which overshoots near the equivalence point, where the curve steepens fast: on the weakest modelled
    eluate it left 1 pool at pH 7.0 in 10 more than 0.05 off (worst 0.47), the curve 1 in 30 (worst 0.13, a pool that
    came out low and could not take less Tris than 8 tips add). Later trims of an acid pool: the secant, applying 0.9.
    Titrant goes in for a shift towards more titrant; eluate for a shift back. A Tris pool that needs less Tris than 8
    tips can add (MIN_TRIM_ML) gets that smallest addition plus the eluate that brings it back to the right ratio, if the
    lane has room. A correction of 60-100% of the smallest 8-tip addition is rounded up to it (it still moves closer).
    v17: with the dilute Tris on the deck, a Tris addition of up to DIL_TRIM_MAX_ML of it comes from it instead (10x the
    volume for the same Tris; 8 tips can add 1/10 of the old smallest addition): returned as dd, the mL of dilute Tris (db, as always,
    in mL 1 M Tris - mmol); with realistic dispensing (1.4 uL per tip) a weak eluate's Tris pool at pH 7.0 could come out
    0.19 low and be left there. Returns (acid, eluate, Tris mmol, how, dilute Tris mL)."""
    T, E = pl["target"], pl["E"]
    room = max(LANE_MAX_ML - pl["A"] - pl["E"] - pl["B"] - pl.get("W", 0.0), 0.0)
    dil = dil_m()
    base = pl["route"] == "base"
    sign = 1.0 if base else -1.0                       # d(ratio)/d(pH): Tris raises the pH, acid lowers it
    hist = pl.get("hist", [])
    how = ""
    if len(hist) >= 2 and base and base_curve is not None:
        u_m, u_t = base_curve.ratio_for(m)[0], base_curve.ratio_for(T)[0]
        if u_m is not None and u_t:
            r_now = hist[-1][0]
            r_aim = r_now + BASE_TRIM_LATER * (u_t - u_m)
            how = f"run's curve from the pool's last reading, {int(round(100 * BASE_TRIM_LATER))}% of the shift"
    if not how and len(hist) >= 2 and abs(hist[-1][1] - hist[-2][1]) >= 0.05 and \
            sign * (hist[-1][0] - hist[-2][0]) / (hist[-1][1] - hist[-2][1]) > 0:
        (r1, m1), (r2, m2) = hist[-2], hist[-1]
        slope = (r2 - r1) / (m2 - m1)
        r_now, r_aim = r2, r2 + 0.9 * slope * (T - m)
        how = "secant through the pool's own readings"
    if not how:
        c = base_curve if base else curve
        if c is None:
            return 0.0, 0.0, 0.0, "no curve", 0.0
        r_now, _, _ = c.ratio_for(m)
        r_t, _, _ = c.ratio_for(T)
        if not r_t:
            return 0.0, 0.0, 0.0, "no ratio for the target", 0.0
        if r_now is None:
            r_now = 0.0 if (m > T) != base else r_t * 1.2
        sd = base_curve.sd_ph(T) if base and base_curve else None
        gain = BASE_TRIM_FIRST if base and m > T and sd is not None and sd > POOL_TOL else TRIM_GAIN
        r_aim = r_now + gain * (r_t - r_now)
        how = f"run's curve, {int(round(100 * gain))}% of the shift"
    if r_aim <= 0:
        return 0.0, 0.0, 0.0, how, 0.0
    da = de = db = 0.0
    if r_aim > r_now:                                   # more titrant
        if base:
            db = (r_aim - r_now) * E
            if dil > 0 and db / dil <= DIL_TRIM_MAX_ML:  # v17: from the dilute Tris - 8 tips can add 1/10 as much Tris
                dd = db / dil
                if dd >= MIN_TRIM_ML:
                    de = 0.0
                else:                                   # the smallest 8-tip addition, and the eluate that brings it back
                    dd = MIN_TRIM_ML                    # to the ratio aimed at (0.4 mL eluate per 8 uL of Tris at
                    db = dd * dil                       # 0.013 mL/mL: room for it, unlike 80 uL of 1 M Tris)
                    de = max((r_now * E + db) / r_aim - E, 0.0)
                if dd + de > room:
                    return 0.0, 0.0, 0.0, how, 0.0      # the lane is full
                return 0.0, de, dd * dil, how + f", from the {dil:g} M Tris" + (" with eluate" if de else ""), dd
            if db < 0.6 * MIN_TRIM_ML:
                db = MIN_TRIM_ML
                de = (r_now * E + db) / r_aim - E
        else:
            da = (r_aim - r_now) * E
    else:                                               # back towards the eluate
        de = (r_now / r_aim - 1) * E
    da, de, db = (max(x, 0.0) for x in (da, de, db))
    total = da + de + db
    if 0.6 * MIN_TRIM_ML <= total < MIN_TRIM_ML:
        da, de, db = (x * MIN_TRIM_ML / total for x in (da, de, db))
        total = MIN_TRIM_ML
    if total > room:
        if base and db and de:
            return 0.0, 0.0, 0.0, how, 0.0              # the eluate to dilute the smallest Tris addition does not fit
        da, de, db = (x * room / total for x in (da, de, db))
    return da, de, db, how, 0.0

def lane_check(pl, m, curve, base_curve):
    """First reading of a pool against what its recipe gives on the run's curve (the wells). 12 Oct: a lane pool read
    7.51 where its recipe gives 7.0 by every model of the run - as if it held 12% more Tris than planned (or was not
    mixed when read). A large gap is flagged in the log and pools.csv; the trims work from the pool's own readings."""
    c = base_curve if pl["route"] == "base" else curve
    r_nom = (pl["B"] if pl["route"] == "base" else pl["A"]) / pl["E"] if pl["E"] else 0.0
    if c is None or r_nom <= 0:
        return
    r_meas, _, _ = c.ratio_for(m)
    planned = pl["target"] + pl.get("aim", 0.0)
    if r_meas is None:
        return
    pl["lane_factor"] = r_meas / r_nom
    if abs(m - planned) > LANE_FLAG:
        pl["note"] = (f"first reading {m:.2f} where the recipe gives {planned:.2f} on the run's curve - as if the lane held "
                      f"{100 * (pl['lane_factor'] - 1):+.0f}% {'Tris' if pl['route'] == 'base' else 'acid'}; check mixing (stir and re-read) and trough levels")
        nlog(f"Pool lane {pl['lane']}: LANE CHECK - {pl['note']}")

def pools_round(k):
    st = load_state()
    order = st.get("measure", {}).get(str(k), [])
    vals = nums((f"pools_measured_{k}.csv")) if order else []
    curve = make_curve(st, with_val=True)[0]
    base_curve = make_base_curve(st)[0] if st.get("base", {}).get("c2") else None
    by_lane = {pl["lane"]: pl for pl in st.get("pools", [])}
    trims, makes, nxt = {}, {}, []
    for lane, m in zip(order, vals):
        pl = by_lane[lane]
        pl["readings"].append(m)
        pl.setdefault("hist", []).append(((pl["B"] if pl["route"] == "base" else pl["A"]) / pl["E"] if pl["E"] else 0.0, m))
        pl.setdefault("vol", []).append((pl["A"] + pl["E"] + pl["B"] + pl.get("W", 0.0)) / pl["E"] if pl["E"] else 0.0)   # v17
        if not (PH_VALID[0] <= m <= PH_VALID[1]):
            pl["status"] = "not measured"
            nlog(f"Pool lane {lane}: no reading - left as it is"); continue
        T = pl["target"]
        if len(pl["readings"]) == 1:
            lane_check(pl, m, curve, base_curve)
        if abs(m - T) <= POOL_TOL + 1e-9 or k >= 3:          # (7.03 - 7.00 is 0.0300000000000002 in floating point)
            pl["status"] = "done"
            nlog(f"Pool lane {lane}: pH {m:.2f} (target {T:g}, error {m - T:+.2f}) - " + ("within tolerance" if abs(m - T) <= POOL_TOL + 1e-9 else "last round"))
            continue
        da, de, db, how, dd = trim_for(pl, m, curve, base_curve)
        if da + de + (dd if dd else db) < MIN_TRIM_ML - 1e-9:
            pl["status"] = "done"
            why = ("the Tris it needs is less than 8 tips can add, and the lane has no room for the eluate to dilute the "
                   "smallest addition") if pl["route"] == "base" and m < T and not dil_m() else "correction below what 8 tips can add, or lane full"
            pl["note"] = why
            nlog(f"Pool lane {lane}: pH {m:.2f} (target {T:g}) - {why}"); continue
        pl["A"] += da; pl["E"] += de; pl["B"] += db; pl["trims"] += 1
        if dd:                                          # v17: the dilute Tris - its water, and the mL added
            pl["W"] = pl.get("W", 0.0) + dd - db
            pl["Bd"] = pl.get("Bd", 0.0) + dd
        trims[lane] = (da, de, dd, 2) if dd else (da, de, db); nxt.append(lane)
        nlog(f"Pool lane {lane}: pH {m:.2f} (target {T:g}, error {m - T:+.2f}) - trim {da:.3f} mL acid + {de:.3f} mL eluate + "
             + (f"{dd:.3f} mL {dil_m():g} M Tris ({db:.4f} mmol)" if dd else f"{db:.3f} mL Tris") + f" ({how})")
    if k == 1:                                        # pools that did not fit in the troughs in the first block
        used = [0.0, 0.0, 0.0]
        for pl in st.get("pools", []):
            if pl["status"] != "deferred":
                continue
            A, E, B = pl["plan"]
            if used[0] + E > TROUGH_ML or used[1] + A > TROUGH_ML or used[2] + B > TROUGH_ML:
                pl["status"], pl["note"] = "refused", "the troughs cannot hold enough for it in this run"
                nlog(f"Pool lane {pl['lane']}: not made - {pl['note']}"); continue
            used = [used[0] + E, used[1] + A, used[2] + B]
            pl.update(A=A, E=E, B=B, status="made")
            makes[pl["lane"]] = (A, E, B); nxt.append(pl["lane"])
            nlog(f"Pool lane {pl['lane']} pH {pl['target']:g}: made now - {A:.2f} mL acid + {E:.2f} mL eluate + {B:.3f} mL Tris")
        write_lines(P("pool_make2.csv"), lane_rows(makes), fmt="{:g}")
    nxt = sorted(set(nxt))
    write_lines(P(f"pool_trim{k}.csv"), lane_rows(trims), fmt="{:g}")
    need = [sum(x[1] for x in list(trims.values()) + list(makes.values())),
            sum(x[0] for x in list(trims.values()) + list(makes.values())),
            sum(x[2] for x in list(trims.values()) + list(makes.values()) if len(x) < 4)]   # v17: 1 M Tris only
    write_lines(P(f"needs_round{k + 1}.csv"), need, fmt="{:.1f}")
    work = 1 if (trims or makes) else 0
    if k >= 3:
        work, nxt = 0, []
    st.setdefault("measure", {})[str(k + 1)] = nxt
    write_measure(nxt)
    write_lines(P(f"round{k + 1}.csv"), [work, len(nxt), len(makes)], fmt="{}")
    pool_tris(st)
    save_state(st)
    write_pool_final(st)

def write_pool_final(st):
    rows = []
    by_lane = {pl["lane"]: pl for pl in st.get("pools", [])}
    for lane in range(1, LANES + 1):
        pl = by_lane.get(lane)
        last = pl["readings"][-1] if pl and pl["readings"] else 0.0
        rows += [last if pl and pl["route"] != "refused" else (-1.0 if pl else 0.0), pl["target"] if pl else 0.0]
    write_lines(P("pool_final.csv"), rows, fmt="{:.2f}")

def summary():
    st = load_state()
    ap = acid_params()
    drift = st.get("acid", {}).get("drift", 0.0) or 0.0
    val = [x - drift for x in nums(("validation_results.csv"))]
    errs = [abs(m - t) for m, t in zip(val, ap["targets"]) if PH_VALID[0] <= m <= PH_VALID[1]]
    corr = nums(("corr_ratios.csv")) or [0.0] * 4
    tv = nums(("tris_val_measured.csv"))
    tt = (nums("tris_params.csv") or [0] * 4)[:4]
    terrs = [abs(m - t) for m, t in zip(tv, tt) if PH_VALID[0] <= m <= PH_VALID[1]]
    tc = nums(("tris_corr.csv")) or [-1.0] * 4                    # -1 = no VI Tris titration this run
    b = st.get("base", {})
    bv = nums("base_val_measured.csv")
    berrs = [abs(m - t) for m, t, r in zip(bv, base_targets(), b.get("ratios", [])) if r and PH_VALID[0] <= m <= PH_VALID[1]]
    bc = (nums("base_corr.csv") or [-1.0] * 4) if b.get("corr") else [-1.0] * 4   # 0 = at or below the eluate's pH
    # v14: the neutralisation numbers in the messages are the stored ones - the ratio map's entries (Tris-equation fit,
    # or the chemistry model where the map takes it: near the equivalence point, or when the VI fit was not usable) -
    # mL Tris per mL of acidified product, from each acid target to nb_proc, and from the process pool to the Tris
    # targets (in place of the Tris refit's own numbers). Without a map: the VI line, as v13.
    neut = [-1.0] * 4
    line = vi_line(st)
    acid_t = ap["targets"]
    try:
        mp = pd.read_csv(P("ratio_map_run.csv"))
        mp = mp[mp.Route == "neutralise"]
    except Exception:
        mp = None
    def stored(start, target):
        if mp is None or mp.empty:
            return None
        d = mp[((mp.Start_pH - start).abs() < 1e-6) & ((mp.Target_pH - target).abs() < 1e-6)]
        return float(d.Per_mL_start.iloc[0]) if len(d) and np.isfinite(d.Per_mL_start.iloc[0]) else None
    if len(corr) >= 4 and tt:
        for k, R_ in enumerate(corr[:4]):
            if R_ > 0:
                x = stored(acid_t[k], tt[3])
                if x is None and line:
                    te, _ = vi_te(line, R_ / (1 + R_))
                    x = te * (1 + 10 ** (tt[3] - line["pk"]))
                if x is not None:
                    neut[k] = x
        if tc and tc[0] != -1.0:
            tc = [stored(acid_t[3], t) or x for t, x in zip(tt, tc[:4])]
    ph0 = st.get("acid", {}).get("ph0") or 0.0
    write_lines(P("run_summary.csv"), [max(errs) if errs else -1, *corr[:4], max(terrs) if terrs else -1, *tc[:4],
                                       max(berrs) if berrs else -1, *bc[:4], *neut, ph0, drift], fmt="{:.4f}")
    write_pool_final(st)
    write_lines(P("pools_any.csv"), [len(st.get("pools", []))], fmt="{}")

# ----------------------------------------------------------------------------- command line
def main():
    mode = sys.argv[1].strip().lower() if len(sys.argv) > 1 else ""
    stage = sys.argv[2].strip().lower() if len(sys.argv) > 2 else ""
    with open(P("algorithm_version.csv"), "w") as f:
        f.write(VERSION + "\n")
    nlog(f"=== {datetime.now():%Y-%m-%d %H:%M:%S} {VERSION} {mode} {stage} ===")
    if mode == "ratio":
        return lookup(sys.argv[2:])
    fn = {("check", ""): check, ("acid", "0"): acid0, ("acid", "1"): acid1, ("acid", "2"): acid2, ("acid", "3"): acid3,
          ("tris", "1"): tris1, ("tris", "2"): tris2, ("summary", ""): summary,
          ("base", "0"): base0, ("base", "1"): base1, ("base", "2"): base2, ("base", "3"): base3,
          ("pools", "range"): ranges, ("pools", "0"): pools0,
          ("pools", "1"): lambda: pools_round(1), ("pools", "2"): lambda: pools_round(2), ("pools", "3"): lambda: pools_round(3)}.get((mode, stage))
    if fn is None:
        return nlog(f"ERROR: unknown command {' '.join(sys.argv[1:])}")
    fn()

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        import traceback
        try:
            nlog(f"ERROR: {e!r}\n{traceback.format_exc()}")
            if len(sys.argv) > 1 and sys.argv[1].lower() == "acid":
                write_lines(P("status.csv"), [0], fmt="{}")
            if len(sys.argv) > 2 and sys.argv[1].lower() == "base" and sys.argv[2] in ("1", "2"):
                write_lines(P("base_status.csv"), [0], fmt="{}")     # the robot skips the rest of the base titration
            if len(sys.argv) > 2 and sys.argv[1].lower() == "pools" and sys.argv[2] == "0":
                write_lines(P("pool_status.csv"), [0, 0, 0, 0, 0, 0], fmt="{}")
            for name in ("round2.csv", "round3.csv"):
                if len(sys.argv) > 1 and sys.argv[1].lower() == "pools":
                    write_lines(P(name), [0, 0, 0], fmt="{}")
        except Exception:
            pass
