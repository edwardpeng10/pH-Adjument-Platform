import os
import sys
import math
from datetime import datetime, timedelta
# EDP_VERSION = "EDP-v17-2026-10"   (EDP_tecan_gp_master.py Check: all six EDP files must carry the same version)

if sys.stdout is None: sys.stdout = open(os.devnull, "w")      # started with pythonw: there is no console
if sys.stderr is None: sys.stderr = open(os.devnull, "w")

# =====================================================================================
# EDP time estimates, v17: how long the robot works until it next needs someone, for the messages.
#   EDP_time_estimate.py plan <block>    count the robot work in <block> from the dispense plans -> time_est.csv
#   EDP_time_estimate.py start <block>   the operator clicked OK: the block starts (logged in timing_log.csv)
#   EDP_time_estimate.py end <block>     the block finished: log the actual time and refine the calibration
# Blocks (v13): acid_c1 (acid Cycle 1 + base Cycle 1), acid_c2 (acid Cycle 2 + base Cycle 2), val (acid validation + VI
#   Tris titration Cycle 1 + base validation), tris_val (VI Tris validation), pools (the custom pools, after the pool
#   request), round2 (trims and deferred pools), round3 (second trims).
# time_est.csv, one number per line, read by EVOware: minutes of robot work; clock time it ends, as HH.MM (always
#   rounded to a minute ending in 5, so it reads the same whether EVOware shows 14.05 or 14.050000); 1 = long enough
#   to step away, 0 = stay nearby; minutes of the short acid Cycle 2 block (shown at plate check 1).
# The work is counted from the same plan files the robot reads: per plate column the portions of each liquid, the
#   washes and the 5 x 900 uL mix; per bulk pool the 8-tip portions, the washes and the 5 + 3 mix cycles. v15 washes:
#   a block's product goes into all its wells (lanes) first, then one wash; then per column (pool) its acid, a wash if
#   Tris follows, its Tris, the mix and one wash. v16: with Use_8_tips (options.csv, third value) every two portions of
#   a liquid are one 8-tip trip (tips 5-8 put the second portion into the same wells), an odd last one a 4-tip trip.
#   v17: the base plans (plan_base_*.csv) have 20 numbers per column - a 4th liquid, the dilute Tris (Product 2W RIGHT
#   well), dispensed after any 1 M Tris with a wash between them.
#   The seconds per operation below are first estimates;
#   every finished block compares its real duration with the count, and the median of the last 12 blocks
#   (timing_calibration.csv, kept between runs) scales all estimates. Only the standard library: starts in < 1 s.
# =====================================================================================
DIR = os.environ.get("TECAN_DATA_DIR", r"C:\ProgramData\Tecan\EVOware\database\tecan_data")
OUT = os.path.join(DIR, "time_est.csv")
LOG = os.path.join(DIR, "timing_log.csv")
CAL = os.path.join(DIR, "timing_calibration.csv")      # kept between runs (not cleared by EDP_tecan_reset.py)

T_PORTION = 22.0      # s, one 4-tip aspirate (liquid detection) + dispense of up to 980 uL
T_WASH = 20.0         # s, one wash (3 mL waste + 4 mL cleaner)
T_MIX4 = 55.0         # s, 5 x 900 uL mix in one plate column (4 tips)
T_COLUMN = 10.0       # s, per plate column: moves, plan imports, loop logic
T_PORTION8 = 25.0     # s, one 8-tip aspirate from a trough + dispense into an LpH lane
T_PORTION2 = 30.0     # s, v16: one 8-tip aspirate + tips 5-8 and tips 1-4 dispensing into the same plate column
T_MIXCYCLE8 = 11.0    # s, one 980 uL mix cycle with 8 tips in an LpH lane
T_MIX8 = 8.0          # s, per 8-tip mix command (moves, liquid detection)
T_MARKER = 2.0        # s, the timing call that ends a block
POOL_UL = 30000.0     # Pool_Volume in the EVOware script (uL per LpH lane)
PORTION_UL = 980.0
AWAY_MIN = 15         # minutes: at least this long and the message says the operator can step away
POOLS = ("P32", "P34", "P36")
POOL_RATIO_ROW = {"P32": 3, "P34": 2, "P36": 1}       # ratios_clean.csv rows: pH 3.6, 3.4, 3.2 -> Pool 3, 2, 1
TYPICAL_ACID_C2 = (2, 3, 0)                           # portions of eluate, acid, Tris in a typical acid Cycle 2 column


def numbers(name):
    try:
        with open(os.path.join(DIR, name)) as f:
            out = []
            for tok in f.read().replace(",", " ").replace(";", " ").split():
                try:
                    out.append(float(tok.strip('"')))
                except ValueError:
                    pass
            return out
    except Exception:
        return []

def eight_tips():
    """v16: Use_8_tips from options.csv (neut, base_on, Use_8_tips); 1 when the file says nothing."""
    o = numbers("options.csv")
    return not (len(o) > 2 and o[2] == 0)

def trips_seconds(n):
    """n portions of one liquid into one plate column: v16 - two per 8-tip trip, an odd last one with tips 1-4."""
    n = int(n)
    return (n // 2) * T_PORTION2 + (n % 2) * T_PORTION if eight_tips() else n * T_PORTION

def column_seconds(ne, na, nb, nd=0):
    """One plate column: its eluate portions (in the block's product part), its acid, Tris and (v17, base columns)
    dilute Tris portions with a wash between any two of them, the mix, one wash. The block adds one wash after all its
    product (block_seconds)."""
    reagents = (na > 0) + (nb > 0) + (nd > 0)
    return (T_COLUMN + trips_seconds(ne) + trips_seconds(na) + trips_seconds(nb) + trips_seconds(nd)
            + max(reagents, 1) * T_WASH + T_MIX4)

def plan_seconds(name):
    """All columns in a plan file: 15 numbers per column (portions, then uL for tips 1-4, per liquid); v17: 20 in the
    base plans (the dilute Tris as a 4th liquid)."""
    rows = numbers(name)
    per = 20 if name.startswith("plan_base_") else 15
    total, cols = 0.0, 0
    for i in range(0, len(rows) - per + 1, per):
        ne, na, nb = int(rows[i]), int(rows[i + 5]), int(rows[i + 10])
        nd = int(rows[i + 15]) if per == 20 else 0
        if ne + na + nb + nd:
            total += column_seconds(ne, na, nb, nd)
            cols += 1
    return total, cols

def portions(ul):
    return int(ul // PORTION_UL) + (1 if ul % PORTION_UL > 1e-6 else 0)

def lane_seconds(acid_ul, eluate_ul, tris_ul, trim):
    """One lane of the 6-well plate, 8 tips (uL per tip). Making (v15): its eluate in the block's product part, then
    acid or Tris, mix, one wash (the block adds one wash after all the eluate: lanes_seconds). Trim: eluate portions
    each followed by a wash, acid, Tris, mix, wash."""
    na, ne, nb = portions(acid_ul), portions(eluate_ul), portions(tris_ul)
    if na + ne + nb == 0:
        return 0.0
    washes = ne + 1 if trim else 1 + (na > 0 and nb > 0)     # v13: no pool holds both acid and Tris
    mix_s = 2 * T_MIX8 + 8 * T_MIXCYCLE8                     # 5 + 3 mix cycles
    if nb:
        mix_s += T_MIX8 + 5 * T_MIXCYCLE8                    # v14: 5 more cycles where Tris went in (denser than eluate)
    return (na + ne + nb) * T_PORTION8 + washes * T_WASH + mix_s

def lanes_seconds(name, trim):
    rows = numbers(name)                              # 6 lanes x [flag, acid, eluate, Tris uL per tip]
    sec = sum(lane_seconds(rows[i + 1], rows[i + 2], rows[i + 3], trim)
              for i in range(0, len(rows) - 3, 4) if rows[i] > 0)
    return sec + (T_WASH if sec and not trim else 0.0)  # v15: one wash after the eluate of all the pools made

TYPICAL_ACID_C2 = (2, 1, 0)                           # portions of eluate, acid, Tris in a typical acid Cycle 2 column
def typical_base():
    """Portions in a base titration column before it is planned: v17 with the dilute Tris (options.csv, 4th value)
    ~2.5 mL eluate and up to 0.98 mL 0.1 M Tris; without it ~2.95 mL eluate and < 0.2 mL 1 M Tris."""
    o = numbers("options.csv")
    return (3, 0, 0, 1) if len(o) > 3 and o[3] > 0 else (4, 0, 1, 0)

def block_seconds(block):
    """Robot seconds of work in a block, before calibration. A plate block (acid_c1, acid_c2, val, tris_val) adds the
    one wash after all its product (v15)."""
    sec = T_WASH if block in ("acid_c1", "acid_c2", "val", "tris_val") else 0.0
    if block == "acid_c1":
        sec += plan_seconds("plan_c1.csv")[0] + plan_seconds("plan_base_c1.csv")[0]
    elif block == "acid_c2":
        s, n = plan_seconds("plan_c2.csv")
        sec += s if n else 2 * column_seconds(*TYPICAL_ACID_C2)
        s, n = plan_seconds("plan_base_c2.csv")              # before it is planned: one typical column if base Cycle 1 ran
        stopped = numbers("base_status.csv")[:1] == [0.0]   # (and the base titration was not stopped at plate check 1)
        sec += s if n else (column_seconds(*typical_base()) if plan_seconds("plan_base_c1.csv")[1] and not stopped else 0.0)
    elif block == "val":
        sec += plan_seconds("plan_val.csv")[0] + plan_seconds("plan_tris_c1.csv")[0] + plan_seconds("plan_base_val.csv")[0]
    elif block == "tris_val":
        sec += plan_seconds("plan_tris_val.csv")[0]
    elif block == "pools":
        sec += lanes_seconds("pool_make1.csv", False)
    elif block == "round2":
        sec += lanes_seconds("pool_trim1.csv", True) + lanes_seconds("pool_make2.csv", False)
    elif block == "round3":
        sec += lanes_seconds("pool_trim2.csv", True)
    return sec + T_MARKER

def scale():
    vals = sorted(numbers("timing_calibration.csv")[-12:])
    if not vals:
        return 1.0
    mid = len(vals) // 2
    med = vals[mid] if len(vals) % 2 else 0.5 * (vals[mid - 1] + vals[mid])
    return min(max(med, 0.5), 2.5)

def shown_minutes(seconds):
    """Whole minutes below 20, then rounded to 5."""
    m = seconds / 60.0
    if m < 20:
        return max(1, int(round(m)))
    return int(5 * round(m / 5.0))

def clock(minutes_from_now):
    """Clock time as HH.MM with the minutes rounded to the nearest value ending in 5 (05, 15 ... 55)."""
    t = datetime.now() + timedelta(minutes=minutes_from_now)
    total = t.hour * 60 + t.minute + t.second / 60.0
    k = int(round((total - 5) / 10.0))
    total = (10 * k + 5) % (24 * 60)
    return f"{total // 60}.{total % 60:02d}"

def log(event, block, raw_min, actual_min="", factor=""):
    try:
        new = not os.path.exists(LOG) or os.path.getsize(LOG) == 0
        with open(LOG, "a") as f:
            if new:
                f.write("timestamp,event,block,estimate_min,actual_min,scale\n")
            f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S},{event},{block},{raw_min},{actual_min},{factor}\n")
    except Exception:
        pass

def last_row(event, block):
    try:
        with open(LOG) as f:
            rows = [l.strip().split(",") for l in f.read().splitlines()[1:]]
    except Exception:
        return None
    for r in reversed(rows):
        if len(r) >= 4 and r[1] == event and r[2] == block:
            return r
    return None

def plan(block):
    raw = block_seconds(block)
    k = scale()
    minutes = shown_minutes(raw * k)
    nxt = shown_minutes(block_seconds("acid_c2") * k) if block == "acid_c1" else 0
    with open(OUT, "w") as f:
        f.write(f"{minutes}\n{clock(raw * k / 60.0)}\n{1 if minutes >= AWAY_MIN else 0}\n{nxt}\n")
    log("plan", block, round(raw / 60.0, 2), "", round(k, 3))

def start(block):
    log("start", block, round(block_seconds(block) / 60.0, 2), "", round(scale(), 3))

def end(block):
    s = last_row("start", block)
    if not s:
        return log("end", block, "", "", "")
    try:
        t0 = datetime.strptime(s[0], "%Y-%m-%d %H:%M:%S")
        raw = float(s[3])
    except Exception:
        return log("end", block, "", "", "")
    actual = (datetime.now() - t0).total_seconds() / 60.0
    log("end", block, round(raw, 2), round(actual, 2), round(scale(), 3))
    if raw >= 1.0 and 1.0 <= actual <= 600.0 and 0.3 <= actual / raw <= 4.0:    # a real block, not a test or a stall
        vals = numbers("timing_calibration.csv")[-11:] + [actual / raw]
        try:
            with open(CAL, "w") as f:
                f.write("\n".join(f"{v:.4f}" for v in vals) + "\n")
        except Exception:
            pass

if __name__ == "__main__":
    try:
        event, block = sys.argv[1].lower(), sys.argv[2].lower()
        {"plan": plan, "start": start, "end": end}[event](block)
    except Exception as e:
        try:
            with open(OUT, "w") as f:                  # never leave the message without numbers
                f.write("0\n0\n0\n0\n")
            with open(os.path.join(DIR, "python_log.txt"), "a") as f:
                f.write(f"Time estimate failed ({' '.join(sys.argv[1:])}): {e!r}\n")
        except Exception:
            pass
