import os
import re
import sys
from datetime import datetime
# EDP_VERSION = "EDP-v17-2026-10"   (EDP_tecan_gp_master.py Check: all six EDP files must carry the same version)

if sys.stdout is None: sys.stdout = open(os.devnull, "w")      # started with pythonw: there is no console
if sys.stderr is None: sys.stderr = open(os.devnull, "w")

# Checks one typed entry of pH readings and writes them in the format the analysis expects.
#   EDP_parse_readings.py <entry file> <output file> <number of readings>
# Accepts values separated by commas, semicolons or spaces. Each must be 0 (well not measured) or between pH 1 and 14.
# Writes entry_ok.csv = 1 if the entry is usable, 0 if it should be typed again.
DIR = os.environ.get("TECAN_DATA_DIR", r"C:\ProgramData\Tecan\EVOware\database\tecan_data")

def read_text(name):
    try:
        with open(os.path.join(DIR, name), "r") as f:
            return f.read()
    except Exception:
        return ""

def run_details(entry):
    """One optional line: product type, protein mg/mL, room temperature C, initials, probe reading in pH 7.00 buffer,
    buffer description (everything after the fifth comma). Writes run_metadata.csv for the archive."""
    text = re.sub(r"^\s*\"?\s*N/?A\s*\"?\s*$", "", read_text(entry).strip(), flags=re.I).replace('"', "")
    parts = [p.strip() for p in text.split(",", 5)] + [""] * 6
    num = lambda x: float(re.sub(r"[^0-9.]", "", x) or 0) if re.search(r"\d", x) else 0.0
    batch = [q.strip().strip('"') for q in read_text("batch_info.csv").strip().split(",")] + ["", ""]
    well = re.findall(r"[\d.]+", read_text("well_limits.csv")) or ["0"]
    q = lambda x: '"' + x.replace('"', "") + '"'
    row = [q(batch[0]), q(batch[1]), f"{num(parts[0]):g}", f"{num(parts[1]):g}", f"{num(parts[2]):g}",
           q(parts[5] or "N/A"), q(parts[3] or "N/A"), f"{num(parts[4]):g}", well[0]]
    with open(os.path.join(DIR, "run_metadata.csv"), "w", newline="") as f:
        f.write(",".join(row) + "\r\n")
    with open(os.path.join(DIR, "entry_ok.csv"), "w") as f:
        f.write("1\n")

def ratio_pairs(entry):
    """Custom ratios requested as pairs of pH values, e.g. 3.2-7.4, 3.4-7.4 (up to 10 pairs; blank = none)."""
    vals = [float(x) for x in re.findall(r"\d+(?:\.\d+)?", read_text(entry).replace('"', " "))]
    vals = vals[:20] if len(vals) % 2 == 0 else vals[:len(vals) - 1][:20]
    n = len(vals) // 2
    with open(os.path.join(DIR, "custom_ratios_input.csv"), "w", newline="") as f:
        f.write(",".join(f"{v:.6f}" for v in vals + [0.0] * (20 - len(vals))) + "\r\n")
    with open(os.path.join(DIR, "custom_count.csv"), "w") as f:
        f.write(f"{n}\n")
    with open(os.path.join(DIR, "entry_ok.csv"), "w") as f:
        f.write("1\n")

def main():
    if len(sys.argv) >= 3 and sys.argv[1] == "meta":
        return run_details(sys.argv[2])
    if len(sys.argv) >= 3 and sys.argv[1] == "pairs":
        return ratio_pairs(sys.argv[2])
    entry, out = sys.argv[1], sys.argv[2]
    if sys.argv[3] == "pools":                       # v11: as many values as pools listed for this round
        try:
            with open(os.path.join(DIR, "pools_expected.csv")) as f:
                expected = int(float(f.read().split()[0]))
        except Exception:
            expected = 0
    else:
        expected = int(sys.argv[3])
    text = ""
    try:
        with open(os.path.join(DIR, entry), "r") as f:
            text = f.read()
    except Exception:
        pass
    text = re.sub(r"^\s*\"?\s*N/?A", "", text, flags=re.I)                  # the box's default text, if left in front
    tokens = [t for t in re.split(r"[,;\s]+", text.replace('"', " ")) if t and t.upper() not in ("N/A", "NA")]
    values, problem = [], ""
    for t in tokens:
        try:
            values.append(float(t))
        except ValueError:
            problem = f"'{t}' is not a number"
            break
    if not problem and len(values) != expected:
        problem = f"{len(values)} values instead of {expected}"
    if not problem:
        bad = [v for v in values if not (v == 0 or 1.0 <= v <= 14.0)]
        if bad:
            problem = f"value(s) {bad} outside pH 1-14"
    accepted_file = os.path.join(DIR, "entries_accepted.csv")
    key = ",".join(f"{v:.2f}" for v in values)
    if not problem and expected >= 4 and not os.path.basename(out).startswith(("validation", "pools_measured")):   # curves only
        try:
            with open(accepted_file) as f:
                for line in f.read().splitlines():
                    name, _, vals = line.partition("|")
                    if name != out and vals == key:
                        problem = f"identical to the readings already entered for {name} - check the right column was measured"
                        break
        except FileNotFoundError:
            pass
    ok = not problem
    if ok:
        with open(accepted_file, "a") as f:
            f.write(f"{out}|{key}\n")
        with open(os.path.join(DIR, out), "w", newline="") as f:
            f.write(",".join(f"{v:.6f}" for v in values) + "\r\n")
    with open(os.path.join(DIR, "entry_ok.csv"), "w") as f:
        f.write("1\n" if ok else "0\n")
    try:                                               # one row per typed entry: used to see where a run's time goes
        with open(os.path.join(DIR, "entries_log.csv"), "a") as f:
            if f.tell() == 0:
                f.write("timestamp,file,expected,accepted,values,problem\n")
            f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S},{out},{expected},{int(ok)},"
                    f"\"{' '.join(f'{v:.2f}' for v in values)}\",\"{problem}\"\n")
    except Exception:
        pass
    try:
        with open(os.path.join(DIR, "python_log.txt"), "a") as f:
            f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S} ENTRY {out}: " + (", ".join(f"{v:.2f}" for v in values) if ok else f"REJECTED ({problem}): {text.strip()}") + "\n")
    except Exception:
        pass

if __name__ == "__main__":
    main()
