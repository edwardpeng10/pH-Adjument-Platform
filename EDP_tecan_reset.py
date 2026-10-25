import os
import sys
import shutil
from datetime import datetime
# EDP_VERSION = "EDP-v17-2026-10"   (EDP_tecan_gp_master.py Check: all six EDP files must carry the same version)

if sys.stdout is None: sys.stdout = open(os.devnull, "w")
if sys.stderr is None: sys.stderr = open(os.devnull, "w")

# =====================================================================================
# EDP_tecan_reset.py, v17 - run at the start of every run.
# Empties every working file a run writes, found by name prefix, so nothing from a previous run can be read by mistake.
# Kept between runs: last_acid_ratios.csv, titrant_config.csv, timing_calibration.csv (the time estimates' calibration).
# python_log.txt starts fresh each run; the previous one is moved to logs\python_log_<date_time>.txt (the archive also
# keeps a copy with every run).
# v14 restores what the v10 reset did and the v11 rewrite dropped: a new run_id.csv every run (RUN_<date_time>; the
# archive ties a run's rows together with it - without it every run would be RUN_unknown), titrant_config.csv created
# if missing, and run_metadata.csv / probe_check.csv / tecan_instructions*.txt cleared (a stale run_metadata.csv would
# give every later run an old batch's names in the archive); val_status.csv is cleared too.
# =====================================================================================
DIR = os.environ.get("TECAN_DATA_DIR", r"C:\ProgramData\Tecan\EVOware\database\tecan_data")
KEEP = {"last_acid_ratios.csv", "titrant_config.csv", "timing_calibration.csv"}
PREFIXES = ("plan_", "needs_", "status", "ratios_", "corr_", "validation_", "cycle2_", "acid_", "current_ph_data",
            "method_used", "pool", "tris_", "setup_", "round", "time_est", "timing_log", "run_summary", "run_state",
            "entry", "entries_", "titration_params", "batch_info", "well_limits", "predicted_ratios", "chemistry_fit",
            "algorithm_version", "options", "custom_", "neutralisation_", "bulk_ph", "start_ph", "recheck", "last_acid_plan",
            "base_", "check_", "ratio_", "run_metadata", "probe_check", "tecan_instructions", "val_status", "drift_",
            "workbook_location")                                      # v17: where the run's workbook went (final message)

def main():
    for name in os.listdir(DIR):
        path = os.path.join(DIR, name)
        if name in KEEP or not os.path.isfile(path) or name.lower().endswith(".py"):
            continue
        if name.startswith(PREFIXES) and name.lower().endswith((".csv", ".txt", ".json")):
            try:
                open(path, "w").close()
            except Exception:
                pass
    try:                                                   # one identifier tying this run's archive rows together
        with open(os.path.join(DIR, "run_id.csv"), "w") as f:
            f.write(f"RUN_{datetime.now():%Y%m%d_%H%M%S}\n")
    except Exception:
        pass
    default = os.path.join(DIR, "titrant_config.csv")      # lab settings, not a working file: only created if missing
    if not os.path.exists(default) or os.path.getsize(default) == 0:
        try:
            with open(default, "w") as f:
                f.write("acid_name,acid_conc_M,base_name,base_conc_M\n")
                f.write("Acetic acid,1.0,Tris,1.0\n")
        except Exception:
            pass
    log = os.path.join(DIR, "python_log.txt")
    try:
        if os.path.exists(log) and os.path.getsize(log) > 0:
            os.makedirs(os.path.join(DIR, "logs"), exist_ok=True)
            shutil.move(log, os.path.join(DIR, "logs", f"python_log_{datetime.now():%Y%m%d_%H%M%S}.txt"))
        with open(log, "w") as f:
            f.write(f"=== {datetime.now():%Y-%m-%d %H:%M:%S} new run (EDP v17) ===\n")
    except Exception:
        pass

if __name__ == "__main__":
    try:
        main()
    except Exception:
        pass
