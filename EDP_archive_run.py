import sys
import os
import shutil
import glob
import pandas as pd
import numpy as np
import datetime
# EDP_VERSION = "EDP-v17-2026-10"   (EDP_tecan_gp_master.py Check: all six EDP files must carry the same version)

if sys.stdout is None: sys.stdout = open(os.devnull, "w")      # started with pythonw: there is no console
if sys.stderr is None: sys.stderr = open(os.devnull, "w")

# Archives one titration into the ML data store, or (with "Session") closes off a run.
#   EDP_archive_run.py Acid                  -> archive the acid titration
#   EDP_archive_run.py Base P32|P34|P36      -> archive one neutralisation pool
#   EDP_archive_run.py Session               -> write the run-level row once the run is finished
#   EDP_archive_run.py Pools                 -> v11: archive the pools made in the 6-well plate (pools.csv), and
#                                               their readings as lane points (measurements.csv, Point_Type "Pool")
#   EDP_archive_run.py Base TR               -> v11: the VI neutralisation data, Tris into the acidified eluate (recipe
#                                               wells at the process pool)
#   EDP_archive_run.py Base EL               -> v13: the base titration, Tris straight into the eluate (Material
#                                               "eluate", Initial_pH = the eluate's own pH)
#   EDP_archive_run.py Map                   -> v14: the run's ratio map (ratio_map.csv: acid, base and neutralisation
#                                               ratios from any start pH) and the chemistry-model fingerprint
#                                               (chemistry_fits.csv, Phase "Model")
#
# Tables written to ML_data_archive:
#   experiment_summary.csv  one row per titration (conditions, status, metadata)
#   measurements.csv        every curve point, with its cycle
#   target_results.csv      target, recommended volume, validated pH, error
#   experiment_warnings.csv anything that needed flagging
#   predictions.csv         what the chemistry predicted from Cycle 1 alone, against what was measured (shadow mode)
#   chemistry_fits.csv      fitted chemistry for each curve: the parameters a predictive model is built on
#   sessions.csv            one row per run (operator, temperature, probe checks, titrants, timings)
#   runs/<Run_ID>/...       a copy of the working files behind every archived titration
# v8: target_results gains Corrected_Ratio and Correction_Pct (the validation refit, EDP_tecan_gp_master.py cycle 3);
#   measurements gains Material (eluate / recipe made in the well / bulk pool) and Excluded (typo exclusions), and
#   holds the validation wells as Cycle 3; chemistry_fits gains the Tris fit (pKa', recipe and real-pool
#   equivalence, their ratio); neutralisation Initial_pH is the bulk pool's measured pH; sessions gains Robot_Min
#   and Stop_Min from timing_log.csv. Older files are migrated automatically (a _before_<timestamp>.csv is kept).

WORKING_DIR = os.environ.get("TECAN_DATA_DIR", r"C:\ProgramData\Tecan\EVOware\database\tecan_data")
ARCHIVE_DIR = os.environ.get("TECAN_ARCHIVE_DIR", r"C:\ProgramData\Tecan\EVOware\database\ML_data_archive")

FILES = {
    "BATCH": os.path.join(WORKING_DIR, "batch_info.csv"),
    "META": os.path.join(WORKING_DIR, "run_metadata.csv"),
    "TITRANT": os.path.join(WORKING_DIR, "titrant_config.csv"),
    "RUNID": os.path.join(WORKING_DIR, "run_id.csv"),
    "PARAMS": os.path.join(WORKING_DIR, "titration_params.csv"),
    "CURVE": os.path.join(WORKING_DIR, "current_ph_data.csv"),
    "RATIOS": os.path.join(WORKING_DIR, "ratios_clean.csv"),
    "VAL": os.path.join(WORKING_DIR, "validation_results.csv"),
    "METHOD": os.path.join(WORKING_DIR, "method_used.txt"),
    "STATUS": os.path.join(WORKING_DIR, "status.csv"),
    "PREDICT": os.path.join(WORKING_DIR, "predicted_ratios.csv"),
    "CHEM": os.path.join(WORKING_DIR, "chemistry_fit.csv"),
    "PLAN": os.path.join(WORKING_DIR, "acid_plan.txt"),
    "VERSION": os.path.join(WORKING_DIR, "algorithm_version.csv"),
    "ENTRIES": os.path.join(WORKING_DIR, "entries_log.csv"),
    "PROBE": os.path.join(WORKING_DIR, "probe_check.csv"),
    "LOG": os.path.join(WORKING_DIR, "python_log.txt"),
    "CORR": os.path.join(WORKING_DIR, "corr_ratios.csv"),
    "TRIS": os.path.join(WORKING_DIR, "tris_fit.csv"),
    "TIMING": os.path.join(WORKING_DIR, "timing_log.csv"),
    "TIMECAL": os.path.join(WORKING_DIR, "timing_calibration.csv"),
    "DRIFT": os.path.join(WORKING_DIR, "drift_check.csv"),
    "STATE": os.path.join(WORKING_DIR, "run_state.json"),
}
POOL_TAGGED = ["PARAMS", "CURVE", "RATIOS", "VAL", "METHOD", "STATUS", "PREDICT", "CHEM", "CORR", "TRIS"]
PRODUCT_TYPES = {1: "Unadjusted ProA eluate", 2: "LpH VI", 3: "nLpH VI", 4: "AEX product", 5: "Adjusted CEX load",
                 6: "CEX eluate", 7: "Other"}

SUMMARY_COLUMNS = ["Experiment_ID", "Run_ID", "Pool", "Batch_ID", "Product_Name", "Product_Type", "Phase", "Timestamp",
                   "Status", "Status_Code", "Algorithm_Version", "Product_Vol_mL", "Max_Titrant_mL", "Max_Well_Fill_mL",
                   "Titrant_Name", "Titrant_Conc_M", "Initial_pH", "Process_Target_pH", "Protein_Conc_mgmL",
                   "Temperature_C", "Buffer", "Operator", "Probe_Check_Start", "Plan", "Bulk_Pool_pH", "Well_Volume_mL",
                   "Real_Recipe_Factor"]

# --- DATA CLEANING HELPERS ---
def clean_str(val):
    if pd.isna(val): return "Unknown"
    s = str(val).strip()
    if s in ["", "nan", "None", "0", "0.0"]:
        return "Unknown"
    return s

def clean_float(val, treat_zero_as_missing=False):
    if pd.isna(val): return np.nan
    try:
        f = float(val)
        if treat_zero_as_missing and f == 0.0:
            return np.nan
        return f
    except:
        return np.nan

def safe_read_csv(filepath, header=None):
    if not os.path.exists(filepath) or os.stat(filepath).st_size == 0:
        return None
    try:
        return pd.read_csv(filepath, header=header)
    except:
        return None

def read_text(path, default=""):
    try:
        with open(path) as f:
            return f.read().strip()
    except Exception:
        return default

def append_table(path, rows, columns=None):
    """Append rows, creating the file with a header. If the file predates a column, rewrite it once (keeping a backup)."""
    df = pd.DataFrame(rows)
    if columns:
        df = df.reindex(columns=columns)
    if os.path.exists(path) and os.stat(path).st_size > 0:
        try:
            old = pd.read_csv(path)
        except Exception:
            old = None
        if old is not None and list(old.columns) != list(df.columns):
            stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            shutil.copyfile(path, path.replace(".csv", f"_before_{stamp}.csv"))
            pd.concat([old.reindex(columns=df.columns), df], ignore_index=True).to_csv(path, index=False)
            return
    df.to_csv(path, mode="a", index=False, header=not os.path.exists(path) or os.stat(path).st_size == 0)

def use_pool_files(pool):
    """Neutralisation pools run side by side, each with its own working files (e.g. status_P32.csv)."""
    for key in POOL_TAGGED:
        root, ext = os.path.splitext(FILES[key])
        tagged = f"{root}_{pool}{ext}"
        if os.path.exists(tagged):
            FILES[key] = tagged
    FILES["PLAN"] = os.path.join(WORKING_DIR, f"neutralisation_plan_{pool}.txt")

def read_metadata():
    """Run-level details entered by the operator, with the titrants taken from the lab's configuration file."""
    meta = dict(batch_id="Unknown", product_name="Unknown", product_type="Unknown", protein=np.nan, temperature=np.nan,
                buffer="Unknown", operator="Unknown", probe_start=np.nan, max_fill=np.nan,
                acid_name="Unknown", acid_conc=np.nan, base_name="Unknown", base_conc=np.nan)
    df = safe_read_csv(FILES["META"])
    if df is not None and not df.empty:
        row = df.iloc[0]
        get = lambda i: row.iloc[i] if len(row) > i else np.nan
        meta.update(batch_id=clean_str(get(0)), product_name=clean_str(get(1)),
                    protein=clean_float(get(3), True), temperature=clean_float(get(4), True),
                    buffer=clean_str(get(5)), operator=clean_str(get(6)),
                    probe_start=clean_float(get(7), True), max_fill=clean_float(get(8), True))
        code = clean_float(get(2))
        meta["product_type"] = PRODUCT_TYPES.get(int(code), clean_str(get(2))) if not pd.isna(code) else clean_str(get(2))
    else:                                                    # older runs: fall back to batch_info.csv
        df_b = safe_read_csv(FILES["BATCH"])
        if df_b is not None and not df_b.empty:
            row = df_b.iloc[0]
            meta.update(batch_id=clean_str(row.iloc[0]) if len(row) > 0 else "Unknown",
                        product_name=clean_str(row.iloc[1]) if len(row) > 1 else "Unknown",
                        product_type=clean_str(row.iloc[2]) if len(row) > 2 else "Unknown")
    t = safe_read_csv(FILES["TITRANT"], header=0)
    if t is not None and not t.empty and len(t.columns) >= 4:
        row = t.iloc[0]
        meta.update(acid_name=clean_str(row.iloc[0]), acid_conc=clean_float(row.iloc[1]),
                    base_name=clean_str(row.iloc[2]), base_conc=clean_float(row.iloc[3]))
    return meta

def copy_working_files(run_id, experiment_id):
    """Keep the working files behind this titration, so any archived row can be traced back to what was measured."""
    dest = os.path.join(ARCHIVE_DIR, "runs", run_id, experiment_id)
    try:
        os.makedirs(dest, exist_ok=True)
        for key in ("PARAMS", "CURVE", "RATIOS", "VAL", "METHOD", "STATUS", "PREDICT", "CHEM", "PLAN", "CORR", "TRIS"):
            src = FILES[key]
            if os.path.exists(src) and os.stat(src).st_size > 0:
                shutil.copyfile(src, os.path.join(dest, os.path.basename(src)))
        tag = experiment_id.split("_")[-1] if experiment_id.split("_")[-1].startswith("P") else ""
        names = ["cycle2_measured.csv", "cycle2_vols.csv", "plan_c1.csv", "plan_c2.csv", "plan_val.csv"]
        candidates = [f"neutralisation_c1_measured_{tag}.csv"] + [n.replace(".csv", f"_{tag}.csv") for n in names] \
            if tag else ["acid_c1_measured.csv"] + names
        if experiment_id.endswith("_EL"):                   # v13: base titration (Tris straight into the eluate)
            candidates = ["base_c1_measured.csv", "base_c2_measured.csv", "base_val_measured.csv", "plan_base_c1.csv",
                          "plan_base_c2.csv", "plan_base_val.csv", "base_corr.csv", "acid_c1_measured.csv"]
        elif experiment_id.endswith("_TR"):                 # VI neutralisation data (recipe wells)
            candidates = ["tris_c1_measured.csv", "tris_val_measured.csv", "plan_tris_c1.csv", "plan_tris_val.csv",
                          "tris_corr.csv", "tris_ratios.csv"]
        for candidate in candidates:
            src = os.path.join(WORKING_DIR, candidate)
            if os.path.exists(src) and os.stat(src).st_size > 0:
                shutil.copyfile(src, os.path.join(dest, candidate))
    except Exception:
        pass

BASE_FILES = dict(FILES)

def main():
    phase_tag = sys.argv[1] if len(sys.argv) > 1 else "Base"
    pool = sys.argv[2].strip() if len(sys.argv) > 2 else ""
    if not os.path.exists(ARCHIVE_DIR):
        os.makedirs(ARCHIVE_DIR)
    run_id = read_text(FILES["RUNID"]) or "RUN_unknown"
    meta = read_metadata()
    if phase_tag.lower() == "session":
        return write_session(run_id, meta)
    if phase_tag.lower() == "pools":
        return archive_pools(run_id, meta)
    if phase_tag.lower() == "map":
        return archive_map(run_id, meta)
    if pool.lower() == "all":                                   # every selected pool, one after another
        flags = [clean_float(x) for x in read_text(os.path.join(WORKING_DIR, "pools_selected.csv")).replace(",", " ").split()]
        for i, tag in enumerate(("P32", "P34", "P36")):
            if i < len(flags) and flags[i] == 1:
                FILES.clear(); FILES.update(BASE_FILES)
                archive_one(phase_tag, tag, run_id, meta)
        return
    archive_one(phase_tag, pool, run_id, meta)

def archive_one(phase_tag, pool, run_id, meta):
    if pool:
        use_pool_files(pool)

    now = datetime.datetime.now()
    experiment_id = f"EXP_{now.strftime('%Y%m%d_%H%M%S')}" + (f"_{pool}" if pool else "")
    timestamp_str = now.strftime("%Y-%m-%d %H:%M:%S")
    warnings_list = []

    def log_warning(warning_type, details):
        warnings_list.append({"Experiment_ID": experiment_id, "Run_ID": run_id, "Timestamp": timestamp_str,
                              "Warning_Type": warning_type, "Details": details})

    if meta["batch_id"] == "Unknown": log_warning("Missing Metadata", "Batch_ID was not provided.")
    if meta["product_name"] == "Unknown": log_warning("Missing Metadata", "Product_Name was not provided.")
    if meta["acid_name"] == "Unknown": log_warning("Missing Metadata", "Titrants were not recorded (titrant_config.csv).")
    if pd.isna(meta["protein"]): log_warning("Missing Metadata", "Protein concentration was not entered.")

    # --- parameters ---
    prod_vol, max_titrant = 1.0, 0.25
    target_phs = [np.nan, np.nan, np.nan, np.nan]
    process_target_ph, deviation_tolerance = np.nan, 0.20
    df_p = safe_read_csv(FILES["PARAMS"])
    if df_p is not None and not df_p.empty:
        row = df_p.iloc[0]
        if len(row) > 0: prod_vol = clean_float(row.iloc[0]) if not pd.isna(clean_float(row.iloc[0])) else 1.0
        if len(row) > 1: max_titrant = clean_float(row.iloc[1]) if not pd.isna(clean_float(row.iloc[1])) else 0.25
        for i in range(4):
            if len(row) > i + 2:
                target_phs[i] = clean_float(row.iloc[i + 2], treat_zero_as_missing=True)
        if len(row) > 5:
            process_target_ph = clean_float(row.iloc[5], treat_zero_as_missing=True)
        if len(row) > 6:
            custom_tol = clean_float(row.iloc[6])
            if not pd.isna(custom_tol):
                deviation_tolerance = custom_tol
    if pd.isna(process_target_ph):
        log_warning("Missing Parameter", "Process_Target_pH was left blank or zero.")

    status_code, status_text = "Unknown", "Unknown"
    if os.path.exists(FILES["STATUS"]):
        status_code = clean_str(read_text(FILES["STATUS"]))
        status_text = {"0": "Cycle2_Needed", "1": "Final", "2": "Final_Cycle2", "3": "Uncertain"}.get(status_code, "Unknown")

    df_curve = safe_read_csv(FILES["CURVE"], header=0)
    initial_ph = np.nan
    is_base = phase_tag.lower().startswith("base")
    bulk_ph = clean_float(read_text(os.path.join(WORKING_DIR, f"bulk_ph_{pool}.csv")), True) if pool else np.nan
    if is_base and pool != "EL":                 # v8: Cycle 1 starts near the equivalence point; the pool's own pH
        initial_ph = bulk_ph                     # is the bulk pool's measured pH
    elif df_curve is not None and 'pH' in df_curve.columns and len(df_curve) > 0:
        initial_ph = clean_float(df_curve['pH'].iloc[0], treat_zero_as_missing=True)
        if pd.isna(initial_ph):
            log_warning("Missing Measurement", "Initial pH reading failed or recorded as 0.0.")
    tris = {}
    df_t = safe_read_csv(FILES["TRIS"], header=0)
    if df_t is not None and not df_t.empty:
        tris = {k: clean_float(v) for k, v in df_t.iloc[0].items()}
    factor = tris.get("real_recipe_factor", np.nan)
    if is_base and not pd.isna(factor) and abs(factor - 1) > 0.05:
        log_warning("Bulk Pool Differs", f"The bulk pool needed {factor:.3f} x the Tris of its recipe (made in the wells) "
                    "- the ratios are for the bulk pool as the robot draws it.")

    # --- TABLE 1: EXPERIMENT SUMMARY ---
    titrant_name = meta["base_name"] if phase_tag.lower().startswith("base") else meta["acid_name"]
    titrant_conc = meta["base_conc"] if phase_tag.lower().startswith("base") else meta["acid_conc"]
    summary = {
        "Experiment_ID": experiment_id, "Run_ID": run_id, "Pool": pool or "", "Batch_ID": meta["batch_id"],
        "Product_Name": meta["product_name"], "Product_Type": meta["product_type"], "Phase": phase_tag,
        "Timestamp": timestamp_str, "Status": status_text, "Status_Code": status_code,
        "Algorithm_Version": read_text(FILES["VERSION"], "Unknown"), "Product_Vol_mL": prod_vol,
        "Max_Titrant_mL": max_titrant, "Max_Well_Fill_mL": meta["max_fill"], "Titrant_Name": titrant_name,
        "Titrant_Conc_M": titrant_conc, "Initial_pH": initial_ph, "Process_Target_pH": process_target_ph,
        "Protein_Conc_mgmL": meta["protein"], "Temperature_C": meta["temperature"], "Buffer": meta["buffer"],
        "Operator": meta["operator"], "Probe_Check_Start": meta["probe_start"],
        "Plan": read_text(FILES["PLAN"])[:200],
        "Bulk_Pool_pH": bulk_ph,
        "Well_Volume_mL": meta["max_fill"],
        "Real_Recipe_Factor": factor,
    }
    append_table(os.path.join(ARCHIVE_DIR, "experiment_summary.csv"), [summary], SUMMARY_COLUMNS)

    # --- TABLE 2: MEASUREMENTS ---
    if df_curve is not None and not df_curve.empty:
        meas_rows = []
        for idx, row in df_curve.iterrows():
            vol = clean_float(row.get('Volume_Added', np.nan))
            if pd.isna(vol): continue
            pt_idx = idx + 1
            cyc = row.get('Cycle', np.nan)
            cycle = int(cyc) if not pd.isna(cyc) else (1 if pt_idx <= 12 else 2)
            measured_ph = clean_float(row.get('pH', np.nan), treat_zero_as_missing=True)
            if pd.isna(measured_ph) and cycle != 3:      # missing validation wells are reported under Targets
                log_warning("Missing Measurement", f"Probe recorded 0.0 or NaN at Point {pt_idx} (Vol {vol} mL).")
            valid = clean_float(row.get('Valid', 1))
            stock = clean_float(row.get('Stock_M', 1.0))          # v17: the Tris stock of a base well (1 M or the dilute Tris)
            meas_rows.append({
                "Experiment_ID": experiment_id, "Run_ID": run_id, "Point_Index": pt_idx, "Cycle": cycle,
                "Volume_Added_mL": vol, "Ratio": round(vol / prod_vol, 5) if prod_vol else np.nan,
                "Measured_pH": measured_ph,
                "GPR_Predicted_pH": clean_float(row.get('GPR_Predicted', np.nan), treat_zero_as_missing=True),
                "GPR_Variance": clean_float(row.get('GPR_Variance', np.nan)),
                "Point_Type": {1: "Initial", 2: "Cycle 2", 3: "Validation"}.get(cycle, "Cycle 2"),
                "Material": ("eluate" if not is_base or pool == "EL" else "recipe made in the well" if (cycle == 1 or pool == "TR") else "bulk pool"),
                "Excluded": 0 if pd.isna(measured_ph) or pd.isna(valid) or valid == 1 else 1,
                "Titrant_Stock_M": stock if not pd.isna(stock) else 1.0})
        if meas_rows:
            append_table(os.path.join(ARCHIVE_DIR, "measurements.csv"), meas_rows,
                         ["Experiment_ID", "Run_ID", "Point_Index", "Cycle", "Volume_Added_mL", "Ratio",
                          "Measured_pH", "GPR_Predicted_pH", "GPR_Variance", "Point_Type", "Material", "Excluded",
                          "Titrant_Stock_M"])

    # --- TABLE 3: TARGET RESULTS (and the shadow prediction for each target) ---
    ratios, validations, methods, corrected = [], [], [], []
    df_r = safe_read_csv(FILES["RATIOS"])
    if df_r is not None: ratios = df_r.iloc[:, 0].tolist()
    df_k = safe_read_csv(FILES["CORR"])
    if df_k is not None: corrected = df_k.iloc[:, 0].tolist()
    df_v = safe_read_csv(FILES["VAL"])
    if df_v is not None and not df_v.empty: validations = df_v.iloc[0].tolist()
    if os.path.exists(FILES["METHOD"]):
        methods = [line.strip() for line in read_text(FILES["METHOD"]).splitlines() if line.strip()]
    predicted = {}
    df_pred = safe_read_csv(FILES["PREDICT"], header=0)
    fit_rms = np.nan
    if df_pred is not None and not df_pred.empty:
        for _, r in df_pred.iterrows():
            predicted[clean_float(r.get("target_pH"))] = clean_float(r.get("predicted_ratio"))
        fit_rms = clean_float(df_pred.iloc[0].get("fit_rms_pH"))

    target_rows, prediction_rows, counter = [], [], 1
    for i, t_ph in enumerate(target_phs):
        if pd.isna(t_ph): continue
        ratio = clean_float(ratios[i]) if i < len(ratios) else np.nan
        rec_vol = ratio * prod_vol if not pd.isna(ratio) else np.nan
        val_ph = clean_float(validations[i], treat_zero_as_missing=True) if i < len(validations) else np.nan
        method = clean_str(methods[i]) if i < len(methods) else "Unknown"
        abs_err = round(abs(val_ph - t_ph), 4) if not pd.isna(val_ph) else np.nan
        if pd.isna(val_ph):
            log_warning("Missing Validation", f"No validation pH for Target {t_ph}.")
        elif abs_err > deviation_tolerance:
            log_warning("High Deviation", f"Target {t_ph} validated at {val_ph} (error {abs_err} > {deviation_tolerance}).")
        corr = clean_float(corrected[i]) if i < len(corrected) else np.nan
        target_rows.append({
            "Experiment_ID": experiment_id, "Run_ID": run_id, "Target_Index": counter, "Target_pH": t_ph,
            "Is_Process_Target": "Yes" if t_ph == process_target_ph else "No",
            "Ratio": ratio, "Recommended_Volume_mL": round(rec_vol, 4) if not pd.isna(rec_vol) else np.nan,
            "Validated_pH": val_ph, "Method": method, "Absolute_Error": abs_err,
            "Corrected_Ratio": round(corr, 4) if not pd.isna(corr) else np.nan,
            "Correction_Pct": round(100 * (corr / ratio - 1), 2) if not pd.isna(corr) and not pd.isna(ratio) and ratio else np.nan})
        pred_ratio = predicted.get(t_ph, np.nan)
        prediction_rows.append({
            "Experiment_ID": experiment_id, "Run_ID": run_id, "Phase": phase_tag, "Pool": pool or "",
            "Target_pH": t_ph, "Predicted_Ratio_Cycle1": pred_ratio, "Fit_RMS_pH": fit_rms,
            "Final_Ratio": ratio, "Validated_pH": val_ph, "Method": method,
            "Predicted_Minus_Final_Pct": round(100 * (pred_ratio - ratio) / ratio, 2)
            if not pd.isna(pred_ratio) and not pd.isna(ratio) and ratio else np.nan})
        counter += 1
    # validation check shown to the operator: signed error per target, worst error, limit, pass (1) / fail (0)
    errs = [r["Validated_pH"] - r["Target_pH"] for r in target_rows]
    measured = [e for e in errs if not pd.isna(e)]
    worst = max((abs(e) for e in measured), default=0.0)
    check = [round(e, 2) if not pd.isna(e) else 0.0 for e in errs] + [0.0] * (4 - len(errs))
    check += [round(worst, 2), deviation_tolerance, int(bool(measured) and worst <= deviation_tolerance), len(measured)]
    try:
        with open(os.path.join(WORKING_DIR, f"validation_check{('_' + pool) if pool else ''}.csv"), "w") as f:
            f.write("\n".join(f"{x}" for x in check[:4] + check[-4:]) + "\n")
    except Exception:
        pass
    if target_rows:
        append_table(os.path.join(ARCHIVE_DIR, "target_results.csv"), target_rows,
                     ["Experiment_ID", "Run_ID", "Target_Index", "Target_pH", "Is_Process_Target", "Ratio",
                      "Recommended_Volume_mL", "Validated_pH", "Method", "Absolute_Error", "Corrected_Ratio",
                      "Correction_Pct"])
    if prediction_rows:
        append_table(os.path.join(ARCHIVE_DIR, "predictions.csv"), prediction_rows,
                     ["Experiment_ID", "Run_ID", "Phase", "Pool", "Target_pH", "Predicted_Ratio_Cycle1", "Fit_RMS_pH",
                      "Final_Ratio", "Validated_pH", "Method", "Predicted_Minus_Final_Pct"])

    # --- TABLE 5: CHEMISTRY FIT ---
    df_c = safe_read_csv(FILES["CHEM"], header=0)
    if df_c is not None and not df_c.empty:
        r = df_c.iloc[0]
        append_table(os.path.join(ARCHIVE_DIR, "chemistry_fits.csv"), [{
            "Experiment_ID": experiment_id, "Run_ID": run_id, "Phase": phase_tag, "Pool": pool or "",
            "Product_Name": meta["product_name"], "Titrant": clean_str(r.get("titrant")),
            "Acetate_M": clean_float(r.get("acetate_M")), "Strong_Ions_M": clean_float(r.get("strong_ions_M")),
            "Fit_RMS_pH": clean_float(r.get("fit_rms_pH")), "N_Points": clean_float(r.get("n_points")),
            "Product_Vol_mL": prod_vol, "Initial_pH": initial_ph,
            "Tris_pKa": tris.get("tris_pka", np.nan), "Te_Recipe": tris.get("te_recipe", np.nan),
            "Te_Real_Pool": tris.get("te_real", np.nan), "Real_Recipe_Factor": factor,
            "Tris_Fit_RMS_pH": tris.get("fit_rms_pH", np.nan)}],
            ["Experiment_ID", "Run_ID", "Phase", "Pool", "Product_Name", "Titrant", "Acetate_M", "Strong_Ions_M",
             "Fit_RMS_pH", "N_Points", "Product_Vol_mL", "Initial_pH", "Tris_pKa", "Te_Recipe", "Te_Real_Pool",
             "Real_Recipe_Factor", "Tris_Fit_RMS_pH"])

    # --- TABLE 4: WARNINGS ---
    if warnings_list:
        append_table(os.path.join(ARCHIVE_DIR, "experiment_warnings.csv"), warnings_list,
                     ["Experiment_ID", "Run_ID", "Timestamp", "Warning_Type", "Details"])
    copy_working_files(run_id, experiment_id)

def number_lines(name):
    """Numbers in a working file, one per line or comma-separated; [] if missing or empty."""
    out = []
    for tok in read_text(os.path.join(WORKING_DIR, name)).replace(",", " ").replace(";", " ").split():
        f = clean_float(tok.strip('"'))
        if not pd.isna(f):
            out.append(f)
    return out

def write_run_summary():
    """v8: the numbers the final message shows (run_summary.csv, one per line): worst validation error for the acid
    and each pool (-1 = not validated), the acid ratios after the validation refit, each pool's Tris ratio at the
    process pH after the refit (0 = none), and each bulk pool's Tris need relative to its recipe (0 = not measured)."""
    rows = []
    for tag in ("", "_P32", "_P34", "_P36"):
        chk = number_lines(f"validation_check{tag}.csv")
        rows.append(f"{chk[4]:.2f}" if len(chk) >= 8 and chk[7] > 0 else "-1")
    acid = number_lines("corr_ratios.csv")
    rows += [f"{acid[i]:.4f}" if len(acid) == 4 else "0" for i in range(4)]
    for tag in ("P32", "P34", "P36"):
        corr = number_lines(f"corr_ratios_{tag}.csv")
        rows.append(f"{corr[3]:.4f}" if len(corr) == 4 else "0")
    for tag in ("P32", "P34", "P36"):
        df_t = safe_read_csv(os.path.join(WORKING_DIR, f"tris_fit_{tag}.csv"), header=0)
        f = clean_float(df_t.iloc[0].get("real_recipe_factor")) if df_t is not None and not df_t.empty else np.nan
        rows.append(f"{f:.3f}" if not pd.isna(f) else "0")
    try:
        with open(os.path.join(WORKING_DIR, "run_summary.csv"), "w") as f:
            f.write("\n".join(rows) + "\n")
    except Exception:
        pass

def archive_pools(run_id, meta):
    """v11: one row per pool made in the 6-well plate (run_state.json, written by the analysis): lane, target, route
    (v13: acid = eluate + acid, Tris = eluate + Tris; Final_Ratio is acid or Tris per mL eluate), planned ratio, volumes, the pH after each round, trims, final pH and error, and - for low-pH pools - the Tris
    needed per mL of pool at each neutralisation target. Every reading also goes to measurements.csv as a lane point
    (Point_Type "Pool") at the pool's nominal composition."""
    import json
    try:
        with open(os.path.join(WORKING_DIR, "run_state.json")) as f:
            st = json.load(f)
    except Exception:
        return
    pools = st.get("pools", [])
    if not pools:
        return
    tt = number_lines("tris_params.csv")
    now = datetime.datetime.now()
    rows, points = [], []
    for pl in pools:
        experiment_id = f"EXP_{now.strftime('%Y%m%d_%H%M%S')}_LANE{pl['lane']}"
        r = pl.get("readings", [])
        final = r[-1] if r else np.nan
        A, E, B = pl.get("A", 0.0), pl.get("E", 0.0), pl.get("B", 0.0)
        base = pl.get("route") == "base"                      # v13: eluate + Tris (Tris per mL eluate)
        ratio = (B / E if E else np.nan) if base else (A / E if E else np.nan)
        row = {"Experiment_ID": experiment_id, "Run_ID": run_id, "Timestamp": now.strftime("%Y-%m-%d %H:%M:%S"),
               "Batch_ID": meta["batch_id"], "Product_Name": meta["product_name"], "Product_Type": meta["product_type"],
               "Lane": pl["lane"], "Target_pH": pl["target"], "Route": {"acid": "acid", "base": "Tris", "basic": "acid then Tris"}.get(pl.get("route"), pl.get("route")),
               "Status": pl.get("status"), "Planned_Ratio": pl.get("planned_ratio", np.nan),
               "Acid_mL": round(A, 3), "Eluate_mL": round(E, 3), "Tris_mL": round(B, 3),
               "Dilute_Tris_mL": round(pl.get("Bd", 0.0), 3),       # v17: trims from the dilute Tris (its mmol are in Tris_mL)
               "Volume_mL": round(A + E + B + pl.get("W", 0.0), 2),
               "Final_Ratio": round(ratio, 5) if not pd.isna(ratio) else np.nan,
               "pH_1": r[0] if len(r) > 0 else np.nan, "pH_2": r[1] if len(r) > 1 else np.nan, "pH_3": r[2] if len(r) > 2 else np.nan,
               "Trims": pl.get("trims", 0), "Final_pH": final,
               "Final_Error": round(final - pl["target"], 3) if not pd.isna(final) and final else np.nan,
               "Note": pl.get("note", ""), "Lane_Factor": round(pl["lane_factor"], 4) if pl.get("lane_factor") else np.nan,
               "Aim_pH": round(pl.get("aim", 0.0), 3), "Algorithm_Version": read_text(FILES["VERSION"], "Unknown")}
        for k, x in enumerate(pl.get("tris_ratios", [])[:4]):
            row[f"Tris_Ratio_pH_{tt[k]:g}" if k < len(tt) else f"Tris_Ratio_{k + 1}"] = round(x, 5)
        rows.append(row)
        hist = pl.get("hist") or []
        vols = pl.get("vol") or []                            # v17: total mL per mL eluate when each reading was taken
        for i, m in enumerate(r):
            ri = hist[i][0] if i < len(hist) else ratio        # v17: the pool's ratio when that reading was taken (trims change it)
            # v17: the Tris stock as mmol per mL of the Tris solutions added so far - below 1 after trims from the 0.1 M Tris
            sol = vols[i] - 1 - (A / E if E else 0.0) if i < len(vols) else 0.0
            stock = round(min(ri / sol, 1.0), 4) if base and sol > 1e-9 and not pd.isna(ri) else 1.0
            points.append({"Experiment_ID": experiment_id, "Run_ID": run_id, "Point_Index": i + 1, "Cycle": i + 1,
                           "Volume_Added_mL": round(ri, 5) if not pd.isna(ri) else np.nan,
                           "Ratio": round(ri, 5) if not pd.isna(ri) else np.nan, "Measured_pH": m,
                           "GPR_Predicted_pH": np.nan, "GPR_Variance": np.nan, "Point_Type": "Pool",
                           "Material": "pool (lane, Tris)" if base else "pool (lane)", "Excluded": 0, "Titrant_Stock_M": stock})
        if not pd.isna(final) and final and abs(final - pl["target"]) > 0.05:
            append_table(os.path.join(ARCHIVE_DIR, "experiment_warnings.csv"), [{
                "Experiment_ID": experiment_id, "Run_ID": run_id, "Timestamp": row["Timestamp"], "Warning_Type": "Pool Off Target",
                "Details": f"Lane {pl['lane']} ended at pH {final:.2f} against target {pl['target']:.2f}."}],
                ["Experiment_ID", "Run_ID", "Timestamp", "Warning_Type", "Details"])
    cols = list(dict.fromkeys(k for row in rows for k in row))
    append_table(os.path.join(ARCHIVE_DIR, "pools.csv"), rows, cols)
    if points:
        append_table(os.path.join(ARCHIVE_DIR, "measurements.csv"), points,
                     ["Experiment_ID", "Run_ID", "Point_Index", "Cycle", "Volume_Added_mL", "Ratio",
                      "Measured_pH", "GPR_Predicted_pH", "GPR_Variance", "Point_Type", "Material", "Excluded", "Titrant_Stock_M"])
    dest = os.path.join(ARCHIVE_DIR, "runs", run_id)
    try:
        os.makedirs(dest, exist_ok=True)
        for name in ("run_state.json", "pool_final.csv", "pools_measured_1.csv", "pools_measured_2.csv", "pools_measured_3.csv"):
            src = os.path.join(WORKING_DIR, name)
            if os.path.exists(src) and os.stat(src).st_size > 0:
                shutil.copyfile(src, os.path.join(dest, name))
    except Exception:
        pass

def archive_map(run_id, meta):
    """v14: the run's ratio map (ratio_map_run.csv, written by the analysis) appended to ratio_map.csv with the run's
    identity, and the chemistry model's fingerprint of the product as a chemistry_fits.csv row (Phase "Model")."""
    import json
    now = datetime.datetime.now()
    try:
        df = pd.read_csv(os.path.join(WORKING_DIR, "ratio_map_run.csv"))
    except Exception:
        df = pd.DataFrame()
    if not df.empty:
        df.insert(0, "Run_ID", run_id)
        df.insert(1, "Batch_ID", meta["batch_id"]); df.insert(2, "Product_Name", meta["product_name"])
        df.insert(3, "Product_Type", meta["product_type"]); df.insert(4, "Timestamp", now.strftime("%Y-%m-%d %H:%M:%S"))
        df["Algorithm_Version"] = read_text(FILES["VERSION"], "Unknown")
        append_table(os.path.join(ARCHIVE_DIR, "ratio_map.csv"), df.to_dict("records"), list(df.columns))
    try:
        with open(os.path.join(WORKING_DIR, "run_state.json")) as f:
            st = json.load(f)
    except Exception:
        st = {}
    ch = st.get("chem")
    if ch:
        q = ch["q"]
        append_table(os.path.join(ARCHIVE_DIR, "chemistry_fits.csv"), [{
            "Experiment_ID": f"EXP_{now.strftime('%Y%m%d_%H%M%S')}_MODEL", "Run_ID": run_id, "Phase": "Model", "Pool": "",
            "Product_Name": meta["product_name"], "Titrant": "acid + Tris", "Acetate_M": round(q[0], 5),
            "Strong_Ions_M": round(q[1], 5), "Fit_RMS_pH": round(ch["rms"], 4), "N_Points": ch["n"], "Product_Vol_mL": 1.0,
            "Initial_pH": st.get("acid", {}).get("ph0"), "Tris_pKa": np.nan, "Te_Recipe": np.nan, "Te_Real_Pool": np.nan,
            "Real_Recipe_Factor": np.nan, "Tris_Fit_RMS_pH": np.nan, "Carboxyl_M": round(q[2], 5), "Histidine_M": round(q[3], 5),
            "Tris_In_Eluate_M": round(q[4], 5), "Model_Temp_C": round(q[5], 2),
            "Tris_Strength": round(q[6], 4) if len(q) > 6 else 1.0,          # v17: 1.0 = the model used 1:1 stocks
            "Tris_Strength_Fitted": 1 if len(q) > 6 else 0,
            "RMS_Acid_pH": ch.get("rms_acid"), "RMS_Base_pH": ch.get("rms_base"), "RMS_Neutralisation_pH": ch.get("rms_neutralisation"),
            "Base_Equivalence": st.get("base", {}).get("e")}],
            ["Experiment_ID", "Run_ID", "Phase", "Pool", "Product_Name", "Titrant", "Acetate_M", "Strong_Ions_M",
             "Fit_RMS_pH", "N_Points", "Product_Vol_mL", "Initial_pH", "Tris_pKa", "Te_Recipe", "Te_Real_Pool",
             "Real_Recipe_Factor", "Tris_Fit_RMS_pH", "Carboxyl_M", "Histidine_M", "Tris_In_Eluate_M", "Model_Temp_C",
             "Tris_Strength", "Tris_Strength_Fitted", "RMS_Acid_pH", "RMS_Base_pH", "RMS_Neutralisation_pH", "Base_Equivalence"])
    dest = os.path.join(ARCHIVE_DIR, "runs", run_id)
    try:
        os.makedirs(dest, exist_ok=True)
        for name in ("ratio_map_run.csv", "run_state.json"):
            src = os.path.join(WORKING_DIR, name)
            if os.path.exists(src) and os.stat(src).st_size > 0:
                shutil.copyfile(src, os.path.join(dest, name))
    except Exception:
        pass

def write_session(run_id, meta):
    """One row per run: conditions, probe checks, how many entries were typed and how long the run took.
    (v11: the final message's numbers are written by EDP_tecan_gp_master.py Summary.)"""
    now = datetime.datetime.now()
    probe_end = np.nan
    vals = read_text(FILES["PROBE"]).replace('"', ' ').replace(",", " ").split()
    for v in vals:
        f = clean_float(v, treat_zero_as_missing=True)
        if not pd.isna(f):
            probe_end = f
    entries = safe_read_csv(FILES["ENTRIES"], header=0)
    typed = int(len(entries)) if entries is not None else 0
    rejected = int((entries["accepted"] == 0).sum()) if entries is not None and "accepted" in entries else 0
    wells = int(entries.loc[entries["accepted"] == 1, "expected"].sum()) if entries is not None and "accepted" in entries else 0
    started = np.nan
    if run_id.startswith("RUN_"):
        try:
            started = datetime.datetime.strptime(run_id[4:], "%Y%m%d_%H%M%S").strftime("%Y-%m-%d %H:%M:%S")
        except Exception:
            pass
    summary_path = os.path.join(ARCHIVE_DIR, "experiment_summary.csv")
    n_exp = 0
    df_s = safe_read_csv(summary_path, header=0)
    if df_s is not None and "Run_ID" in df_s:
        n_exp = int((df_s["Run_ID"] == run_id).sum())
    minutes = np.nan
    if isinstance(started, str):
        minutes = round((now - datetime.datetime.strptime(started, "%Y-%m-%d %H:%M:%S")).total_seconds() / 60, 1)
    drift = {}                                     # v16: the drift check of plate check 3, in pH and in % of the reading
    try:
        import json
        with open(FILES["STATE"]) as f:
            drift = (json.load(f).get("acid") or {}).get("drift_check") or {}
    except Exception:
        pass
    robot = np.nan                                 # v8: minutes the robot actually worked (timing_log.csv)
    df_tl = safe_read_csv(FILES["TIMING"], header=0)
    if df_tl is not None and "event" in df_tl and "actual_min" in df_tl:
        done = pd.to_numeric(df_tl.loc[df_tl["event"] == "end", "actual_min"], errors="coerce").dropna()
        if len(done):
            robot = round(float(done.sum()), 1)
    append_table(os.path.join(ARCHIVE_DIR, "sessions.csv"), [{
        "Run_ID": run_id, "Started": started, "Finished": now.strftime("%Y-%m-%d %H:%M:%S"), "Minutes": minutes,
        "Operator": meta["operator"], "Batch_ID": meta["batch_id"], "Product_Name": meta["product_name"],
        "Product_Type": meta["product_type"], "Protein_Conc_mgmL": meta["protein"], "Buffer": meta["buffer"],
        "Temperature_C": meta["temperature"], "Acid_Name": meta["acid_name"], "Acid_Conc_M": meta["acid_conc"],
        "Base_Name": meta["base_name"], "Base_Conc_M": meta["base_conc"], "Probe_Check_Start": meta["probe_start"],
        "Probe_Check_End": probe_end,
        "Probe_Drift": round(probe_end - meta["probe_start"], 3) if not pd.isna(probe_end) and not pd.isna(meta["probe_start"]) else np.nan,
        "Drift_Check_pH": drift.get("mean_ph", np.nan), "Drift_Check_Pct": drift.get("mean_pct", np.nan),
        "Drift_Taken_Off_pH": drift.get("applied_ph", np.nan),
        "Titrations": n_exp, "Entries_Typed": typed, "Entries_Rejected": rejected, "Wells_Measured": wells,
        "Algorithm_Version": read_text(FILES["VERSION"], "Unknown"), "Robot_Min": robot,
        "Stop_Min": round(minutes - robot, 1) if not pd.isna(minutes) and not pd.isna(robot) else np.nan}],
        ["Run_ID", "Started", "Finished", "Minutes", "Operator", "Batch_ID", "Product_Name", "Product_Type",
         "Protein_Conc_mgmL", "Buffer", "Temperature_C", "Acid_Name", "Acid_Conc_M", "Base_Name", "Base_Conc_M",
         "Probe_Check_Start", "Probe_Check_End", "Probe_Drift", "Drift_Check_pH", "Drift_Check_Pct", "Drift_Taken_Off_pH",
         "Titrations", "Entries_Typed", "Entries_Rejected",
         "Wells_Measured", "Algorithm_Version", "Robot_Min", "Stop_Min"])
    dest = os.path.join(ARCHIVE_DIR, "runs", run_id)
    try:
        os.makedirs(dest, exist_ok=True)
        for key in ("ENTRIES", "META", "TITRANT", "LOG", "TIMING", "TIMECAL", "DRIFT", "STATE"):
            if os.path.exists(FILES[key]) and os.stat(FILES[key]).st_size > 0:
                shutil.copyfile(FILES[key], os.path.join(dest, os.path.basename(FILES[key])))
    except Exception:
        pass

if __name__ == "__main__":
    main()
