# EDP pH adjustment platform

Automated titration and pH adjustment of bioprocess intermediates on a Tecan Freedom EVO (Freedom EVOware 2.8), with a
Python analysis that turns each run into validated acid, base and neutralisation ratios. This is the data-collection
stage of a predictive, machine-learning-integrated pH adjustment platform (AstraZeneca placement project, Purification
Process Sciences, BPD).

> **Internal.** The scripts, documentation and examples contain AstraZeneca product codes, batch IDs, run data and
> network paths. Keep this repository private or internal unless these have been removed and sharing has been approved.

**Status:** v17 (22 October 2026) is built and tested in the simulator, and has not yet run on the robot. v16 ran end to
end on 22 October 2026.

Start with [`EDP_platform_handoff.md`](EDP_platform_handoff.md). It explains the platform, the analysis maths, every
robot run so far and the simulated results.

## What a run does

- **One run, up to three titrations:** the acid titration (1 M acetic acid), plus two optional Tris titrations chosen at
  the start: VI neutralisation (Tris into acidified product) and the base titration of the eluate (Tris straight into
  the eluate). 24-well plates, 3 mL per well, read with a macro pH probe, one typed entry per plate.
- **Analysis:** each titration is fitted with physically shaped equations and an error model; typos are dropped;
  validation wells are made at the calculated ratios and the fit is redone with them; probe drift is measured and
  corrected; and a charge-balance chemistry model is fitted to every reading of the run.
- **Output:** a ratio map (acid, base and neutralisation ratios from any pH to any other), and optional custom pools at
  any pH in the measured range, trimmed until they read within 0.03 pH.
- **Data for machine learning:** every run is archived as tables (`wells.csv`, `measurements.csv`,
  `chemistry_fits.csv`, `ratio_map.csv` and others) and as an Excel workbook per run.

## Repository layout

| Path | Contents |
|---|---|
| `esc/` | `EDP_Titration_python_ML_hybrid.esc`, the EVOware script (v17) |
| `python/` | the six EDP Python files that run on the robot PC (`...\database\tecan_data`), all tagged `EDP-v17-2026-10` |
| `tools/` | the script generator (`gen17.py`), the envelope check, the EVOware simulator, scenarios, benches and replays of real runs |
| `reference/` | scripts that have run on the robot (the envelope references, including the original GCB script) and earlier versions in `v10_previous/` to `v16_previous/` |
| `EDP_platform_handoff.md` | full documentation |

## Installing on the robot PC

1. Back up the current script and Python files (v16 is also in `reference/v16_previous/`).
2. Copy all six `python/*.py` files into `C:\ProgramData\Tecan\EVOware\database\tecan_data\`. The script checks at the
   start that all six are the same version.
3. Copy `esc/EDP_Titration_python_ML_hybrid.esc` into `C:\ProgramData\Tecan\EVOware\database\Scripts\` as a file.
   Never open and save an `.esc` in a text editor: it damages the well selections. Delete any
   `EDP_Titration_python_ML_hybrid.rec` next to it.
4. Make the 0.1 M Tris: 10.00 mL of the same 1 M Tris stock into a 100 mL volumetric flask, water to the mark (Class A
   pipette and flask: its strength enters every base ratio).
5. Set up the deck as in the handoff (section 3), plus the 0.1 M Tris in the Product 2W RIGHT well (about 50 mL, at least
   25 mL). Nothing in the door row.
6. Open the script in EVOware, check the worktable and save.
7. On the first run, check that the final message says the workbook went to folder 1 (the H: drive).

## Development loop

Run from the repository root, with Python 3 and numpy, scipy, pandas and openpyxl:

```
python tools/gen17.py reference/EDP_v12_12Oct_2026.esc esc/EDP_Titration_python_ML_hybrid.esc
python tools/envelope.py esc/EDP_Titration_python_ML_hybrid.esc
python tools/scen17.py
```

- `gen17.py` builds the script, `envelope.py` checks that it stays within what has already run on the robot, and
  `scen17.py` runs the 23 simulated scenarios.
- Benches: `run_lab17.py` (the three titrations and the ratio map), `trim_lab17.py` (Tris pools) and
  `option_minutes17.py` (robot time per start-prompt combination).
- Replays of real runs through the analysis: `replay_21oct.py` and `replay_22oct.py`.
- Section 10 of the handoff lists every tool and the EVOware rules any script change must follow.

## v17 in brief (22 October 2026)

- **0.1 M Tris for the base wells.** On 22 October the base wells near pH 7 got 47-65 µL of 1 M Tris, and 1 µL moved
  the pH by 0.05-0.09: the robot's small dispenses, not the calculation, limited the base validation. Every base well
  that needs at most one 980 µL portion of the 0.1 M Tris now uses it, and small Tris pool trims do too. Ratios stay in
  mL of 1 M Tris (mmol) per mL of eluate. Simulated base validation: 0.075 to 0.037 pH RMS on the 22 October eluate,
  0.27-0.29 to 0.04-0.06 on weakly buffered eluates.
- **Results workbook to the H: drive** (by H:'s network name if Python cannot see the letter; Public Documents if H:
  cannot be reached). The final message says which folder.
- **For machine learning:** `ML_data_archive\wells.csv` holds every well and pool reading of every run, with its
  composition, Tris stock, total volume, pH, whether the fit left it out, and when it was typed.
- **Probe reminders:** plate checks with base wells say to rinse the probe and wait for a steady reading.
- **Analysis:** per-well Tris stock in the fit's weights and in the chemistry model, a 3-sigma outlier rule for base
  readings, and a fallback for the VI fit that the simulator found could fail on weakly buffered eluates.

## A note on `.esc` files

`.gitattributes` marks `.esc` files so that Git never changes their line endings or characters. The first line of an
`.esc` is a checksum of the rest of the file, and the well selections use special characters that a text editor or a
line-ending conversion would damage.
