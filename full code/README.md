# Analysis code: individual rates of change in plasma Alzheimer's disease biomarkers

R code for the data preparation and every analysis, table and figure in "Individual rates of
change in plasma Alzheimer's disease biomarkers: reliability, minimal detectable change and
concurrent validity in two cohorts" (Yang et al., submitted to Alzheimer's & Dementia, 2026).
Pipeline version 0.34.0 (`PIPELINE_VERSION` in `00_setup.R`).

## Data

No data are included. Both cohorts are available to qualified investigators on application:

- ADNI, including the FNIH Biomarkers Consortium plasma measures
  (`FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv`): https://adni.loni.usc.edu
- A4/LEARN: https://www.a4studydata.org

The authors are not permitted to redistribute either dataset. The scripts read the released
files unchanged (file names are listed in `00_setup.R`, `RAW_FILES` and `A4_FILES`) and never
write to the raw folders. `make_test_data.R` builds small synthetic files in the same layouts so
the code can be exercised without the data; its output is random and not an estimate of anything. To use
it: `Sys.setenv(FAKE_JM_DIR = "<empty folder>")`, `source("make_test_data.R")` with the working directory set
to the code folder, then point `ADNI_JM_DIR` at `<empty folder>/FAKE_TEST`, `A4_DIR` at `/tmp/fake_a4` (where the A4/LEARN
section writes) and `ADNI_OUT_DIR` at an output folder, and run `run_all.R`.

## Requirements

R 4.5 with `dplyr`, `tidyr`, `readr`, `stringr`, `lubridate`, `purrr`, `tibble`, `lme4`,
`ggplot2`, `patchwork`, `ggrepel`, `scales` and `lcmm` (one section of `08_sensitivity_analyses.R`
only); `rstudioapi` is optional.

## Paths

`00_setup.R` holds three paths, each with a default and an environment variable that overrides it:

| Setting | Meaning | Environment variable |
|---|---|---|
| `BASE_DIR` | folder with the raw ADNI files (read only) | `ADNI_JM_DIR` |
| `OUT_ROOT` | folder under which `pipeline_out/`, `a4learn_out/` and `combined/` are written | `ADNI_OUT_DIR` |
| `A4_DIR` | folder with the A4/LEARN release (read only) | `A4_DIR` |

Set them with `Sys.setenv()` before sourcing, or edit the defaults in `00_setup.R`.

## How to run

The scripts live in one flat folder and source `00_setup.R` by file name, so keep them together.
`run_all.R` runs everything in order for both cohorts; open it in RStudio and press Source, or
source it by its full path. `STEPS` at the top selects files by their two-digit prefix. The data
build (01 to 03) takes a few minutes; the analyses (04 to 11) several hours, most of it the
A4/LEARN bootstrap in 07 and the power simulation in 08. The ADNI bootstraps in 09 save every draw
and resume if interrupted; the A4/LEARN bootstrap in 07 does not. `N_BOOT` sets the number of
participant resamples for both (default 1,000, as in the paper). `N_BOOT = 0` skips the resampling in
07 and in the first section of 09; the later sections of 09, and 10 and 11, need the saved draws, so a
run without a bootstrap should set `STEPS` to `"01"` to `"08"`.

Each file re-sources `00_setup.R` at the start of its sections and sets the cohort where it matters,
so any single file can also be sourced on its own once the files before it have run, with the working
directory set to the code folder (`setwd()`); `03_analysis_frames.R` takes the cohort from `COHORT` and is
run twice by the driver, once per cohort.

## Files, in run order

| File | What it does |
|---|---|
| `00_setup.R` | Settings (paths, cohort switch, raw file names, constants, helper functions), the figure theme, `define_units()` (the analysis units by cohort), `load_change_objects()` (the frames and helpers shared by the change-model sections) and the driver helper `run_steps()`. Sourced at the start of every section. |
| `01_data_adni.R` | ADNI data build: the plasma table (one row per person, visit and measure), participants, reference measures (amyloid PET, CSF, cognition, MRI, diagnosis) matched to the nearest plasma visit, modifiers (creatinine, eGFR, BMI), and the check tables. Writes `pipeline_out/`. |
| `02_data_a4learn.R` | The A4/LEARN release in the same layout: plasma pTau217 series, participants, reference measures, modifiers. Writes `a4learn_out/`. |
| `03_analysis_frames.R` | The master table (one row per person per plasma visit) and the analysis frames: persons x platforms x occasions per biomarker on the log, z and baseline-z scales, per-person slopes, per-person criterion slopes over the plasma span; for A4/LEARN also the frame check against ADNI. Run once per cohort. |
| `04_level_and_change_models.R` | Level variance components and agreement between platforms; the multi-platform change model, reliability of a one-platform slope and of the mean of platforms, overall and by amyloid status; minimal detectable change on the z scale; slope agreement; validity of single-platform and pooled slopes against criterion change; sampling requirement (T1 to T6). Leave-one-platform-out and the split of the occasion-level variance into the part shared by every platform at a blood draw and the assay-specific part (T19, T20). The upper and lower bounds on the gain from averaging (R1), each measure fitted alone (R4), the pTau217 slope ICC on the three immunoassays (R6). |
| `05_minimal_detectable_change.R` | Minimal detectable change in each assay's own units: within-person CV, mean rate of change, the 95% minimal detectable change between two measurements and for a three-sample slope as fold change and per cent, as a rise and as the equivalent fall, and the years of mean amyloid-positive change needed to exceed it (T14, R7). |
| `06_sampling_design_and_trials.R` | Trial sample sizes per platform and for the mean of platforms (T15); reliability under alternative sampling schedules and the span needed for target reliabilities (T16, T16b, F6); validity by observed spread of sampling times (T17); the mean of platforms under both bounds (R2); trial sizes on amyloid-positive variance components (R3); validity by SST fifth allowing for the precision of both the plasma and the Centiloid slope (R14b). |
| `07_a4learn.R` | A4/LEARN: the random-slope change model per unit and slope reliability (T2); the baseline to week-12 test-retest and its comparison with the long change model's residual (T7, F1); the within-person design test (T8, T8b, F8); minimal detectable change (T3, T14); validity (T5); sampling requirement (T6, T16b, F3); the participant bootstrap (T9). The fourth (extension) sample: four slope designs from each participant's own samples, with participant-bootstrap intervals for the three design contrasts (`combined/T9`, `T9b`; Figure 4C). The line through baseline, week 240 and the extension sample against twelve-week variance, and the curvature of the mean trajectory (R8a, R8b). Plasma samples per participant in the placebo arm and LEARN (R14a). |
| `08_sensitivity_analyses.R` | C2N pTau217 floor (T10, T10b, R5a to R5c), covariate adjustment (T11), scale robustness (T12), level agreement per platform pair (T13), latent-process harmonisation with `lcmm::multlcmm` (T8, T8b, T8c, F5), multitrait-multimethod analysis of the Roche and Quanterix block (T9, T9b, T9c), simulated power for a platform-specific slope variance (T21, F13), sample storage time (R9a to R9c), hippocampal volume raw against ICV-adjusted (R16). |
| `09_bootstrap_intervals.R` | The likelihood-ratio test of the platform-specific slope variance and residual diagnostics (T2d, T2e, T2f, F4); the participant bootstrap, 1,000 resamples saved draw by draw (T7); and every interval and paired comparison derived from the same resamples: minimal detectable change and trial sizes (T3b, T15b), the averaging lower bound, fitted-alone reliabilities, amyloid-positive trial sizes and immunoassay ICC (R10), span needed (R11), paired differences and ratios (R12, R12b), the amyloid-positive minus amyloid-negative validity difference (printed), every paired comparison quoted in the text and in Table S26 (R13a, R13b, R13c), the reduction in trial size from five rather than three samples (R15a), the pooled minus single-platform pTau217 validity against Centiloid change (R15b), and intervals for the whole-cohort trial sizes of the pTau217 mean of platforms under the lower bound (R17). Also the observed-data validity by amyloid group and the Pearson versions of the design-test correlations for both cohorts (T7c, T17p, T8p, T8bp). |
| `10_tables.R` | Manuscript tables as CSV and markdown (Table 1 participant characteristics, Table 2 variance components and reliabilities, Table 3 minimal detectable change, the supplementary assay tables) and the cross-cohort tables C1 to C7 with figure F30. No model fits. |
| `11_figures.R` | All figures: working versions, journal versions with legends, the manuscript composites Figure 1 to Figure 4 and Supplementary Figures S1 to S3, S5 and S6, with `manuscript_figure_legends.md`. No model fits. |
| `run_all.R` | The driver. |
| `make_test_data.R` | Synthetic test files in the layouts of the real releases. Not part of the pipeline. |

Each file is divided into sections; each section carries the header of the script it was
developed as, so the provenance of every table is visible.

## Outputs

All under `OUT_ROOT`: `pipeline_out/` (ADNI frames and checks) and `pipeline_out/paper2/`
(ADNI tables T1 to T21 and R1 to R12b, figures, bootstrap draws); `a4learn_out/` and
`a4learn_out/paper2/` (the same for A4/LEARN); `combined/` (C1 to C7, the fourth-sample table,
F30). `run_all.R` writes a console log to `pipeline_out/`.

## Checksums

`MANIFEST.sha256` lists the SHA-256 checksum of every script and this README as archived.
