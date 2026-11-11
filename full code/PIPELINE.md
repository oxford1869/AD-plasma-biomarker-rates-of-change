# The analysis pipeline, step by step

A specification of what the code in this folder does, written for someone checking the code
against it. Each step gives its inputs, the method or model as the code implements it, the
files it writes, and where the result appears in the manuscript. File names in the left margin
are the scripts in this folder; the "section" numbers are the numbered sections inside each
script (every section begins with a banner naming the script it was developed as). Table names
such as T2 or R10 are the CSV files in `pipeline_out/paper2/` (ADNI) or `a4learn_out/paper2/`
(A4/LEARN); the supplementary table numbers (S1, S9b, ...) are those of the manuscript.

Notation used throughout: a participant *i*, a platform (assay) *p*, a blood draw (occasion)
*o*, time *t* in years since the participant's first plasma sample. All models are fitted on
`z_log_aligned`: the natural log of the delivered value, standardised within each measure
(mean 0, SD 1 over all observations of that measure), with the sign set so that higher means
more abnormal (the Aβ42/40 ratio is reversed). Mixed models are fitted with `lme4::lmer`,
REML unless stated, optimiser bobyqa for the change, level and sample models and lme4's default
for the small ICC and test-retest models. "SST" is the spread of a participant's sampling times,
Σ(t − mean t)², and `median_sst` is its median over participants with at least two draws
(12.855 year² in ADNI); in stage 8 sections 1 and 3 and stage 9 sections 1, 4 and 7 the median is
taken over all participants in the frame, which is the same value when every participant in the
frame has two or more draws (stage 9 section 4 checks this by requiring its spans to equal T16b).

---

## Stage 0. Settings (`00_setup.R`)

Paths (`BASE_DIR`, `OUT_ROOT`, `A4_DIR`), the cohort switch `COHORT` ("adni" or "a4learn",
which chooses between `pipeline_out/` and `a4learn_out/`), the raw file names, constants
(`VISIT_GAP_DAYS = 7`, `CENTILOID_CUTS` 20 and 37, `MATCH_WINDOW_DAYS` criteria 180, modifiers
365, diagnosis 180, `CRITERION_SLOPE_PAD_DAYS = 180`), helper functions (date parsing, column
picking, `slope_by()` for per-participant least-squares slopes, `say()`), the measure lookup
table (which column of the FNIH file is which biomarker on which platform, and its direction),
the figure theme, and two functions used by later steps:

- `define_units(long)`: the analysis units by cohort. ADNI: the five pTau217 measures (ALZpath,
  C2N concentration, C2N %pTau217, Fujirebio, Janssen) on the same participants. A4/LEARN: one
  assay (Lilly MSD), so the units are subsamples: A4 placebo, LEARN, A4 solanezumab, placebo
  plus LEARN (primary), all three arms. Also the floor rule (`always_floor`: participants whose
  C2N pTau217 is at the lower limit of quantification at every draw).
- `load_change_objects()`: the frames with the factor columns the change models need
  (`RID_f`, `arm_f`, `batch_f`, `sample_f` = participant × draw), `median_sst`, the three trial
  schedules, the mean amyloid-positive slope per measure, the two model formulas of step 4
  (reduced model and sample model), and three helpers (`n_slope`, `fit_info`, `vget`).

---

## Stage 1. ADNI data build (`01_data_adni.R`, `COHORT = adni`)

Inputs: the FNIH Biomarkers Consortium plasma file (`FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv`),
ADNIMERGE, the Berkeley amyloid PET file (Centiloids), the UPenn Roche Elecsys CSF file,
LABDATA, PTDEMOG, MMSE, VITALS. Nothing in the raw folder is written to.

| Section | What it does | Writes |
|---|---|---|
| 1 (plasma) | One row per participant, visit and measure. Dates within 7 days are one visit; `visit_no` numbers every sample a participant has (shared across platforms), `occasion_no` the samples each platform measured. ng/mL converted to pg/mL. Lower and upper limit flags from TESTVALUE and COMMENTS. Fujirebio, Roche and Quanterix Aβ42/40 ratios derived from the components. Quanterix-only samples flagged. | `plasma_long.rds`, `plasma_wide.rds/.csv` |
| 2 (participants) | One row per participant: number of plasma visits, first and last plasma date, span; `in_core_5lab` (at least one plasma visit measured by any of the five main, non-Quanterix, platforms); baseline characteristics from ADNIMERGE and PTDEMOG; age at first plasma sample; diagnosis at first plasma. | `participants.rds/.csv`, `checks/flow_counts.csv` |
| 3 (reference measures) | Amyloid PET (Centiloid, summary SUVR, status), CSF (Aβ42, Aβ40, tau, pTau181 and ratios), cognition (MMSE, ADAS-Cog13, mPACC, CDR-SB), MRI (hippocampus, ICV, adjusted hippocampal volume) and diagnosis in one tall table, each assessment given its nearest plasma visit and the gap in days. | `criteria_long.rds/.csv`, `dx_long.rds` |
| 4 (modifiers) | Creatinine, eGFR (CKD-EPI 2021 and 2009), weight, height, BMI at every visit. | `modifiers_long.rds/.csv` |
| 5 (checks) | Assertions and coverage tables on the plasma table. | `checks/*.csv` |

## Stage 2. A4/LEARN data build (`02_data_a4learn.R`, `COHORT = a4learn`)

Inputs: the A4/LEARN release (external data: pTau217, amyloid and tau SUVR, labs; derived data:
SUBJINFO, SV, PACC; raw data: ptdemog, vitals, cdr, mmse; the visits data dictionary). Days since
T0 are mapped onto a nominal calendar; only differences are used. Centiloid = 183 × SUVR − 177
(Avid florbetapir pipeline, whole-cerebellum reference).

| Section | What it does | Writes |
|---|---|---|
| 1 | The Lilly MSD pTau217 series in the `plasma_long` layout: visit code and kind, substudy, arm, on-treatment flag, the baseline to week-12 short-interval pair flag; raw instrument values kept for out-of-range results with flags. | `a4learn_out/plasma_long.rds` |
| 2 | Participants: age at first plasma draw, sex, education, APOE, arm, treatment dates, BMI; diagnosis CN for all. | `a4learn_out/participants.rds` |
| 3 | Reference measures in the `criteria_long` layout: florbetapir SUVR and Centiloids, amyloid status, tau PET, CSF, MMSE, PACC, CDR. | `a4learn_out/criteria_long.rds` |
| 4 | Modifiers: creatinine and eGFR, HbA1c, CRP, haemoglobin, vitals, BMI; diagnosis (CN throughout). | `a4learn_out/modifiers_long.rds`, `a4learn_out/dx_long.rds` |

## Stage 3. Master table and analysis frames (`03_analysis_frames.R`, run once per cohort)

| Section | What it does | Writes |
|---|---|---|
| 1 (master) | `master_long`: every measurement from every source in one tall table. `master_visit`: one row per participant per plasma visit with the plasma measures as columns, the nearest reference measure within 180 days, the nearest modifier within 365 days, diagnosis within 180 days, and per-participant anchors: baseline Centiloid, `amyloid_bl_cl20` (baseline Centiloid above 20, a strict inequality, the amyloid status used throughout; the manuscript says "above 20 Centiloids"), `amyloid_bl_cl37`, `amyloid_ever_cl20`. | `master_long.rds`, `master_visit.rds/.csv` |
| 2 (analysis frame) | Sample: participants with `in_core_5lab`, Quanterix-only samples excluded. `frame_paper2_long`: participant × platform × occasion per biomarker with `log_value`, `z_log` (standardised within measure), `z_log_aligned` (× direction), `z_bl` (standardised on occasion 1 only), `batch` = platform:RUN_NO, `complete_3x` (three samples on every platform), covariates and amyloid status. `frame_paper2_slopes` (and `_wide`): per participant, per measure ordinary least-squares slope of `z_log_aligned` on `plasma_time` (`slope_by`: two or more draws with positive time spread; `n_obs`, `span`, first and last value kept). `frame_paper2_criteria`: per participant least-squares slope of each reference measure over the plasma span padded by 180 days at each end (Centiloid, CSF pTau181/Aβ42, CSF Aβ42, CSF pTau181, MMSE, ADAS-Cog13, mPACC in its digit and Trails B versions, CDR-SB, adjusted hippocampal volume, and, where present, PACC and tau PET). | `frame_paper2_long.rds`, `frame_paper2_slopes.rds`, `frame_paper2_slopes_wide.rds/.csv`, `frame_paper2_criteria.rds`, `checks/paper2_design.csv` |
| 3 | `frame_paper3_long` and `frame_paper3_person`: participant × visit × measure with every exposure scale as columns, and person-level outcomes; used by the fourth-sample and R8 analyses (stage 7, sections 2 and 3) and by the A4 flow counts (section 4). | `frame_paper3_long.rds`, `frame_paper3_person.rds/.csv`, `frame_paper3_specs.csv`, `checks/paper3_outcomes.csv` |
| 4 (A4 only) | Column-by-column comparison of the A4/LEARN frames with the ADNI frames; flow counts by substudy and arm; the analysis samples (primary placebo plus LEARN; sensitivity all arms). | console and `checks/` |

---

## Stage 4. Level and change (`04_level_and_change_models.R`, ADNI)

**Section 1 (T1 to T6).**

*Level model, per biomarker (T1, Figure 2A, Table 2).* Crossed random effects on
`z_log_aligned`, all platforms and all visits:

    y ~ 1 + (1|participant) + (1|platform) + (1|participant:platform) + (1|sample)
          + (1|platform:visit) + (1|batch)

where `sample` is participant × draw (the aliquot shared by every platform) and `batch` is the
analytical run within platform. Each component is reported as a share of the total. The
participant component is between-participant variance in true level; the relative G coefficient
(T1b) is participant / (participant + participant × platform + sample + residual); an absolute
G (participant / total) is written beside it.

*Change model, per biomarker (T2, Table 2, Figure 2B).*

    y ~ platform * t + (1 + t | participant) + (1 + t || participant:platform) + (1|batch)

Components: var(true slope) = the participant random-slope variance; var(platform-specific
slope) = the participant:platform random-slope variance; residual. Reliability of a
participant's observed three-sample slope on one platform

    rel_one = var_true / (var_true + var_platform_slope + var_residual / median_sst)

and for the mean of the *k* platforms measuring the biomarker

    rel_mean = var_true / (var_true + (var_platform_slope + var_residual / median_sst) / k).

The same model is refitted in amyloid-negative and amyloid-positive participants (T2b) and,
for pTau217, without the participants always at the C2N floor (T2c). These joint-model
reliabilities count variation shared across platforms at a draw as true change and are the
paper's *upper bounds* (see section 2 and R1).

*Minimal detectable change on the z scale (T3).* Per platform, a model with participant
intercepts and slopes, `y ~ t + (1 + t | participant)`; the residual SD is the within-person,
within-platform measurement SD (`sem`). MDC95 between two measurements = 1.96 × √2 × sem;
MDC95 of a three-sample slope = 1.96 × sem / √`median_sst` (until November 2026 the SST of three evenly
spaced draws over the median span was used here, 12.5 against 12.855; the bootstrap in stage 9 section 2
recomputes the median SST on each resample). The share of
participants whose own slope exceeds the slope MDC is reported in each direction.

*Agreement of slopes between platforms (T4, Figure S2, Table 2 ICC).* Pearson and Spearman
correlations and tertile kappa for each pair of measures; the ICC of the per-participant
slopes across measures from `slope ~ 1 + (1|participant) + (1|measure)`, single measure
(ICC(2,1)) and mean of measures (ICC(2,k)). The pooled participant slope (T4b) is the fixed
plus random slope for the participant from the change model.

*Validity (T5, Figure S4).* Spearman correlation of each single-measure slope and of the
pooled slope with the participant's concurrent slope in Centiloid, CSF pTau181/Aβ42, MMSE,
ADAS-Cog13, mPACC and adjusted hippocampal volume; the sign of hippocampal volume, mPACC and
MMSE is reversed so that positive always means the plasma slope tracks worsening.

*Sampling requirement (T6).* For *k* equally spaced draws over *S* years,
SST = Σ(t − mean t)² and reliability as above with that SST in place of `median_sst`; a grid of
k = 2, 3, 4, 5, 6 draws and S = 2, 4, 6, 8, 10 years.

**Section 2 (T19, T20).** *Leave one platform out:* the reduced change model (section 1's
model without the platform-specific slope term, see stage 9 section 1) refitted, for each
biomarker with three or more platforms (pTau217), with each platform dropped; reported: correlation of the resulting pooled slope with the
all-platform pooled slope, the implied single-assay reliability, validity against Centiloid
change. *Occasion-variance split:* the reduced model with a participant × draw random effect
added,

    y ~ platform * t + (1 + t | participant) + (1 | participant:platform) + (1 | sample) + (1|batch)

(the "sample model"). The former residual splits into the part shared by every platform at a
draw (`var_shared`) and the assay-specific part (`var_assay`). Reported per biomarker with an
identifiability flag (needs three or more platforms and a non-singular fit): the pooled-slope
reliability with the shared part retained, `var_true / (var_true + (var_shared +
var_assay / k) / median_sst)`, and with all variation averaged.

**Section 3 (R1, Table S10b inputs, Results).** The upper bound (T2) and the lower bound
(T20: one platform `var_true / (var_true + (var_shared + var_assay) / median_sst)`; mean of
platforms with the shared part retained) on the reliability gain from averaging, per
biomarker, with the identifiability note.

**Section 4 (R4, Results, Figure 2B grey bands).** Each of the 15 measures fitted alone with
`y ~ t + (1 + t | participant)` (the T3 model): reliability
`var_slope / (var_slope + sd_residual² / SST)` at the ADNI median SST and at the A4 placebo
SST (13.533), beside the joint-model value. This is the like-for-like comparison with the
single-assay A4/LEARN cohort and is not a bound.

**Section 5 (R6).** The pTau217 slope ICC recomputed on the three immunoassays only (ALZpath,
Fujirebio, Janssen), beside the five-measure value, because the two C2N readouts share one
assay.

## Stage 5. Minimal detectable change in assay units (`05_minimal_detectable_change.R`)

**Section 1 (T14, Table 3, Figure 3).** For each measure, `sd_log` = SD of the log value (what
`z_log` divided by), so a z quantity *q* is a log change of `q × sd_log` and a fold change of
`exp(q × sd_log)`. Reported per platform: within-person CV = √(exp((sem × sd_log)²) − 1); mean
change per year in all participants (the fixed-effect slope of the per-platform T3 model) and in
amyloid-positive participants (the mean least-squares slope of amyloid-positive participants with
three or more samples), undoing the sign alignment so per
cent changes read in the delivered direction; MDC95 between two measurements as fold change
and per cent; MDC95 of a three-sample slope per year; years of mean amyloid-positive change
needed before two measurements differ by more than the two-point MDC95
(`mdc95_two_points / |mean amyloid-positive slope|`); the between-participant CV of
participant means.

**Section 2 (R7, Table 3, Abstract).** The two-point MDC95 as a rise
(`mdc95_two_points_pct`) and as the equivalent fall, `100 × (1 − 1 / fold)`.

## Stage 6. Sampling design and trials (`06_sampling_design_and_trials.R`)

**Section 1 (T15 to T17; Figure 4A; Tables S9, S10).**

*Trial sample size (T15).* Two-arm trial, plasma slope as endpoint, treatment slows the mean
amyloid-positive rate by 25% or 50%; three schedules (18 months with 4 draws; 24 months with 3
draws; 24 months with 5 draws), each with its SST. Per platform: var(observed slope) =
var_true (T2) + sd_residual² (T3) / SST; for the mean of *k* platforms the error term is the
mean of the platforms' sd_residual² divided by *k*. Participants per arm

    n = ceiling( 2 (z_0.975 + z_0.80)² × var(observed slope) / (effect × |mean amyloid-positive slope|)² ),

and the same for change from baseline over the schedule. The components come from the whole
ADNI cohort, so the numbers describe an ADNI-like population.

*Sampling schedules (T16, T16b, F6).* Reliability depends on the design only through SST, so
evenly spaced, front-loaded, back-loaded and replicate-draw schedules and the observed ADNI
design are compared on the one curve `var_true / (var_true + var_residual / SST)` (Figure 4A).
The span needed for reliability 0.7, 0.8 and 0.9 with *k* evenly spaced draws:
`sst_needed = ratio × target / (1 − target)` with `ratio = var_residual / var_true` (divided by
the number of platforms for the mean), and `span = √(sst_needed × 12 (k − 1) / (k (k + 1)))`.

*Validity by observed SST (T17, Table S14).* Participants' single-measure slopes against
concurrent Centiloid change, by quintile of observed SST, with the reliability the components
predict at each quintile's median SST.

**Section 2 (R2; Table S9b, S10).** The same trial-size and span calculations for the mean of
platforms under both bounds: the upper bound exactly as T15; the lower bound from the sample
model refitted unrounded, with error `var_shared + var_assay / k`. The upper rows are checked
to reproduce T15.

**Section 3 (R3; Figure 4B, Table S10, Results).** Trial sizes with the variance components
refitted in amyloid-positive participants only: the reduced model and the sample model per
biomarker, and the per-platform T3 model, in amyloid-positive participants; fit diagnostics
(singular, converged) are written (R3b). Per platform and for the mean of platforms under
each bound, the same *n* formula with the amyloid-positive components, beside the
whole-cohort T15 value.

## Stage 7. A4/LEARN (`07_a4learn.R`)

**Section 1 (`COHORT = a4learn`).** Per unit (A4 placebo, LEARN, A4 solanezumab, placebo plus
LEARN, all arms) and floor rule (as delivered; excluding always-at-floor): the random-slope
change model on the Lilly pTau217 series, `y ~ t + (1 + t | participant)` for single-arm units and
`y ~ arm * t + (1 + t | participant)` for the combined units (so the slope variance is within arm), var(true slope),
residual, reliability at the unit's observed median SST (T2; Table S12). Baseline to week-12
test-retest in A4 (the short-interval pair): within-person SD and CV on the log scale, ICC(2,1),
Bland-Altman limits, two-point MDC95, share of pairs beyond it, and the ratio of the long
change model's residual variance to the 12-week residual variance (T7, F1; Results). Design
test within participants with baseline, week-12 and week-240 draws: the least-squares slope
from each subset of their own draws (baseline + week 12; baseline + week 240; week 12 + week
240; all three), with the subset's median SST, predicted reliability, correlation with the
three-draw slope and validity against Centiloid change (T8, F8; Table S13a); validity by
observed SST band (T8b; Table S13b). MDC per unit on the z and fold scales (T3, T14); validity
of least-squares and BLUP slopes against Centiloid, PACC, MMSE, CDR-SB, tau PET and CSF change
(T5); sampling requirement and span needed (T6, T16b, F3); participant bootstrap of the
slope reliability per unit and of the 12-week test-retest ICC(2,1) and within-person CV, `N_BOOT`
draws (default 1,000), seed 20260906, no cache (T9).

**Section 2 (`COHORT = adni`, reads `a4learn_out/`).** Participants with a fourth (extension)
plasma sample. Four designs from each participant's own samples: A first, second, third
(the standard window); B first, second, fourth (the longer window); C first, third, fourth; D all four.
For each: median SST, the reliability predicted from the A4 placebo anchor (reliability 0.601
at SST 13.533, odds linear in SST), the reliability refitted within these participants with
the random-slope model, and validity against Centiloid change; the three contrasts A vs B (window), B vs D (fourth
sample), B vs C (position of the middle sample) as differences in predicted reliability, refitted
reliability and validity, each with a 95% percentile interval from a participant bootstrap within
the same participants: draw *b* resamples participants with replacement (seed `P2_SEED4` + *b*,
default 20260908), each copy with its own id, and recomputes the three quantities per design on the
resample; `P2_BOOT4` draws (default 1,000). Writes `combined/T9_fourth_sample.csv`,
`T9b`, `F9` (Figure 4C, Table S13c).

**Section 3 (R8; Results, Methods).** In the four-sample participants whose second sample is the week-12 draw (before 0.5 years), the signed
residual of the line through samples 1, 3 and 4 (week 12 left out), scaled so its square is
the residual sum of squares, and the baseline to week-12 difference, both centred within arm;
the ratio of the residual variance of the line to the 12-week variance (per-measurement), with
a participant bootstrap interval (`P2_BOOT_R8` draws, default 2,000, seed 20261106), for all
arms (primary) and the placebo arm (R8a). Curvature of the mean trajectory: linear against
quadratic fixed effect of time, random intercept and slope, ML, likelihood-ratio test and the
quadratic coefficient with its interval (R8b).

## Stage 8. Sensitivity analyses (`08_sensitivity_analyses.R`, ADNI)

| Section | Analysis | Writes |
|---|---|---|
| 1 | All fits in this section use the reduced change model (no platform-specific slope). *C2N floor* (the pTau217 model under five treatments of the lower-limit values: as delivered; excluding always-at-floor participants; excluding floor observations; half-floor substitution; cross-platform multiple imputation, each floor value drawn from below the limit using a regression on the same draw's ALZpath, Fujirebio and Janssen values, `N_IMP` = 10 imputations pooled). *Covariate adjustment*: the reduced model with age at first plasma (centred at 72), sex, APOE ε4 and diagnosis, each also interacting with time. *Scale robustness*: the change decomposition on the raw-scale z and on the baseline-standardised log scale. *Level agreement per platform pair* at each visit: mean log difference and limits of agreement, fold ratio, centred ICC (Figure S1). | T10, T10b, T11, T12, T13 |
| 2 | *C2N floor by amyloid status*: observations and participants at the floor by amyloid group (R5a); the two C2N readouts (concentration and %pTau217) refitted alone in three samples (all data; excluding floor draws; participants never at the floor) with reliability, MDC rise and fall, years to exceed, trial size (R5b); agreement of the two C2N slopes in the same three samples (R5c). | R5a, R5b, R5c |
| 3 | *Latent-process harmonisation*, `lcmm::multlcmm`: per biomarker the platforms are measurements of one latent trajectory with platform-specific linear links and errors, participant latent intercept and slope; per-platform error SD on the latent scale (T8), latent participant slopes against the lme4 pooled and single-platform slopes (T8b), latent-scale reliability (T8c), link functions (F5). | T8, T8b, T8c, F5 |
| 4 | *Multitrait-multimethod*, Roche and Quanterix block (pTau181, GFAP, NfL, Aβ42/40 on the same samples): level model `y ~ analyte * vendor + (1|participant) + (1|participant:analyte) + (1|participant:vendor) + (1|participant:analyte:vendor) + (1|participant:visit) + (1|participant:visit:vendor) + (1|batch)`; slope model `y ~ analyte * vendor * t + (1 + t || participant) + (1 + t || participant:analyte) + (1 + t || participant:vendor) + (1|participant:analyte:vendor) + (1|participant:visit) + (1|batch)`; the participant:vendor terms are the method factor; likelihood-ratio tests of the vendor terms. | T9, T9b, T9c |
| 5 | *Power for a platform-specific slope variance*: data simulated from the reduced model's components for pTau217 on the real participants, times and platform structure, a platform-specific slope SD injected at 0, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6 and 0.8 of the true between-participant slope SD (0 to 0.64 of the variance), full and reduced models refitted, boundary-corrected likelihood-ratio test; power = share of `POWER_NREP` (200) replications rejecting at 0.05; the smallest ratio detectable with 80% power (Table S17, Figure S8). | T21, F13 |
| 6 | *Storage time*: years from collection to assay (all assays January to March 2024) per sample (R9a); the standardised log value of the first sample regressed on storage years, adjusted for age at the sample, amyloid status and diagnosis, the coefficient being −d in SD units and the implied slope bias +d per year as a per cent of the mean annual rate (all and amyloid-positive), with 95% CIs for the coefficient and for the share of the amyloid-positive rate (R9b); slopes regressed on storage (R9c). Storage years is collinear with enrolment era and batch, so the coefficient bounds a bias rather than attributing it (Table S26). | R9a, R9b, R9c |

## Stage 9. The bootstrap and every interval (`09_bootstrap_intervals.R`, ADNI)

**Section 1 (T2d, T2e, T2f, F4, T7; Table 2 intervals, Figure 2B, Table S25; the residual plots F4 become Figure S3 when stage 11 copies them).**
*Likelihood-ratio test* of the platform-specific slope variance: the change model with and
without `(t || participant:platform)`'s slope term, ML fits, statistic 2Δlog-likelihood,
boundary-corrected p = 0.5 × P(χ²₁ > statistic); if p > 0.05 for every biomarker the reduced
model (no platform-specific slope) is used in the bootstrap. *Residual diagnostics*: residual
SD by platform, fitted tertile and floor status, QQ summary. *Participant bootstrap*: draw *b*
resamples participants with replacement with `set.seed(20260904 + b)`, each copy given a new
id; the change model (reduced) is refitted per biomarker in the full resample and in its
amyloid-negative and amyloid-positive subsets; saved per draw: G coefficients, `rel_one` and
`rel_mean` (all and by amyloid status), slope ICCs, Spearman validity of the pooled slope
against Centiloid, CSF pTau181/Aβ42, MMSE, ADAS-Cog13, mPACC and adjusted hippocampal volume
change, with the signs of hippocampal volume, mPACC and MMSE reversed as in T5 so that positive
always means the plasma slope tracks worsening. Every draw is appended to `bootstrap_draws.rds` and the run
resumes from it; T7 is the 2.5th and 97.5th percentiles over all saved draws (`n_draws`
reported). `N_BOOT` default 1,000.

**Section 2 (T3b, T15b; Table 3 and Table S10 intervals).** The same resamples and seeds.
Within each draw: the reduced change model per biomarker (var_true, residual), the per-platform
T3 model (sem, mean rate), per-participant least-squares slopes, `sd_log` and the median span,
all recomputed on the resample; then the two-point and slope MDC95 as per cent, years to exceed at
the amyloid-positive rate, within-person CV, amyloid-positive rate, share with a reliable increase,
and the T15 slope-endpoint trial sizes for the 25% effect only. Cached in
`mdc_bootstrap_draws.rds`; `P2_MDC_BOOT` default 1,000.

**Section 3 (R10; Tables S9b, S10, S16c; Figures 2B and 4B).** The same resamples and seeds.
Within each draw: B1 the pTau217 sample model (lower bound of averaging: one platform, mean of
platforms with the shared part retained, and the spans needed for 0.7, 0.8, 0.9); B2 each of
the 15 measures fitted alone at the ADNI and A4 SSTs; B3 trial sizes on amyloid-positive
components for every platform and for the mean of platforms under the upper bound, and for
pTau217 also under the lower bound, three schedules, two effects; B4 the immunoassay ICC. Draw 0 is the full data and must reproduce R1, R3a, R4 and
R6. Singular draws are kept and their proportion reported. Cached in `R10_bootstrap_draws.rds`.

**Section 4 (R11; Table S9b).** Intervals for the span needed, derived from the saved T7 draws
without refitting: `ratio = (1 / rel_one − 1) × median_sst` (one platform) or
`(1 / rel_mean − 1) × median_sst` (mean), then the span formula of stage 6; the point values
must equal T16b, which is checked.

**Section 5 (R12, R12b; Tables S10b, S16c; every head-to-head statement in the Results).**
Paired differences and ratios computed within each draw, so that comparisons have valid
intervals: D1 trial sizes (mean of platforms under each bound minus Fujirebio; Roche minus
Quanterix pTau181; ratios C2N %pTau217 / Fujirebio and Roche / Quanterix; GFAP and NfL pairs);
D2 reliabilities (each pTau217 measure fitted alone minus the joint value; the lower bounds
minus the joint values; GFAP and NfL); D3 the immunoassay ICC minus the five-measure ICC;
D4 the variance components of the singular amyloid-positive fits (GFAP, NfL, Aβ42/40), to
state which component sat at the boundary.

**Section 6 (printed only).** The amyloid-positive minus amyloid-negative difference in pTau217
validity against Centiloid change within each saved draw, 2.5th and 97.5th percentiles
(Results; the construction of Table S25).

**Section 7 (T7c, T17p; A4 T8p, T8bp; Tables S13a, S13b, S14).** T7c: observed-data (not
bootstrap) reliabilities and Spearman and Pearson validity of the pooled slope by amyloid group,
the reduced change model refitted within each group as in the bootstrap, same six criteria and
sign orientation; T17p: the ADNI validity
by SST quintile with Pearson beside Spearman; T8p: the A4 within-participant design subsets with
Pearson; T8bp: the A4 validity by SST band with Pearson.

**Section 8 (R13a, R13b, R13c; Table 2 and Table S25 intervals, the paired comparisons in the
Results).** From the saved draws, without refitting: R13a from the main bootstrap, per biomarker,
the gain from averaging (mean of platforms minus one platform), each biomarker's reliability minus
pTau217's and Aβ42/40 minus NfL, and amyloid-positive minus amyloid-negative differences in the
two reliabilities and in validity against Centiloid change; R13b from the MDC bootstrap, platform
differences in within-person CV, two-point and slope MDC and years to exceed, for every pTau217
pair and for Roche against Quanterix; R13c the share of singular sample-model fits in the R10
bootstrap. Percentile intervals over the 1,000 draws.

## Stage 10. Tables (`10_tables.R`)

Section 1 (ADNI): Table 1 participant characteristics overall and by amyloid status at the
first plasma sample, including, for baseline-negative participants, conversion within the sampling
period (first positive scan dated on or before the last five-platform plasma sample) and after it; Table 2 level components, G coefficient, slope components, slope
reliability with intervals (all and by amyloid status), slope ICC, per biomarker; Table 3
minimal detectable change per platform on the delivered scale; Table S1 assays (units,
replicates, observations, % at the lower limit, batches, external QC CVs); Table S2 the
analytical share of the occasion-level variance (QC CV² / within-person CV²). Section 2:
cross-cohort tables C1 design, C2 slope reliability with intervals, C3 the A4 12-week
test-retest beside each cohort's change-model residual, C4 validity, C5 minimal detectable
change, C6 span needed, C7 the ADNI SST quintiles beside the A4 within-person design subsets,
and figure F30. No model fits.

## Stage 11. Figures (`11_figures.R`)

From the finished tables and the two analysis frames (`frame_paper2_long`, `frame_paper2_slopes`); no model fits. Working figures F7 to F12; the journal versions with legends;
the manuscript composites: Figure 1 (design schematic, pTau217 trajectories), Figure 2 (level
components; slope reliability with upper-bound symbols, pTau217 lower-bound open symbols from
R10, fitted-alone grey bands from R10), Figure 3 (slopes against the MDC band; years to exceed
it), Figure 4 (reliability against SST with schedules; amyloid-positive trial sizes with
intervals from R10, both bounds for the averaged pTau217 measure; the A4 four-sample test from
`combined/T9`); Supplementary Figures S1 to S5; `manuscript_figure_legends.md`.

---

## Order of execution and dependencies

`run_all.R` runs: 01 (adni), 02 (a4learn), 03 (adni), 03 (a4learn), 04, 05, 06, 07, 08, 09,
10, 11. The dependencies that fix this order: 04 needs the frames (03); 05 and 06 need T2 and
T3 (04); 06 section 2 needs T2, T3 (04) and T15 (06 section 1) and refits the sample model from the
frames; 07 needs the A4 frames (02, 03) and, for sections 2 and 3, `frame_paper3_long` and
`frame_paper3_person` (03); 08 needs T3, T4b (04) and T14 (05); 09 needs the frames of both
cohorts, T2 and T16b (04, 06), T14 (05, for a check), R1, R3a, R4, R6 (04, 06), its own T2d and
the A4 T2 (07); 10 needs T7 (09) and both cohorts' tables; 11 needs R10 (09), `combined/T9` (07)
and the frames. Every section
re-sources `00_setup.R` and, where it matters, sets `COHORT`, so `OUT_DIR` always points at
the right cohort's folder.

## Reproducibility notes

- Seeds: ADNI participant bootstraps `20260904 + draw` in stage 9 sections 1 to 3, so the
  three bootstraps use the same resamples and paired comparisons within a draw are valid;
  A4/LEARN bootstrap 20260906 (stage 7 section 1); R8 bootstrap 20261106 (stage 7 section 3); the cross-platform multiple imputation 20260904
  (stage 8 section 1); the power simulation 20260908 (stage 8 section 5); the fourth-sample contrast bootstrap `20260908 + draw` (stage 7
  section 2, `P2_SEED4`). Stage 4 section 2 sets 20260907 (`P2_SEED`) but draws no random numbers as the code stands.
- The cached draws (`bootstrap_draws.rds`, `mdc_bootstrap_draws.rds`, `R10_bootstrap_draws.rds`
  in `pipeline_out/paper2/`) are resumed, not recomputed, if present; delete them to force a
  full recomputation.
- Rounding: tables are written at full precision except where a script rounds for display
  (T3 percentages to one decimal, T4 and T5 correlations to three decimals in both cohorts, T19
  to three, T20 to four).
