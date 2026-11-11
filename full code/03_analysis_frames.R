# =============================================================================================
# 03_analysis_frames.R
# Master table and analysis frames, for whichever cohort COHORT names.
#
# Run twice by run_all.R, with COHORT = adni and COHORT = a4learn. Sections 1 to 3 (formerly 06, 07b, 07c)
# join the data and derive the frames; section 4 (formerly 09e) runs only for A4/LEARN and compares its
# frames with the ADNI frames.
# =============================================================================================


# =============================================================================================
# Section 1. Master table: every measurement in one tall table, one row per plasma visit
# Original script: 06_master.R
# =============================================================================================
# =============================================================================
# 06_master.R
# Joins everything into two master tables.
#
#   master_long.rds   Every measurement from every source, one row each, in
#                     one tall table with the same columns: domain, source,
#                     measure, exam_date, value, flags, nearest plasma visit.
#                     Use this when a model needs its own matching rule.
#
#   master_visit.rds  One row per person per plasma visit. Participant
#                     characteristics, the 15 plasma measures as columns, and
#                     for every criterion and modifier the nearest assessment
#                     within the window set in 00_setup.R (MATCH_WINDOW_DAYS),
#                     with its gap in days. Baseline anchors (amyloid status
#                     and CSF at the first plasma visit) are added per person.
#                     This is the frame most models start from.
#
# Both are saved as .rds and .csv. checks/master_visit_missingness.csv shows
# how complete each column is at each plasma visit.
# =============================================================================

# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
source("00_setup.R")

pl   <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
crit <- readRDS(file.path(OUT_DIR, "criteria_long.rds"))
mods <- readRDS(file.path(OUT_DIR, "modifiers_long.rds"))
dx   <- if (file.exists(file.path(OUT_DIR, "dx_long.rds"))) readRDS(file.path(OUT_DIR, "dx_long.rds")) else NULL

# ---- 1. master_long --------------------------------------------------------------
plasma_rows <- pl %>%
  filter(role == "measure") %>%
  transmute(RID, domain = "plasma", source = platform, measure = PLASMA_BIOMARKER,
            biomarker, platform,
            exam_date = visit_date, VISCODE2, value, value_raw = TESTVALUE_raw,
            comparator, lloq_flag, uloq_flag, imputed_flag, derived,
            visit_no, occasion_no, visit_date, gap_days = 0L, plasma_time,
            RUN_NO, RUNDATE, quanterix_only_visit)
crit_rows <- crit %>%
  transmute(RID, domain = "criterion", source, measure, biomarker = NA_character_, platform = NA_character_,
            exam_date, VISCODE2, value, value_raw, comparator, lloq_flag, uloq_flag,
            imputed_flag = FALSE, derived, visit_no, occasion_no = NA_integer_, visit_date, gap_days, plasma_time,
            RUN_NO = NA_character_, RUNDATE = as.Date(NA), quanterix_only_visit = NA, tracer)
mod_rows <- mods %>%
  transmute(RID, domain = "modifier", source, measure, biomarker = NA_character_, platform = NA_character_,
            exam_date, VISCODE2, value, value_raw = NA_character_, comparator = NA_character_,
            lloq_flag = FALSE, uloq_flag = FALSE, imputed_flag = FALSE, derived,
            visit_no, occasion_no = NA_integer_, visit_date, gap_days, plasma_time,
            RUN_NO = NA_character_, RUNDATE = as.Date(NA), quanterix_only_visit = NA)
dx_rows <- if (!is.null(dx)) dx %>%
  transmute(RID, domain = "diagnosis", source = "adnimerge", measure = "dx",
            biomarker = NA_character_, platform = NA_character_,
            exam_date, VISCODE2, value = NA_real_, value_raw = dx, comparator = NA_character_,
            lloq_flag = FALSE, uloq_flag = FALSE, imputed_flag = FALSE, derived = FALSE,
            visit_no, occasion_no = NA_integer_, visit_date, gap_days, plasma_time,
            RUN_NO = NA_character_, RUNDATE = as.Date(NA), quanterix_only_visit = NA) else NULL

master_long <- bind_rows(plasma_rows, crit_rows, mod_rows, dx_rows) %>%
  left_join(part %>% select(RID, in_core_5lab, in_quanterix, quanterix_only_person), by = "RID") %>%
  arrange(RID, domain, source, measure, exam_date)
save_rds(master_long, "master_long")
write_csv(master_long, file.path(OUT_DIR, "master_long.csv"))
say("master_long: %d rows, %d people; rows by domain:", nrow(master_long), n_distinct(master_long$RID))
print(master_long %>% count(domain))

# ---- 2. master_visit: one row per person per plasma visit --------------------------
visits <- pl %>%
  distinct(RID, visit_no, visit_date, plasma_time, VISCODE2,
           n_platforms_at_visit, quanterix_only_visit,
           # A4/LEARN only (absent in ADNI): design columns carried through to the frames
           across(any_of(c("substudy", "arm", "visit_kind", "visit_week", "on_treatment", "short_interval_pair"))))

plasma_wide <- pl %>%
  filter(role == "measure") %>%
  select(RID, visit_no, PLASMA_BIOMARKER, value) %>%
  pivot_wider(names_from = PLASMA_BIOMARKER, values_from = value,
              values_fn = function(x) if (length(x) > 1) NA_real_ else x)
# Floor flags for the measures where it matters
floor_wide <- pl %>%
  filter(role == "measure", PLASMA_BIOMARKER %in% c("C2N_plasma_ptau217", "Fuji_plasma_ptau217")) %>%
  transmute(RID, visit_no, name = paste0(PLASMA_BIOMARKER, "_floor"), lloq_flag) %>%
  pivot_wider(names_from = name, values_from = lloq_flag, values_fn = any)
occ_wide <- pl %>%
  filter(role == "measure", !is.na(occasion_no)) %>%
  distinct(RID, visit_no, platform, occasion_no) %>%
  mutate(name = paste0("occasion_", tolower(platform))) %>%
  select(-platform) %>%
  pivot_wider(names_from = name, values_from = occasion_no)

# Nearest assessment of each measure within the window, as value + gap columns
nearest_wide <- function(df, window) {
  near <- df %>%
    filter(!is.na(visit_no), abs(gap_days) <= window) %>%
    group_by(RID, visit_no, measure) %>%
    slice_min(abs(gap_days), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(RID, visit_no, measure, value, gap_days)
  vals <- near %>% select(-gap_days) %>% pivot_wider(names_from = measure, values_from = value)
  gaps <- near %>% select(-value) %>%
    mutate(measure = paste0(measure, "_gap")) %>%
    pivot_wider(names_from = measure, values_from = gap_days)
  left_join(vals, gaps, by = c("RID", "visit_no"))
}
crit_wide <- nearest_wide(crit, MATCH_WINDOW_DAYS[["criteria"]])
mod_wide  <- nearest_wide(mods, MATCH_WINDOW_DAYS[["modifiers"]])
dx_wide <- if (!is.null(dx)) dx %>%
  filter(!is.na(visit_no), abs(gap_days) <= MATCH_WINDOW_DAYS[["dx"]]) %>%
  group_by(RID, visit_no) %>% slice_min(abs(gap_days), n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, visit_no, dx, dx_gap = gap_days) else NULL

# Per-person anchors at the first plasma visit
anchor <- function(df, m, new, window) {
  df %>% filter(measure == m, visit_no == 1, abs(gap_days) <= window) %>%
    group_by(RID) %>% slice_min(abs(gap_days), n = 1, with_ties = FALSE) %>% ungroup() %>%
    transmute(RID, !!new := value, !!paste0(new, "_gap") := gap_days)
}
anchors <- part %>% select(RID) %>%
  left_join(anchor(crit, "centiloid", "centiloid_bl", MATCH_WINDOW_DAYS[["criteria"]]), by = "RID") %>%
  left_join(anchor(crit, "csf_ptau_abeta42", "csf_ptau_abeta42_bl", MATCH_WINDOW_DAYS[["criteria"]]), by = "RID") %>%
  left_join(anchor(mods, "egfr_2021", "egfr_2021_bl", MATCH_WINDOW_DAYS[["modifiers"]]), by = "RID") %>%
  mutate(amyloid_bl_cl20 = centiloid_bl > CENTILOID_CUTS[["cl20"]],
         amyloid_bl_cl37 = centiloid_bl > CENTILOID_CUTS[["cl37"]])
if ("bmi" %in% mods$measure)
  anchors <- anchors %>% left_join(anchor(mods, "bmi", "bmi_bl", MATCH_WINDOW_DAYS[["modifiers"]]), by = "RID")

# Ever amyloid positive on any scan, and the first positive scan date
amy_ever <- crit %>%
  filter(measure == "centiloid") %>%
  group_by(RID) %>%
  summarise(centiloid_max = max(value),
            amyloid_ever_cl20 = any(value > CENTILOID_CUTS[["cl20"]]),
            first_positive_cl20_date = suppressWarnings(min(exam_date[value > CENTILOID_CUTS[["cl20"]]])),
            n_pet_scans = n(), .groups = "drop") %>%
  mutate(first_positive_cl20_date = as.Date(ifelse(is.finite(first_positive_cl20_date),
                                                   first_positive_cl20_date, NA), origin = "1970-01-01"))

part_cols <- part %>%
  select(RID, PTID, in_core_5lab, in_quanterix, quanterix_only_person,
         n_plasma_visits, first_plasma, last_plasma, plasma_span_years,
         age_first_plasma, PTGENDER, PTEDUCAT, APOE4, APOE4_carrier, DX_bl, DX_at_first_plasma,
         COLPROT, ORIGPROT, any_of(c("treatment_start", "treatment_end", "completed")))

master_visit <- visits %>%
  left_join(part_cols, by = "RID") %>%
  mutate(age_at_visit = age_first_plasma + plasma_time) %>%
  left_join(anchors, by = "RID") %>%
  left_join(amy_ever, by = "RID") %>%
  left_join(occ_wide, by = c("RID", "visit_no")) %>%
  left_join(plasma_wide, by = c("RID", "visit_no")) %>%
  left_join(floor_wide, by = c("RID", "visit_no")) %>%
  left_join(crit_wide, by = c("RID", "visit_no")) %>%
  left_join(mod_wide, by = c("RID", "visit_no")) %>%
  { if (!is.null(dx_wide)) left_join(., dx_wide, by = c("RID", "visit_no")) else . } %>%
  arrange(RID, visit_no)

# Columns the 07 frame scripts refer to by name. A cohort without a measure (for
# example A4/LEARN has no ADAS13 or MRI) gets the column as NA so the scripts run unchanged.
expected_cols <- c("adas13", "cdrsb", "hcv_adj", "hippocampus", "icv", "mmse", "mpacc_digit", "mpacc_trailsb",
                   "amyloid_status", "centiloid", "summary_suvr", "csf_abeta40", "csf_abeta42", "csf_abeta42_40",
                   "csf_ptau181", "csf_ptau_abeta42", "csf_ptau_tau", "csf_tau", "bmi", "creatinine", "egfr_2009",
                   "egfr_2021", "height_cm", "weight_kg", "dx", "dx_gap", "bmi_bl", "bmi_bl_gap",
                   "csf_ptau_abeta42_bl", "csf_ptau_abeta42_bl_gap", "egfr_2021_bl", "egfr_2021_bl_gap",
                   "occasion_c2n", "occasion_fujirebio", "occasion_quanterix", "occasion_roche", "occasion_alzpath", "occasion_janssen")
for (cc in setdiff(expected_cols, names(master_visit))) master_visit[[cc]] <- if (cc %in% c("dx")) NA_character_ else NA_real_
for (cc in setdiff(paste0(expected_cols[!str_detect(expected_cols, "_gap$|_bl$|^occasion_")], "_gap"), names(master_visit)))
  master_visit[[cc]] <- NA_integer_
save_rds(master_visit, "master_visit")
write_csv(master_visit, file.path(OUT_DIR, "master_visit.csv"))
say("master_visit: %d rows (person-visits), %d people, %d columns",
    nrow(master_visit), n_distinct(master_visit$RID), ncol(master_visit))

# ---- 3. Completeness by visit, core sample only ------------------------------------
core <- master_visit %>% filter(in_core_5lab, !quanterix_only_visit)
n_by_visit <- core %>% count(visit_no)
miss <- core %>%
  group_by(visit_no) %>%
  summarise(across(-RID, ~ round(100 * mean(!is.na(.x)), 1)), .groups = "drop") %>%
  pivot_longer(-visit_no, names_to = "column", values_to = "pct_present") %>%
  pivot_wider(names_from = visit_no, values_from = pct_present, names_prefix = "visit_")
say("Core plasma visits by visit number: %s",
    paste(sprintf("%d: %d", n_by_visit$visit_no, n_by_visit$n), collapse = ", "))
write_csv(miss, file.path(CHECK_DIR, "master_visit_completeness.csv"))
key_cols <- c("centiloid", "csf_ptau_abeta42", "mmse", "adas13", "mpacc_digit", "hcv_adj",
              "egfr_2021", "bmi", "dx")
say("Per cent of core plasma visits with each criterion/modifier within window:")
print(as.data.frame(miss %>% filter(column %in% key_cols)))

# Anchors summary
say("Baseline anchors (core sample): centiloid at visit 1 for %d, CSF pTau/Ab42 for %d, eGFR for %d",
    sum(!is.na(anchors$centiloid_bl[anchors$RID %in% core$RID])),
    sum(!is.na(anchors$csf_ptau_abeta42_bl[anchors$RID %in% core$RID])),
    sum(!is.na(anchors$egfr_2021_bl[anchors$RID %in% core$RID])))
print(core %>% distinct(RID, amyloid_bl_cl20, amyloid_ever_cl20) %>% count(amyloid_bl_cl20, amyloid_ever_cl20))


# =============================================================================================
# Section 2. Analysis frame: persons x platforms x occasions per biomarker, per-person slopes, criterion slopes
# Original script: 07b_frame_paper2_reliability.R
# =============================================================================================
# =============================================================================
# 07b_frame_paper2_reliability.R
# Analysis frame for Paper 2: reliability and validity of individual change.
#
# The core object is persons x platforms x occasions for each biomarker, on a
# common (log, standardised) scale, with the analytical run kept as a batch
# facet. From it the script also builds per-person, per-platform slopes and
# per-person criterion slopes (amyloid, CSF, cognition, MRI) for validity.
#
# Sample   core 5-lab people, Quanterix-only samples excluded. Quanterix rows
#          are included where the same sample was measured (in_quanterix).
#          occasion_no is the within-platform occasion (1, 2, 3); visit_no is
#          the shared sample index and is what aligns platforms.
#
# Output   frame_paper2_long.rds     person x platform x occasion x biomarker
#          frame_paper2_slopes.rds   person x platform x biomarker slopes
#          frame_paper2_slopes_wide.rds  slopes side by side across platforms
#          frame_paper2_criteria.rds person-level criterion change
#          checks/paper2_design.csv  how balanced the design is
# =============================================================================

source("00_setup.R")

ml   <- readRDS(file.path(OUT_DIR, "master_long.rds"))
mv   <- readRDS(file.path(OUT_DIR, "master_visit.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
crit <- readRDS(file.path(OUT_DIR, "criteria_long.rds"))

measures <- MEASURE_LOOKUP %>% filter(role == "measure")

# ---- 1. Persons x platforms x occasions ----------------------------------------
long <- ml %>%
  filter(domain == "plasma", in_core_5lab, !quanterix_only_visit, !is.na(value)) %>%
  select(RID, biomarker, platform, measure, visit_no, occasion_no, plasma_time, visit_date,
         value, lloq_flag, imputed_flag, RUN_NO, RUNDATE) %>%
  left_join(measures %>% select(measure = PLASMA_BIOMARKER, direction), by = "measure") %>%
  mutate(log_value = safe_log(value)) %>%
  group_by(measure) %>%
  mutate(z_log = as.numeric(scale(log_value)),
         z_log_aligned = z_log * direction,
         # also a within-measure baseline-standardised scale: centred and
         # scaled on occasion 1 only, so change is in baseline-SD units
         z_bl = (log_value - mean(log_value[occasion_no == 1], na.rm = TRUE)) /
                sd(log_value[occasion_no == 1], na.rm = TRUE),
         z_bl_aligned = z_bl * direction) %>%
  ungroup() %>%
  # time relative to each person's first plasma sample is plasma_time already;
  # add time centred at the person's mean so slopes and intercepts are less correlated
  group_by(RID) %>%
  mutate(time_c = plasma_time - mean(unique(plasma_time))) %>%
  ungroup() %>%
  left_join(part %>% select(RID, age_first_plasma, PTGENDER, PTEDUCAT, APOE4_carrier,
                            DX_at_first_plasma, in_quanterix), by = "RID") %>%
  left_join(mv %>% distinct(RID, amyloid_bl_cl20, amyloid_bl_cl37, amyloid_ever_cl20, centiloid_bl),
            by = "RID") %>%
  # batch is nested in platform; make an unambiguous label
  mutate(batch = paste(platform, RUN_NO, sep = ":")) %>%
  arrange(RID, biomarker, platform, occasion_no)

# Design balance per biomarker: people with every platform at every occasion
platforms_per_bm <- long %>% distinct(biomarker, platform) %>% count(biomarker, name = "n_platforms")
complete <- long %>%
  distinct(RID, biomarker, platform, visit_no) %>%
  count(RID, biomarker, name = "n_cells") %>%
  left_join(platforms_per_bm, by = "biomarker") %>%
  mutate(complete_3x = n_cells == 3 * n_platforms)
long <- long %>% left_join(complete %>% select(RID, biomarker, complete_3x), by = c("RID", "biomarker"))

design <- long %>%
  group_by(biomarker) %>%
  summarise(n_platforms = n_distinct(platform), n_people = n_distinct(RID),
            n_people_complete = n_distinct(RID[complete_3x]),
            n_cells = n(), n_batches = n_distinct(batch),
            pct_floor = round(100 * mean(lloq_flag), 1), .groups = "drop")
print(as.data.frame(design))
write_csv(design, file.path(CHECK_DIR, "paper2_design.csv"))

# ---- 2. Per-person, per-platform slopes ------------------------------------------
slopes <- slope_by(long, "z_log_aligned", "plasma_time", RID, biomarker, platform, measure) %>%
  left_join(long %>% group_by(RID, measure) %>%
              summarise(n_floor = sum(lloq_flag), any_floor = any(lloq_flag), .groups = "drop"),
            by = c("RID", "measure")) %>%
  left_join(part %>% select(RID, age_first_plasma, PTGENDER, APOE4_carrier, DX_at_first_plasma), by = "RID") %>%
  left_join(mv %>% distinct(RID, amyloid_bl_cl20, amyloid_ever_cl20, centiloid_bl), by = "RID")

# Slopes side by side across measures, one row per person x biomarker.
# Keyed on measure rather than platform because C2N contributes two pTau217
# measures (absolute and %pTau217). Short labels: platform, plus "_ratio"
# for the C2N ratio.
short <- measures %>%
  transmute(measure = PLASMA_BIOMARKER,
            label = case_when(PLASMA_BIOMARKER == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio",
                              TRUE ~ platform_from_prefix(PLASMA_BIOMARKER)))
slopes_wide <- slopes %>%
  left_join(short, by = "measure") %>%
  select(RID, biomarker, label, slope) %>%
  pivot_wider(names_from = label, values_from = slope, names_prefix = "slope_") %>%
  left_join(slopes %>% left_join(short, by = "measure") %>% select(RID, biomarker, label, y_first) %>%
              pivot_wider(names_from = label, values_from = y_first, names_prefix = "bl_"),
            by = c("RID", "biomarker"))

# ---- 3. Criterion change per person (validity targets) ---------------------------
# Use every assessment from CRITERION_SLOPE_PAD_DAYS before the first plasma
# sample to the same after the last, so the criterion slope covers the same
# years as the plasma slope.
span <- part %>% select(RID, first_plasma, last_plasma)
crit_in <- crit %>%
  filter(measure %in% c("centiloid", "csf_ptau_abeta42", "csf_abeta42", "csf_ptau181",
                        "mmse", "adas13", "mpacc_digit", "mpacc_trailsb", "cdrsb", "hcv_adj",
                        "pacc", "tau_pet")) %>%
  inner_join(span, by = "RID") %>%
  filter(exam_date >= first_plasma - CRITERION_SLOPE_PAD_DAYS,
         exam_date <= last_plasma  + CRITERION_SLOPE_PAD_DAYS) %>%
  mutate(t = as.numeric(exam_date - first_plasma) / YEAR_DAYS)
crit_slopes <- slope_by(crit_in, "value", "t", RID, measure) %>%
  rename(crit_n = n_obs, crit_span = span, crit_first = y_first, crit_last = y_last,
         crit_slope = slope, crit_slope_se = slope_se) %>%
  select(-intercept, -t_first, -t_last)
crit_wide <- crit_slopes %>%
  select(RID, measure, crit_slope, crit_n, crit_span) %>%
  pivot_wider(names_from = measure, values_from = c(crit_slope, crit_n, crit_span),
              names_glue = "{measure}_{.value}")
names(crit_wide) <- str_replace(names(crit_wide), "_crit_", "_")

criteria <- part %>% filter(in_core_5lab) %>% select(RID) %>%
  left_join(crit_wide, by = "RID") %>%
  left_join(mv %>% distinct(RID, centiloid_bl, amyloid_bl_cl20, amyloid_bl_cl37, amyloid_ever_cl20,
                            centiloid_max, first_positive_cl20_date, n_pet_scans), by = "RID")

# ---- 4. Save ------------------------------------------------------------------------
save_rds(long,        "frame_paper2_long")
save_rds(slopes,      "frame_paper2_slopes")
save_rds(slopes_wide, "frame_paper2_slopes_wide")
save_rds(criteria,    "frame_paper2_criteria")
write_csv(slopes_wide, file.path(OUT_DIR, "frame_paper2_slopes_wide.csv"))

say("frame_paper2_long: %d cells, %d people, %d biomarkers, %d measures",
    nrow(long), n_distinct(long$RID), n_distinct(long$biomarker), n_distinct(long$measure))
say("Slopes: %d person x measure; with 3 points: %d", nrow(slopes), sum(slopes$n_obs == 3))
say("Criterion slopes with 2+ points: centiloid %d, csf pTau/Ab42 %d, mmse %d, hcv_adj %d",
    sum(criteria$centiloid_n >= 2, na.rm = TRUE), sum(criteria$csf_ptau_abeta42_n >= 2, na.rm = TRUE),
    sum(criteria$mmse_n >= 2, na.rm = TRUE), if ("hcv_adj_n" %in% names(criteria)) sum(criteria$hcv_adj_n >= 2, na.rm = TRUE) else 0L)
# Quick look: cross-platform slope correlations per biomarker
cors <- slopes_wide %>%
  group_by(biomarker) %>%
  group_modify(~ {
    m <- .x %>% select(starts_with("slope_")) %>% select(where(~ !all(is.na(.x))))
    if (ncol(m) < 2) return(tibble(pair = NA_character_, r = NA_real_, n = NA_integer_))
    cm <- cor(m, use = "pairwise.complete.obs")
    p <- which(upper.tri(cm), arr.ind = TRUE)
    tibble(pair = paste(str_remove(colnames(cm)[p[, 1]], "slope_"),
                        str_remove(colnames(cm)[p[, 2]], "slope_"), sep = " vs "),
           r = round(cm[p], 3),
           n = sapply(seq_len(nrow(p)), function(i) sum(complete.cases(m[, c(p[i, 1], p[i, 2])]))))
  }) %>% ungroup()
print(as.data.frame(cors))
write_csv(cors, file.path(CHECK_DIR, "paper2_slope_correlations.csv"))


# =============================================================================================
# Section 3. Person x visit x measure frame with every exposure scale (used by the fourth-sample test)
# Original script: 07c_frame_paper3_multiverse.R
# =============================================================================================
# =============================================================================
# 07c_frame_paper3_multiverse.R
# Analysis frame for Paper 3: how modelling choices change what a plasma
# trajectory is said to predict (multiverse / specification curve), using
# cross-platform consistency as the criterion for an estimator.
#
# The frame is deliberately generous: every scale of the exposure, every
# concurrent outcome at every visit, and every person-level outcome, so that
# each specification is a filter and a formula rather than new data work.
#
# Sample   core 5-lab people, Quanterix-only samples excluded.
#
# Output   frame_paper3_long.rds    person x visit x measure, exposure scales +
#                                   concurrent outcomes at that visit
#          frame_paper3_person.rds  person-level outcomes and subgroup flags
#          frame_paper3_specs.csv   the specification grid (choices and defaults)
#          checks/paper3_outcomes.csv
# =============================================================================

source("00_setup.R")

ml   <- readRDS(file.path(OUT_DIR, "master_long.rds"))
mv   <- readRDS(file.path(OUT_DIR, "master_visit.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
crit <- readRDS(file.path(OUT_DIR, "criteria_long.rds"))
p2c  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))   # criterion slopes from 07b
measures <- MEASURE_LOOKUP %>% filter(role == "measure")

# ---- 1. Person-level outcomes and subgroup flags ------------------------------------
last_visit <- mv %>%
  filter(in_core_5lab, !quanterix_only_visit) %>%
  group_by(RID) %>% slice_max(visit_no, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(RID, last_plasma_time = plasma_time,
            centiloid_last = centiloid, amyloid_last_cl20 = centiloid > CENTILOID_CUTS[["cl20"]],
            amyloid_last_cl37 = centiloid > CENTILOID_CUTS[["cl37"]],
            mmse_last = mmse, adas13_last = adas13, hcv_adj_last = hcv_adj,
            csf_ptau_abeta42_last = csf_ptau_abeta42, dx_last = dx)

dx_rank <- function(x) case_when(x == "CN" ~ 1L, x == "MCI" ~ 2L, x == "Dementia" ~ 3L, TRUE ~ NA_integer_)
# Diagnosis progression: any move from CN to MCI/Dementia or MCI to Dementia
# between the first and last plasma visit (using dx at those visits)
dx_prog <- mv %>%
  filter(in_core_5lab, !quanterix_only_visit, !is.na(dx)) %>%
  group_by(RID) %>%
  summarise(dx_first = dx[which.min(visit_no)], dx_last = dx[which.max(visit_no)],
            dx_ever_dementia = any(dx == "Dementia"), .groups = "drop") %>%
  mutate(dx_progressed = dx_rank(dx_last) > dx_rank(dx_first))

person <- part %>% filter(in_core_5lab) %>%
  select(RID, PTID, age_first_plasma, PTGENDER, PTEDUCAT, APOE4, APOE4_carrier,
         DX_bl, DX_at_first_plasma, n_plasma_visits, plasma_span_years, in_quanterix) %>%
  left_join(mv %>% distinct(RID, centiloid_bl, centiloid_bl_gap, amyloid_bl_cl20, amyloid_bl_cl37,
                            amyloid_ever_cl20, centiloid_max, first_positive_cl20_date, n_pet_scans,
                            csf_ptau_abeta42_bl, bmi_bl, egfr_2021_bl), by = "RID") %>%
  left_join(last_visit, by = "RID") %>%
  left_join(dx_prog %>% select(RID, dx_first, dx_ever_dementia, dx_progressed), by = "RID") %>%
  left_join(p2c %>% select(RID, ends_with("_slope"), ends_with("_n")), by = "RID") %>%
  mutate(
    # Concurrent-outcome design flags. "Prevalent" = positive already at the
    # first plasma visit; these people contaminate any incident analysis.
    prevalent_amyloid_cl20 = amyloid_bl_cl20,
    incident_amyloid_cl20  = !amyloid_bl_cl20 & amyloid_last_cl20,
    baseline_negative_cl20 = !amyloid_bl_cl20,
    baseline_negative_cl37 = !amyloid_bl_cl37,
    baseline_cn            = DX_at_first_plasma == "CN",
    # First positive scan relative to the plasma series, in years
    first_positive_time = as.numeric(first_positive_cl20_date - part$first_plasma[match(RID, part$RID)]) / YEAR_DAYS,
    first_positive_after_last_plasma = !is.na(first_positive_time) & first_positive_time > last_plasma_time
  )

# ---- 2. Long frame: exposure scales + concurrent outcomes at each visit ------------
long <- ml %>%
  filter(domain == "plasma", in_core_5lab, !quanterix_only_visit, !is.na(value)) %>%
  select(RID, biomarker, platform, measure, visit_no, occasion_no, plasma_time, value,
         lloq_flag, imputed_flag, RUN_NO) %>%
  left_join(measures %>% select(measure = PLASMA_BIOMARKER, direction), by = "measure") %>%
  group_by(measure) %>%
  mutate(
    log_value   = safe_log(value),
    log1p_value = log1p(value),                         # what the double-cohort paper used
    z_raw       = as.numeric(scale(value)),
    z_log       = as.numeric(scale(log_value)),
    z_log_aligned = z_log * direction,
    # floor handling options as columns: as delivered, set to NA, half the floor
    value_floor_na   = if_else(lloq_flag, NA_real_, value),
    value_floor_half = if_else(lloq_flag, value / 2, value),
    log_floor_na     = safe_log(value_floor_na),
    log_floor_half   = safe_log(value_floor_half)
  ) %>%
  ungroup() %>%
  group_by(RID) %>%
  mutate(time_c = plasma_time - mean(unique(plasma_time))) %>%
  ungroup() %>%
  # concurrent outcomes at the same plasma visit
  left_join(mv %>% select(RID, visit_no, age_at_visit,
                          any_of(c("centiloid", "centiloid_gap", "amyloid_status", "summary_suvr",
                                   "mmse", "adas13", "mpacc_digit", "mpacc_trailsb", "cdrsb",
                                   "hippocampus", "icv", "hcv_adj", "csf_ptau_abeta42", "csf_abeta42",
                                   "csf_ptau181", "bmi", "dx", "dx_gap", "pacc", "tau_pet",
                                   "substudy", "arm", "on_treatment", "short_interval_pair"))),
            by = c("RID", "visit_no")) %>%
  mutate(amyloid_cl20 = centiloid > CENTILOID_CUTS[["cl20"]],
         amyloid_cl37 = centiloid > CENTILOID_CUTS[["cl37"]]) %>%
  left_join(person %>% select(RID, PTGENDER, PTEDUCAT, APOE4_carrier, age_first_plasma,
                              DX_at_first_plasma, amyloid_bl_cl20, baseline_negative_cl20,
                              prevalent_amyloid_cl20, incident_amyloid_cl20, baseline_cn),
            by = "RID") %>%
  arrange(RID, biomarker, platform, visit_no)

# ---- 3. Specification grid -----------------------------------------------------------
specs <- tribble(
  ~dimension,          ~options,                                                        ~default,
  "exposure_measure",  "each of the 15 measures; primary = the 5 pTau217 measures",     "C2N_plasma_ptau217",
  "exposure_scale",    "value | log_value | log1p_value | z_log | z_log_aligned",        "z_log",
  "floor_handling",    "as_delivered | floor_na | floor_half | exclude_always_floor | censored", "as_delivered",
  "trajectory_model",  "lme_slope | hlme_2class | hlme_3class | hlme_4class | gamm",     "lme_slope",
  "class_summary",     "modal_class | posterior_probs | BCH | random_slope_draws",       "modal_class",
  "linkage",           "two_stage | joint (multlcmm/Jointlcmm)",                          "two_stage",
  "outcome_type",      "concurrent_level | concurrent_rate | status_at_last | incident (baseline-negative only)", "concurrent_rate",
  "outcome",           "centiloid | amyloid_cl20 | amyloid_cl37 | csf_ptau_abeta42 | mmse | adas13 | mpacc_digit | hcv_adj | dx_progressed", "centiloid",
  "sample",            "all_core | baseline_negative_cl20 | baseline_negative_cl37 | baseline_cn | complete_3_visits", "all_core",
  "covariates",        "none | age_sex | age_sex_apoe | age_sex_apoe_educ",               "age_sex_apoe",
  "time_scale",        "plasma_time | age_at_visit",                                     "plasma_time",
  "platform_criterion","agreement of the estimate across the 5 pTau217 platforms (rank correlation, kappa of class labels)", "kappa_and_rank_r"
)
write_csv(specs, file.path(OUT_DIR, "frame_paper3_specs.csv"))

# ---- 4. Save and summarise -------------------------------------------------------------
save_rds(person, "frame_paper3_person")
save_rds(long,   "frame_paper3_long")
write_csv(person, file.path(OUT_DIR, "frame_paper3_person.csv"))

outc <- tibble(
  outcome = c("amyloid positive at first plasma (CL20)", "baseline negative, positive at last plasma (incident)",
              "baseline negative, negative at last plasma", "ever positive on any scan",
              "first positive scan after last plasma", "DX progressed first to last plasma visit",
              "concurrent centiloid slope (2+ scans)", "concurrent CSF pTau/Ab42 slope (2+)",
              "concurrent MMSE slope (2+)", "concurrent hcv_adj slope (2+)"),
  n = c(sum(person$prevalent_amyloid_cl20, na.rm = TRUE), sum(person$incident_amyloid_cl20, na.rm = TRUE),
        sum(person$baseline_negative_cl20 & !person$amyloid_last_cl20, na.rm = TRUE),
        sum(person$amyloid_ever_cl20, na.rm = TRUE), sum(person$first_positive_after_last_plasma, na.rm = TRUE),
        sum(person$dx_progressed, na.rm = TRUE),
        sum(person$centiloid_n >= 2, na.rm = TRUE), sum(person$csf_ptau_abeta42_n >= 2, na.rm = TRUE),
        sum(person$mmse_n >= 2, na.rm = TRUE), if ("hcv_adj_n" %in% names(person)) sum(person$hcv_adj_n >= 2, na.rm = TRUE) else 0L)
)
print(as.data.frame(outc))
write_csv(outc, file.path(CHECK_DIR, "paper3_outcomes.csv"))
say("frame_paper3_long: %d rows, %d people, %d measures; person frame %d rows",
    nrow(long), n_distinct(long$RID), n_distinct(long$measure), nrow(person))


# =============================================================================================
# Section 4. A4/LEARN frame check against the ADNI frames (A4/LEARN only)
# Original script: 09e_a4_frame_check.R
# =============================================================================================
if (Sys.getenv("COHORT", "adni") == "a4learn") {
  # =============================================================================
  # 09e_a4_frame_check.R
  # After 09a to 09d, 06_master.R and 07a to 07c have been run with COHORT = "a4learn":
  # compare the A4/LEARN frames with the ADNI frames column by column, and write
  # the flow and design tables the notes need.
  #
  # Output   a4learn_out/checks/a4_frame_columns.csv   column present in ADNI / A4LEARN / both
  #          a4learn_out/checks/a4_flow_counts.csv     people by substudy and arm at each stage
  #          a4learn_out/checks/a4_analysis_samples.csv  the samples the analysis scripts will use
  # =============================================================================

  source("00_setup.R")
  if (COHORT != "a4learn") stop("Set Sys.setenv(COHORT = 'a4learn') before running the 09 scripts.")
  ADNI_OUT <- file.path(OUT_ROOT, "pipeline_out")

  frames <- c("plasma_long", "participants", "criteria_long", "modifiers_long", "master_long", "master_visit",
              "frame_paper2_long", "frame_paper2_slopes",
              "frame_paper2_criteria", "frame_paper3_long", "frame_paper3_person")
  cols <- bind_rows(lapply(frames, function(f) {
    a <- file.path(OUT_DIR, paste0(f, ".rds")); b <- file.path(ADNI_OUT, paste0(f, ".rds"))
    na <- if (file.exists(a)) names(readRDS(a)) else character()
    nb <- if (file.exists(b)) names(readRDS(b)) else character()
    tibble(frame = f, column = union(nb, na), in_adni = column %in% nb, in_a4learn = column %in% na)
  })) %>% mutate(status = case_when(in_adni & in_a4learn ~ "both", in_adni ~ "ADNI only", TRUE ~ "A4/LEARN only"))
  write_csv(cols, file.path(CHECK_DIR, "a4_frame_columns.csv"))
  say("Columns by frame:"); print(as.data.frame(cols %>% count(frame, status) %>% pivot_wider(names_from = status, values_from = n, values_fill = 0)))
  if (!dir.exists(ADNI_OUT)) say("(ADNI pipeline_out not found at %s; the ADNI side of the comparison is empty)", ADNI_OUT)

  part <- readRDS(file.path(OUT_DIR, "participants.rds"))
  p3   <- readRDS(file.path(OUT_DIR, "frame_paper3_person.rds")) %>% left_join(part %>% select(RID, substudy, arm), by = "RID")
  flow <- bind_rows(
    part %>% count(substudy, arm, name = "n") %>% mutate(stage = "any plasma draw"),
    part %>% filter(n_plasma_visits >= 2) %>% count(substudy, arm, name = "n") %>% mutate(stage = "2+ draws"),
    part %>% filter(n_plasma_visits >= 3) %>% count(substudy, arm, name = "n") %>% mutate(stage = "3+ draws"),
    part %>% filter(has_short_interval_pair) %>% count(substudy, arm, name = "n") %>% mutate(stage = "baseline and week-12 pair"),
    p3 %>% filter(n_plasma_visits >= 2, centiloid_n >= 2) %>% count(substudy, arm, name = "n") %>% mutate(stage = "2+ draws and 2+ amyloid scans"),
    p3 %>% filter(n_plasma_visits >= 3, centiloid_n >= 2) %>% count(substudy, arm, name = "n") %>% mutate(stage = "3+ draws and 2+ amyloid scans"),
    p3 %>% filter(n_plasma_visits >= 2, csf_ptau_abeta42_n >= 2) %>% count(substudy, arm, name = "n") %>% mutate(stage = "2+ draws and 2 CSF")
  ) %>% select(stage, substudy, arm, n) %>% pivot_wider(names_from = c(substudy, arm), values_from = n, values_fill = 0)
  write_csv(flow, file.path(CHECK_DIR, "a4_flow_counts.csv"))
  print(as.data.frame(flow))

  samples <- tibble(
    sample = c("primary: placebo + LEARN", "sensitivity: all three arms", "within A4 placebo", "within LEARN", "within A4 solanezumab"),
    n_people = c(sum(part$arm %in% c("Placebo", "LEARN")), nrow(part), sum(part$arm == "Placebo"), sum(part$arm == "LEARN"), sum(part$arm == "Solanezumab")),
    n_with_outcome = c(sum(p3$arm %in% c("Placebo", "LEARN") & p3$centiloid_n >= 2, na.rm = TRUE), sum(p3$centiloid_n >= 2, na.rm = TRUE),
                       sum(p3$arm == "Placebo" & p3$centiloid_n >= 2, na.rm = TRUE), sum(p3$arm == "LEARN" & p3$centiloid_n >= 2, na.rm = TRUE),
                       sum(p3$arm == "Solanezumab" & p3$centiloid_n >= 2, na.rm = TRUE)))
  write_csv(samples, file.path(CHECK_DIR, "a4_analysis_samples.csv"))
  print(as.data.frame(samples))
}
