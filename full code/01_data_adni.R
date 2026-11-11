# =============================================================================================
# 01_data_adni.R
# ADNI data build: plasma table, participants, reference measures, modifiers, checks.
#
# Reads the raw ADNI files (read only) and writes pipeline_out/: plasma_long, participants, criteria_long,
# modifiers_long and the check tables. Five sections, formerly scripts 01, 02, 03, 04 and 99.
# Run by run_all.R with COHORT = adni; takes a few minutes and fits nothing.
# =============================================================================================


# =============================================================================================
# Section 1. Plasma measurements, one row per person, visit and measure
# Original script: 01_plasma_long.R
# =============================================================================================
# =============================================================================
# 01_plasma_long.R
# Build the tall plasma master from FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv.
#
# Output: one row per person, per visit, per measure (platform x biomarker).
#         Saved as plasma_long.rds. A convenience wide frame, one row per
#         person per visit with one column per measure, is saved alongside as
#         plasma_wide.rds. The wide frame is derived from the tall one and is
#         never edited directly.
#
# A visit is one plasma sample. Dates within VISIT_GAP_DAYS (7) of each other
# for one person are treated as the same sample; on the real file this never
# happens, so visit = distinct date = distinct SAMPLEID prefix. Two indices
# are carried: visit_no (shared across platforms, numbers every sample) and
# occasion_no (numbers the samples each platform actually measured, 1..k).
# Quanterix measured a partly different sample set from the other five labs,
# so quanterix_only flags samples and people nobody else measured.
# =============================================================================


# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
Sys.setenv(COHORT = "adni"); source("00_setup.R")

say("Reading %s", RAW_FILES$plasma)
raw <- read_csv(RAW_FILES$plasma, col_types = cols(.default = col_character()),
                show_col_types = FALSE)
say("Raw rows: %d  columns: %d", nrow(raw), ncol(raw))

# ---- 1. Keep the columns we use, parse the simple ones -----------------------
plasma <- raw %>%
  transmute(
    RID            = as.integer(RID),
    PTID,
    SAMPLEID,
    sample_prefix  = if_else(!is.na(SAMPLEID) & SAMPLEID != "",
                             str_replace(SAMPLEID, "-[^-]*$", ""), NA_character_),
    EXAMDATE_row   = parse_date_dmy(EXAMDATE),
    VISCODE2_row   = recode_viscode(VISCODE2),
    PLASMA_BIOMARKER,
    TESTNAME,
    ASSAYPLATFORM,
    ASSAYVERSION,
    PERFORMINGLAB,
    INSTRUMENTID,
    MATRIX,
    TESTVALUE_raw  = TESTVALUE,
    unit_raw       = UNITS,
    CV             = suppressWarnings(as.numeric(CV)),
    RUN_NO         = as.character(RUN_NO),
    SEQUENCE_NO    = as.character(SEQUENCE_NO),
    RUNDATE        = parse_date_dmy(RUNDATE),
    COMMENTS
  )

# ---- 2. Visit index: cluster each person's dates ------------------------------
visit_index <- plasma %>%
  filter(!is.na(EXAMDATE_row)) %>%
  distinct(RID, EXAMDATE_row) %>%
  arrange(RID, EXAMDATE_row) %>%
  group_by(RID) %>%
  mutate(
    gap       = as.integer(EXAMDATE_row - lag(EXAMDATE_row)),
    new_visit = is.na(gap) | gap > VISIT_GAP_DAYS,
    visit_no  = cumsum(new_visit)
  ) %>%
  group_by(RID, visit_no) %>%
  mutate(
    visit_date     = min(EXAMDATE_row),
    date_span_days = as.integer(max(EXAMDATE_row) - min(EXAMDATE_row))
  ) %>%
  ungroup() %>%
  select(RID, EXAMDATE_row, visit_no, visit_date, date_span_days)

plasma <- plasma %>%
  left_join(visit_index, by = c("RID", "EXAMDATE_row")) %>%
  group_by(RID) %>%
  mutate(plasma_time = as.numeric(visit_date - min(visit_date, na.rm = TRUE)) / YEAR_DAYS) %>%
  group_by(RID, visit_no) %>%
  mutate(VISCODE2 = first(na.omit(VISCODE2_row))) %>%
  ungroup()

# ---- 3. Value and detection-limit flags ---------------------------------------
# The value column may carry a comparator such as "<0.1". The COMMENTS column
# may say "LOD", "LLOQ", or that a value was imputed. Capture all of it before
# coercing to numeric so nothing is lost.
plasma <- plasma %>%
  mutate(
    comparator   = str_trim(str_extract(TESTVALUE_raw, "^\\s*[<>]")),
    value_raw    = if_else(str_detect(coalesce(TESTVALUE_raw, ""), "[0-9]"),
                           suppressWarnings(readr::parse_number(TESTVALUE_raw)),
                           NA_real_),
    comment_lc   = str_to_lower(coalesce(COMMENTS, "")),
    testval_lc   = str_to_lower(coalesce(TESTVALUE_raw, "")),
    lloq_flag    = (!is.na(comparator) & comparator == "<") |
                   str_detect(comment_lc, "lod|lloq|below|detection") |
                   str_detect(testval_lc, "lod|lloq"),
    uloq_flag    = !is.na(comparator) & comparator == ">",
    imputed_flag = str_detect(comment_lc, "imput")
  ) %>%
  select(-comment_lc, -testval_lc)

# ---- 4. Platform, biomarker, role ----------------------------------------------
plasma <- plasma %>%
  mutate(platform = platform_from_prefix(PLASMA_BIOMARKER)) %>%
  left_join(MEASURE_LOOKUP, by = "PLASMA_BIOMARKER") %>%
  mutate(role = coalesce(role, "unmapped"))

unmapped <- plasma %>% filter(role == "unmapped") %>% distinct(PLASMA_BIOMARKER)
if (nrow(unmapped) > 0) {
  say("WARNING: %d PLASMA_BIOMARKER values are not in MEASURE_LOOKUP:", nrow(unmapped))
  print(unmapped)
}

# ---- 5. Standardise units to pg/mL where the unit is a concentration -----------
plasma <- plasma %>%
  mutate(
    unit_lc  = str_to_lower(str_trim(coalesce(unit_raw, ""))),
    is_ng    = unit_lc == "ng/ml",
    unit_std = case_when(
      is_ng               ~ "pg/mL",
      unit_lc == "pg/ml"  ~ "pg/mL",
      unit_lc == "%"      ~ "%",
      unit_lc == ""       ~ NA_character_,
      TRUE                ~ unit_raw
    ),
    value    = if_else(is_ng, value_raw * 1000, value_raw),
    derived  = FALSE
  ) %>%
  select(-unit_lc, -is_ng)

# ---- 6. Derive Ab42/40 for platforms that do not supply it --------------------
components <- plasma %>%
  filter(role == "component", biomarker == "Ab42_40",
         platform %in% DERIVE_RATIO_FOR)

derived_ratio <- components %>%
  group_by(RID, visit_no, platform) %>%
  summarise(
    PTID           = first(PTID),
    SAMPLEID       = first(SAMPLEID),
    sample_prefix  = first(sample_prefix),
    EXAMDATE_row   = first(EXAMDATE_row),
    VISCODE2_row   = first(VISCODE2_row),
    VISCODE2       = first(VISCODE2),
    visit_date     = first(visit_date),
    date_span_days = first(date_span_days),
    plasma_time    = first(plasma_time),
    ASSAYPLATFORM  = first(ASSAYPLATFORM),
    ASSAYVERSION   = first(ASSAYVERSION),
    PERFORMINGLAB  = first(PERFORMINGLAB),
    INSTRUMENTID   = first(INSTRUMENTID),
    MATRIX         = first(MATRIX),
    RUN_NO         = first(RUN_NO),
    SEQUENCE_NO    = first(SEQUENCE_NO),
    RUNDATE        = first(RUNDATE),
    Ab40           = value[component == "Ab40"][1],
    Ab42           = value[component == "Ab42"][1],
    lloq_flag      = any(lloq_flag,    na.rm = TRUE),
    imputed_flag   = any(imputed_flag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(Ab40), !is.na(Ab42), Ab40 > 0) %>%
  mutate(
    PLASMA_BIOMARKER = case_when(
      platform == "Fujirebio" ~ "Fuji_plasma_Ab42_Ab40",
      platform == "Roche"     ~ "Roche_plasma_Ab42_Ab40",
      platform == "Quanterix" ~ "QX_plasma_Ab42_Ab40"
    ),
    TESTNAME      = "Abeta_ratio_derived",
    TESTVALUE_raw = NA_character_,
    value_raw     = Ab42 / Ab40,
    value         = Ab42 / Ab40,
    unit_raw      = NA_character_,
    unit_std      = NA_character_,
    CV            = NA_real_,
    COMMENTS      = "derived in pipeline from Ab42 and Ab40",
    comparator    = NA_character_,
    uloq_flag     = FALSE,
    biomarker     = "Ab42_40",
    role          = "measure",
    component     = NA_character_,
    direction     = -1,
    derived       = TRUE
  ) %>%
  select(-Ab40, -Ab42)

plasma <- bind_rows(plasma, derived_ratio) %>%
  arrange(RID, visit_no, biomarker, platform, PLASMA_BIOMARKER)

# ---- 6b. Occasion index per platform, and the Quanterix-only flags ----------
# occasion_no: within person and platform, the rank of the visit among the
# visits that platform has any value for. For the five main labs this is
# 1, 2, 3 for everyone; visit_no can be 1, 2, 4 when Quanterix measured an
# extra sample in between.
occasion_index <- plasma %>%
  filter(!is.na(value)) %>%
  distinct(RID, platform, visit_no) %>%
  arrange(RID, platform, visit_no) %>%
  group_by(RID, platform) %>%
  mutate(occasion_no = row_number()) %>%
  ungroup()

# A sample is Quanterix-only if no other platform has a value for it; a
# person is Quanterix-only if that is true of all their samples.
sample_platforms <- plasma %>%
  filter(!is.na(value)) %>%
  group_by(RID, visit_no) %>%
  summarise(n_platforms_at_visit = n_distinct(platform),
            quanterix_only_visit = all(platform == "Quanterix"),
            .groups = "drop")
person_flags <- sample_platforms %>%
  group_by(RID) %>%
  summarise(quanterix_only_person = all(quanterix_only_visit), .groups = "drop")

plasma <- plasma %>%
  left_join(occasion_index,  by = c("RID", "platform", "visit_no")) %>%
  left_join(sample_platforms, by = c("RID", "visit_no")) %>%
  left_join(person_flags,    by = "RID")
say("Quanterix-only people: %d; Quanterix-only samples of other people: %d",
    sum(person_flags$quanterix_only_person),
    sample_platforms %>% filter(quanterix_only_visit) %>%
      anti_join(person_flags %>% filter(quanterix_only_person), by = "RID") %>% nrow())

# ---- 7. Final column order and save ---------------------------------------------
# visit_date is the visit-level date shared by every platform; EXAMDATE_row is
# what the lab reported for that row. Use visit_no and visit_date for analysis.
plasma_long <- plasma %>%
  select(
    RID, PTID, visit_no, occasion_no, visit_date, plasma_time, VISCODE2,
    n_platforms_at_visit, quanterix_only_visit, quanterix_only_person,
    SAMPLEID, sample_prefix, EXAMDATE_row, VISCODE2_row, date_span_days,
    biomarker, platform, PLASMA_BIOMARKER, TESTNAME, role, component, direction,
    value, value_raw, TESTVALUE_raw, unit_std, unit_raw,
    lloq_flag, uloq_flag, imputed_flag, comparator, derived,
    ASSAYPLATFORM, ASSAYVERSION, PERFORMINGLAB, INSTRUMENTID, MATRIX,
    CV, RUN_NO, SEQUENCE_NO, RUNDATE, COMMENTS
  )

save_rds(plasma_long, "plasma_long")
say("plasma_long: %d rows, %d people, %d visits, %d measures",
    nrow(plasma_long), n_distinct(plasma_long$RID),
    nrow(distinct(plasma_long, RID, visit_no)),
    n_distinct(plasma_long$PLASMA_BIOMARKER[plasma_long$role == "measure"]))

# ---- 8. Convenience wide frame (derived, never edited) --------------------------
plasma_wide <- plasma_long %>%
  filter(role == "measure") %>%
  select(RID, visit_no, visit_date, VISCODE2, plasma_time,
         quanterix_only_visit, quanterix_only_person, PLASMA_BIOMARKER, value) %>%
  pivot_wider(names_from = PLASMA_BIOMARKER, values_from = value,
              values_fn = function(x) if (length(x) > 1) NA_real_ else x) %>%
  arrange(RID, visit_no)
# If two samples were merged into one visit the pivot above would have had
# duplicates; they are set to NA rather than silently taking one. Duplicates
# are listed by 99_checks.R (duplicate_rows.csv) and should be zero.

save_rds(plasma_wide, "plasma_wide")
write_csv(plasma_wide, file.path(OUT_DIR, "plasma_wide.csv"))
say("plasma_wide: %d rows, %d columns", nrow(plasma_wide), ncol(plasma_wide))


# =============================================================================================
# Section 2. Participants: sampling summary, analytic flags, baseline characteristics
# Original script: 02_participants.R
# =============================================================================================
# =============================================================================
# 02_participants.R
# One row per person: plasma sampling summary, analytic-sample flags, and
# baseline characteristics from ADNIMERGE and PTDEMOG.
#
# Output: participants.rds (and .csv), plus checks/flow_counts.csv giving
#         the numbers for a participant flow diagram.
#
# Age at first plasma sample is ADNIMERGE baseline AGE plus the years between
# the ADNI baseline exam and the first plasma sample. PTDEMOG birth year and
# month are kept alongside so the two can be compared.
# =============================================================================

# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
Sys.setenv(COHORT = "adni"); source("00_setup.R")

pl <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))

# ---- 1. Plasma sampling summary per person -----------------------------------
plasma_person <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  group_by(RID) %>%
  summarise(
    PTID                 = first(PTID),
    n_plasma_visits      = n_distinct(visit_no),
    first_plasma         = min(visit_date),
    last_plasma          = max(visit_date),
    plasma_span_years    = as.numeric(max(visit_date) - min(visit_date)) / YEAR_DAYS,
    n_visits_5lab        = n_distinct(visit_no[platform != "Quanterix"]),
    n_visits_quanterix   = n_distinct(visit_no[platform == "Quanterix"]),
    n_visits_all6        = n_distinct(visit_no[n_platforms_at_visit == 6]),
    quanterix_only_person = first(quanterix_only_person),
    .groups = "drop"
  ) %>%
  mutate(
    # Core sample: measured on the five main platforms. Quanterix sample: also
    # has Quanterix on at least two of those samples.
    in_core_5lab   = n_visits_5lab >= 1,
    in_quanterix   = n_visits_quanterix >= 1 & !quanterix_only_person
  )
say("People in plasma file: %d; core (5-lab) sample: %d; with Quanterix too: %d; Quanterix-only: %d",
    nrow(plasma_person), sum(plasma_person$in_core_5lab),
    sum(plasma_person$in_quanterix), sum(plasma_person$quanterix_only_person))

# ---- 2. ADNIMERGE baseline row -----------------------------------------------
say("Reading %s", RAW_FILES$merge)
merge_raw <- read_adni(RAW_FILES$merge)
m_vis  <- pick_col(merge_raw, c("VISCODE", "VISCODE2"))
m_date <- pick_col(merge_raw, c("EXAMDATE"))
m_dx   <- pick_col(merge_raw, c("DX"), required = FALSE)
m_dxbl <- pick_col(merge_raw, c("DX_BL"), required = FALSE)
m_col  <- pick_col(merge_raw, c("COLPROT"), required = FALSE)
m_orig <- pick_col(merge_raw, c("ORIGPROT"), required = FALSE)
m_cols <- c(AGE = "AGE", PTGENDER = "PTGENDER", PTEDUCAT = "PTEDUCAT", APOE4 = "APOE4",
            PTETHCAT = "PTETHCAT", PTRACCAT = "PTRACCAT", PTMARRY = "PTMARRY")
m_cols <- m_cols[m_cols %in% names(merge_raw)]

merge_all <- merge_raw %>%
  mutate(RID = as.integer(RID),
         exam_date = parse_date_ymd(.data[[m_date]]),
         viscode   = recode_viscode(.data[[m_vis]])) %>%
  filter(RID %in% plasma_person$RID)

# Baseline = the "bl" row; if a person has none, the earliest dated row.
merge_bl <- merge_all %>%
  arrange(RID, viscode != "bl", exam_date) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  transmute(
    RID,
    adni_bl_date    = exam_date,
    adni_bl_viscode = viscode,
    across(all_of(unname(m_cols)), ~ .x),
    DX_bl    = if (!is.na(m_dxbl)) .data[[m_dxbl]] else NA_character_,
    COLPROT  = if (!is.na(m_col))  .data[[m_col]]  else NA_character_,
    ORIGPROT = if (!is.na(m_orig)) .data[[m_orig]] else NA_character_
  ) %>%
  mutate(AGE = suppressWarnings(as.numeric(AGE)),
         PTEDUCAT = suppressWarnings(as.numeric(PTEDUCAT)),
         APOE4 = suppressWarnings(as.integer(APOE4)))
say("ADNIMERGE rows for plasma people: %d; people with a baseline row: %d (of which not coded bl: %d)",
    nrow(merge_all), nrow(merge_bl), sum(merge_bl$adni_bl_viscode != "bl"))

# Diagnosis at the ADNIMERGE visit nearest the first plasma sample
dx_first <- NULL
if (!is.na(m_dx)) {
  dx_first <- merge_all %>%
    filter(!is.na(.data[[m_dx]]), .data[[m_dx]] != "") %>%
    transmute(RID, dx_date = exam_date, DX_at_first_plasma = .data[[m_dx]]) %>%
    inner_join(plasma_person %>% select(RID, first_plasma), by = "RID") %>%
    mutate(dx_gap_days = as.integer(dx_date - first_plasma)) %>%
    group_by(RID) %>% slice_min(abs(dx_gap_days), n = 1, with_ties = FALSE) %>% ungroup() %>%
    select(RID, DX_at_first_plasma, dx_gap_days)
}

# ---- 3. PTDEMOG: sex, birth month/year, education as a cross-check -----------
say("Reading %s", RAW_FILES$demog)
demog_raw <- read_adni(RAW_FILES$demog)
d_yy <- pick_col(demog_raw, c("PTDOBYY"), required = FALSE)
d_mm <- pick_col(demog_raw, c("PTDOBMM"), required = FALSE)
d_dob <- pick_col(demog_raw, c("PTDOB"), required = FALSE)
demog <- demog_raw %>%
  mutate(RID = as.integer(RID)) %>%
  filter(RID %in% plasma_person$RID) %>%
  transmute(
    RID,
    PTGENDER_demog = if ("PTGENDER" %in% names(demog_raw)) PTGENDER else NA_character_,
    PTDOBYY = if (!is.na(d_yy)) suppressWarnings(as.integer(.data[[d_yy]])) else NA_integer_,
    PTDOBMM = if (!is.na(d_mm)) suppressWarnings(as.integer(.data[[d_mm]])) else NA_integer_,
    PTDOB_raw = if (!is.na(d_dob)) .data[[d_dob]] else NA_character_
  ) %>%
  group_by(RID) %>%
  summarise(across(everything(), ~ first(na.omit(.x))), .groups = "drop") %>%
  mutate(PTGENDER_demog = case_when(PTGENDER_demog %in% c("1", "Male")   ~ "Male",
                                    PTGENDER_demog %in% c("2", "Female") ~ "Female",
                                    TRUE ~ PTGENDER_demog))
# If PTDOB is "MM/YYYY" style, fill month and year from it
if (any(!is.na(demog$PTDOB_raw)) && all(is.na(demog$PTDOBYY))) {
  demog <- demog %>%
    mutate(PTDOBMM = suppressWarnings(as.integer(str_extract(PTDOB_raw, "^[0-9]{1,2}"))),
           PTDOBYY = suppressWarnings(as.integer(str_extract(PTDOB_raw, "[0-9]{4}"))))
}

# ---- 4. Assemble --------------------------------------------------------------
participants <- plasma_person %>%
  left_join(merge_bl, by = "RID") %>%
  left_join(demog, by = "RID") %>%
  mutate(
    age_first_plasma = AGE + as.numeric(first_plasma - adni_bl_date) / YEAR_DAYS,
    # Cross-check from birth year/month, mid-month
    # PTDEMOG codes missing month or year as -4 or similar; only build a date
    # from plausible values, and parse strictly so a bad row gives NA.
    dob_ok = !is.na(PTDOBYY) & !is.na(PTDOBMM) & PTDOBMM >= 1 & PTDOBMM <= 12 &
             PTDOBYY >= 1900 & PTDOBYY <= 2000,
    dob_mid = as.Date(ifelse(dob_ok, sprintf("%04d-%02d-15", PTDOBYY, PTDOBMM), NA_character_),
                      format = "%Y-%m-%d"),
    age_first_plasma_dob = as.numeric(first_plasma - dob_mid) / YEAR_DAYS,
    age_diff_check = age_first_plasma - age_first_plasma_dob,
    sex_mismatch   = !is.na(PTGENDER_demog) & !is.na(PTGENDER) & PTGENDER_demog != PTGENDER,
    APOE4_carrier  = case_when(is.na(APOE4) ~ NA, APOE4 >= 1 ~ TRUE, TRUE ~ FALSE),
    in_merge       = !is.na(adni_bl_date)
  )
if (!is.null(dx_first)) participants <- participants %>% left_join(dx_first, by = "RID")

say("Age at first plasma: median %.1f (range %.1f to %.1f); missing %d",
    median(participants$age_first_plasma, na.rm = TRUE),
    min(participants$age_first_plasma, na.rm = TRUE),
    max(participants$age_first_plasma, na.rm = TRUE),
    sum(is.na(participants$age_first_plasma)))
say("Age check vs PTDEMOG birth date: n=%d, max |diff| %.2f years; sex mismatches: %d",
    sum(!is.na(participants$age_diff_check)),
    max(abs(participants$age_diff_check), na.rm = TRUE),
    sum(participants$sex_mismatch, na.rm = TRUE))

save_rds(participants, "participants")
write_csv(participants, file.path(OUT_DIR, "participants.csv"))

# ---- 5. Flow counts -----------------------------------------------------------
flow <- tribble(
  ~step, ~n,
  "People in FNIH plasma file",                     nrow(participants),
  "  Quanterix-only (excluded from cross-platform)", sum(participants$quanterix_only_person),
  "Core sample: 5 main platforms",                  sum(participants$in_core_5lab),
  "  with 3 samples on the 5 platforms",            sum(participants$n_visits_5lab == 3),
  "  also measured by Quanterix",                   sum(participants$in_quanterix),
  "  with all 6 platforms on all 3 samples",        sum(participants$n_visits_all6 == 3),
  "Core sample matched to ADNIMERGE",               sum(participants$in_core_5lab & participants$in_merge),
  "  with APOE4",                                   sum(participants$in_core_5lab & !is.na(participants$APOE4)),
  "  with age at first plasma",                     sum(participants$in_core_5lab & !is.na(participants$age_first_plasma))
)
print(as.data.frame(flow))
write_csv(flow, file.path(CHECK_DIR, "flow_counts.csv"))

# Baseline table for the core sample
bl_tab <- participants %>%
  filter(in_core_5lab) %>%
  summarise(
    n = n(),
    age_mean = mean(age_first_plasma, na.rm = TRUE), age_sd = sd(age_first_plasma, na.rm = TRUE),
    female_n = sum(PTGENDER == "Female", na.rm = TRUE),
    educ_mean = mean(PTEDUCAT, na.rm = TRUE),
    apoe4_carrier_n = sum(APOE4_carrier, na.rm = TRUE), apoe4_missing = sum(is.na(APOE4)),
    span_median = median(plasma_span_years), span_min = min(plasma_span_years), span_max = max(plasma_span_years)
  )
print(as.data.frame(bl_tab))
write_csv(bl_tab, file.path(CHECK_DIR, "baseline_core_sample.csv"))
if ("DX_bl" %in% names(participants)) {
  dx_tab <- participants %>% filter(in_core_5lab) %>% count(DX_bl, DX_at_first_plasma)
  print(as.data.frame(dx_tab))
  write_csv(dx_tab, file.path(CHECK_DIR, "diagnosis_core_sample.csv"))
}
say("participants: %d rows written", nrow(participants))


# =============================================================================================
# Section 3. Reference measures: amyloid PET, CSF, cognition, MRI, diagnosis
# Original script: 03_criteria_long.R
# =============================================================================================
# =============================================================================
# 03_criteria_long.R
# The outcome and criterion measures in the same tall shape as the plasma
# file: one row per person per assessment date per measure, with the nearest
# plasma visit attached.
#
# Sources and measures
#   amyloid   UCBERKELEY_AMY_6MM   centiloid, summary_suvr, amyloid_status (+ tracer)
#   csf       UPENN Roche Elecsys  csf_abeta42, csf_abeta40, csf_tau, csf_ptau181,
#                                  csf_ptau_abeta42, csf_abeta42_40, csf_ptau_tau (derived)
#   cognition ADNIMERGE            mmse, adas13, mpacc_digit, mpacc_trailsb, cdrsb
#   mri       ADNIMERGE            hippocampus, icv, hcv_adj (derived)
#   diagnosis ADNIMERGE            dx (character, kept in dx_long.rds separately)
#
# Matching: every assessment is linked to the plasma visit of the same person
# closest in time. gap_days is assessment date minus plasma date, so negative
# means the assessment came before the blood draw. No tolerance is applied
# here; each paper sets its own window when it builds its frame.
#
# Output: criteria_long.rds (numeric measures), dx_long.rds (diagnosis),
#         checks/criteria_coverage.csv
# =============================================================================

# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
Sys.setenv(COHORT = "adni"); source("00_setup.R")

pl <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
plasma_visits <- pl %>%
  distinct(RID, visit_no, visit_date) %>%
  group_by(RID) %>% mutate(first_plasma = min(visit_date)) %>% ungroup()
rids <- unique(plasma_visits$RID)

# A helper to bring one source into the tall shape
to_long <- function(df, date_col, viscode_col, measure_cols, source) {
  df %>%
    transmute(RID = as.integer(RID),
              exam_date = .data[[date_col]],
              VISCODE2  = if (!is.na(viscode_col)) recode_viscode(.data[[viscode_col]]) else NA_character_,
              across(all_of(measure_cols))) %>%
    filter(RID %in% rids, !is.na(exam_date)) %>%
    pivot_longer(all_of(measure_cols), names_to = "measure", values_to = "value_raw") %>%
    mutate(source = source)
}

# ---- 1. Amyloid PET ------------------------------------------------------------
say("Reading %s", RAW_FILES$amyloid)
amy_raw <- read_adni(RAW_FILES$amyloid)
a_date <- pick_col(amy_raw, c("SCANDATE", "EXAMDATE"))
a_vis  <- pick_col(amy_raw, c("VISCODE2", "VISCODE"), required = FALSE)
a_cl   <- pick_col(amy_raw, c("CENTILOIDS", "CENTILOID"))
a_suvr <- pick_col(amy_raw, c("SUMMARY_SUVR", "SUMMARYSUVR"), required = FALSE)
a_stat <- pick_col(amy_raw, c("AMYLOID_STATUS", "AMYLOID_STATUS_COMPOSITE_REF"), required = FALSE)
a_trac <- pick_col(amy_raw, c("TRACER"), required = FALSE)

amy <- amy_raw %>%
  mutate(scan_date = parse_date_ymd(.data[[a_date]]),
         centiloid = .data[[a_cl]],
         summary_suvr   = if (!is.na(a_suvr)) .data[[a_suvr]] else NA_character_,
         amyloid_status = if (!is.na(a_stat)) .data[[a_stat]] else NA_character_,
         tracer         = if (!is.na(a_trac)) .data[[a_trac]] else NA_character_)
amy_long <- to_long(amy, "scan_date", a_vis, c("centiloid", "summary_suvr", "amyloid_status"), "amyloid") %>%
  left_join(amy %>% transmute(RID = as.integer(RID), exam_date = scan_date, tracer) %>% distinct(),
            by = c("RID", "exam_date"))

# ---- 2. CSF Elecsys ------------------------------------------------------------
say("Reading %s", RAW_FILES$csf)
csf_raw <- read_adni(RAW_FILES$csf)
c_date <- pick_col(csf_raw, c("EXAMDATE", "DRAWDATE"))
c_vis  <- pick_col(csf_raw, c("VISCODE2", "VISCODE"), required = FALSE)
c_map  <- c(csf_abeta42 = "ABETA42", csf_abeta40 = "ABETA40", csf_tau = "TAU", csf_ptau181 = "PTAU")
c_map  <- c_map[c_map %in% names(csf_raw)]
if (!"ABETA42" %in% c_map) stop("CSF file has no ABETA42 column. Columns: ", paste(names(csf_raw), collapse = ", "))

csf <- csf_raw %>% mutate(draw_date = parse_date_ymd(.data[[c_date]]))
for (nm in names(c_map)) csf[[nm]] <- csf[[c_map[[nm]]]]
csf_long <- to_long(csf, "draw_date", c_vis, names(c_map), "csf")

# ---- 3. Cognition, MRI and diagnosis from ADNIMERGE ----------------------------
say("Reading %s", RAW_FILES$merge)
merge_raw <- read_adni(RAW_FILES$merge)
m_date <- pick_col(merge_raw, c("EXAMDATE"))
m_vis  <- pick_col(merge_raw, c("VISCODE2", "VISCODE"))
m_map  <- c(mmse = "MMSE", adas13 = "ADAS13", mpacc_digit = "MPACCDIGIT",
            mpacc_trailsb = "MPACCTRAILSB", cdrsb = "CDRSB",
            hippocampus = "HIPPOCAMPUS", icv = "ICV")
present <- m_map[m_map %in% names(merge_raw)]
absent  <- m_map[!m_map %in% names(merge_raw)]
if (length(absent) > 0) say("ADNIMERGE columns not found (skipped): %s", paste(absent, collapse = ", "))

merge_df <- merge_raw %>% mutate(exam_date = parse_date_ymd(.data[[m_date]]))
for (nm in names(present)) merge_df[[nm]] <- merge_df[[present[[nm]]]]
cog_cols <- intersect(names(present), c("mmse", "adas13", "mpacc_digit", "mpacc_trailsb", "cdrsb"))
mri_cols <- intersect(names(present), c("hippocampus", "icv"))
cog_long <- to_long(merge_df, "exam_date", m_vis, cog_cols, "cognition")
mri_long <- if (length(mri_cols) > 0) to_long(merge_df, "exam_date", m_vis, mri_cols, "mri") else NULL

dx_long <- NULL
if ("DX" %in% names(merge_raw)) {
  dx_long <- merge_df %>%
    transmute(RID = as.integer(RID), exam_date, VISCODE2 = recode_viscode(.data[[m_vis]]), dx = DX) %>%
    filter(RID %in% rids, !is.na(exam_date), !is.na(dx), dx != "")
}

# ---- 4. Combine, parse values, derive ratios ----------------------------------
crit <- bind_rows(amy_long, csf_long, cog_long, mri_long) %>%
  filter(!is.na(value_raw), value_raw != "")
pc <- parse_censored(crit$value_raw)
crit <- crit %>%
  mutate(value = pc$value, comparator = pc$comparator,
         lloq_flag = !is.na(comparator) & comparator == "<",
         uloq_flag = !is.na(comparator) & comparator == ">",
         derived = FALSE) %>%
  filter(!is.na(value)) %>%
  distinct(RID, exam_date, measure, .keep_all = TRUE)

# Derived CSF ratios and adjusted hippocampal volume, computed within a date
derive <- function(df, new, num, den, f = function(a, b) a / b, src) {
  df %>% filter(measure %in% c(num, den)) %>%
    select(RID, exam_date, VISCODE2, measure, value, lloq_flag, uloq_flag) %>%
    pivot_wider(names_from = measure, values_from = c(value, lloq_flag, uloq_flag)) %>%
    filter(!is.na(.data[[paste0("value_", num)]]), !is.na(.data[[paste0("value_", den)]])) %>%
    transmute(RID, exam_date, VISCODE2, measure = new, source = src,
              value_raw = NA_character_,
              value = f(.data[[paste0("value_", num)]], .data[[paste0("value_", den)]]),
              comparator = NA_character_,
              lloq_flag = .data[[paste0("lloq_flag_", num)]] | .data[[paste0("lloq_flag_", den)]],
              uloq_flag = .data[[paste0("uloq_flag_", num)]] | .data[[paste0("uloq_flag_", den)]],
              derived = TRUE)
}
derived_rows <- list()
if (all(c("csf_ptau181", "csf_abeta42") %in% crit$measure))
  derived_rows$a <- derive(crit, "csf_ptau_abeta42", "csf_ptau181", "csf_abeta42", src = "csf")
if (all(c("csf_abeta42", "csf_abeta40") %in% crit$measure))
  derived_rows$b <- derive(crit, "csf_abeta42_40", "csf_abeta42", "csf_abeta40", src = "csf")
if (all(c("csf_ptau181", "csf_tau") %in% crit$measure))
  derived_rows$c <- derive(crit, "csf_ptau_tau", "csf_ptau181", "csf_tau", src = "csf")
if (all(c("hippocampus", "icv") %in% crit$measure)) {
  # Head-size adjusted hippocampal volume. The slope 0.0017 and reference ICV
  # 1,498,858 mm3 are the values used in the earlier double-cohort analysis;
  # re-estimate the slope on the current sample before using this in a paper.
  derived_rows$d <- derive(crit, "hcv_adj", "hippocampus", "icv",
                           f = function(h, i) h - 0.0017 * (i - 1498858), src = "mri")
}
crit <- bind_rows(crit, bind_rows(derived_rows))

# ---- 5. Match every assessment to the nearest plasma visit -------------------
crit <- match_to_plasma(crit, "exam_date", plasma_visits) %>%
  arrange(RID, source, measure, exam_date) %>%
  group_by(RID, measure) %>%
  mutate(assessment_no = row_number(),
         n_assessments = n()) %>%
  ungroup() %>%
  select(RID, source, measure, exam_date, VISCODE2, assessment_no, n_assessments,
         value, value_raw, comparator, lloq_flag, uloq_flag, derived,
         visit_no, visit_date, gap_days, plasma_time, tracer)

save_rds(crit, "criteria_long")
write_csv(crit, file.path(OUT_DIR, "criteria_long.csv"))
say("criteria_long: %d rows, %d people, %d measures",
    nrow(crit), n_distinct(crit$RID), n_distinct(crit$measure))

if (!is.null(dx_long)) {
  dx_long <- match_to_plasma(dx_long, "exam_date", plasma_visits) %>%
    arrange(RID, exam_date)
  save_rds(dx_long, "dx_long")
  say("dx_long: %d rows, %d people", nrow(dx_long), n_distinct(dx_long$RID))
}

# ---- 6. Coverage: how much of each criterion sits near a plasma visit --------
cov <- crit %>%
  group_by(source, measure) %>%
  summarise(
    n_rows        = n(),
    n_people      = n_distinct(RID),
    n_people_ge2  = n_distinct(RID[n_assessments >= 2]),
    n_people_ge3  = n_distinct(RID[n_assessments >= 3]),
    n_within_90d  = sum(abs(gap_days) <= 90),
    n_within_180d = sum(abs(gap_days) <= 180),
    n_within_365d = sum(abs(gap_days) <= 365),
    people_within_180d = n_distinct(RID[abs(gap_days) <= 180]),
    median_abs_gap = median(abs(gap_days)),
    n_censored    = sum(lloq_flag | uloq_flag),
    first_date    = min(exam_date), last_date = max(exam_date),
    .groups = "drop"
  )
print(as.data.frame(cov))
write_csv(cov, file.path(CHECK_DIR, "criteria_coverage.csv"))

# Per plasma visit: is there an assessment of each criterion within 180 days?
near <- crit %>%
  filter(abs(gap_days) <= 180) %>%
  distinct(RID, visit_no, measure) %>%
  count(measure, visit_no, name = "n_plasma_visits_with_match") %>%
  pivot_wider(names_from = visit_no, values_from = n_plasma_visits_with_match,
              names_prefix = "visit_", values_fill = 0)
print(as.data.frame(near))
write_csv(near, file.path(CHECK_DIR, "criteria_match_by_plasma_visit.csv"))


# =============================================================================================
# Section 4. Modifiers: creatinine, eGFR, weight, height, BMI
# Original script: 04_modifiers.R
# =============================================================================================
# =============================================================================
# 04_modifiers.R
# Physiological modifiers in the same tall shape: creatinine and eGFR at every
# LABDATA visit, and weight, height and BMI at every VITALS visit if that file
# is present. Each assessment is matched to the nearest plasma visit.
#
#   lab     LABDATA.CSV   creatinine (mg/dL, RCT392), egfr_2021, egfr_2009
#   vitals  VITALS.csv    weight_kg, height_cm, bmi   (optional file)
#
# eGFR needs age at the lab visit and sex. Age is ADNIMERGE baseline age plus
# the years since the ADNI baseline exam, both taken from participants.rds.
# Height is missing at many ADNI visits (it is measured at screening), so BMI
# uses the person's median measured height if the visit has none.
#
# Output: modifiers_long.rds (and .csv), checks/modifiers_coverage.csv
# =============================================================================

# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
Sys.setenv(COHORT = "adni"); source("00_setup.R")

pl   <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
plasma_visits <- pl %>%
  distinct(RID, visit_no, visit_date) %>%
  group_by(RID) %>% mutate(first_plasma = min(visit_date)) %>% ungroup()
rids <- unique(plasma_visits$RID)
person_age <- part %>% select(RID, AGE, adni_bl_date, PTGENDER)

# ---- 1. Creatinine and eGFR from LABDATA ---------------------------------------
say("Reading %s", RAW_FILES$lab)
lab_raw <- read_adni(RAW_FILES$lab)
l_date <- pick_col(lab_raw, c("EXAMDATE", "VISDATE"))
l_vis  <- pick_col(lab_raw, c("VISCODE2", "VISCODE"), required = FALSE)
l_cr   <- pick_col(lab_raw, c("RCT392"))

lab <- lab_raw %>%
  transmute(RID = as.integer(RID),
            exam_date = parse_date_ymd(.data[[l_date]]),
            VISCODE2  = if (!is.na(l_vis)) recode_viscode(.data[[l_vis]]) else NA_character_,
            cr_raw    = .data[[l_cr]]) %>%
  filter(RID %in% rids, !is.na(exam_date)) %>%
  # Non-numeric codes (-1, SCC09, ADC259 and the like) mean not done
  mutate(creatinine = if_else(str_detect(coalesce(cr_raw, ""), "^\\s*[0-9]*\\.?[0-9]+\\s*$"),
                              suppressWarnings(as.numeric(cr_raw)), NA_real_),
         creatinine = if_else(!is.na(creatinine) & creatinine > 0 & creatinine < 20,
                              creatinine, NA_real_)) %>%
  filter(!is.na(creatinine)) %>%
  distinct(RID, exam_date, .keep_all = TRUE) %>%
  left_join(person_age, by = "RID") %>%
  mutate(age_at_exam = AGE + as.numeric(exam_date - adni_bl_date) / YEAR_DAYS,
         egfr_2021 = egfr_ckdepi_2021(creatinine, age_at_exam, PTGENDER),
         egfr_2009 = egfr_ckdepi_2009(creatinine, age_at_exam, PTGENDER))
say("LABDATA: %d creatinine values for %d people; %d rows dropped as non-numeric or implausible",
    nrow(lab), n_distinct(lab$RID),
    sum(lab_raw$RID %in% as.character(rids)) - nrow(lab))
say("Creatinine (mg/dL): median %.2f, range %.2f to %.2f; eGFR 2021 median %.0f",
    median(lab$creatinine), min(lab$creatinine), max(lab$creatinine), median(lab$egfr_2021, na.rm = TRUE))

lab_long <- lab %>%
  select(RID, exam_date, VISCODE2, creatinine, egfr_2021, egfr_2009, age_at_exam) %>%
  pivot_longer(c(creatinine, egfr_2021, egfr_2009), names_to = "measure", values_to = "value") %>%
  mutate(source = "lab", derived = measure != "creatinine")

# ---- 2. Weight, height, BMI from VITALS (optional) -------------------------------
vit_long <- NULL
if (file.exists(RAW_FILES$vitals)) {
  say("Reading %s", RAW_FILES$vitals)
  vit_raw <- read_adni(RAW_FILES$vitals)
  v_date <- pick_col(vit_raw, c("EXAMDATE", "VISDATE"))
  v_vis  <- pick_col(vit_raw, c("VISCODE2", "VISCODE"), required = FALSE)
  v_wt   <- pick_col(vit_raw, c("VSWEIGHT"))
  v_wtu  <- pick_col(vit_raw, c("VSWTUNIT"), required = FALSE)
  v_ht   <- pick_col(vit_raw, c("VSHEIGHT"), required = FALSE)
  v_htu  <- pick_col(vit_raw, c("VSHTUNIT"), required = FALSE)
  num <- function(x) { x <- suppressWarnings(as.numeric(x)); if_else(!is.na(x) & x > 0, x, NA_real_) }

  vit <- vit_raw %>%
    transmute(RID = as.integer(RID),
              exam_date = parse_date_ymd(.data[[v_date]]),
              VISCODE2  = if (!is.na(v_vis)) recode_viscode(.data[[v_vis]]) else NA_character_,
              wt = num(.data[[v_wt]]),
              wtu = if (!is.na(v_wtu)) .data[[v_wtu]] else NA_character_,
              ht = if (!is.na(v_ht)) num(.data[[v_ht]]) else NA_real_,
              htu = if (!is.na(v_htu)) .data[[v_htu]] else NA_character_) %>%
    filter(RID %in% rids, !is.na(exam_date)) %>%
    # ADNI unit codes: weight 1 = pounds, 2 = kilograms; height 1 = inches, 2 = cm.
    # If the unit column is absent, guess from the magnitude.
    mutate(weight_kg = case_when(wtu == "1" ~ wt * 0.45359237,
                                 wtu == "2" ~ wt,
                                 is.na(wtu) & wt > 130 ~ wt * 0.45359237,
                                 TRUE ~ wt),
           height_cm = case_when(htu == "1" ~ ht * 2.54,
                                 htu == "2" ~ ht,
                                 is.na(htu) & ht < 100 ~ ht * 2.54,
                                 TRUE ~ ht),
           weight_kg = if_else(weight_kg >= 30 & weight_kg <= 250, weight_kg, NA_real_),
           height_cm = if_else(height_cm >= 120 & height_cm <= 220, height_cm, NA_real_)) %>%
    filter(!is.na(weight_kg) | !is.na(height_cm)) %>%
    distinct(RID, exam_date, .keep_all = TRUE) %>%
    group_by(RID) %>%
    mutate(height_person = median(height_cm, na.rm = TRUE),
           height_used   = coalesce(height_cm, height_person),
           height_source = case_when(!is.na(height_cm) ~ "visit",
                                     !is.na(height_person) ~ "person median",
                                     TRUE ~ NA_character_)) %>%
    ungroup() %>%
    mutate(bmi = weight_kg / (height_used / 100)^2)
  say("VITALS: %d rows for %d people; weight at %d, height measured at %d, BMI computable at %d",
      nrow(vit), n_distinct(vit$RID), sum(!is.na(vit$weight_kg)), sum(!is.na(vit$height_cm)),
      sum(!is.na(vit$bmi)))
  vit_long <- vit %>%
    select(RID, exam_date, VISCODE2, weight_kg, height_cm, bmi, height_source) %>%
    pivot_longer(c(weight_kg, height_cm, bmi), names_to = "measure", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(source = "vitals", derived = measure == "bmi")
} else {
  say("VITALS.csv not found at %s; BMI skipped. Download VITALS.csv from ADNI to add it.",
      RAW_FILES$vitals)
}

# ---- 3. Combine and match to plasma visits --------------------------------------
mods <- bind_rows(lab_long, vit_long) %>%
  match_to_plasma("exam_date", plasma_visits) %>%
  arrange(RID, source, measure, exam_date) %>%
  group_by(RID, measure) %>%
  mutate(assessment_no = row_number(), n_assessments = n()) %>%
  ungroup() %>%
  select(RID, source, measure, exam_date, VISCODE2, assessment_no, n_assessments,
         value, derived, age_at_exam, height_source,
         visit_no, visit_date, gap_days, plasma_time)

save_rds(mods, "modifiers_long")
write_csv(mods, file.path(OUT_DIR, "modifiers_long.csv"))
say("modifiers_long: %d rows, %d people, %d measures", nrow(mods), n_distinct(mods$RID), n_distinct(mods$measure))

# ---- 4. Coverage -----------------------------------------------------------------
cov <- mods %>%
  group_by(source, measure) %>%
  summarise(n_rows = n(), n_people = n_distinct(RID),
            n_people_ge2 = n_distinct(RID[n_assessments >= 2]),
            n_people_ge3 = n_distinct(RID[n_assessments >= 3]),
            n_within_180d = sum(abs(gap_days) <= 180),
            n_within_365d = sum(abs(gap_days) <= 365),
            people_within_365d = n_distinct(RID[abs(gap_days) <= 365]),
            people_before_or_at_first_plasma = n_distinct(RID[plasma_time <= 0.5]),
            median_abs_gap = median(abs(gap_days)),
            median_value = median(value), p05 = quantile(value, .05), p95 = quantile(value, .95),
            .groups = "drop")
print(as.data.frame(cov))
write_csv(cov, file.path(CHECK_DIR, "modifiers_coverage.csv"))

near <- mods %>%
  filter(abs(gap_days) <= 365) %>%
  distinct(RID, visit_no, measure) %>%
  count(measure, visit_no, name = "n") %>%
  pivot_wider(names_from = visit_no, values_from = n, names_prefix = "visit_", values_fill = 0)
print(as.data.frame(near))
write_csv(near, file.path(CHECK_DIR, "modifiers_match_by_plasma_visit.csv"))


# =============================================================================================
# Section 5. Checks and coverage tables on the plasma table
# Original script: 99_checks.R
# =============================================================================================
# =============================================================================
# 99_checks.R
# Assertions and coverage tables for plasma_long. Run after 01_plasma_long.R.
# Prints a summary and writes every table to pipeline_out/checks/.
#
# The tables that decide the plan:
#   coverage_by_measure.csv   how many people have 1, 2, 3+ visits per measure
#   lod_by_measure_visit.csv  how often each measure sits at the detection floor
#   units_by_measure.csv      what unit each measure was reported in
# =============================================================================


# Always reload the shared config when it is in the working directory, so an
# old copy left in the R session cannot be used by mistake. If it is not
# there, fall back to whatever is already loaded, or stop.
Sys.setenv(COHORT = "adni"); source("00_setup.R")

pl <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
say("Loaded plasma_long: %d rows, %d people, pipeline %s",
    nrow(pl), n_distinct(pl$RID), attr(pl, "pipeline_version"))

write_check <- function(df, name) {
  write_csv(df, file.path(CHECK_DIR, paste0(name, ".csv")))
  invisible(df)
}

# ---- 1. Basic shape and unmapped measures -----------------------------------
say("\n== 1. Shape ==")
say("Date range: %s to %s",
    format(min(pl$visit_date, na.rm = TRUE)), format(max(pl$visit_date, na.rm = TRUE)))
say("Rows with unparsed EXAMDATE: %d", sum(is.na(pl$EXAMDATE_row)))
say("Rows where TESTVALUE was present but did not parse: %d",
    sum(is.na(pl$value_raw) & !is.na(pl$TESTVALUE_raw) & !pl$derived))

unmapped <- pl %>% filter(role == "unmapped") %>% count(PLASMA_BIOMARKER)
say("Unmapped PLASMA_BIOMARKER values: %d", nrow(unmapped))
if (nrow(unmapped) > 0) print(unmapped)
write_check(unmapped, "unmapped_measures")

# ---- 2. Uniqueness -----------------------------------------------------------
say("\n== 2. Uniqueness and visit integrity ==")
dups <- pl %>%
  count(RID, visit_no, PLASMA_BIOMARKER) %>%
  filter(n > 1)
say("Duplicate RID x visit x measure rows: %d", nrow(dups))
if (nrow(dups) > 0) print(as.data.frame(head(dups, 20)))
write_check(dups, "duplicate_rows")

# Split visits: one visit whose rows carry more than one EXAMDATE across
# platforms (dates within VISIT_GAP_DAYS were merged into the visit). Each
# row here is a lab-to-lab date disagreement that a date-keyed index would
# have split into two visits.
split_visits <- pl %>%
  filter(date_span_days > 0) %>%
  group_by(RID, visit_no, visit_date, date_span_days) %>%
  summarise(
    dates     = paste(sort(unique(format(EXAMDATE_row))), collapse = " | "),
    platforms = paste(sort(unique(platform)), collapse = ", "),
    viscodes  = paste(sort(unique(coalesce(VISCODE2_row, "NA"))), collapse = ", "),
    prefixes  = paste(sort(unique(coalesce(sample_prefix, "NA"))), collapse = ", "),
    .groups = "drop"
  )
say("Visits with more than one reported date (merged within %d days): %d",
    VISIT_GAP_DAYS, nrow(split_visits))
if (nrow(split_visits) > 0) print(as.data.frame(head(split_visits, 20)))
write_check(split_visits, "split_visits")

# Gap between consecutive visits per person. Anything under ~150 days is
# suspicious for ADNI and would suggest VISIT_GAP_DAYS is too small.
visit_gaps <- pl %>%
  distinct(RID, visit_no, visit_date) %>%
  arrange(RID, visit_no) %>%
  group_by(RID) %>%
  mutate(gap_days = as.integer(visit_date - lag(visit_date))) %>%
  ungroup() %>%
  filter(!is.na(gap_days))
say("Gap between consecutive visits (days): min %d, p05 %d, median %d, max %d",
    min(visit_gaps$gap_days), as.integer(quantile(visit_gaps$gap_days, .05)),
    as.integer(median(visit_gaps$gap_days)), max(visit_gaps$gap_days))
say("Visit gaps under 150 days: %d", sum(visit_gaps$gap_days < 150))
write_check(visit_gaps, "visit_gaps")

# Visits per person, with distinct dates and distinct SAMPLEID prefixes
# alongside, so the sample-ID semantics can be read off rather than assumed.
visits_per_person <- pl %>%
  group_by(RID) %>%
  summarise(n_visits   = n_distinct(visit_no),
            n_dates    = n_distinct(EXAMDATE_row),
            n_prefixes = n_distinct(sample_prefix, na.rm = TRUE),
            .groups = "drop")
say("Visits per person:")
print(table(visits_per_person$n_visits))
say("People with more distinct dates than visits (merged split dates): %d",
    sum(visits_per_person$n_dates > visits_per_person$n_visits))
say("People where SAMPLEID prefixes = visits: %d; prefixes = 1: %d; other: %d",
    sum(visits_per_person$n_prefixes == visits_per_person$n_visits),
    sum(visits_per_person$n_prefixes == 1 & visits_per_person$n_visits > 1),
    sum(visits_per_person$n_prefixes != visits_per_person$n_visits &
        !(visits_per_person$n_prefixes == 1 & visits_per_person$n_visits > 1)))
write_check(visits_per_person, "visits_per_person")

# VISCODE2 consistency: one code per visit, one visit per code
vis_multi <- pl %>%
  filter(!is.na(VISCODE2_row)) %>%
  distinct(RID, visit_no, VISCODE2_row) %>%
  count(RID, visit_no) %>% filter(n > 1)
code_multi <- pl %>%
  filter(!is.na(VISCODE2)) %>%
  distinct(RID, VISCODE2, visit_no) %>%
  count(RID, VISCODE2) %>% filter(n > 1)
say("Visits carrying more than one VISCODE2: %d", nrow(vis_multi))
say("RID x VISCODE2 spanning more than one visit: %d", nrow(code_multi))
say("Visits with no VISCODE2 from any platform: %d",
    pl %>% distinct(RID, visit_no, VISCODE2) %>% filter(is.na(VISCODE2)) %>% nrow())
write_check(vis_multi, "visit_multiple_viscodes")
write_check(code_multi, "viscode_multiple_visits")

# ---- 3. Coverage by measure: the scoping question ---------------------------
say("\n== 3. Coverage by measure ==")
coverage <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  distinct(RID, PLASMA_BIOMARKER, biomarker, platform, visit_no) %>%
  count(biomarker, platform, PLASMA_BIOMARKER, RID, name = "n_visits") %>%
  group_by(biomarker, platform, PLASMA_BIOMARKER) %>%
  summarise(
    n_people      = n(),
    n_ge2_visits  = sum(n_visits >= 2),
    n_ge3_visits  = sum(n_visits >= 3),
    n_ge4_visits  = sum(n_visits >= 4),
    median_visits = median(n_visits),
    .groups = "drop"
  ) %>%
  arrange(biomarker, desc(n_ge2_visits))
print(as.data.frame(coverage))
write_check(coverage, "coverage_by_measure")

# Rows by visit number per measure
rows_by_visit <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  count(biomarker, platform, PLASMA_BIOMARKER, visit_no) %>%
  pivot_wider(names_from = visit_no, values_from = n, names_prefix = "visit_",
              values_fill = 0)
write_check(rows_by_visit, "rows_by_measure_visit")

# ---- 4. Detection floor and imputation --------------------------------------
say("\n== 4. Detection floor ==")
lod <- pl %>%
  filter(role == "measure") %>%
  group_by(biomarker, platform, PLASMA_BIOMARKER, visit_no) %>%
  summarise(
    n         = n(),
    n_lloq    = sum(lloq_flag, na.rm = TRUE),
    pct_lloq  = round(100 * mean(lloq_flag, na.rm = TRUE), 1),
    n_imputed = sum(imputed_flag, na.rm = TRUE),
    n_uloq    = sum(uloq_flag, na.rm = TRUE),
    .groups = "drop"
  )
lod_summary <- lod %>%
  group_by(biomarker, platform, PLASMA_BIOMARKER) %>%
  summarise(n = sum(n), n_lloq = sum(n_lloq),
            pct_lloq = round(100 * n_lloq / n, 1),
            n_imputed = sum(n_imputed), .groups = "drop") %>%
  arrange(desc(pct_lloq))
print(as.data.frame(lod_summary))
write_check(lod, "lod_by_measure_visit")
write_check(lod_summary, "lod_by_measure")

# Distinct floor values, to confirm the substituted constant per measure
floor_values <- pl %>%
  filter(role == "measure", lloq_flag) %>%
  count(PLASMA_BIOMARKER, value_raw, name = "n_rows") %>%
  arrange(PLASMA_BIOMARKER, desc(n_rows))
write_check(floor_values, "lod_floor_values")

# People by floor exposure on each measure: never, sometimes, always
floor_people <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  group_by(PLASMA_BIOMARKER, RID) %>%
  summarise(n_obs = n(), n_floor = sum(lloq_flag), .groups = "drop") %>%
  mutate(floor_group = case_when(
    n_floor == 0      ~ "never",
    n_floor == n_obs  ~ "always",
    TRUE              ~ "sometimes"
  )) %>%
  group_by(PLASMA_BIOMARKER) %>%
  summarise(
    n_people      = n(),
    n_never       = sum(floor_group == "never"),
    n_sometimes   = sum(floor_group == "sometimes"),
    n_always      = sum(floor_group == "always"),
    pct_any_floor = round(100 * mean(floor_group != "never"), 1),
    pct_always    = round(100 * mean(floor_group == "always"), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(pct_any_floor))
print(as.data.frame(floor_people))
write_check(floor_people, "people_by_floor_exposure")

# Per-person floor pattern for C2N pTau217, for joining to other analyses
c2n_floor_by_person <- pl %>%
  filter(PLASMA_BIOMARKER == "C2N_plasma_ptau217", !is.na(value)) %>%
  arrange(RID, visit_no) %>%
  group_by(RID) %>%
  summarise(n_obs = n(), n_floor = sum(lloq_flag),
            pattern = paste(as.integer(lloq_flag), collapse = ""),
            .groups = "drop") %>%
  mutate(floor_group = case_when(n_floor == 0 ~ "never",
                                 n_floor == n_obs ~ "always",
                                 TRUE ~ "sometimes"))
say("C2N pTau217 floor patterns (1 = at floor, ordered by visit):")
print(c2n_floor_by_person %>% count(pattern, sort = TRUE) %>% head(10))
write_check(c2n_floor_by_person, "c2n_ptau217_floor_by_person")

# ---- 5. Units and value ranges -----------------------------------------------
say("\n== 5. Units ==")
units <- pl %>%
  count(biomarker, platform, PLASMA_BIOMARKER, role, unit_raw, unit_std) %>%
  arrange(biomarker, platform)
print(as.data.frame(units))
write_check(units, "units_by_measure")

ranges <- pl %>%
  filter(!is.na(value)) %>%
  group_by(biomarker, platform, PLASMA_BIOMARKER, role, unit_std) %>%
  summarise(n = n(),
            min = min(value), p05 = quantile(value, .05), median = median(value),
            p95 = quantile(value, .95), max = max(value), .groups = "drop") %>%
  mutate(across(c(min, p05, median, p95, max), ~ signif(.x, 4)))
print(as.data.frame(ranges))
write_check(ranges, "value_ranges_by_measure")

# C2N supplies its own Ab42/40. Recompute from components and compare.
c2n_ratio_check <- pl %>%
  filter(platform == "C2N", biomarker == "Ab42_40") %>%
  select(RID, visit_no, PLASMA_BIOMARKER, value) %>%
  pivot_wider(names_from = PLASMA_BIOMARKER, values_from = value,
              values_fn = first) %>%
  mutate(recomputed = C2N_plasma_Abeta42 / C2N_plasma_Abeta40,
         diff = C2N_plasma_Abeta42_Abeta40 - recomputed) %>%
  filter(!is.na(diff))
say("C2N Ab42/40 supplied vs recomputed: n=%d, max |diff|=%s",
    nrow(c2n_ratio_check), format(max(abs(c2n_ratio_check$diff)), digits = 3))
write_check(c2n_ratio_check, "c2n_ratio_recompute")

# ---- 6. Batch, CV, sample identity --------------------------------------------
say("\n== 6. Batch and analytical precision ==")
run_by_platform <- pl %>%
  filter(role == "measure") %>%
  distinct(platform, RID, visit_no, RUN_NO, RUNDATE) %>%
  count(platform, RUN_NO, RUNDATE) %>%
  arrange(platform, RUNDATE, RUN_NO)
say("Distinct runs per platform:")
print(run_by_platform %>% count(platform, name = "n_runs"))
write_check(run_by_platform, "runs_by_platform")

# Does batch vary within a visit number? If yes, batch is separable from occasion.
run_by_visit <- pl %>%
  filter(role == "measure", !is.na(RUN_NO)) %>%
  distinct(platform, visit_no, RUN_NO) %>%
  count(platform, visit_no, name = "n_runs_at_visit")
write_check(run_by_visit, "runs_by_platform_visit")

cv_pop <- pl %>%
  group_by(platform) %>%
  summarise(n = n(), n_cv = sum(!is.na(CV)), pct_cv = round(100 * mean(!is.na(CV)), 1),
            .groups = "drop")
print(as.data.frame(cv_pop))
write_check(cv_pop, "cv_populated_by_platform")

# Sample identity: platforms per visit, and whether every platform at a visit
# reports the same SAMPLEID prefix (same aliquot) or not.
platforms_per_visit <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  group_by(RID, visit_no) %>%
  summarise(n_platforms = n_distinct(platform),
            n_prefixes  = n_distinct(sample_prefix, na.rm = TRUE),
            n_sampleids = n_distinct(SAMPLEID, na.rm = TRUE),
            .groups = "drop")
say("Platforms per visit:")
print(table(platforms_per_visit$n_platforms))
say("Visits with more than one SAMPLEID prefix across platforms: %d",
    sum(platforms_per_visit$n_prefixes > 1))
say("Rows with no SAMPLEID: %d", sum(is.na(pl$SAMPLEID) | pl$SAMPLEID == ""))
write_check(platforms_per_visit, "platforms_per_visit")

# Visits where only one platform has values: which platform, and are these
# extra visits (4th, 5th) or people who otherwise have no measures at all?
single_platform_visits <- pl %>%
  filter(role == "measure", !is.na(value)) %>%
  group_by(RID, visit_no, visit_date, VISCODE2) %>%
  filter(n_distinct(platform) == 1) %>%
  summarise(platform = first(platform),
            measures = paste(sort(unique(PLASMA_BIOMARKER)), collapse = ", "),
            n_rows   = n(), .groups = "drop") %>%
  left_join(visits_per_person %>% select(RID, n_visits), by = "RID")
say("Visits with values from one platform only: %d", nrow(single_platform_visits))
if (nrow(single_platform_visits) > 0) {
  print(single_platform_visits %>% count(platform, visit_no, n_visits))
  print(as.data.frame(head(single_platform_visits, 10)))
}
write_check(single_platform_visits, "single_platform_visits")

# People with no measure values at all: what rows do they have?
no_measure_people <- pl %>%
  group_by(RID) %>%
  filter(!any(role == "measure" & !is.na(value))) %>%
  count(RID, platform, role, has_value = !is.na(value), name = "n_rows") %>%
  ungroup()
say("People with no measure values: %d", n_distinct(no_measure_people$RID))
if (nrow(no_measure_people) > 0) print(as.data.frame(head(no_measure_people, 20)))
write_check(no_measure_people, "no_measure_people")

# Rows per person per platform, to see whether the 4th and 5th visits carry
# any values on any platform.
rows_by_person_visit <- pl %>%
  filter(role == "measure") %>%
  group_by(RID, visit_no, visit_date, VISCODE2) %>%
  summarise(n_rows = n(), n_values = sum(!is.na(value)),
            platforms = paste(sort(unique(platform[!is.na(value)])), collapse = ", "),
            .groups = "drop")
say("Rows and values at visit 4+:")
print(as.data.frame(rows_by_person_visit %>% filter(visit_no >= 4)))
write_check(rows_by_person_visit, "rows_by_person_visit")

# ---- 7. Missingness in the wide frame ------------------------------------------
say("\n== 7. Wide-frame missingness ==")
pw <- readRDS(file.path(OUT_DIR, "plasma_wide.rds"))
miss <- tibble(column = names(pw),
               n_missing = sapply(pw, function(x) sum(is.na(x))),
               pct_missing = round(100 * sapply(pw, function(x) mean(is.na(x))), 1)) %>%
  arrange(desc(pct_missing))
print(as.data.frame(miss))
write_check(miss, "wide_missingness")

say("\nAll check tables written to %s", CHECK_DIR)
