# =============================================================================================
# 02_data_a4learn.R
# A4/LEARN data build in the same layout as ADNI.
#
# Reads the A4/LEARN release (read only) and writes a4learn_out/: plasma_long, participants, criteria_long,
# modifiers_long. Four sections, formerly scripts 09a to 09d. Run by run_all.R with COHORT = a4learn.
# =============================================================================================


# =============================================================================================
# Section 1. Plasma pTau217 series (Lilly MSD)
# Original script: 09a_a4_plasma_long.R
# =============================================================================================
# =============================================================================
# 09a_a4_plasma_long.R
# A4/LEARN: the plasma pTau217 series in the tall format the ADNI pipeline uses.
#
# Source   External Data/biomarker_pTau217.csv (Lilly MSD electrochemiluminescence
#          pTau217, U/mL), Derived Data/SUBJINFO.csv (arm, treatment dates),
#          Documents/Data Dictionaries/visits_datadic.csv (visit names).
#
# One row per draw. Columns mirror ADNI's plasma_long so that 06_master.R and the
# 07 frame scripts run unchanged with COHORT = "a4learn". Extra columns:
#   substudy, arm            A4 or LEARN; Placebo, Solanezumab or LEARN
#   visit_week, visit_kind   scheduled week; screening / baseline / follow-up /
#                            early_termination / open_label
#   on_treatment             draw taken between treatment start and end (solanezumab arm)
#   short_interval_pair      A4 baseline and week-12 draws (about 12 weeks apart):
#                            a short-interval test-retest pair
#   days_t0                  days since T0 (A4 randomisation; LEARN baseline exam)
#
# Values: the laboratory reports "<LLOQ" and ">ULOQ" as text and keeps the raw
# instrument reading in ORRESRAW. value = the raw reading (as delivered by the
# instrument) with lloq_flag / uloq_flag set; value_raw is the same; the text is
# kept in TESTVALUE_raw. Rows marked not done (haemolysed, insufficient volume,
# freeze-thaw limit) are dropped and counted.
#
# Dates: A4/LEARN release days since T0; the pipeline's Date columns are those
# days on a nominal calendar (A4_EPOCH + days). Only differences are used.
#
# Output   a4learn_out/plasma_long.rds, plasma_wide.csv
#          checks/a4_plasma_flow.csv, a4_draws_by_visit.csv, a4_draws_per_person.csv,
#          a4_floor_by_substudy.csv, a4_design_sst.csv
# =============================================================================

Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
if (COHORT != "a4learn") stop("Set Sys.setenv(COHORT = 'a4learn') before running the 09 scripts.")

raw  <- read_adni(A4_FILES$ptau217)
subj <- read_adni(A4_FILES$subjinfo)
vis  <- read_adni(A4_FILES$visits)

# ---- 1. Rows and exclusions --------------------------------------------------------
flow <- tibble(step = "rows in file", n_rows = nrow(raw), n_people = n_distinct(raw$BID))
r <- raw %>%
  filter(SUBSTUDY %in% c("A4", "LEARN")) %>%
  mutate(not_done = !is.na(STAT) & STAT != "", reason = REASND)
flow <- bind_rows(flow, tibble(step = "A4 or LEARN", n_rows = nrow(r), n_people = n_distinct(r$BID)))
not_done <- r %>% filter(not_done) %>% count(SUBSTUDY, reason, name = "n")
r <- r %>% filter(!not_done)
flow <- bind_rows(flow, tibble(step = "result present", n_rows = nrow(r), n_people = n_distinct(r$BID)))

# ---- 2. Values and flags -----------------------------------------------------------
pc <- parse_censored(r$ORRES)
r <- r %>%
  mutate(
    value_reported = pc$value,
    comparator     = case_when(str_detect(ORRES, "^<") ~ "<", str_detect(ORRES, "^>") ~ ">", TRUE ~ NA_character_),
    lloq_flag      = !is.na(comparator) & comparator == "<",
    uloq_flag      = !is.na(comparator) & comparator == ">",
    value_raw      = suppressWarnings(as.numeric(ORRESRAW)),
    value          = coalesce(if_else(lloq_flag | uloq_flag, value_raw, value_reported), value_raw),
    days_t0        = as.numeric(COLLECTION_DATE_DAYS_T0),
    VISCODE        = as.integer(VISCODE)
  )
say("Below LLOQ %d, above ULOQ %d, values missing after parsing %d",
    sum(r$lloq_flag), sum(r$uloq_flag), sum(is.na(r$value)))
r <- r %>% filter(!is.na(value), !is.na(days_t0))
flow <- bind_rows(flow, tibble(step = "numeric value and date", n_rows = nrow(r), n_people = n_distinct(r$BID)))

# ---- 3. Visit labels ---------------------------------------------------------------
vmap <- vis %>%
  transmute(SUBSTUDY, VISCODE = as.integer(VISCODE), VISNAME, ATTRIBS) %>%
  mutate(visit_week = suppressWarnings(as.integer(str_match(VISNAME, "wk(\\d+)")[, 2])),
         visit_kind = case_when(
           str_detect(ATTRIBS, "^sc") ~ "screening",
           str_detect(ATTRIBS, "^bl") ~ "baseline",
           ATTRIBS %in% c("et1") ~ "early_termination",
           ATTRIBS %in% c("et2") ~ "open_label_final",
           str_detect(ATTRIBS, "ole") ~ "open_label",
           str_detect(ATTRIBS, "^w") ~ "follow_up",
           TRUE ~ "other"),
         visit_week = if_else(visit_kind == "baseline", 0L, visit_week)) %>%
  distinct(SUBSTUDY, VISCODE, .keep_all = TRUE)
r <- r %>% left_join(vmap %>% select(SUBSTUDY, VISCODE, VISNAME, visit_week, visit_kind), by = c("SUBSTUDY", "VISCODE"))

# ---- 4. Arm and treatment exposure -------------------------------------------------
arm <- subj %>%
  transmute(BID, arm = case_when(SUBSTUDY == "LEARN" ~ "LEARN", TX == "Solanezumab" ~ "Solanezumab",
                                 TX == "Placebo" ~ "Placebo", TRUE ~ NA_character_),
            trt_start = as.numeric(TRTSTDT_DAYS_T0), trt_end = as.numeric(TRTENDT_DAYS_T0))
r <- r %>% left_join(arm, by = "BID") %>%
  mutate(on_treatment = arm == "Solanezumab" & !is.na(trt_start) & days_t0 >= trt_start &
                        (is.na(trt_end) | days_t0 <= trt_end))

# ---- 5. Visits: order draws, collapse near-duplicates, time since first draw --------
r <- r %>% arrange(BID, days_t0) %>%
  group_by(BID) %>%
  mutate(visit_no = cumsum(c(TRUE, diff(days_t0) > VISIT_GAP_DAYS))) %>%
  ungroup()
dups <- r %>% count(BID, visit_no) %>% filter(n > 1)
if (nrow(dups)) {
  say("%d person-visits have more than one draw within %d days; keeping the first", nrow(dups), VISIT_GAP_DAYS)
  r <- r %>% group_by(BID, visit_no) %>% slice_min(days_t0, n = 1, with_ties = FALSE) %>% ungroup()
}
r <- r %>% group_by(BID) %>%
  mutate(first_days = min(days_t0), plasma_time = (days_t0 - first_days) / YEAR_DAYS,
         n_draws = n(),
         short_interval_pair = SUBSTUDY == "A4" & visit_kind %in% c("baseline", "follow_up") &
                               VISCODE %in% c(6L, 9L) & sum(VISCODE %in% c(6L, 9L)) == 2) %>%
  ungroup()

# ---- 6. Assemble in the ADNI plasma_long layout -------------------------------------
plasma_long <- r %>%
  transmute(
    RID = a4_rid(BID), PTID = BID,
    visit_no, occasion_no = visit_no,
    visit_date = a4_date(days_t0), plasma_time,
    VISCODE2 = paste0("v", VISCODE),
    n_platforms_at_visit = 1L, quanterix_only_visit = FALSE, quanterix_only_person = FALSE,
    SAMPLEID = SAMPLE_ID, sample_prefix = SAMPLE_ID,
    EXAMDATE_row = a4_date(days_t0), VISCODE2_row = paste0("v", VISCODE), date_span_days = 0L,
    biomarker = "pTau217", platform = "Lilly MSD", PLASMA_BIOMARKER = "Lilly_plasma_ptau217",
    TESTNAME = "ptau217", role = "measure", component = NA_character_, direction = 1,
    value, value_raw, TESTVALUE_raw = ORRES, unit_std = "U/mL", unit_raw = ORRESU,
    lloq_flag, uloq_flag, imputed_flag = FALSE, comparator, derived = FALSE,
    ASSAYPLATFORM = "MSD", ASSAYVERSION = "Lilly pTau217 ECL", PERFORMINGLAB = NAM,
    INSTRUMENTID = "MSD Sector S 600MM", MATRIX = SPEC, CV = NA_real_,
    RUN_NO = NA_character_, SEQUENCE_NO = NA_character_, RUNDATE = as.Date(NA), COMMENTS = COMMENT2,
    # A4/LEARN design columns
    substudy = SUBSTUDY, arm, days_t0, visit_code = VISCODE, visit_name = VISNAME, visit_week, visit_kind,
    on_treatment, short_interval_pair
  ) %>% arrange(RID, visit_no)
stopifnot(!any(duplicated(plasma_long$RID[!duplicated(plasma_long$PTID)])))   # BID -> RID one to one

save_rds(plasma_long, "plasma_long")
wide <- plasma_long %>% select(RID, PTID, substudy, arm, visit_no, visit_code, days_t0, plasma_time, value, lloq_flag, uloq_flag)
write_csv(wide, file.path(OUT_DIR, "plasma_wide.csv"))

# ---- 7. Checks -----------------------------------------------------------------------
flow <- bind_rows(flow, tibble(step = "draws kept", n_rows = nrow(plasma_long), n_people = n_distinct(plasma_long$RID)))
write_csv(flow, file.path(CHECK_DIR, "a4_plasma_flow.csv"))
write_csv(not_done, file.path(CHECK_DIR, "a4_plasma_not_done.csv"))
by_visit <- plasma_long %>% group_by(substudy, arm, visit_code, visit_name, visit_kind) %>%
  summarise(n = n(), days_median = median(days_t0), days_min = min(days_t0), days_max = max(days_t0),
            pct_floor = round(100 * mean(lloq_flag), 1), .groups = "drop") %>% arrange(substudy, arm, visit_code)
write_csv(by_visit, file.path(CHECK_DIR, "a4_draws_by_visit.csv"))
per_person <- plasma_long %>% group_by(RID, substudy, arm) %>%
  summarise(k = n(), span = max(plasma_time), sst = sum((plasma_time - mean(plasma_time))^2),
            pair = any(short_interval_pair), .groups = "drop")
write_csv(per_person %>% count(substudy, arm, k), file.path(CHECK_DIR, "a4_draws_per_person.csv"))
design <- per_person %>% filter(k >= 2) %>% group_by(substudy, arm) %>%
  summarise(people = n(), with_3_plus = sum(k >= 3), with_pair = sum(pair),
            span_median = round(median(span), 2), sst_median = round(median(sst), 2),
            sst_q25 = round(quantile(sst, .25), 2), sst_q75 = round(quantile(sst, .75), 2), .groups = "drop")
write_csv(design, file.path(CHECK_DIR, "a4_design_sst.csv"))
floor <- plasma_long %>% group_by(substudy, arm) %>%
  summarise(draws = n(), pct_lloq = round(100 * mean(lloq_flag), 1), pct_uloq = round(100 * mean(uloq_flag), 2),
            people_always_floor = n_distinct(RID[ave(lloq_flag, RID, FUN = all) == 1]), .groups = "drop")
write_csv(floor, file.path(CHECK_DIR, "a4_floor_by_substudy.csv"))

say("plasma_long: %d draws, %d people (A4 placebo %d, solanezumab %d, LEARN %d)",
    nrow(plasma_long), n_distinct(plasma_long$RID),
    n_distinct(plasma_long$RID[plasma_long$arm == "Placebo"]), n_distinct(plasma_long$RID[plasma_long$arm == "Solanezumab"]),
    n_distinct(plasma_long$RID[plasma_long$arm == "LEARN"]))
print(as.data.frame(by_visit %>% select(-visit_name)))
print(as.data.frame(design))
print(as.data.frame(floor))


# =============================================================================================
# Section 2. Participants
# Original script: 09b_a4_participants.R
# =============================================================================================
# =============================================================================
# 09b_a4_participants.R
# A4/LEARN: one row per person, in the layout of ADNI's participants.rds.
#
# Source   Derived Data/SUBJINFO.csv (arm, age at consent, sex, education, race,
#          ethnicity, APOE genotype, baseline BMI, enrolment and treatment days,
#          completion), Raw Data/ptdemog.csv (screening age and sex as a check),
#          a4learn_out/plasma_long.rds (plasma series per person).
#
# Conventions carried over from ADNI: RID integer id; PTID the BID; PTGENDER
# "Female"/"Male"; PTEDUCAT years; APOE4 = number of e4 alleles; APOE4_carrier;
# DX_bl and DX_at_first_plasma "CN" (both substudies enrolled cognitively
# unimpaired people; CDR global is kept separately in the criteria); in_core_5lab
# TRUE for everyone with a plasma draw (the analysis scripts filter on it);
# in_quanterix FALSE; COLPROT/ORIGPROT hold substudy and arm.
#
# Age: SUBJINFO gives age at consent (AGEYR) and enrolment days since consent
# and since T0, so age at any day t (since T0) = AGEYR + (t - consent_day_t0)/365.25.
#
# Output   a4learn_out/participants.rds, participants.csv
#          checks/a4_participants_by_arm.csv, a4_participants_missing.csv, a4_age_sex_check.csv
# =============================================================================

Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
if (COHORT != "a4learn") stop("Set Sys.setenv(COHORT = 'a4learn') before running the 09 scripts.")

subj <- read_adni(A4_FILES$subjinfo)
dem  <- read_adni(A4_FILES$ptdemog)
pl   <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))

num <- function(x) suppressWarnings(as.numeric(x))

# ---- 1. Plasma series per person ----------------------------------------------------
series <- pl %>% group_by(RID, PTID) %>%
  summarise(n_plasma_visits = n_distinct(visit_no), first_plasma = min(visit_date), last_plasma = max(visit_date),
            plasma_span_years = as.numeric(max(visit_date) - min(visit_date)) / YEAR_DAYS,
            first_plasma_days = min(days_t0), substudy = first(substudy), arm = first(arm),
            has_short_interval_pair = any(short_interval_pair), .groups = "drop")

# ---- 2. SUBJINFO --------------------------------------------------------------------
s <- subj %>%
  filter(SUBSTUDY %in% c("A4", "LEARN")) %>%
  transmute(
    PTID = BID, substudy_info = SUBSTUDY,
    arm_info = case_when(SUBSTUDY == "LEARN" ~ "LEARN", TX == "Solanezumab" ~ "Solanezumab", TX == "Placebo" ~ "Placebo"),
    AGEYR = num(AGEYR),
    PTGENDER = case_when(SEX == "1" ~ "Female", SEX == "2" ~ "Male", TRUE ~ NA_character_),
    PTEDUCAT = num(EDCCNTU),
    race = case_when(RACE == "1" ~ "White", RACE == "2" ~ "Black or African American", RACE == "58" ~ "Asian",
                     RACE == "79" ~ "Native Hawaiian or Other Pacific Islander", RACE == "84" ~ "American Indian or Alaskan Native",
                     RACE == "100" ~ "More than one race", TRUE ~ "Unknown or not reported"),
    ethnicity = case_when(ETHNIC == "50" ~ "Hispanic or Latino", ETHNIC == "56" ~ "Not Hispanic or Latino", TRUE ~ "Unknown or not reported"),
    APOE_genotype = APOEGN,
    APOE4 = str_count(coalesce(APOEGN, ""), "E4"),
    APOE4_carrier = case_when(is.na(APOEGN) | APOEGN == "" ~ NA, TRUE ~ str_detect(APOEGN, "E4")),
    bmi_bl_info = num(BMIBL), height_cm_bl = num(HGTCMBL), weight_kg_bl = num(WGTKGBL),
    enrol_days_consent = num(ENRLDT_DAYS_CONSENT), enrol_days_t0 = num(ENRLDT_DAYS_T0),
    consent_day_t0 = enrol_days_t0 - enrol_days_consent,
    randomised_day_t0 = num(RNDMDT_DAYS_T0),
    treatment_start = num(TRTSTDT_DAYS_T0), treatment_end = num(TRTENDT_DAYS_T0),
    discontinued_day_t0 = num(DISCDTC_DAYS_T0),
    completed = SUBJCOMPTR == "1",
    gds_over_5_bl = GDSSEVBL == "1"
  )

# screening demographics: age at the screening exam and the exam day are present for both
# substudies (SUBJINFO's consent and enrolment days are A4 only), so they anchor age
d <- dem %>% transmute(PTID = BID, sex_demog = case_when(PTGENDER == "1" ~ "Male", PTGENDER == "2" ~ "Female"),
                       age_screen = num(PTAGE), screen_day_t0 = num(EXAMDATE_DAYS_T0), educ_demog = num(PTEDUCAT)) %>%
  filter(!is.na(age_screen)) %>% distinct(PTID, .keep_all = TRUE)

# ---- 3. Assemble --------------------------------------------------------------------
participants <- series %>%
  left_join(s, by = "PTID") %>%
  left_join(d, by = "PTID") %>%
  mutate(
    RID, PTID,
    age_first_plasma = coalesce(age_screen + (first_plasma_days - screen_day_t0) / YEAR_DAYS,
                                AGEYR + (first_plasma_days - consent_day_t0) / YEAR_DAYS),
    AGE = coalesce(age_screen, AGEYR),
    completed = if_else(substudy == "LEARN", NA, completed),
    DX_bl = "CN", DX_at_first_plasma = "CN", dx_gap_days = NA_integer_,
    COLPROT = substudy, ORIGPROT = arm,
    n_visits_5lab = n_plasma_visits, n_visits_quanterix = 0L, n_visits_all6 = 0L,
    quanterix_only_person = FALSE, in_core_5lab = n_plasma_visits >= 1, in_quanterix = FALSE,
    in_merge = !is.na(AGEYR),
    adni_bl_date = a4_date(0), adni_bl_viscode = "bl",
    PTGENDER_demog = sex_demog, sex_mismatch = !is.na(sex_demog) & !is.na(PTGENDER) & sex_demog != PTGENDER,
    age_diff_check = age_screen - (AGEYR + (screen_day_t0 - consent_day_t0) / YEAR_DAYS)   # screening age vs consent age (A4 only)
  ) %>%
  select(RID, PTID, substudy, arm, n_plasma_visits, first_plasma, last_plasma, plasma_span_years,
         n_visits_5lab, n_visits_quanterix, n_visits_all6, quanterix_only_person, in_core_5lab, in_quanterix,
         adni_bl_date, adni_bl_viscode, AGE, PTGENDER, PTEDUCAT, APOE4, APOE4_carrier, APOE_genotype,
         DX_bl, COLPROT, ORIGPROT, PTGENDER_demog, age_first_plasma, sex_mismatch, age_diff_check,
         in_merge, DX_at_first_plasma, dx_gap_days,
         race, ethnicity, bmi_bl_info, height_cm_bl, weight_kg_bl, age_screen, screen_day_t0, consent_day_t0, randomised_day_t0,
         treatment_start, treatment_end, discontinued_day_t0, completed, gds_over_5_bl, has_short_interval_pair)

save_rds(participants, "participants")
write_csv(participants, file.path(OUT_DIR, "participants.csv"))

# ---- 4. Checks ------------------------------------------------------------------------
by_arm <- participants %>% group_by(substudy, arm) %>%
  summarise(n = n(), age_mean = round(mean(age_first_plasma, na.rm = TRUE), 1), pct_female = round(100 * mean(PTGENDER == "Female", na.rm = TRUE)),
            educ_mean = round(mean(PTEDUCAT, na.rm = TRUE), 1), pct_apoe4 = round(100 * mean(APOE4_carrier, na.rm = TRUE)),
            draws_mean = round(mean(n_plasma_visits), 2), span_median = round(median(plasma_span_years), 2),
            with_pair = sum(has_short_interval_pair), completed = sum(completed, na.rm = TRUE), .groups = "drop")
write_csv(by_arm, file.path(CHECK_DIR, "a4_participants_by_arm.csv"))
missing <- participants %>% summarise(across(c(age_first_plasma, PTGENDER, PTEDUCAT, APOE4_carrier, bmi_bl_info, treatment_start),
                                             ~ sum(is.na(.x)))) %>% pivot_longer(everything(), names_to = "column", values_to = "n_missing")
write_csv(missing, file.path(CHECK_DIR, "a4_participants_missing.csv"))
chk <- participants %>% summarise(n = n(), sex_mismatch = sum(sex_mismatch, na.rm = TRUE),
                                  age_diff_abs_median = round(median(abs(age_diff_check), na.rm = TRUE), 2),
                                  age_diff_abs_max = round(max(abs(age_diff_check), na.rm = TRUE), 2))
write_csv(chk, file.path(CHECK_DIR, "a4_age_sex_check.csv"))
say("participants: %d people", nrow(participants)); print(as.data.frame(by_arm)); print(as.data.frame(chk))


# =============================================================================================
# Section 3. Reference measures
# Original script: 09c_a4_criteria_long.R
# =============================================================================================
# =============================================================================
# 09c_a4_criteria_long.R
# A4/LEARN: reference measures (the "criteria") in the layout of ADNI's criteria_long.
#
#   centiloid, summary_suvr, amyloid_status   florbetapir composite SUVR (whole cerebellum
#                                             reference, Avid pipeline) and its Centiloid value
#                                             (Navitsky 2018: 183 x SUVR - 177); status = CL > 20
#   tau_pet                                   flortaucipir MUBADA mask SUVR (cerebellar crus reference)
#   csf_abeta42, csf_abeta40, csf_abeta42_40, csf_tau, csf_ptau181, csf_ptau_abeta42, csf_ptau_tau
#                                             from the central laboratory file (screening LP and week 240)
#   mmse, pacc, cdrsb, cdr_global             cognition and function (dates from the visit file where
#                                             the assessment file has none)
#
# Every assessment is matched to the person's nearest plasma draw (visit_no,
# gap_days, plasma_time), as in 03_criteria_long.R. Only people with a plasma
# draw are kept.
#
# Output   a4learn_out/criteria_long.rds, criteria_long.csv
#          checks/a4_criteria_coverage.csv, a4_criteria_match_by_plasma_visit.csv
# =============================================================================

Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
if (COHORT != "a4learn") stop("Set Sys.setenv(COHORT = 'a4learn') before running the 09 scripts.")

pl   <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
num  <- function(x) suppressWarnings(as.numeric(x))
ids  <- part %>% select(PTID, RID)
plasma_visits <- pl %>% distinct(RID, visit_no, visit_date) %>%
  group_by(RID) %>% mutate(first_plasma = min(visit_date)) %>% ungroup()
sv <- read_adni(A4_FILES$sv) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISITCD), sv_days = num(SVSTDTC_DAYS_T0)) %>%
  filter(!is.na(sv_days)) %>% distinct(PTID, VISCODE, .keep_all = TRUE)

# ---- 1. Amyloid PET -----------------------------------------------------------------
# Some scans (LEARN week 240 and early termination; about a tenth of tau PET) have no scan
# date in the imaging file; the visit file's date for the same visit is used instead.
sv_day <- function(ptid, viscode) sv$sv_days[match(paste(ptid, viscode), paste(sv$PTID, sv$VISCODE))]
amy <- read_adni(A4_FILES$amyloid) %>%
  filter(BRAIN_REGION == "Composite_Summary", SCAN_ANALYZED %in% c("Yes", "YES", "Y", NA)) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), VISCODE2 = paste0("v", VISCODE), days = num(SCAN_DATE_DAYS_T0),
            summary_suvr = num(SUVR_CER), tracer = LIGAND) %>%
  mutate(date_source = if_else(is.na(days), "visit file", "scan file"), days = coalesce(days, sv_day(PTID, VISCODE))) %>%
  filter(!is.na(summary_suvr), !is.na(days)) %>% select(-VISCODE) %>%
  mutate(centiloid = a4_centiloid(summary_suvr), amyloid_status = as.numeric(centiloid > CENTILOID_CUTS[["cl20"]])) %>%
  pivot_longer(c(summary_suvr, centiloid, amyloid_status), names_to = "measure", values_to = "value") %>%
  mutate(source = "florbetapir_avid", derived = measure != "summary_suvr")

# ---- 2. Tau PET (subset) ---------------------------------------------------------------
tau <- read_adni(A4_FILES$tau)
tau_region <- if ("MUBADA Mask" %in% tau$BRAIN_REGION) "MUBADA Mask" else NA_character_
tau <- if (is.na(tau_region)) tibble() else tau %>%
  filter(BRAIN_REGION == tau_region) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), VISCODE2 = paste0("v", VISCODE), days = num(SCAN_DATE_DAYS_T0),
            value = coalesce(num(SUVR_CRUS), num(SUVR_PERSI)), tracer = LIGAND) %>%
  mutate(days = coalesce(days, sv_day(PTID, VISCODE))) %>%
  filter(!is.na(value), !is.na(days)) %>% select(-VISCODE) %>%
  mutate(measure = "tau_pet", source = "flortaucipir_mubada", derived = FALSE)

# ---- 3. CSF from the central laboratory file ---------------------------------------------
lab <- read_adni(A4_FILES$lab)
csf_map <- c("A-Beta 1-42, CSF Free-221" = "csf_abeta42", "CSF Free A-Beta1-40-221" = "csf_abeta40",
             "P-Tau-181, CSF-221" = "csf_ptau181", "Total Tau, CSF-221" = "csf_tau")
csf <- lab %>% filter(LBTEST %in% names(csf_map)) %>%
  transmute(PTID = BID, VISCODE2 = paste0("v", as.integer(VISCODE)), days = num(LBDTM_DAYS_T0),
            measure = unname(csf_map[LBTEST]), value = coalesce(num(SIRESN), num(RPTRESN)), unit = coalesce(SIU, RPTU)) %>%
  filter(!is.na(value), !is.na(days)) %>%
  group_by(PTID, VISCODE2, days, measure) %>% summarise(value = mean(value), .groups = "drop")
csf_ratios <- csf %>% select(-any_of("unit")) %>%
  pivot_wider(names_from = measure, values_from = value) %>%
  mutate(csf_abeta42_40  = if ("csf_abeta40" %in% names(.)) csf_abeta42 / csf_abeta40 else NA_real_,
         csf_ptau_abeta42 = if ("csf_ptau181" %in% names(.)) csf_ptau181 / csf_abeta42 else NA_real_,
         csf_ptau_tau     = if (all(c("csf_ptau181", "csf_tau") %in% names(.))) csf_ptau181 / csf_tau else NA_real_) %>%
  select(PTID, VISCODE2, days, any_of(c("csf_abeta42_40", "csf_ptau_abeta42", "csf_ptau_tau"))) %>%
  pivot_longer(-c(PTID, VISCODE2, days), names_to = "measure", values_to = "value") %>% filter(!is.na(value))
csf_all <- bind_rows(csf %>% mutate(derived = FALSE), csf_ratios %>% mutate(derived = TRUE)) %>%
  mutate(source = "central_lab_csf", tracer = NA_character_)

# ---- 4. Cognition and function --------------------------------------------------------------
mmse <- read_adni(A4_FILES$mmse) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), value = num(MMSCORE)) %>% filter(!is.na(value)) %>%
  left_join(sv, by = c("PTID", "VISCODE")) %>%
  transmute(PTID, VISCODE2 = paste0("v", VISCODE), days = sv_days, measure = "mmse", value, source = "mmse", derived = FALSE)
pacc <- read_adni(A4_FILES$pacc) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), value = num(PACC)) %>% filter(!is.na(value)) %>%
  left_join(sv, by = c("PTID", "VISCODE")) %>%
  transmute(PTID, VISCODE2 = paste0("v", VISCODE), days = sv_days, measure = "pacc", value, source = "pacc_derived", derived = TRUE)
cdr <- read_adni(A4_FILES$cdr) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), days = num(CDADTC_DAYS_T0),
            cdrsb = num(CDRSB), cdr_global = num(CDGLOBAL)) %>%
  left_join(sv, by = c("PTID", "VISCODE")) %>% mutate(days = coalesce(days, sv_days)) %>%
  pivot_longer(c(cdrsb, cdr_global), names_to = "measure", values_to = "value") %>% filter(!is.na(value)) %>%
  transmute(PTID, VISCODE2 = paste0("v", VISCODE), days, measure, value, source = "cdr", derived = FALSE)
cog <- bind_rows(mmse, pacc, cdr) %>% filter(!is.na(days)) %>% mutate(tracer = NA_character_)

say("Amyloid scans dated from the visit file: %d of %d", sum(amy$date_source == "visit file") / 3, nrow(amy) / 3)
amy <- amy %>% select(-date_source)
# ---- 5. Combine, match to plasma draws -----------------------------------------------------
crit <- bind_rows(amy, tau, csf_all, cog) %>%
  inner_join(ids, by = "PTID") %>%
  mutate(exam_date = a4_date(days), value_raw = as.character(value), comparator = NA_character_,
         lloq_flag = FALSE, uloq_flag = FALSE) %>%
  select(RID, source, measure, exam_date, VISCODE2, value, value_raw, comparator, lloq_flag, uloq_flag, derived, tracer) %>%
  distinct() %>%
  arrange(RID, measure, exam_date) %>%
  group_by(RID, measure) %>% mutate(assessment_no = row_number(), n_assessments = n()) %>% ungroup()
crit <- match_to_plasma(crit, "exam_date", plasma_visits) %>%
  select(RID, source, measure, exam_date, VISCODE2, assessment_no, n_assessments, value, value_raw, comparator,
         lloq_flag, uloq_flag, derived, visit_no, visit_date, gap_days, plasma_time, tracer)

save_rds(crit, "criteria_long")
write_csv(crit, file.path(OUT_DIR, "criteria_long.csv"))

# ---- 6. Checks ------------------------------------------------------------------------------
cov <- crit %>% group_by(measure) %>%
  summarise(people = n_distinct(RID), assessments = n(),
            people_2plus = n_distinct(RID[n_assessments >= 2]),
            within_90d = sum(abs(gap_days) <= 90), within_180d = sum(abs(gap_days) <= 180), within_365d = sum(abs(gap_days) <= 365),
            .groups = "drop")
write_csv(cov, file.path(CHECK_DIR, "a4_criteria_coverage.csv"))
byv <- crit %>% filter(abs(gap_days) <= MATCH_WINDOW_DAYS[["criteria"]]) %>%
  distinct(RID, visit_no, measure) %>% count(measure, visit_no) %>%
  pivot_wider(names_from = visit_no, values_from = n, names_prefix = "plasma_visit_")
write_csv(byv, file.path(CHECK_DIR, "a4_criteria_match_by_plasma_visit.csv"))
say("criteria_long: %d rows, %d people, %d measures", nrow(crit), n_distinct(crit$RID), n_distinct(crit$measure))
print(as.data.frame(cov)); print(as.data.frame(byv))


# =============================================================================================
# Section 4. Modifiers
# Original script: 09d_a4_modifiers_long.R
# =============================================================================================
# =============================================================================
# 09d_a4_modifiers_long.R
# A4/LEARN: physiological modifiers in the layout of ADNI's modifiers_long, plus a
# diagnosis file (dx_long) so that 06_master.R runs unchanged.
#
#   creatinine (mg/dL; the laboratory reports umol/L, divided by 88.42), egfr_2021,
#   egfr_2009 (CKD-EPI, age at the draw, sex), hba1c, crp, haemoglobin
#                                  central laboratory file: screening, baseline and the
#                                  scheduled clinic visits, A4 arm only (LEARN: screening only)
#   weight_kg, height_cm, bmi, sbp, dbp, pulse
#                                  vitals file, every visit, both substudies. Height is measured
#                                  at screening only and carried forward for BMI.
#   dx                             "CN" for everyone at every visit (both substudies enrolled
#                                  cognitively unimpaired people); CDR global is in the criteria
#
# Output   a4learn_out/modifiers_long.rds, dx_long.rds, modifiers_long.csv
#          checks/a4_modifiers_coverage.csv, a4_modifiers_match_by_plasma_visit.csv
# =============================================================================

Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
if (COHORT != "a4learn") stop("Set Sys.setenv(COHORT = 'a4learn') before running the 09 scripts.")

pl   <- readRDS(file.path(OUT_DIR, "plasma_long.rds"))
part <- readRDS(file.path(OUT_DIR, "participants.rds"))
num  <- function(x) suppressWarnings(as.numeric(x))
ids  <- part %>% select(PTID, RID, PTGENDER, AGE, screen_day_t0)
plasma_visits <- pl %>% distinct(RID, visit_no, visit_date) %>%
  group_by(RID) %>% mutate(first_plasma = min(visit_date)) %>% ungroup()
sv <- read_adni(A4_FILES$sv) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISITCD), sv_days = num(SVSTDTC_DAYS_T0)) %>%
  filter(!is.na(sv_days)) %>% distinct(PTID, VISCODE, .keep_all = TRUE)
age_at <- function(AGE, screen_day_t0, days) AGE + (days - screen_day_t0) / YEAR_DAYS   # AGE = age at the screening exam
# last observation carried forward within a person (height measured at screening only)
carry_forward <- function(x) { for (i in seq_along(x)) if (is.na(x[i]) && i > 1) x[i] <- x[i - 1]; x }

# ---- 1. Central laboratory ------------------------------------------------------------
lab <- read_adni(A4_FILES$lab)
lab_map <- c("Creatinine (Rate Blanked)" = "creatinine_umol", "Hemoglobin A1c,Variant(-70)" = "hba1c",
             "C-Reactive Protein" = "crp", "Hemoglobin" = "haemoglobin")
lb <- lab %>% filter(LBTEST %in% names(lab_map)) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), days = num(LBDTM_DAYS_T0),
            measure = unname(lab_map[LBTEST]), value = coalesce(num(SIRESN), num(RPTRESN)), unit = coalesce(SIU, RPTU)) %>%
  filter(!is.na(value), !is.na(days)) %>%
  group_by(PTID, VISCODE, days, measure) %>% summarise(value = mean(value), unit = first(unit), .groups = "drop")
creat <- lb %>% filter(measure == "creatinine_umol") %>%
  inner_join(ids, by = "PTID") %>%
  mutate(creatinine = if_else(str_detect(tolower(coalesce(unit, "")), "umol"), value / 88.42, value),
         age_at_exam = age_at(AGE, screen_day_t0, days),
         egfr_2021 = egfr_ckdepi_2021(creatinine, age_at_exam, PTGENDER),
         egfr_2009 = egfr_ckdepi_2009(creatinine, age_at_exam, PTGENDER)) %>%
  select(PTID, VISCODE, days, creatinine, egfr_2021, egfr_2009, age_at_exam) %>%
  pivot_longer(c(creatinine, egfr_2021, egfr_2009), names_to = "measure", values_to = "value") %>%
  mutate(derived = measure != "creatinine", source = "central_lab")
other_lab <- lb %>% filter(measure != "creatinine_umol") %>% inner_join(ids %>% select(PTID, AGE, screen_day_t0), by = "PTID") %>%
  mutate(age_at_exam = age_at(AGE, screen_day_t0, days), derived = FALSE, source = "central_lab") %>%
  select(PTID, VISCODE, days, measure, value, age_at_exam, derived, source)

# ---- 2. Vitals -----------------------------------------------------------------------
vit <- read_adni(A4_FILES$vitals) %>%
  transmute(PTID = BID, VISCODE = as.integer(VISCODE), weight_kg = num(STDWT), height_cm = num(STDHT),
            sbp = num(VSBPSYS), dbp = num(VSBPDIA), pulse = num(VSPULSE)) %>%
  left_join(sv, by = c("PTID", "VISCODE")) %>% rename(days = sv_days) %>% filter(!is.na(days)) %>%
  arrange(PTID, days) %>% group_by(PTID) %>%
  mutate(height_source = if_else(!is.na(height_cm), "measured", "carried"),
         height_cm = carry_forward(height_cm)) %>%
  ungroup() %>%
  mutate(bmi = if_else(!is.na(weight_kg) & !is.na(height_cm), weight_kg / (height_cm / 100)^2, NA_real_)) %>%
  inner_join(ids %>% select(PTID, AGE, screen_day_t0), by = "PTID") %>%
  mutate(age_at_exam = age_at(AGE, screen_day_t0, days)) %>%
  pivot_longer(c(weight_kg, height_cm, bmi, sbp, dbp, pulse), names_to = "measure", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(derived = measure == "bmi", source = "vitals",
         height_source = if_else(measure %in% c("height_cm", "bmi"), height_source, NA_character_)) %>%
  select(PTID, VISCODE, days, measure, value, age_at_exam, derived, source, height_source)

# ---- 3. Combine and match to plasma draws -------------------------------------------------
mods <- bind_rows(creat, other_lab, vit) %>%
  inner_join(ids %>% select(PTID, RID), by = "PTID") %>%
  mutate(exam_date = a4_date(days), VISCODE2 = paste0("v", VISCODE)) %>%
  select(RID, source, measure, exam_date, VISCODE2, value, derived, age_at_exam, any_of("height_source")) %>%
  distinct(RID, measure, exam_date, .keep_all = TRUE) %>%
  arrange(RID, measure, exam_date) %>%
  group_by(RID, measure) %>% mutate(assessment_no = row_number(), n_assessments = n()) %>% ungroup()
if (!"height_source" %in% names(mods)) mods$height_source <- NA_character_
mods <- match_to_plasma(mods, "exam_date", plasma_visits) %>%
  select(RID, source, measure, exam_date, VISCODE2, assessment_no, n_assessments, value, derived,
         age_at_exam, height_source, visit_no, visit_date, gap_days, plasma_time)
save_rds(mods, "modifiers_long")
write_csv(mods, file.path(OUT_DIR, "modifiers_long.csv"))

# ---- 4. Diagnosis: cognitively unimpaired at enrolment; one row per plasma draw ------------
dx <- pl %>% distinct(RID, visit_no, visit_date, VISCODE2, plasma_time) %>%
  transmute(RID, dx = "CN", exam_date = visit_date, VISCODE2, visit_no, visit_date, gap_days = 0L, plasma_time)
save_rds(dx, "dx_long")

# ---- 5. Checks ----------------------------------------------------------------------------
cov <- mods %>% group_by(measure) %>%
  summarise(people = n_distinct(RID), assessments = n(), people_2plus = n_distinct(RID[n_assessments >= 2]),
            within_365d_of_a_draw = sum(abs(gap_days) <= 365), .groups = "drop")
write_csv(cov, file.path(CHECK_DIR, "a4_modifiers_coverage.csv"))
byv <- mods %>% filter(abs(gap_days) <= MATCH_WINDOW_DAYS[["modifiers"]]) %>%
  distinct(RID, visit_no, measure) %>% count(measure, visit_no) %>%
  pivot_wider(names_from = visit_no, values_from = n, names_prefix = "plasma_visit_")
write_csv(byv, file.path(CHECK_DIR, "a4_modifiers_match_by_plasma_visit.csv"))
say("modifiers_long: %d rows, %d people, %d measures", nrow(mods), n_distinct(mods$RID), n_distinct(mods$measure))
print(as.data.frame(cov)); print(as.data.frame(byv))
