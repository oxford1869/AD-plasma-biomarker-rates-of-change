# =============================================================================================
# make_test_data.R
# Synthetic test files in the layouts of the real releases, for exercising the code without the data.
#
# Three sections, formerly zz_make_fake_fnih.R, zz_make_fake_others.R and zz_make_fake_a4learn.R. Not part of
# run_all.R. The output is random and is not an estimate of anything; never point these at the real data folders.
# =============================================================================================


# =============================================================================================
# Section 1. Fake FNIH plasma file
# Original script: zz_make_fake_fnih.R
# =============================================================================================
# Builds a small fake FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv with the same
# columns and quirks as the real one, for testing the pipeline only.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
set.seed(1)

# TEST HELPER ONLY. Do not run this against your real data folder.
# It writes into a subfolder called FAKE_TEST and refuses to overwrite anything.
base <- Sys.getenv("FAKE_JM_DIR", unset = tempdir())
out_dir <- file.path(base, "FAKE_TEST")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
target <- file.path(out_dir, "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv")
if (file.exists(target)) stop("Refusing to overwrite existing file: ", target)
cat("Fake data will be written to:", out_dir, "\n")

n_people <- 40
people <- tibble(RID = sample(1:6000, n_people)) %>%
  mutate(PTID = sprintf("%03d_S_%04d", sample(1:200, n_people, TRUE), RID),
         n_visits = sample(1:3, n_people, replace = TRUE, prob = c(.15, .35, .5)),
         first_date = as.Date("2010-01-01") + sample(0:1500, n_people, TRUE))

visits <- people %>%
  mutate(visit = map(n_visits, seq_len)) %>%
  unnest(visit) %>%
  mutate(EXAMDATE = first_date + (visit - 1) * round(runif(n(), 700, 1500)),
         VISCODE2 = c("bl", "m24", "m48")[visit],
         SAMPLEID = sprintf("AA8%05d-%02d", RID, visit * 2))

# measure spec: PLASMA_BIOMARKER, TESTNAME, lab fields, unit, typical scale
spec <- tribble(
  ~PLASMA_BIOMARKER, ~TESTNAME, ~ASSAYPLATFORM, ~ASSAYVERSION, ~PERFORMINGLAB, ~INSTRUMENTID, ~UNITS, ~mu, ~sd,
  "C2N_plasma_Abeta40",         "Abeta40",       "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "pg/mL", 450, 60,
  "C2N_plasma_Abeta42",         "Abeta42",       "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "pg/mL", 45,  8,
  "C2N_plasma_Abeta42_Abeta40", "Abeta_ratio",   "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "",      0.10, 0.015,
  "C2N_plasma_ptau217",         "pTau217",       "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "pg/mL", 2.5, 2.5,
  "C2N_plasma_nptau217",        "npTau217",      "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "pg/mL", 55,  15,
  "C2N_plasma_ptau217_ratio",   "pTau217_ratio", "Precivity", "PrecivityAD2", "C2N Diagnostics", "LC-MS/MS", "%",     4,   3,
  "Fuji_plasma_ptau217",        "pTau217",       "Lumipulse", "G1200",        "Indiana U",       "Lumipulse G", "pg/mL", 0.15, 0.1,
  "Fuji_plasma_Ab40",           "Ab40",          "Lumipulse", "G1200",        "Indiana U",       "Lumipulse G", "pg/mL", 350, 50,
  "Fuji_plasma_Ab42",           "Ab42",          "Lumipulse", "G1200",        "Indiana U",       "Lumipulse G", "pg/mL", 35,  6,
  "AlzPath_plasma_ptau217",     "pTau217",       "Simoa",     "ALZpath",      "Quanterix",       "HD-X",   "pg/mL", 0.4, 0.3,
  "Janssen_plasma_ptau217",     "pTau217",       "Simoa",     "Janssen",      "Quanterix",       "HD-X",   "pg/mL", 0.05, 0.04,
  "Roche_plasma_Ab40",          "Ab40",          "Elecsys",   "cobas",        "U Gothenburg",    "cobas e", "ng/mL", 0.35, 0.05,
  "Roche_plasma_Ab42",          "Ab42",          "Elecsys",   "cobas",        "U Gothenburg",    "cobas e", "pg/mL", 45,  8,
  "Roche_plasma_GFAP",          "GFAP",          "Elecsys",   "cobas",        "U Gothenburg",    "cobas e", "ng/mL", 0.15, 0.05,
  "Roche_plasma_NfL",           "NfL",           "Elecsys",   "cobas",        "U Gothenburg",    "cobas e", "pg/mL", 4,   1.5,
  "Roche_plasma_ptau181",       "pTau181",       "Elecsys",   "cobas",        "U Gothenburg",    "cobas e", "pg/mL", 1.2, 0.5,
  "QX_plasma_ptau181",          "pTau181",       "Simoa",     "pTau181 v2",   "Quanterix",       "HD-X",   "pg/mL", 20,  8,
  "QX_plasma_Ab40",             "Ab40",          "Simoa",     "N4PE",         "Quanterix",       "HD-X",   "pg/mL", 130, 25,
  "QX_plasma_Ab42",             "Ab42",          "Simoa",     "N4PE",         "Quanterix",       "HD-X",   "pg/mL", 7.5, 1.5,
  "QX_plasma_GFAP",             "GFAP",          "Simoa",     "N4PE",         "Quanterix",       "HD-X",   "pg/mL", 220, 70,
  "QX_plasma_NfL",              "NfL",           "Simoa",     "N4PE",         "Quanterix",       "HD-X",   "pg/mL", 25,  10
)

rows <- tidyr::crossing(visits, spec) %>%
  mutate(
    val = pmax(rnorm(n(), mu, sd), 0.01),
    # Roche skipped some early visits for some people
    drop_roche = grepl("^Roche_", PLASMA_BIOMARKER) & visit == 1 & (RID %% 3 == 0),
    # Quanterix Ab42 missing for some
    drop_qx42 = PLASMA_BIOMARKER == "QX_plasma_Ab42" & (RID %% 4 == 0),
    # C2N pTau217 floor
    at_floor = PLASMA_BIOMARKER == "C2N_plasma_ptau217" & val < 1.2,
    TESTVALUE = case_when(
      at_floor ~ "0.65",
      PLASMA_BIOMARKER == "Roche_plasma_ptau181" & val < 0.6 ~ "<0.6",
      TRUE ~ formatC(val, format = "f", digits = 3)
    ),
    COMMENTS = case_when(
      at_floor ~ "LOD",
      TRUE ~ ""
    )
  ) %>%
  filter(!drop_roche, !drop_qx42)

# C2N ratio imputed flag follows the pTau217 floor
floor_keys <- rows %>% filter(at_floor) %>% distinct(RID, EXAMDATE)
rows <- rows %>%
  left_join(floor_keys %>% mutate(ratio_imp = TRUE), by = c("RID", "EXAMDATE")) %>%
  mutate(COMMENTS = if_else(PLASMA_BIOMARKER == "C2N_plasma_ptau217_ratio" & coalesce(ratio_imp, FALSE),
                            "based on imputed p-tau217 concentration", COMMENTS))

lab_run <- rows %>% distinct(PERFORMINGLAB, RID, EXAMDATE) %>%
  group_by(PERFORMINGLAB) %>%
  mutate(RUN_NO = as.character(sample(1:6, n(), TRUE)),
         RUNDATE = format(as.Date("2024-02-20") + as.integer(RUN_NO), "%d/%m/%Y")) %>%
  ungroup()

# Split-visit quirk: for three people, Roche reports visit 2 one day later
# than the other labs. Same SAMPLEID, different EXAMDATE, missing VISCODE2.
split_rids <- head(sort(unique(rows$RID[rows$visit == 2])), 3)
rows <- rows %>%
  mutate(is_split = RID %in% split_rids & visit == 2 & PERFORMINGLAB == "U Gothenburg",
         EXAMDATE = if_else(is_split, EXAMDATE + 1, EXAMDATE),
         VISCODE2 = if_else(is_split, NA_character_, VISCODE2))
lab_run <- lab_run %>%
  bind_rows(rows %>% filter(is_split) %>% distinct(PERFORMINGLAB, RID, EXAMDATE) %>%
              mutate(RUN_NO = "9", RUNDATE = "29/02/2024"))

fake <- rows %>%
  left_join(lab_run, by = c("PERFORMINGLAB", "RID", "EXAMDATE")) %>%
  transmute(
    STUDY_ID = "FNIH BC ADNI 00002", SAMPLEID, PTID, RID,
    EXAMDATE = format(EXAMDATE, "%d/%m/%Y"), VISCODE2,
    ASSAYPLATFORM, ASSAYVERSION, PERFORMINGLAB, INSTRUMENTID,
    MATRIX = "K2 EDTA plasma", TESTNAME, PLASMA_BIOMARKER, TESTVALUE, UNITS,
    CV = if_else(PERFORMINGLAB == "Quanterix", formatC(runif(n(), 2, 8), digits = 2, format = "f"), ""),
    RUN_NO, SEQUENCE_NO = "", RUNDATE, COMMENTS
  ) %>%
  arrange(RID, EXAMDATE, PERFORMINGLAB, TESTNAME)

write_csv(fake, target, na = "")
cat("Fake file written:", nrow(fake), "rows,", n_distinct(fake$RID), "people\n")


# =============================================================================================
# Section 2. Fake ADNIMERGE, PET, CSF, lab and demographic files
# Original script: zz_make_fake_others.R
# =============================================================================================
# TEST HELPER ONLY. Builds small fake versions of ADNIMERGE_DATA.csv,
# UCBERKELEY_AMY_6MM.CSV, UPENNBIOMK_ROCHE_ELECSYS.csv, LABDATA.CSV,
# PTDEMOG.csv and MMSE.CSV next to the fake plasma file, using the same
# people and dates. Never run this against your real data folder. It writes
# only into FAKE_JM_DIR/FAKE_TEST and refuses to overwrite.
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
set.seed(2)

base <- Sys.getenv("FAKE_JM_DIR", unset = tempdir())
out_dir <- file.path(base, "FAKE_TEST")
plasma_file <- file.path(out_dir, "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv")
if (!file.exists(plasma_file)) stop("Run zz_make_fake_fnih.R first: ", plasma_file)
targets <- file.path(out_dir, c("ADNIMERGE_DATA.csv", "UCBERKELEY_AMY_6MM.CSV",
                                "UPENNBIOMK_ROCHE_ELECSYS.csv", "LABDATA.CSV",
                                "PTDEMOG.csv", "MMSE.CSV", "VITALS.csv"))
if (any(file.exists(targets))) stop("Refusing to overwrite existing fake files in ", out_dir)

pl <- read_csv(plasma_file, col_types = cols(.default = col_character())) %>%
  transmute(RID = as.integer(RID), PTID, date = as.Date(EXAMDATE, "%d/%m/%Y"), VISCODE2) %>%
  distinct() %>% arrange(RID, date)

people <- pl %>% group_by(RID, PTID) %>%
  summarise(first = min(date), .groups = "drop") %>%
  mutate(bl_date  = first - sample(0:1200, n(), TRUE),
         AGE      = round(runif(n(), 60, 85), 1),
         PTGENDER = sample(c("Male", "Female"), n(), TRUE),
         PTEDUCAT = sample(12:20, n(), TRUE),
         APOE4    = sample(0:2, n(), TRUE, prob = c(.6, .3, .1)),
         DX_bl    = sample(c("CN", "SMC", "EMCI", "LMCI", "AD"), n(), TRUE),
         COLPROT  = sample(c("ADNI2", "ADNI3", "ADNIGO"), n(), TRUE),
         ORIGPROT = sample(c("ADNI1", "ADNI2", "ADNIGO"), n(), TRUE))

# ADNIMERGE: baseline plus a visit every ~12 months, many near plasma dates
visits <- people %>%
  mutate(k = sample(4:9, n(), TRUE)) %>%
  mutate(v = map(k, ~ seq_len(.x))) %>% unnest(v) %>%
  mutate(EXAMDATE = bl_date + round((v - 1) * runif(n(), 330, 400)),
         VISCODE  = if_else(v == 1, "bl", paste0("m", (v - 1) * 12)),
         Years_bl = as.numeric(EXAMDATE - bl_date) / 365.25,
         MMSE     = pmin(30, round(29 - 0.4 * (v - 1) + rnorm(n(), 0, 1.5))),
         ADAS13   = pmax(0, round(12 + 1.2 * (v - 1) + rnorm(n(), 0, 4), 1)),
         mPACCdigit = round(-0.5 * (v - 1) + rnorm(n(), 0, 2), 2),
         CDRSB    = pmax(0, round(0.5 * (v - 1) + rnorm(n(), 0, 1), 1)),
         DX       = sample(c("CN", "MCI", "Dementia", ""), n(), TRUE, prob = c(.5, .35, .1, .05)),
         Hippocampus = round(7000 - 60 * (v - 1) + rnorm(n(), 0, 400)),
         ICV      = round(1500000 + rnorm(n(), 0, 120000)),
         ABETA    = if_else(runif(n()) < .3, as.character(round(runif(n(), 500, 2000))), ""),
         PTAU     = if_else(ABETA != "", as.character(round(runif(n(), 10, 60), 1)), "")) %>%
  select(RID, PTID, VISCODE, COLPROT, ORIGPROT, EXAMDATE, DX_bl, AGE, PTGENDER, PTEDUCAT,
         APOE4, MMSE, ADAS13, mPACCdigit, CDRSB, DX, Hippocampus, ICV, ABETA, PTAU, Years_bl)
# one RID with a missing DX_bl row order quirk: drop the baseline row for one person
visits <- visits %>% filter(!(RID == min(RID) & VISCODE == "bl"))
write_csv(visits %>% mutate(EXAMDATE = format(EXAMDATE, "%Y-%m-%d")),
          file.path(out_dir, "ADNIMERGE_DATA.csv"), na = "")

# Amyloid PET: 1 to 3 scans per person, most within 6 months of a plasma date
amy <- pl %>% group_by(RID) %>% slice_sample(n = 2) %>% ungroup() %>%
  mutate(SCANDATE = date + sample(-200:200, n(), TRUE),
         VISCODE  = VISCODE2,
         TRACER   = sample(c("FBP", "FBB"), n(), TRUE),
         CENTILOIDS = round(rnorm(n(), 25, 30), 2),
         SUMMARY_SUVR = round(1 + CENTILOIDS / 150, 4),
         AMYLOID_STATUS = as.integer(CENTILOIDS > 20)) %>%
  select(RID, VISCODE, SCANDATE, TRACER, SUMMARY_SUVR, CENTILOIDS, AMYLOID_STATUS)
write_csv(amy %>% mutate(SCANDATE = format(SCANDATE, "%Y-%m-%d")),
          file.path(out_dir, "UCBERKELEY_AMY_6MM.CSV"), na = "")

# CSF Elecsys: about half the people, 1 to 2 draws, some with comparators
csf <- pl %>% filter(RID %% 2 == 0) %>% group_by(RID) %>% slice_sample(n = 2) %>% ungroup() %>%
  mutate(EXAMDATE = date + sample(-60:60, n(), TRUE),
         ABETA42  = round(runif(n(), 400, 1900)),
         ABETA40  = round(runif(n(), 8000, 25000)),
         TAU      = round(runif(n(), 150, 500)),
         PTAU     = round(runif(n(), 12, 50), 1),
         ABETA42  = if_else(ABETA42 > 1700, ">1700", as.character(ABETA42)),
         PTAU     = if_else(PTAU < 13, "<8", as.character(PTAU))) %>%
  select(RID, VISCODE2, EXAMDATE, ABETA42, ABETA40, TAU, PTAU)
write_csv(csf %>% mutate(EXAMDATE = format(EXAMDATE, "%Y-%m-%d")),
          file.path(out_dir, "UPENNBIOMK_ROCHE_ELECSYS.csv"), na = "")

# LABDATA: creatinine (RCT392) at several visits, with the odd non-numeric code
lab <- visits %>% select(RID, VISCODE, EXAMDATE) %>%
  mutate(VISCODE2 = VISCODE,
         RCT392 = as.character(round(runif(n(), 0.6, 1.6), 2)),
         RCT392 = if_else(runif(n()) < .03, sample(c("-1", "SCC09", ""), n(), TRUE), RCT392)) %>%
  select(RID, VISCODE, VISCODE2, EXAMDATE, RCT392)
write_csv(lab %>% mutate(EXAMDATE = format(EXAMDATE, "%Y-%m-%d")),
          file.path(out_dir, "LABDATA.CSV"), na = "")

# PTDEMOG: sex as 1/2, birth month and year
demog <- people %>%
  transmute(RID, PTID, VISCODE = "sc",
            PTGENDER = if_else(PTGENDER == "Male", 1L, 2L),
            PTDOBMM = sample(1:12, n(), TRUE),
            PTDOBYY = as.integer(format(bl_date, "%Y")) - floor(AGE),
            PTEDUCAT)
write_csv(demog, file.path(out_dir, "PTDEMOG.csv"), na = "")

# MMSE.CSV: same visits, MMSCORE
mmse <- visits %>% transmute(RID, VISCODE, VISCODE2 = VISCODE, VISDATE = EXAMDATE, MMSCORE = MMSE)
write_csv(mmse %>% mutate(VISDATE = format(VISDATE, "%Y-%m-%d")),
          file.path(out_dir, "MMSE.CSV"), na = "")

# VITALS: weight and height with mixed units (1 = lbs/in, 2 = kg/cm), some junk
vit <- visits %>% transmute(RID, VISCODE, VISCODE2 = VISCODE, EXAMDATE,
                            VSWTUNIT = sample(1:2, n(), TRUE),
                            VSHTUNIT = sample(1:2, n(), TRUE)) %>%
  mutate(wt_kg = rnorm(n(), 75, 12), ht_cm = rnorm(n(), 168, 9),
         VSWEIGHT = if_else(VSWTUNIT == 1, round(wt_kg * 2.20462, 1), round(wt_kg, 1)),
         VSHEIGHT = if_else(VSHTUNIT == 1, round(ht_cm / 2.54, 1), round(ht_cm, 1)),
         VSHEIGHT = if_else(runif(n()) < .3, NA_real_, VSHEIGHT),
         VSWEIGHT = if_else(runif(n()) < .02, -4, VSWEIGHT)) %>%
  select(RID, VISCODE, VISCODE2, EXAMDATE, VSWEIGHT, VSWTUNIT, VSHEIGHT, VSHTUNIT)
write_csv(vit %>% mutate(EXAMDATE = format(EXAMDATE, "%Y-%m-%d")),
          file.path(out_dir, "VITALS.csv"), na = "")

cat("Fake companion files written to", out_dir, "\n")


# =============================================================================================
# Section 3. Fake A4/LEARN release
# Original script: zz_make_fake_a4learn.R
# =============================================================================================
# =============================================================================
# zz_make_fake_a4learn.R
# Writes a small fake A4/LEARN release with the file layouts of the real one, so
# the 09 scripts can be tested without the data. Never touches real files.
#   Rscript zz_make_fake_a4learn.R /tmp/fake_a4        (60 people)
#   Rscript zz_make_fake_a4learn.R /tmp/fake_a4_big 8  (8 times larger, for testing the analysis scripts)
# Then: A4_DIR=/tmp/fake_a4 ADNI_OUT_DIR=/tmp/out_a4 COHORT=a4learn Rscript 09a_a4_plasma_long.R
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(readr); library(purrr) })
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "/tmp/fake_a4"
SCALE <- if (length(args) > 1) as.numeric(args[2]) else 1   # 2nd argument multiplies the sample sizes
for (d in c("External Data", "Derived Data", "Raw Data", "Documents/Data Dictionaries")) dir.create(file.path(root, d), recursive = TRUE, showWarnings = FALSE)
set.seed(4)

n_pl <- round(25 * SCALE); n_so <- round(20 * SCALE); n_le <- round(15 * SCALE); N <- n_pl + n_so + n_le
subj <- tibble(SUBSTUDY = c(rep("A4", n_pl + n_so), rep("LEARN", n_le)),
               BID = sprintf("B%08d", 10000000 + seq_len(N)),
               TX = c(rep("Placebo", n_pl), rep("Solanezumab", n_so), rep(NA, n_le)),
               AGEYR = round(runif(N, 65, 85), 2), SEX = sample(1:2, N, TRUE),
               RACE = sample(c(1, 1, 1, 2, 58), N, TRUE), EDCCNTU = sample(12:20, N, TRUE), ETHNIC = sample(c(50, 56, 56, 56), N, TRUE),
               APOEGN = sample(c("E3/E3", "E3/E4", "E4/E4", "E2/E3"), N, TRUE, prob = c(.5, .35, .1, .05)),
               ENRLDT_DAYS_CONSENT = sample(60:120, N, TRUE), ENRLDT_DAYS_T0 = 0,
               RNDMDT_DAYS_T0 = 0, TRTSTDT_DAYS_T0 = 0, TRTENDT_DAYS_T0 = sample(1500:1800, N, TRUE),
               DISCDTC_DAYS_T0 = NA, HGTCMBL = round(runif(N, 150, 185), 1), WGTKGBL = round(runif(N, 55, 100), 1),
               SUBJCOMPTR = sample(c(0, 1), N, TRUE, prob = c(.2, .8)), GDSSEVBL = 0) %>%
  mutate(BMIBL = WGTKGBL / (HGTCMBL / 100)^2, APOEGNPRSNFLG = as.integer(grepl("E4", APOEGN)))
subj$TRTENDT_DAYS_T0[subj$SUBSTUDY == "LEARN"] <- NA; subj$TRTSTDT_DAYS_T0[subj$SUBSTUDY == "LEARN"] <- NA
# a few screen failures
sf <- tibble(SUBSTUDY = "SF", BID = sprintf("B%08d", 10000900 + 1:5), TX = NA, AGEYR = 70, SEX = 1, RACE = 1, EDCCNTU = 16, ETHNIC = 56,
             APOEGN = "E3/E3", ENRLDT_DAYS_CONSENT = 0, ENRLDT_DAYS_T0 = NA, RNDMDT_DAYS_T0 = NA, TRTSTDT_DAYS_T0 = NA, TRTENDT_DAYS_T0 = NA,
             DISCDTC_DAYS_T0 = NA, HGTCMBL = 170, WGTKGBL = 70, SUBJCOMPTR = 0, GDSSEVBL = 0, BMIBL = 24.2, APOEGNPRSNFLG = 0)
write_csv(bind_rows(subj, sf), file.path(root, "Derived Data", "SUBJINFO.csv"), na = "NA")

# visits dictionary (subset)
vis <- tribble(~SUBSTUDY, ~VISCODE, ~VISNAME, ~VISORDER, ~ATTRIBS, ~VISITGROUP,
  "A4", 1, "Visit 1 (Screening)", 1, "sc1", "Screening", "A4", 2, "Visit 2 (Screening PET)", 2, "sc2", "Screening",
  "A4", 5, "Visit 5 (Screening LP)", 5, "sc5", "Screening", "A4", 6, "Visit 6 (Baseline)", 6, "bl6", "Blinded",
  "A4", 9, "Visit 9 (wk12 Clinic)", 9, "w12", "Blinded", "A4", 12, "Visit 12 (wk24 Clinic)", 12, "w24", "Blinded",
  "A4", 18, "Visit 18 (wk48 Clinic)", 18, "w48", "Blinded", "A4", 33, "Visit 33 (wk108 Clinic)", 33, "w108", "Blinded",
  "A4", 48, "Visit 48 (wk168 Clinic)", 48, "w168", "Blinded", "A4", 66, "Visit 66 (wk240 End DB/Start OLE)", 66, "w240", "Blinded",
  "A4", 997, "Final OLE [ET or Study Close]", 997, "et2", "Open-label", "A4", 999, "Early Termination [DOUBLE BLIND]", 999, "et1", "Blinded",
  "LEARN", 1, "A4 Visit 1 (Screening)", 1, "sc1", "Screening", "LEARN", 2, "A4 Visit 2 (Screening PET)", 2, "sc2", "Screening",
  "LEARN", 6, "Visit 1 (Baseline: A4 v4-6)", 6, "bl", "Blinded", "LEARN", 12, "Visit 2 (wk24: A4 v12)", 12, "w24", "Blinded",
  "LEARN", 24, "Visit 4 (wk72: A4 v24)", 24, "w72", "Blinded", "LEARN", 48, "Visit 8 (wk168: A4 v48)", 48, "w168", "Blinded",
  "LEARN", 66, "Visit 11 (wk240: A4v66)", 66, "w240", "Blinded", "LEARN", 999, "Early Termination (Protocol)", 999, "et1", "Unscheduled")
write_csv(vis, file.path(root, "Documents", "Data Dictionaries", "visits_datadic.csv"))

# visit days per person
sched <- list(A4 = c(`1` = -70, `2` = -45, `5` = -20, `6` = 0, `9` = 85, `12` = 170, `18` = 340, `33` = 760, `48` = 1180, `66` = 1700),
              LEARN = c(`1` = -90, `2` = -60, `6` = 0, `12` = 170, `24` = 505, `48` = 1180, `66` = 1750))
sv <- subj %>% select(SUBSTUDY, BID) %>% rowwise() %>%
  do({ s <- sched[[.$SUBSTUDY]]; tibble(SUBSTUDY = .$SUBSTUDY, BID = .$BID, VISITCD = as.integer(names(s)),
                                          VISIT = paste("Visit", names(s)), SVSTDTC_DAYS_T0 = s + sample(-5:5, length(s), TRUE)) }) %>% ungroup() %>%
  mutate(SVSTDTC_DAYS_CONSENT = SVSTDTC_DAYS_T0 + 90, SVTYPE = "Standard", SVUSEDTC_DAYS_CONSENT = SVSTDTC_DAYS_CONSENT, SVUSEDTC_DAYS_T0 = SVSTDTC_DAYS_T0)
write_csv(sv %>% mutate(VISITCD = sprintf("%03d", VISITCD)), file.path(root, "Derived Data", "SV.csv"))
day_of <- function(bid, vc) sv$SVSTDTC_DAYS_T0[match(paste(bid, vc), paste(sv$BID, sv$VISITCD))]

# plasma pTau217: A4 v6, v9, v66 (+ some v999); LEARN v1, v24, v66
truth <- subj %>% transmute(BID, SUBSTUDY, lev = if_else(SUBSTUDY == "A4", rnorm(n(), 0.3, 0.08), rnorm(n(), 0.18, 0.04)), sl = rnorm(n(), 0.02, 0.01))
pt <- truth %>% rowwise() %>% do({
  vcs <- if (.$SUBSTUDY == "A4") c(6, 9, 66) else c(1, 24, 66)
  if (runif(1) < 0.15) vcs <- c(vcs[1:2], 999)
  tibble(SUBSTUDY = .$SUBSTUDY, BID = .$BID, VISCODE = vcs,
         COLLECTION_DATE_DAYS_T0 = ifelse(vcs == 999, sample(300:1200, 1), day_of(.$BID, vcs)),
         val = .$lev + .$sl * COLLECTION_DATE_DAYS_T0 / 365.25 + rnorm(length(vcs), 0, 0.03)) }) %>% ungroup() %>%
  mutate(SAMPLE_ID = sprintf("%032x", row_number()), TESTCD = "PTAU217", TEST = "Phosphorylated-tau-217",
         ORRESRAW = round(val, 3), ORRES = case_when(val < 0.15 ~ "<LLOQ", val > 1.2 ~ ">ULOQ", TRUE ~ as.character(round(val, 3))),
         ORRESU = "U/mL", STAT = if_else(runif(n()) < 0.03, "Not Done", NA_character_), REASND = if_else(!is.na(STAT), "Hemolyzed", NA_character_),
         NAM = "Lilly Clinical Diagnostics Laboratory", SPEC = "Plasma", METHOD = "P-tau217 ECL MSD",
         COLLECTION_DATE_DAYS_CONSENT = COLLECTION_DATE_DAYS_T0 + 90, COMMENT = "not for diagnosis",
         COMMENT2 = if_else(ORRES %in% c("<LLOQ", ">ULOQ"), "Raw results in ORRESRAW fall outside the analytical measuring range", NA_character_)) %>%
  select(SUBSTUDY, BID, VISCODE, SAMPLE_ID, TESTCD, TEST, ORRES, ORRESU, STAT, REASND, NAM, SPEC, METHOD,
         COLLECTION_DATE_DAYS_CONSENT, COLLECTION_DATE_DAYS_T0, COMMENT, ORRESRAW, COMMENT2)
write_csv(pt, file.path(root, "External Data", "biomarker_pTau217.csv"), na = "")

# amyloid PET: v2 and v66 (80%), composite + one region
amy <- subj %>% rowwise() %>% do({
  vcs <- if (runif(1) < 0.8) c(2, 66) else 2
  base <- if (.$SUBSTUDY == "A4") runif(1, 1.15, 1.6) else runif(1, 0.9, 1.1)
  tibble(SUBSTUDY = .$SUBSTUDY, BID = .$BID, VISCODE = vcs, scan_number = seq_along(vcs),
         scan_date_DAYS_T0 = day_of(.$BID, vcs), suvr = base + 0.015 * scan_date_DAYS_T0 / 365.25 + rnorm(length(vcs), 0, 0.02)) }) %>% ungroup() %>%
  mutate(scan_date_DAYS_CONSENT = scan_date_DAYS_T0 + 90, ligand = "Florbetapir", scan_analyzed = "Yes") %>%
  tidyr::crossing(brain_region = c("Composite_Summary", "lprecuneus_gm")) %>%
  mutate(suvr_cer = round(suvr + if_else(brain_region == "Composite_Summary", 0, 0.05), 3), suvr_persi = NA, suvr_crus = NA) %>%
  select(SUBSTUDY, BID, VISCODE, scan_number, scan_date_DAYS_CONSENT, scan_date_DAYS_T0, ligand, scan_analyzed, brain_region, suvr_cer, suvr_persi, suvr_crus)
write_csv(amy, file.path(root, "External Data", "imaging_SUVR_amyloid.csv"), na = "NA")
tau <- amy %>% filter(brain_region == "Composite_Summary", SUBSTUDY == "A4") %>% slice_sample(prop = 0.4) %>%
  mutate(ligand = "Flortaucipir", brain_region = "MUBADA Mask", suvr_cer = NA, suvr_persi = round(runif(n(), 1, 1.4), 3), suvr_crus = round(runif(n(), 1, 1.4), 3))
write_csv(tau, file.path(root, "External Data", "imaging_SUVR_tau.csv"), na = "NA")

# central lab: creatinine etc. at v1, v6, v9, v12, v18, v33, v48, v66 (A4); v1 (LEARN); CSF at v5 and v66 for some A4
labvis <- function(sub) if (sub == "A4") c(1, 6, 9, 12, 18, 33, 48, 66) else 1
lab <- subj %>% rowwise() %>% do({
  vcs <- labvis(.$SUBSTUDY)
  tibble(SUBSTUDY = .$SUBSTUDY, BID = .$BID, VISCODE = vcs, LBDTM_DAYS_T0 = day_of(.$BID, vcs)) }) %>% ungroup() %>%
  tidyr::crossing(LBTEST = c("Creatinine (Rate Blanked)", "Hemoglobin A1c,Variant(-70)", "C-Reactive Protein", "Hemoglobin")) %>%
  mutate(BATTRNAM = "CHEMISTRY", LBTESTCD = substr(LBTEST, 1, 6), LBSPEC = "SERUM",
         SIRESN = case_when(LBTEST == "Creatinine (Rate Blanked)" ~ round(runif(n(), 55, 110), 0), LBTEST == "Hemoglobin A1c,Variant(-70)" ~ round(runif(n(), 5, 6.5), 1),
                            LBTEST == "C-Reactive Protein" ~ round(runif(n(), 0.3, 5), 2), TRUE ~ round(runif(n(), 120, 160), 0)),
         SIU = case_when(LBTEST == "Creatinine (Rate Blanked)" ~ "umol/L", LBTEST == "Hemoglobin A1c,Variant(-70)" ~ "%", LBTEST == "C-Reactive Protein" ~ "mg/L", TRUE ~ "g/L"),
         RPTRESN = SIRESN, RPTRESC = as.character(SIRESN), RPTU = SIU)
csf <- subj %>% filter(SUBSTUDY == "A4") %>% slice_sample(prop = 0.4) %>% rowwise() %>% do({
  vcs <- if (runif(1) < 0.6) c(5, 66) else 5
  tibble(SUBSTUDY = .$SUBSTUDY, BID = .$BID, VISCODE = vcs, LBDTM_DAYS_T0 = day_of(.$BID, vcs)) }) %>% ungroup() %>%
  tidyr::crossing(LBTEST = c("A-Beta 1-42, CSF Free-221", "CSF Free A-Beta1-40-221", "P-Tau-181, CSF-221", "Total Tau, CSF-221")) %>%
  mutate(BATTRNAM = "BIOMARKER", LBTESTCD = "CSF", LBSPEC = "CSF",
         SIRESN = case_when(grepl("42", LBTEST) ~ round(runif(n(), 400, 1200)), grepl("40", LBTEST) ~ round(runif(n(), 5000, 12000)),
                            grepl("P-Tau", LBTEST) ~ round(runif(n(), 15, 60)), TRUE ~ round(runif(n(), 150, 500))),
         SIU = "pg/mL", RPTRESN = SIRESN, RPTRESC = as.character(SIRESN), RPTU = SIU)
labfile <- bind_rows(lab, csf) %>% mutate(VISIT = paste("Visit", VISCODE), VISITTYP = "Scheduled", VISITMOD = NA, LBNAM = "Central", SAMPLE_ID = NA,
  LBDTM_DAYS_CONSENT = LBDTM_DAYS_T0 + 90, RCVDTM_DAYS_CONSENT = NA, RCVDTM_DAYS_T0 = NA, SPECMNUM = NA, FASTSTAT = "N", BATTRID = NA,
  TSTSTAT = "DONE", TSTDTM_DAYS_CONSENT = NA, TSTDTM_DAYS_T0 = NA, TSTTYP = NA, RPTRESNP = NA, RPTNRLO = NA, RPTNRHI = NA, CNVRESC = NA, CNVRESN = NA,
  CNVRESNP = NA, CNVNRLO = NA, CNVNRHI = NA, CNVU = NA, SIRESC = NA, SIRESNP = NA, SINRLO = NA, SINRHI = NA, RPTRTYP = NA, ALRTFL = NA, DELTFL = NA)
write_csv(labfile, file.path(root, "External Data", "clrm_lab.csv"), na = "")

# vitals, mmse, pacc, cdr at every visit in SV (excluding PET/LP visits)
vv <- sv %>% filter(!VISITCD %in% c(2, 5))
vitals <- vv %>% transmute(SUBSTUDY, BID, VISCODE = VISITCD, DONE = 1, NDREASON = NA, VSWEIGHT = round(runif(n(), 55, 100), 1), VSWTUNIT = 2,
                           STDWT = VSWEIGHT, VSHEIGHT = if_else(VISITCD == 1, round(runif(n(), 150, 185), 1), NA_real_), VSHTUNIT = 2, STDHT = VSHEIGHT,
                           VSBPSYS = round(runif(n(), 110, 160)), VSBPDIA = round(runif(n(), 65, 95)), VSPULSE = round(runif(n(), 55, 90)), VSRESP = 16, VSTEMP = 36.6, VSTEMPUN = 2, VSTMPSRC = 1, STDTEMP = 36.6)
write_csv(vitals, file.path(root, "Raw Data", "vitals.csv"), na = "NA")
mmse <- vv %>% transmute(SUBSTUDY, BID, VISCODE = sprintf("%03d", VISITCD), DONE = "Yes", NDREASON = NA, MMSCORE = sample(26:30, n(), TRUE), MMSE = NA)
write_csv(mmse, file.path(root, "Raw Data", "mmse.csv"), na = "NA")
pacc <- vv %>% transmute(SUBSTUDY, BID, VISCODE = VISITCD, PACC.raw = rnorm(n(), 0, 2), PACC = rnorm(n(), 0, 1), PACCIMPFLG = 0)
write_csv(pacc, file.path(root, "Derived Data", "PACC.csv"), na = "NA")
cdr <- vv %>% transmute(SUBSTUDY, BID, VISCODE = sprintf("%03d", VISITCD), DONE = "Yes", CDGLOBAL = sample(c(0, 0, 0, 0.5), n(), TRUE),
                        CDSOB = CDGLOBAL, CDADTC_DAYS_CONSENT = SVSTDTC_DAYS_T0 + 90, CDADTC_DAYS_T0 = SVSTDTC_DAYS_T0, CDRSB = CDGLOBAL)
write_csv(cdr, file.path(root, "Raw Data", "cdr.csv"), na = "NA")
dem <- subj %>% transmute(SUBSTUDY, BID, VISCODE = "001", EXAMDATE_DAYS_CONSENT = 0, EXAMDATE_DAYS_T0 = -ENRLDT_DAYS_CONSENT,
                          PTGENDER = if_else(SEX == 1, 2, 1), PTAGE = AGEYR, PTETHNIC = 2, PTRACE = "5", PTLANG = 1, PTPLANG = 1, PTEDUCAT = EDCCNTU, PTMARRY = 3, PTNOTRT = 1, PTHOME = 1)
write_csv(dem, file.path(root, "Raw Data", "ptdemog.csv"))
ds <- subj %>% transmute(SUBSTUDY, BID, DSSEQ = 3, DSDECOD = "RANDOMIZED", DSCAT = "PROTOCOL MILESTONE", EPOCH = NA, DSSTDTC_DAYS_CONSENT = 69, DSSTDTC_DAYS_T0 = 0, DSSTDY = "0")
write_csv(ds, file.path(root, "Derived Data", "DS.csv"), na = "NA")
cat("fake A4/LEARN written to", root, "\n")
