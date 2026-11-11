# =============================================================================================
# 10_tables.R
# Manuscript, supplementary and cross-cohort tables.
#
# Sections:
#   1. ADNI manuscript tables as CSV and markdown: participant characteristics (Table 1), variance components
#      and reliabilities (Table 2), minimal detectable change (Table 3), assay characteristics and the analytical
#      share of the occasion-level variance (supplementary) (formerly 08b8)
#   2. cross-cohort tables C1 to C7 and figure F30 from the finished tables of both cohorts (formerly 10_combine)
# No model fits.
# =============================================================================================


# =============================================================================================
# Section 1. ADNI manuscript tables
# Original script: 08b8_paper2_tables.R
# =============================================================================================
# =============================================================================
# 08b8_paper2_tables.R
# Paper 2, task 15: manuscript tables, built from the frames and the T-tables so
# that every number in the text can be checked against a file. No model fits.
#
#   Table1_participants.csv        participant characteristics, overall and by
#                                  amyloid PET status at first plasma sample
#   Table2_reliability.csv         per biomarker: level variance components,
#                                  G coefficient, slope variance components,
#                                  slope reliability (one platform, mean of
#                                  platforms; all, amyloid-, amyloid+) with
#                                  bootstrap intervals, slope ICC
#   Table3_mdc.csv                 per platform: units, median level, within-
#                                  person CV, mean change per year (all, amyloid+),
#                                  mdc95 two points and slope, years to exceed,
#                                  % with reliable change
#   TableS1_assays.csv             per measure: platform, technology, units, n
#                                  observations, n people, % at LLOQ, batches
#   tables_paper2.md               all four as markdown, for pasting or pandoc
#
# Cells are formatted as text ready for the manuscript (e.g. "0.40 (0.34 to 0.46)").
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
P2_DIR  <- file.path(OUT_DIR, "paper2")
TAB_DIR <- file.path(P2_DIR, "tables"); dir.create(TAB_DIR, showWarnings = FALSE)
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(TAB_DIR, paste0(name, ".csv"))); invisible(df) }
f2  <- function(x, d = 2) ifelse(is.na(x), "", formatC(x, format = "f", digits = d))
f0  <- function(x) ifelse(is.na(x), "", formatC(x, format = "f", digits = 0))
ci  <- function(est, lo, hi, d = 2) ifelse(is.na(lo), f2(est, d), sprintf("%s (%s to %s)", f2(est, d), f2(lo, d), f2(hi, d)))
n_pct <- function(x) sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
mean_sd <- function(x, d = 1) sprintf("%s (%s)", f2(mean(x, na.rm = TRUE), d), f2(sd(x, na.rm = TRUE), d))
med_iqr <- function(x, d = 1) if (all(is.na(x))) "" else sprintf("%s (%s to %s)", f2(median(x, na.rm = TRUE), d),
                                      f2(quantile(x, 0.25, na.rm = TRUE), d), f2(quantile(x, 0.75, na.rm = TRUE), d))

biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
bm_label   <- c(pTau217 = "pTau217", pTau181 = "pTau181", GFAP = "GFAP", NfL = "NfL", Ab42_40 = "Ab42/40")
arm_of <- function(measure, platform) case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N %pTau217", TRUE ~ platform)
arm_levels <- c("ALZpath", "C2N", "C2N %pTau217", "Fujirebio", "Janssen", "Roche", "Quanterix")

long  <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>% mutate(arm = arm_of(measure, platform))
crit  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
part  <- readRDS(file.path(OUT_DIR, "participants.rds"))
# Analysis sample = everyone in the Paper 2 frame (393 on the real data: all
# with plasma on the five full-set platforms; 392 of them have every cell of the
# 3 x 5 design, one has a partial series). This matches the models in 08b.
core_ids <- unique(long$RID)

# ---- Table 1 ------------------------------------------------------------------------------------
say("Table 1")
visits <- long %>% distinct(RID, visit_no, plasma_time) %>% arrange(RID, plasma_time) %>%
  group_by(RID) %>% summarise(n_visits = n(), span = max(plasma_time) - min(plasma_time),
                              gap = if (n() > 1) mean(diff(plasma_time)) else NA_real_, .groups = "drop")
pp <- long %>% distinct(RID, age_first_plasma, PTGENDER, PTEDUCAT, APOE4_carrier, DX_at_first_plasma,
                        amyloid_bl_cl20, centiloid_bl, in_quanterix) %>%
  group_by(RID) %>% slice(1) %>% ungroup() %>%
  filter(RID %in% core_ids) %>%
  left_join(visits, by = "RID") %>%
  left_join(crit %>% select(RID, centiloid_n, csf_ptau_abeta42_n, hcv_adj_n, adas13_n, mpacc_digit_n, mmse_n,
                            amyloid_ever_cl20), by = "RID") %>%
  # conversion among baseline-negative participants, split by whether the first positive scan is dated on or
  # before the last plasma sample on the five full-set platforms (within the sampling period the analyses use;
  # Quanterix-only samples excluded) or after it (added November 2026)
  left_join(readRDS(file.path(OUT_DIR, "master_visit.rds")) %>% filter(RID %in% core_ids, !quanterix_only_visit %in% TRUE) %>%
              group_by(RID) %>% summarise(last_plasma_date = max(visit_date, na.rm = TRUE),
                                          first_positive_date = suppressWarnings(first(na.omit(first_positive_cl20_date))),
                                          .groups = "drop") %>%
              mutate(conv_within = !is.na(first_positive_date) & first_positive_date <= last_plasma_date,
                     conv_after  = !is.na(first_positive_date) & first_positive_date >  last_plasma_date) %>%
              select(RID, conv_within, conv_after), by = "RID") %>%
  mutate(group = case_when(amyloid_bl_cl20 %in% TRUE ~ "Amyloid positive", amyloid_bl_cl20 %in% FALSE ~ "Amyloid negative",
                           TRUE ~ "No PET at first sample"))
row1 <- function(d) tibble(
  `N` = as.character(nrow(d)),
  `Age at first plasma sample, years, mean (SD)` = mean_sd(d$age_first_plasma),
  `Women, n (%)` = n_pct(d$PTGENDER == "Female"),
  `Education, years, mean (SD)` = mean_sd(d$PTEDUCAT),
  `APOE e4 carrier, n (%)` = n_pct(d$APOE4_carrier),
  `Cognitively normal, n (%)` = n_pct(d$DX_at_first_plasma == "CN"),
  `Mild cognitive impairment, n (%)` = n_pct(d$DX_at_first_plasma == "MCI"),
  `Dementia, n (%)` = n_pct(d$DX_at_first_plasma == "Dementia"),
  `Amyloid PET Centiloid at first sample, median (IQR)` = med_iqr(d$centiloid_bl, 0),
  `Became amyloid positive during follow-up, n (%)` = if (all(d$amyloid_bl_cl20 %in% FALSE)) n_pct(d$amyloid_ever_cl20) else "",
  `Crossed the amyloid threshold within the sampling period, n (%)` = if (all(d$amyloid_bl_cl20 %in% FALSE)) n_pct(d$conv_within) else "",
  `Crossed at a later scan, after sampling ended, n (%)` = if (all(d$amyloid_bl_cl20 %in% FALSE)) n_pct(d$conv_after) else "",
  `Plasma samples per participant, n` = med_iqr(d$n_visits, 0),
  `Plasma sampling span, years, median (IQR)` = med_iqr(d$span),
  `Interval between samples, years, median (IQR)` = med_iqr(d$gap),
  `Also measured on Quanterix, n (%)` = n_pct(d$in_quanterix),
  `Rate of change available: amyloid PET, n (%)` = n_pct(coalesce(d$centiloid_n >= 2, FALSE)),
  `Rate of change available: CSF pTau181/Ab42, n (%)` = n_pct(coalesce(d$csf_ptau_abeta42_n >= 2, FALSE)),
  `Rate of change available: hippocampal volume, n (%)` = n_pct(coalesce(d$hcv_adj_n >= 2, FALSE)),
  `Rate of change available: cognition (ADAS-Cog13), n (%)` = n_pct(coalesce(d$adas13_n >= 2, FALSE)))
t1 <- bind_rows(row1(pp) %>% mutate(group = "All participants"),
                pp %>% group_by(group) %>% group_modify(~ row1(.x)) %>% ungroup()) %>%
  pivot_longer(-group, names_to = "Characteristic", values_to = "value") %>%
  pivot_wider(names_from = group, values_from = value) %>%
  select(Characteristic, `All participants`, any_of(c("Amyloid negative", "Amyloid positive", "No PET at first sample")))
print(as.data.frame(t1)); out(t1, "Table1_participants")

# ---- Table 2 ------------------------------------------------------------------------------------
say("Table 2")
T1 <- rd("T1_level_variance_components"); T1b <- rd("T1b_G_coefficients_level")
T2 <- rd("T2_change_variance_components"); T2b <- rd("T2b_change_by_amyloid_status")
T4 <- rd("T4_slope_agreement"); T7 <- rd("T7_bootstrap_intervals")
T2d <- if (file.exists(file.path(P2_DIR, "T2d_platform_slope_LRT.csv"))) rd("T2d_platform_slope_LRT") else NULL
pct_of <- function(bm, comp) { v <- T1$pct[T1$biomarker == bm & T1$component == comp]; if (length(v)) v[1] else NA_real_ }
booted <- function(bm, subset, q) {
  r <- T7 %>% filter(biomarker == bm, subset == !!subset, quantity == q)
  if (nrow(r)) c(r$ci_low[1], r$ci_high[1]) else c(NA_real_, NA_real_)
}
t2 <- map_dfr(biomarkers, function(bm) {
  g  <- T1b %>% filter(biomarker == bm)
  ch <- T2 %>% filter(biomarker == bm)
  st <- T2b %>% filter(str_starts(model, bm))
  neg <- st %>% filter(str_detect(model, "amyloid-")); pos <- st %>% filter(str_detect(model, "amyloid\\+"))
  ic <- T4 %>% filter(biomarker == bm) %>% slice(1)
  b <- function(subset, q) booted(bm, subset, q)
  tibble(
    Biomarker = bm_label[[bm]],
    `Platforms, n` = as.character(ch$n_arms),
    `Level: person, % of variance` = f0(pct_of(bm, "person")),
    `Level: person x platform, %` = f0(pct_of(bm, "person x platform")),
    `Level: sample (shared by platforms), %` = f0(pct_of(bm, "sample (person x occasion)")),
    `Level: batch, %` = f0(pct_of(bm, "batch within platform")),
    `Level: residual, %` = f0(pct_of(bm, "residual")),
    `G coefficient, one measurement (95% CI)` = ci(g$G_relative_one_platform, b("all", "G_relative_one_platform")[1], b("all", "G_relative_one_platform")[2]),
    `True slope SD, SD units per year` = f2(sqrt(ch$var_person_slope), 3),
    `Platform-specific slope SD` = f2(sqrt(ch$var_person_arm_slope), 3),
    `Residual SD` = f2(sqrt(ch$var_residual), 2),
    `Slope reliability, one platform (95% CI)` = ci(ch$rel_slope_one_platform, b("all", "rel_slope_one_platform")[1], b("all", "rel_slope_one_platform")[2]),
    `Slope reliability, mean of platforms (95% CI)` = ci(ch$rel_slope_mean_of_arms, b("all", "rel_slope_mean_of_arms")[1], b("all", "rel_slope_mean_of_arms")[2]),
    `One platform, amyloid negative (95% CI)` = ci(neg$rel_slope_one_platform, b("amyloid-", "rel_slope_one_platform")[1], b("amyloid-", "rel_slope_one_platform")[2]),
    `One platform, amyloid positive (95% CI)` = ci(pos$rel_slope_one_platform, b("amyloid+", "rel_slope_one_platform")[1], b("amyloid+", "rel_slope_one_platform")[2]),
    `Slope ICC, single platform / mean (95% CI)` = sprintf("%s / %s", ci(ic$icc_single_arm, b("all", "icc_single_arm")[1], b("all", "icc_single_arm")[2]),
                                                          ci(ic$icc_mean_of_arms, b("all", "icc_mean_of_arms")[1], b("all", "icc_mean_of_arms")[2])))
})
t2w <- t2 %>% pivot_longer(-Biomarker, names_to = "Quantity", values_to = "v") %>%
  pivot_wider(names_from = Biomarker, values_from = v)
print(as.data.frame(t2w)); out(t2w, "Table2_reliability")

# ---- Table 3 ------------------------------------------------------------------------------------
say("Table 3")
T14 <- rd("T14_mdc_fold_change")
units <- long %>% distinct(measure, platform, arm, biomarker) %>%
  left_join(readRDS(file.path(OUT_DIR, "plasma_long.rds")) %>% distinct(PLASMA_BIOMARKER, unit_std) %>%
              rename(measure = PLASMA_BIOMARKER), by = "measure")
t3 <- T14 %>% mutate(arm = arm_of(measure, arm)) %>%
  left_join(units %>% select(measure, unit_std), by = "measure") %>%
  mutate(Biomarker = bm_label[biomarker], arm = factor(str_replace(arm, "C2N_ratio", "C2N %pTau217"), arm_levels)) %>%
  arrange(factor(biomarker, biomarkers), arm) %>%
  transmute(Biomarker, Platform = as.character(arm), Units = coalesce(unit_std, "ratio"),
            `Median level` = mapply(function(v) f2(v, if (is.na(v)) 2 else if (v < 1) 3 else if (v < 10) 2 else 0), median_value),
            `Within-person CV, %` = f0(cv_within_pct),
            `Mean change per year, %, all` = f2(mean_pct_per_year, 1),
            `Mean change per year, %, amyloid positive` = f2(mean_pct_per_year_amy_pos, 1),
            `MDC95 between two measurements, %` = f0(mdc95_two_points_pct),
            `MDC95 of a three-sample slope, % per year` = f0(mdc95_slope_pct_per_year),
            `Years to exceed MDC95 at the amyloid-positive rate` = ifelse(years_to_exceed_mdc_amy_pos > 100, "> 100", f0(years_to_exceed_mdc_amy_pos)),
            `Participants with slope above MDC95, %` = f0(pct_reliable_increase),
            `Participants with slope below minus MDC95, %` = f0(pct_reliable_decrease))
print(as.data.frame(t3)); out(t3, "Table3_mdc")

# ---- Analytical precision (external QC, Schindler et al. 2024, Alzheimer's & Dementia, Appendix C) --
# CV of pooled plasma controls (high and low) run with every batch of the FNIH
# project. Where a platform reports Ab42 and Ab40 separately, the CV of the
# derived ratio is approximated as sqrt(CV42^2 + CV40^2) (independent errors;
# an upper bound). Fujirebio Ab42/40 has no QC figure in the appendix.
# Replication: C2N, Fujirebio and Roche ran samples in singlicate; ALZpath,
# Janssen and Quanterix in duplicate (Schindler et al. 2024, main text).
QC_CV <- tribble(
  ~measure,                       ~cv_high, ~cv_low, ~replicates,
  "AlzPath_plasma_ptau217",         8.02,    8.93,   "duplicate",
  "C2N_plasma_ptau217",             7.50,    NA,     "singlicate",
  "C2N_plasma_ptau217_ratio",       5.24,    NA,     "singlicate",
  "Fuji_plasma_ptau217",           11.47,   33.62,   "singlicate",
  "Janssen_plasma_ptau217",         5.63,    7.65,   "duplicate",
  "Roche_plasma_ptau181",           4.50,    6.84,   "singlicate",
  "QX_plasma_ptau181",              3.16,    5.97,   "duplicate",
  "Roche_plasma_GFAP",              1.95,    3.38,   "singlicate",
  "QX_plasma_GFAP",                 6.51,    4.89,   "duplicate",
  "Roche_plasma_NfL",               2.09,    3.12,   "singlicate",
  "QX_plasma_NfL",                  5.79,    4.89,   "duplicate",
  "C2N_plasma_Abeta42_Abeta40",     7.33,    7.86,   "singlicate",
  "Fuji_plasma_Ab42_Ab40",          NA,      NA,     "singlicate",
  "Roche_plasma_Ab42_Ab40",  sqrt(4.77^2 + 5.50^2), sqrt(3.24^2 + 7.57^2), "singlicate",
  "QX_plasma_Ab42_Ab40",     sqrt(4.08^2 + 4.20^2), sqrt(5.05^2 + 4.07^2), "duplicate") %>%
  mutate(cv_analytical = rowMeans(cbind(cv_high, cv_low), na.rm = TRUE),
         cv_analytical = if_else(is.nan(cv_analytical), NA_real_, cv_analytical),
         # Fujirebio pTau217: the low control sat near the assay's lower limit
         # (CV 34%); the high control (11%) is used as the analytical CV for the
         # split, and the low-control figure is reported alongside
         cv_analytical = if_else(measure == "Fuji_plasma_ptau217", cv_high, cv_analytical),
         cv_reported = case_when(is.na(cv_high) ~ "not reported",
                                 is.na(cv_low) ~ sprintf("%.1f (high control only)", cv_high),
                                 TRUE ~ sprintf("%.1f / %.1f", cv_high, cv_low)))

# ---- Table S2: analytical share of the occasion-level variance ------------------------------------
# Within-person CV (Table 3) is everything about one measurement that is not the
# participant's linear trajectory. On the log scale variances add, so the share
# attributable to the assay itself is (analytical CV / within-person CV)^2, with
# the remainder pre-analytical, short-term biological and non-linearity.
say("Table S2")
T14 <- rd("T14_mdc_fold_change")
ts2 <- T14 %>% select(biomarker, arm, measure, cv_within_pct) %>%
  left_join(QC_CV, by = "measure") %>%
  mutate(arm = arm_of(measure, arm), Biomarker = bm_label[biomarker],
         arm = factor(str_replace(arm, "C2N_ratio", "C2N %pTau217"), arm_levels),
         share_analytical = pmin(100, 100 * (cv_analytical / cv_within_pct)^2),
         cv_other = sqrt(pmax(0, cv_within_pct^2 - cv_analytical^2))) %>%
  arrange(factor(biomarker, biomarkers), arm) %>%
  transmute(Biomarker, Platform = as.character(arm), `Replicates` = replicates,
            `QC CV, % (high / low control)` = cv_reported,
            `Analytical CV used, %` = f2(cv_analytical, 1),
            `Within-person CV, %` = f0(cv_within_pct),
            `Non-analytical CV, %` = f0(cv_other),
            `Analytical share of occasion variance, %` = f0(share_analytical))
print(as.data.frame(ts2)); out(ts2, "TableS2_analytical_share")

# ---- Supplementary Table S1: assays --------------------------------------------------------------
say("Table S1")
pl <- readRDS(file.path(OUT_DIR, "plasma_long.rds")) %>% filter(role == "measure", RID %in% core_ids)
tech <- c(C2N = "Liquid chromatography tandem mass spectrometry", Fujirebio = "Chemiluminescent enzyme immunoassay (Lumipulse G)",
          ALZpath = "Single molecule array (Simoa), ALZpath antibody", Janssen = "Single molecule array (Simoa)",
          Roche = "Electrochemiluminescence immunoassay (Elecsys)", Quanterix = "Single molecule array (Simoa)")
ts1 <- pl %>% group_by(biomarker, platform, PLASMA_BIOMARKER, unit_std) %>%
  summarise(n_obs = n(), n_people = n_distinct(RID), pct_lloq = 100 * mean(lloq_flag, na.rm = TRUE),
            n_batches = n_distinct(RUN_NO), .groups = "drop") %>%
  mutate(arm = arm_of(PLASMA_BIOMARKER, platform), Technology = tech[platform]) %>%
  left_join(QC_CV %>% select(measure, replicates, cv_reported), by = c("PLASMA_BIOMARKER" = "measure")) %>%
  arrange(factor(biomarker, biomarkers), factor(arm, arm_levels)) %>%
  transmute(Biomarker = bm_label[biomarker], Platform = arm, Technology, Units = coalesce(unit_std, "ratio"),
            Replicates = replicates,
            `Observations, n` = n_obs, `Participants, n` = n_people,
            `At lower limit of quantification, %` = f2(pct_lloq, 1), `Analytical batches, n` = n_batches,
            `QC CV, % (high / low control)` = cv_reported)
print(as.data.frame(ts1)); out(ts1, "TableS1_assays")

# ---- markdown -----------------------------------------------------------------------------------
md_table <- function(df, title) {
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  c(paste0("## ", title), "", hdr, sep, rows, "")
}
writeLines(c(sprintf("# Paper 2 tables (pipeline %s)", PIPELINE_VERSION), "",
             md_table(t1, "Table 1. Participant characteristics"),
             md_table(t2w, "Table 2. Sources of variance and reliability per biomarker"),
             md_table(t3, "Table 3. Minimal detectable change per platform on the delivered scale"),
             md_table(ts1, "Table S1. Assays"),
             md_table(ts2, "Table S2. Analytical share of the occasion-level variance")),
           file.path(TAB_DIR, "tables_paper2.md"))
say("\nTables written to %s", TAB_DIR)


# =============================================================================================
# Section 2. Cross-cohort tables (C1 to C7, F30)
# Original script: 10_combine.R
# =============================================================================================
# =============================================================================
# 10_combine.R
# Cross-cohort tables and figures. Reads the finished outputs of BOTH cohorts,
# pipeline_out/ (ADNI) and a4learn_out/ (A4/LEARN), and writes combined/ .
# It fits no models: everything here is assembled from the tables the analysis
# scripts have already written, so it is quick and can be re-run at any time.
#
# What "unit" means differs by cohort, and every combined table carries the
# distinction in its `units_are` column:
#   ADNI       one sample, five pTau217 assays on the same blood draws, so a
#              comparison across units is a comparison across assays.
#   A4/LEARN   one assay (Lilly MSD), several subsamples that differ by design in
#              amyloid burden (A4 placebo, LEARN, A4 solanezumab and the combined
#              placebo + LEARN and all-arms samples), so a comparison across units
#              is replication across populations.
# Nothing here averages the two cohorts; they are placed side by side.
#
# Paper 2 outputs (combined/)
#   C1_design.csv              sample, draws, span and SST per cohort and unit
#   C2_reliability.csv         one-platform slope reliability with bootstrap intervals
#   C3_measurement_error.csv   A4 12-week test-retest beside each cohort's model residual
#   C4_validity.csv            slope against concurrent criterion change, both cohorts
#   C5_mdc.csv                 minimal detectable change on the delivered scale
#   C6_span_needed.csv         span needed for reliability 0.7, 0.8, 0.9
#   C7_validity_by_sst.csv     ADNI SST quintiles (T17) and A4 design subsets (T8, T8b)
#   F30_reliability_vs_sst.pdf reliability against SST, both cohorts, observed designs marked
#
# Paper 2 reads no Paper 3 output. The Paper 3 tables (C10 to C14) and figures (F31 to F33)
# are built by Paper 3's own copy of this script (removed from Paper 2 on 27 September 2026).
#
# Run     source("10_combine.R") after the earlier files have run for
#         both cohorts (it is the last step of that driver). Missing inputs are skipped with
#         a message rather than stopping the run, so a partial set still works.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(ggplot2))

ADNI_OUT <- file.path(OUT_ROOT, "pipeline_out")
A4_OUT   <- file.path(OUT_ROOT, "a4learn_out")
COMB_DIR <- file.path(OUT_ROOT, "combined"); dir.create(COMB_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(COMB_DIR, paste0(name, ".csv"))); invisible(df) }

# read a table from one cohort's paper folder; NULL (with a note) if it is not there
rd <- function(cohort, paper, name) {
  f <- file.path(if (cohort == "ADNI") ADNI_OUT else A4_OUT, paper, paste0(name, ".csv"))
  if (!file.exists(f)) { say("  missing (skipped): %s %s/%s.csv", cohort, paper, name); return(NULL) }
  read_csv(f, show_col_types = FALSE)
}
UNITS <- c(ADNI = "five pTau217 assays, same people",
           `A4/LEARN` = "subsamples of one assay, different people")
say("Combining ADNI (%s) and A4/LEARN (%s)", ADNI_OUT, A4_OUT)

# =============================== Paper 2 =====================================
say("\n== Paper 2 ==")
a_T2  <- rd("ADNI", "paper2", "T2_change_variance_components")
b_T2  <- rd("A4",   "paper2", "T2_change_variance_components")
a_T7  <- rd("ADNI", "paper2", "T7_bootstrap_intervals")
b_T9  <- rd("A4",   "paper2", "T9_bootstrap_ci")
b_T7  <- rd("A4",   "paper2", "T7_test_retest")

# ---- C1 design ---------------------------------------------------------------------
design <- bind_rows(
  if (!is.null(a_T2)) a_T2 %>% transmute(cohort = "ADNI", unit = biomarker, floor = "as_delivered",
    n_people, n_platforms = n_arms, sst_median, median_span = NA_real_),
  if (!is.null(b_T2)) b_T2 %>% filter(floor == "as_delivered") %>%
    transmute(cohort = "A4/LEARN", unit, floor, n_people, n_platforms = 1L, sst_median,
              median_span = median_span_3plus))
if (nrow(design)) {
  design <- design %>% mutate(units_are = UNITS[cohort])
  out(design, "C1_design"); print(as.data.frame(design %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
}

# ---- C2 reliability ------------------------------------------------------------------
ci_adni <- if (!is.null(a_T7)) a_T7 %>% filter(quantity == "rel_slope_one_platform", subset == "all") %>%
  transmute(unit = biomarker, lo95 = ci_low, hi95 = ci_high, n_draws) else NULL
ci_a4 <- if (!is.null(b_T9)) b_T9 %>% filter(quantity == "rel_slope_one_platform") %>%
  transmute(unit, lo95, hi95, n_draws) else NULL
rel <- bind_rows(
  if (!is.null(a_T2)) a_T2 %>% transmute(cohort = "ADNI", unit = biomarker, sst_median,
      var_true_slope = var_person_slope, var_residual,
      reliability_one_platform = rel_slope_one_platform,
      reliability_mean_of_platforms = rel_slope_mean_of_arms) %>%
    { if (is.null(ci_adni)) . else left_join(., ci_adni, by = "unit") },
  if (!is.null(b_T2)) b_T2 %>% filter(floor == "as_delivered") %>%
    transmute(cohort = "A4/LEARN", unit, sst_median, var_true_slope = var_person_slope, var_residual,
              reliability_one_platform = rel_slope_one_platform,
              reliability_mean_of_platforms = NA_real_) %>%
    { if (is.null(ci_a4)) . else left_join(., ci_a4, by = "unit") })
if (nrow(rel)) {
  rel <- rel %>% mutate(units_are = UNITS[cohort])
  out(rel, "C2_reliability")
  say("Slope reliability at the observed design:")
  print(as.data.frame(rel %>% select(cohort, unit, sst_median, var_true_slope, var_residual,
                                     reliability_one_platform, any_of(c("lo95", "hi95"))) %>%
                        mutate(across(where(is.numeric), ~ round(.x, 3)))))
}

# ---- C3 measurement error -------------------------------------------------------------
# The A4 12-week pair measures within-person variation directly; in both cohorts the change
# model's residual estimates the same quantity from years of follow-up. Placing them in one
# table is the check that the decomposition's residual is measurement plus short-term
# biological variation rather than misfit of the linear trajectory.
err <- bind_rows(
  if (!is.null(b_T7)) b_T7 %>% transmute(cohort = "A4/LEARN", unit = sample,
      source = sprintf("repeat plasma sample %.0f days apart", gap_days_median), n = n_pairs,
      sd_within_z = sd_within_z, cv_within_pct, icc = icc_2_1,
      mdc95_two_points_pct, pct_beyond_mdc, ratio_to_change_model = 1 / ratio_change_model_to_test_retest),
  if (!is.null(b_T2)) b_T2 %>% filter(floor == "as_delivered") %>%
    transmute(cohort = "A4/LEARN", unit, source = "change model residual", n = n_people,
              sd_within_z = sd_residual, cv_within_pct = NA_real_, icc = NA_real_,
              mdc95_two_points_pct = NA_real_, pct_beyond_mdc = NA_real_, ratio_to_change_model = NA_real_),
  if (!is.null(a_T2)) a_T2 %>% transmute(cohort = "ADNI", unit = biomarker, source = "change model residual",
      n = n_people, sd_within_z = sqrt(var_residual), cv_within_pct = NA_real_, icc = NA_real_,
      mdc95_two_points_pct = NA_real_, pct_beyond_mdc = NA_real_, ratio_to_change_model = NA_real_))
if (nrow(err)) { out(err, "C3_measurement_error")
  print(as.data.frame(err %>% filter(source != "change model residual" | unit %in% c("pTau217", "A4 placebo")) %>%
                        mutate(across(where(is.numeric), ~ round(.x, 3))))) }

# ---- C4 validity ------------------------------------------------------------------------
a_T5 <- rd("ADNI", "paper2", "T5_validity_vs_criterion_change")
b_T5 <- rd("A4",   "paper2", "T5_validity_vs_criterion_change")
COMMON_CRIT <- c("centiloid_slope", "mmse_slope", "csf_ptau_abeta42_slope")
val <- bind_rows(
  if (!is.null(a_T5)) a_T5 %>% filter(biomarker == "pTau217") %>%
    transmute(cohort = "ADNI", unit = arm, criterion, estimate, n, r),
  if (!is.null(b_T5)) b_T5 %>% transmute(cohort = "A4/LEARN", unit, criterion, estimate, n, r))
if (nrow(val)) {
  val <- val %>% mutate(in_both_cohorts = criterion %in% COMMON_CRIT, units_are = UNITS[cohort]) %>%
    arrange(criterion, cohort, unit)
  out(val, "C4_validity")
  say("pTau217 slope against concurrent Centiloid change:")
  print(as.data.frame(val %>% filter(criterion == "centiloid_slope") %>% select(cohort, unit, estimate, n, r)))
}

# ---- C5 minimal detectable change ---------------------------------------------------------
a_T14 <- rd("ADNI", "paper2", "T14_mdc_fold_change")
b_T14 <- rd("A4",   "paper2", "T14_mdc_fold_change")
mdc <- bind_rows(
  if (!is.null(a_T14)) a_T14 %>% filter(biomarker == "pTau217") %>%
    transmute(cohort = "ADNI", unit = arm, cv_within_pct, mean_pct_per_year, mean_pct_per_year_amy_pos,
              mdc95_two_points_pct, mdc95_slope_pct_per_year, years_to_exceed_mdc_amy_pos,
              pct_reliable_increase),
  if (!is.null(b_T14)) b_T14 %>% transmute(cohort = "A4/LEARN", unit, cv_within_pct, mean_pct_per_year,
              mean_pct_per_year_amy_pos, mdc95_two_points_pct, mdc95_slope_pct_per_year,
              years_to_exceed_mdc_amy_pos, pct_reliable_increase))
if (nrow(mdc)) { mdc <- mdc %>% mutate(units_are = UNITS[cohort]); out(mdc, "C5_mdc")
  print(as.data.frame(mdc %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))) }

# ---- C6 span needed ------------------------------------------------------------------------
a_T16b <- rd("ADNI", "paper2", "T16b_span_needed")
b_T16b <- rd("A4",   "paper2", "T16b_span_needed")
span <- bind_rows(
  if (!is.null(a_T16b)) a_T16b %>% transmute(cohort = "ADNI", unit = biomarker, target, k,
      span_needed_one_platform, span_needed_mean_of_platforms),
  if (!is.null(b_T16b)) b_T16b %>% transmute(cohort = "A4/LEARN", unit, target, k,
      span_needed_one_platform = span_needed_years, span_needed_mean_of_platforms = NA_real_))
if (nrow(span)) { out(span, "C6_span_needed")
  print(as.data.frame(span %>% filter(k == 3, target == 0.7) %>%
                        mutate(across(where(is.numeric), ~ round(.x, 1))))) }

# ---- C7 validity by SST -------------------------------------------------------------------
# ADNI varies SST between people (quintiles); A4 varies it within the same people (subsets of
# their own draws). The two together are the between- and within-person versions of one test.
a_T17 <- rd("ADNI", "paper2", "T17_validity_by_sst")
b_T8  <- rd("A4",   "paper2", "T8_design_subsets")
b_T8b <- rd("A4",   "paper2", "T8b_validity_by_sst")
vs <- bind_rows(
  if (!is.null(a_T17)) a_T17 %>% filter(biomarker == "pTau217") %>%
    transmute(cohort = "ADNI", design = "between people (SST quintile)", unit = biomarker,
              level = paste("quintile", band), n = n_people, sst_median, reliability_predicted,
              r_centiloid_spearman, validity_ratio = validity_ratio_to_top_band,
              sqrt_reliability_ratio = sqrt_reliability_ratio_to_top_band),
  if (!is.null(b_T8)) b_T8 %>% transmute(cohort = "A4/LEARN", design = "within person (own plasma samples)",
              unit = sample, level = as.character(subset), n, sst_median, reliability_predicted,
              r_centiloid_spearman, validity_ratio = validity_ratio_to_full,
              sqrt_reliability_ratio = sqrt_reliability_ratio_to_full),
  if (!is.null(b_T8b)) b_T8b %>% transmute(cohort = "A4/LEARN", design = "between people (SST band)",
              unit, level = as.character(band), n, sst_median, reliability_predicted,
              r_centiloid_spearman, validity_ratio = NA_real_, sqrt_reliability_ratio = NA_real_))
if (nrow(vs)) { out(vs, "C7_validity_by_sst")
  print(as.data.frame(vs %>% filter(design != "between people (SST band)") %>%
                        mutate(across(where(is.numeric), ~ round(.x, 2))))) }

# ---- F30 reliability against SST ------------------------------------------------------------
INK <- "#1c1c1c"; INK2 <- "#5a5751"; MUTED <- "#8d8a83"; GRID <- "#e6e3dc"
if (nrow(rel)) {
  curves <- rel %>% filter(!is.na(var_true_slope), var_true_slope > 0) %>%
    select(cohort, unit, var_true_slope, var_residual) %>% crossing(sst = seq(0.5, 60, by = 0.5)) %>%
    mutate(reliability = var_true_slope / (var_true_slope + var_residual / sst))
  pts <- rel %>% filter(!is.na(var_true_slope), var_true_slope > 0) %>%
    transmute(cohort, unit, sst = sst_median, reliability = reliability_one_platform)
  sub_pts <- if (!is.null(b_T8)) b_T8 %>% filter(sample == "A4 placebo") %>%
    transmute(cohort = "A4/LEARN", unit = "design subsets (A4 placebo)", sst = sst_median,
              reliability = reliability_predicted, label = subset) else NULL
  ends <- curves %>% group_by(cohort, unit) %>% filter(sst == max(sst)) %>% ungroup()
  # internal biomarker keys are fine in the CSVs; the figure gets the printed names
  pretty_unit <- function(x) sub("^Ab42_40$", "Ab42/40", as.character(x))
  ends <- ends %>% mutate(unit_lab = pretty_unit(unit))
  F30 <- ggplot(curves, aes(sst, reliability, colour = unit)) +
    geom_hline(yintercept = c(0.7, 0.9), linetype = "dotted", colour = MUTED) +
    geom_line(linewidth = 0.5) +
    geom_point(data = pts, size = 1.8) +
    { if (requireNamespace("ggrepel", quietly = TRUE))
        ggrepel::geom_text_repel(data = ends, aes(label = unit_lab), size = 2.5, hjust = 0, direction = "y",
                                 nudge_x = 3, segment.colour = MUTED, show.legend = FALSE, seed = 1,
                                 bg.colour = "white", bg.r = 0.14, box.padding = 0.2)
      else geom_text(data = ends, aes(label = unit_lab), size = 2.5, hjust = 1, vjust = -0.6, show.legend = FALSE) } +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.28))) +
    guides(colour = "none") +
    facet_wrap(~ cohort) +
    { if (!is.null(sub_pts)) geom_point(data = sub_pts %>% mutate(unit = "A4 placebo"), shape = 4, size = 2, colour = INK) } +
    labs(x = "Spread of sampling times, sum of squared centred years (SST)",
         y = "Reliability of a one-platform slope", colour = NULL,
         title = "Reliability of an individual slope depends on the design only through SST",
         subtitle = "Points: observed design. Crosses (A4/LEARN): subsets of the same people's own plasma samples.") +
    theme(plot.subtitle = element_text(colour = INK2, size = 7.5))
  ggsave(file.path(COMB_DIR, "F30_reliability_vs_sst.pdf"), F30, width = 9, height = 4.5, device = "pdf", useDingbats = FALSE, bg = "white")
  ggsave(file.path(COMB_DIR, "F30_reliability_vs_sst.png"), F30, width = 9, height = 4.5, dpi = 300, bg = "white")
}

say("\nCombined tables and figures written to %s", COMB_DIR)
