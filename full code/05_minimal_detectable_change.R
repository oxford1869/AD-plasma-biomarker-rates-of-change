# =============================================================================================
# 05_minimal_detectable_change.R
# Minimal detectable change in each assay's own units.
#
# ADNI. Sections:
#   1. the within-person CV, the mean rate of change (all participants and amyloid-positive), the 95% minimal
#      detectable change between two measurements and for a three-sample slope as fold change and per cent, and
#      the years of mean amyloid-positive change needed to exceed it (T14; formerly 08b6 section 11)
#   2. the minimal detectable change as a rise and as the equivalent fall (R7; formerly 08b17 R7)
# =============================================================================================


# =============================================================================================
# Section 1. Minimal detectable change on the delivered scale (T14)
# Original script: 08b6_paper2_mdc_design.R, setup and section 11
# =============================================================================================
# =============================================================================
# 08b6_paper2_mdc_design.R
# Paper 2, tasks 11 to 13. No new model fits: everything is computed from the
# variance components already written by 08b (T2, T3) and from the frame.
#
#   11. Minimal detectable change converted back to the delivered scale.
#       The change models work on within-measure z of the log value, so a
#       z-scale quantity q corresponds to a log change of q * sd_log(measure)
#       and a fold change of exp(q * sd_log). Reported per platform:
#         mdc95 between two single measurements (fold, and per cent)
#         mdc95 of a three-draw slope at the median span (fold per year)
#         mean annual change (per cent per year) in all people and in amyloid-positive people
#         years of mean change needed to exceed the two-point mdc95
#   12. Trial sample-size illustration. Two-arm trial, plasma slope (or change
#       from baseline) as the endpoint, treatment slows the mean amyloid-positive
#       rate of change by 25% or 50%; alpha 0.05 two-sided, power 80%. Per
#       platform and for the mean of platforms, under 18-month (4 draws) and
#       24-month (3 draws) schedules. The between-person true-slope variance and
#       the residual come from the ADNI change model, so the numbers describe an
#       ADNI-like trial population and are an illustration, not a design.
#   13. Sampling schedules with unequal spacing. Reliability of an individual
#       slope depends on the design only through the spread of sampling times
#       (SST = sum of squared centred times), so schedules can be compared
#       directly: evenly spaced, front-loaded, back-loaded, replicate draws at
#       each end, and the observed ADNI design. Also the span needed to reach
#       reliability 0.7 and 0.9 with 2, 3 or 4 evenly spaced draws.
#   14. Validity by observed SST. People differ in span, so SST varies between
#       them; validity of the plasma slope against concurrent Centiloid change
#       is reported by SST quintile, against the reliability the variance
#       components predict. The criterion is windowed on each participant's own plasma span, so the
#       number of reference scans and the reference window are reported per band as well: if they
#       fall together with SST, the between-person gradient is partly a criterion artefact and the
#       within-participant test in the second cohort is the stronger evidence.
#       components predict at that SST. This is the between-person counterpart
#       of the within-person design test in the A4/LEARN script (08b9).
#
# Outputs T14_mdc_fold_change.csv, T15_trial_sample_size.csv,
#         T16_sampling_schedules.csv, T16b_span_needed.csv, T17_validity_by_sst.csv,
#         F6_sampling_schedules.pdf
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
if (COHORT != "adni") stop("08b6 is the ADNI Paper 2 script; the A4/LEARN equivalents are in 08b9. Set Sys.setenv(COHORT = 'adni').")
suppressPackageStartupMessages(library(ggplot2))
if (!exists("theme_paper")) theme_paper <- function(base_size = 10) ggplot2::theme_minimal(base_size = base_size)
P2_DIR <- file.path(OUT_DIR, "paper2")
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
T2 <- rd("T2_change_variance_components")
T3 <- rd("T3_minimal_detectable_change")
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")

# Scale factor per measure: SD of the log value, which is what z_log divided by.
scale_tab <- long %>% group_by(biomarker, arm, measure, direction) %>%
  summarise(sd_log = sd(log_value, na.rm = TRUE), mean_log = mean(log_value, na.rm = TRUE),
            median_value = median(value, na.rm = TRUE), .groups = "drop")

# Between-person spread on the delivered scale: the coefficient of variation, across
# participants, of a participant's mean value. It is the companion to the within-person
# CV: an assay can agree closely between platforms in absolute terms and still have a low
# intraclass correlation when this is small, which is the Ab42/40 case.
between_tab <- long %>%
  group_by(biomarker, arm, measure, RID) %>%
  summarise(person_mean = mean(value, na.rm = TRUE), .groups = "drop_last") %>%
  summarise(between_cv_pct = 100 * sd(person_mean, na.rm = TRUE) / mean(person_mean, na.rm = TRUE),
            n_people_cv = sum(!is.na(person_mean)), .groups = "drop")

# Observed design: SST per person (same for every platform on the core people)
sst_person <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), span = max(plasma_time) - min(plasma_time),
            sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
median_sst  <- median(sst_person$sst[sst_person$k >= 2])
median_span <- median(sst_person$span[sst_person$k == 3])
say("Observed design: median SST %.1f y^2, median span %.1f y among 3-draw people", median_sst, median_span)

# ---- 11. Minimal detectable change on the delivered scale -------------------------------
say("\n== 11. MDC as fold change ==")
# Mean OLS slope in amyloid-positive people at first plasma (z-log units, aligned)
amy_slope <- slopes %>% filter(n_obs >= 3, amyloid_bl_cl20 %in% TRUE) %>%
  group_by(biomarker, arm) %>% summarise(mean_slope_amy_pos_z = mean(slope, na.rm = TRUE),
                                         n_amy_pos = n(), .groups = "drop")
mdc <- T3 %>% select(biomarker, arm, measure, sd_residual, sd_person_slope, mean_slope,
                     mdc95_two_points, slope_se_3pts, mdc95_slope_per_year,
                     pct_reliable_increase, pct_reliable_decrease) %>%
  left_join(scale_tab, by = c("biomarker", "arm", "measure")) %>%
  left_join(between_tab, by = c("biomarker", "arm", "measure")) %>%
  left_join(amy_slope, by = c("biomarker", "arm")) %>%
  mutate(
    # z-log quantities are direction-aligned (rise with pathology); undo the
    # alignment so per-cent changes read in the delivered direction
    mean_pct_per_year         = 100 * (exp(mean_slope * sd_log * direction) - 1),
    mean_pct_per_year_amy_pos = 100 * (exp(mean_slope_amy_pos_z * sd_log * direction) - 1),
    cv_within_pct             = 100 * sqrt(exp((sd_residual * sd_log)^2) - 1),  # within-person, within-platform CV
    mdc95_two_points_fold     = exp(mdc95_two_points * sd_log),
    mdc95_two_points_pct      = 100 * (mdc95_two_points_fold - 1),
    mdc95_slope_fold_per_year = exp(mdc95_slope_per_year * sd_log),
    mdc95_slope_pct_per_year  = 100 * (mdc95_slope_fold_per_year - 1),
    years_to_exceed_mdc_all     = mdc95_two_points / abs(mean_slope),
    years_to_exceed_mdc_amy_pos = mdc95_two_points / abs(mean_slope_amy_pos_z)) %>%
  select(biomarker, arm, measure, direction, median_value, sd_log, cv_within_pct, between_cv_pct,
         mean_pct_per_year, mean_pct_per_year_amy_pos, n_amy_pos,
         mdc95_two_points_fold, mdc95_two_points_pct, mdc95_slope_fold_per_year, mdc95_slope_pct_per_year,
         years_to_exceed_mdc_all, years_to_exceed_mdc_amy_pos, pct_reliable_increase, pct_reliable_decrease) %>%
  arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(mdc %>% select(biomarker, arm, cv_within_pct, mean_pct_per_year, mean_pct_per_year_amy_pos,
                                   mdc95_two_points_pct, mdc95_slope_pct_per_year, years_to_exceed_mdc_amy_pos) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(mdc, "T14_mdc_fold_change")



# =============================================================================================
# Section 2. Minimal detectable change in both directions (R7)
# Original script: 08b17_paper2_revision_nov2026.R, R7
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
load_change_objects()
T14 <- rd("T14_mdc_fold_change")
# =============================================================================
# R7. MDC in both directions
# =============================================================================
say("\n== R7. MDC95 as a rise and as the equivalent fall ==")
R7 <- T14 %>% transmute(biomarker, arm, measure,
                        mdc95_rise_pct = mdc95_two_points_pct,
                        mdc95_fall_pct = 100 * (1 - 1 / mdc95_two_points_fold)) %>%
  arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(R7 %>% mutate(across(where(is.numeric), ~ round(.x, 1)))))
out(R7, "R7_mdc_both_directions")

