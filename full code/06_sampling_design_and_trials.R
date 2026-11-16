# =============================================================================================
# 06_sampling_design_and_trials.R
# Consequences for sampling design and for trials.
#
# ADNI. Sections:
#   1. trial sample sizes per platform and for the mean of platforms (T15), reliability under alternative
#      sampling schedules and the span needed for target reliabilities (T16, T16b, F6), validity by observed
#      spread of sampling times (T17) (formerly 08b6 sections 12 to 14)
#   2. trial sizes and spans for the mean of platforms under the upper and lower bound (R2; formerly 08b17 R2)
#   3. trial sizes on amyloid-positive variance components (R3; formerly 08b17 R3)
#   4. validity by SST fifth allowing for the precision of both the plasma and the Centiloid slope (R14b)
# =============================================================================================


# =============================================================================================
# Section 1. Trial sample size, sampling schedules, validity by SST (T15 to T17)
# Original script: 08b6_paper2_mdc_design.R, setup and sections 12 to 14
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

# Mean OLS slope in amyloid-positive people at first plasma (z-log units, aligned)
amy_slope <- slopes %>% filter(n_obs >= 3, amyloid_bl_cl20 %in% TRUE) %>%
  group_by(biomarker, arm) %>% summarise(mean_slope_amy_pos_z = mean(slope, na.rm = TRUE),
                                         n_amy_pos = n(), .groups = "drop")

# ---- 12. Trial sample size ---------------------------------------------------------------
say("\n== 12. Trial sample-size illustration ==")
# Per-person slope over a trial schedule with k draws over S years on one platform:
#   var(observed slope) = var_true_slope + sd_residual_arm^2 / SST
# Change from baseline at S years (two draws):
#   var(change) = var_true_slope * S^2 + 2 * sd_residual_arm^2
# Treatment effect = fraction f of the mean amyloid-positive slope.
# n per arm = 2 * (z_{0.975} + z_{0.80})^2 * var / delta^2
zsum2 <- (qnorm(0.975) + qnorm(0.80))^2
designs <- tribble(~schedule, ~times,
                   "18 months, 4 draws (0, 6, 12, 18 m)", c(0, 0.5, 1, 1.5),
                   "24 months, 3 draws (0, 12, 24 m)",    c(0, 1, 2),
                   "24 months, 5 draws (0, 6, 12, 18, 24 m)", c(0, 0.5, 1, 1.5, 2)) %>%
  mutate(k = map_int(times, length), span = map_dbl(times, max),
         sst = map_dbl(times, ~ sum((.x - mean(.x))^2)))
arm_comp <- T3 %>% select(biomarker, arm, sd_residual) %>%
  left_join(T2 %>% select(biomarker, var_person_slope, n_arms), by = "biomarker") %>%
  left_join(amy_slope, by = c("biomarker", "arm"))
# mean-of-platforms row per biomarker: residual variance averaged across arms, divided by k arms
mean_comp <- arm_comp %>% group_by(biomarker) %>%
  summarise(arm = "mean of platforms", sd_residual = sqrt(mean(sd_residual^2) / first(n_arms)),
            var_person_slope = first(var_person_slope), n_arms = first(n_arms),
            mean_slope_amy_pos_z = mean(mean_slope_amy_pos_z, na.rm = TRUE), n_amy_pos = first(n_amy_pos),
            .groups = "drop")
trial <- bind_rows(arm_comp, mean_comp) %>%
  crossing(designs %>% select(-times)) %>%
  crossing(effect = c(0.25, 0.50)) %>%
  mutate(var_slope_obs   = var_person_slope + sd_residual^2 / sst,
         delta_slope     = effect * abs(mean_slope_amy_pos_z),
         n_per_arm_slope = ceiling(2 * zsum2 * var_slope_obs / delta_slope^2),
         var_change      = var_person_slope * span^2 + 2 * sd_residual^2,
         delta_change    = effect * abs(mean_slope_amy_pos_z) * span,
         n_per_arm_change_from_baseline = ceiling(2 * zsum2 * var_change / delta_change^2),
         reliability_of_trial_slope = var_person_slope / var_slope_obs) %>%
  select(biomarker, arm, schedule, k, span, effect, mean_slope_amy_pos_z, reliability_of_trial_slope,
         n_per_arm_slope, n_per_arm_change_from_baseline) %>%
  arrange(factor(biomarker, biomarkers), schedule, effect, arm)
print(as.data.frame(trial %>% filter(biomarker == "pTau217", effect == 0.25) %>%
                      select(arm, schedule, reliability_of_trial_slope, n_per_arm_slope, n_per_arm_change_from_baseline) %>%
                      mutate(reliability_of_trial_slope = round(reliability_of_trial_slope, 2))))
out(trial, "T15_trial_sample_size")

# ---- 13. Sampling schedules ------------------------------------------------------------
say("\n== 13. Sampling schedules ==")
schedules <- tribble(~schedule, ~times,
  "2 draws, 2 y apart",                 c(0, 2),
  "2 draws, 4 y apart",                 c(0, 4),
  "2 draws, 6 y apart",                 c(0, 6),
  "3 draws evenly over 4 y",            c(0, 2, 4),
  "3 draws front-loaded (0, 0.5, 4 y)", c(0, 0.5, 4),
  "3 draws back-loaded (0, 3.5, 4 y)",  c(0, 3.5, 4),
  "3 draws evenly over 6 y",            c(0, 3, 6),
  "4 draws evenly over 4 y",            c(0, 4/3, 8/3, 4),
  "4 draws: replicate at each end, 4 y",c(0, 0, 4, 4),
  "5 draws annually over 4 y",          c(0, 1, 2, 3, 4),
  "6 draws: three at each end, 4 y",    c(0, 0, 0, 4, 4, 4),
  "ADNI observed (median SST)",         NA) %>%
  mutate(k = map_int(times, ~ if (all(is.na(.x))) 3L else length(.x)),
         span = map_dbl(times, ~ if (all(is.na(.x))) median_span else max(.x)),
         sst = map_dbl(times, ~ if (all(is.na(.x))) median_sst else sum((.x - mean(.x))^2)))
sched <- T2 %>% select(biomarker, var_person_slope, var_residual, n_arms) %>%
  crossing(schedules %>% select(-times)) %>%
  mutate(rel_one_platform = var_person_slope / (var_person_slope + var_residual / sst),
         rel_mean_of_platforms = var_person_slope / (var_person_slope + var_residual / (sst * n_arms))) %>%
  select(biomarker, schedule, k, span, sst, rel_one_platform, rel_mean_of_platforms) %>%
  arrange(factor(biomarker, biomarkers), sst)
sched_wide <- sched %>% select(biomarker, schedule, k, span, sst, rel_one_platform) %>%
  pivot_wider(names_from = biomarker, values_from = rel_one_platform) %>% arrange(sst)
print(as.data.frame(sched_wide %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(sched, "T16_sampling_schedules")

# Span needed for a target reliability with k evenly spaced draws on one platform.
# SST for k even draws over span S is S^2 * k (k + 1) / (12 (k - 1)).
# Required SST = (var_residual / var_true) * rel / (1 - rel).
need <- T2 %>% select(biomarker, var_person_slope, var_residual, n_arms) %>%
  crossing(target = c(0.7, 0.8, 0.9), k = c(2, 3, 4, 6)) %>%
  mutate(sst_needed_one = (var_residual / var_person_slope) * target / (1 - target),
         sst_needed_mean = sst_needed_one / n_arms,
         span_needed_one_platform = sqrt(sst_needed_one * 12 * (k - 1) / (k * (k + 1))),
         span_needed_mean_of_platforms = sqrt(sst_needed_mean * 12 * (k - 1) / (k * (k + 1)))) %>%
  select(biomarker, target, k, span_needed_one_platform, span_needed_mean_of_platforms) %>%
  arrange(factor(biomarker, biomarkers), target, k)
print(as.data.frame(need %>% filter(k == 3) %>% mutate(across(starts_with("span"), ~ round(.x, 1)))))
out(need, "T16b_span_needed")

# Figure: reliability of a one-platform slope against SST, schedules marked
curve_df <- T2 %>% select(biomarker, var_person_slope, var_residual) %>%
  crossing(sst = seq(0.5, 60, by = 0.5)) %>%
  mutate(rel = var_person_slope / (var_person_slope + var_residual / sst),
         biomarker = factor(biomarker, biomarkers))
pts <- sched %>% filter(biomarker == "pTau217") %>% mutate(biomarker = factor(biomarker, biomarkers))
F6 <- ggplot(curve_df, aes(sst, rel, colour = biomarker)) +
  geom_line() +
  geom_point(data = pts, aes(sst, rel_one_platform), colour = "black", size = 1.8) +
  geom_text(data = pts, aes(sst, rel_one_platform, label = schedule), colour = "black",
            size = 2.4, hjust = -0.05, vjust = -0.4, check_overlap = TRUE) +
  geom_hline(yintercept = c(0.7, 0.9), linetype = "dotted") +
  labs(x = "Spread of sampling times, sum of squared centred years (SST)",
       y = "Reliability of a one-platform slope", colour = NULL,
       title = "Reliability depends on the design only through SST (points: pTau217)") +
  theme_paper(base_size = 10)
ggsave(file.path(P2_DIR, "F6_sampling_schedules.pdf"), F6, width = 8, height = 5.5)

# ---- 14. Validity by observed spread of sampling times -------------------------------------
# Between-person counterpart of the A4/LEARN within-person design test: ADNI spans vary
# between people (median 5.0 y, IQR 4.0 to 6.5), so SST varies too. If the reliability
# formula describes these data, validity against concurrent Centiloid change should rise
# with observed SST, and roughly in proportion to the square root of predicted reliability.
# Weaker than the A4 test, because span is not randomly assigned: people with longer spans
# stayed in the study longer, so this is reported as a supporting check, not a design test.
say("\n== 14. Validity by observed SST ==")
critr <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
sst_bands <- slopes %>% filter(n_obs >= 2) %>%
  select(RID, biomarker, arm, slope) %>%
  left_join(sst_person %>% select(RID, k, span, sst), by = "RID") %>%
  left_join(critr %>% select(RID, centiloid_slope, centiloid_n, centiloid_span), by = "RID") %>%
  filter(!is.na(slope), !is.na(centiloid_slope), !is.na(sst)) %>%
  group_by(biomarker) %>% mutate(band = ntile(sst, 5)) %>% ungroup()
T17 <- sst_bands %>%
  group_by(biomarker, band) %>%
  summarise(n_person_platform = n(), n_people = n_distinct(RID),
            sst_median = median(sst), span_median = median(span), k_median = median(k),
            # the criterion is measured over the participant's own plasma window plus 180 days at each end,
            # so a participant sampled over a short span also has a short window for the reference slope and
            # a noisier criterion. That pushes the observed validity down in the low bands for a reason that
            # has nothing to do with the plasma design, which is what these two columns expose.
            criterion_n_median = median(centiloid_n, na.rm = TRUE),
            criterion_span_median = median(centiloid_span, na.rm = TRUE),
            r_centiloid_spearman = cor(slope, centiloid_slope, method = "spearman"), .groups = "drop") %>%
  left_join(T2 %>% select(biomarker, var_person_slope, var_residual), by = "biomarker") %>%
  mutate(reliability_predicted = var_person_slope / (var_person_slope + var_residual / sst_median)) %>%
  group_by(biomarker) %>%
  mutate(validity_ratio_to_top_band = r_centiloid_spearman / r_centiloid_spearman[band == max(band)],
         sqrt_reliability_ratio_to_top_band = sqrt(reliability_predicted / reliability_predicted[band == max(band)])) %>%
  ungroup() %>%
  select(biomarker, band, n_person_platform, n_people, k_median, span_median, sst_median,
         criterion_n_median, criterion_span_median,
         reliability_predicted, r_centiloid_spearman, validity_ratio_to_top_band, sqrt_reliability_ratio_to_top_band)
print(as.data.frame(T17 %>% filter(biomarker == "pTau217") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(T17, "T17_validity_by_sst")
say("Criterion measurement by SST band (pTau217): if the low bands also have fewer scans over a shorter")
say("window, part of the fall in validity is the criterion, not the plasma design.")
print(as.data.frame(T17 %>% filter(biomarker == "pTau217") %>%
        select(band, sst_median, span_median, criterion_n_median, criterion_span_median,
               r_centiloid_spearman, validity_ratio_to_top_band, sqrt_reliability_ratio_to_top_band) %>%
        mutate(across(where(is.numeric), ~ round(.x, 2)))))

say("\nTables T14 to T17 and F6 written to %s", P2_DIR)


# =============================================================================================
# Section 2. Mean of platforms under both bounds: trial sizes and spans (R2)
# Original script: 08b17_paper2_revision_nov2026.R, R2
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
load_change_objects()
T2  <- rd("T2_change_variance_components")
T3  <- rd("T3_minimal_detectable_change")
T15 <- rd("T15_trial_sample_size")
# =============================================================================
# R2. Trial sizes and spans for the mean of platforms, upper and lower bound
# =============================================================================
say("\n== R2. Mean-of-platforms trial sizes and spans, upper and lower bound ==")
# Upper bound exactly as 08b6: T2 true-slope variance, mean of per-platform T3 residual
# variances divided by the number of platforms.
upper_comp <- T3 %>% select(biomarker, arm, sd_residual) %>%
  left_join(T2 %>% select(biomarker, var_person_slope, n_arms), by = "biomarker") %>%
  group_by(biomarker) %>%
  summarise(var_true = first(var_person_slope),
            var_err_mean = mean(sd_residual^2) / first(n_arms), .groups = "drop") %>%
  mutate(bound = "upper (paper model)", identifiable = TRUE)
# Lower bound: the model with a person-by-sample effect, refitted here unrounded (T20 is
# rounded to 4 decimals); shared part kept, assay-specific part averaged.
lower_full <- map_dfr(biomarkers, function(bm) {
  d  <- long %>% filter(biomarker == bm) %>% mutate(arm_f = droplevels(factor(arm)))
  fs <- lmer(f_sample, data = d, REML = TRUE, control = ctrl)
  tibble(biomarker = bm, n_platforms = n_distinct(d$arm),
         var_true = vget(fs, "RID_f", "plasma_time"), var_shared = vget(fs, "sample_f"),
         var_assay = vget(fs, "Residual"),
         identifiable = n_distinct(d$arm) >= 3 && !isSingular(fs))
})
say("Person-by-sample model refitted (should match T20 to 4 decimals):")
print(as.data.frame(lower_full))
out(lower_full, "R2c_sample_model_components_unrounded")
lower_comp <- lower_full %>%
  transmute(biomarker, var_true, var_err_mean = var_shared + var_assay / n_platforms,
            bound = "lower (shared variation not averaged)", identifiable)
R2_trial <- bind_rows(upper_comp, lower_comp) %>%
  left_join(amy_mean_bm, by = "biomarker") %>%
  crossing(designs) %>% crossing(effect = c(0.25, 0.50)) %>%
  mutate(var_slope_obs = var_true + var_err_mean / sst,
         reliability_of_trial_slope = var_true / var_slope_obs,
         n_per_arm_slope = n_slope(var_slope_obs, effect * abs(mean_slope_amy_pos_z))) %>%
  mutate(across(c(reliability_of_trial_slope, n_per_arm_slope), ~ if_else(identifiable, .x, NA_real_))) %>%
  select(biomarker, bound, identifiable, schedule, effect, reliability_of_trial_slope, n_per_arm_slope) %>%
  arrange(factor(biomarker, biomarkers), schedule, effect, bound)
# Check: the upper rows must reproduce T15 "mean of platforms"
chk <- R2_trial %>% filter(bound == "upper (paper model)") %>%
  left_join(T15 %>% filter(arm == "mean of platforms") %>%
              select(biomarker, schedule, effect, n_T15 = n_per_arm_slope),
            by = c("biomarker", "schedule", "effect"))
say("Check against T15 (mean of platforms): %d of %d rows identical",
    sum(chk$n_per_arm_slope == chk$n_T15, na.rm = TRUE), nrow(chk))
print(as.data.frame(R2_trial %>% filter(biomarker == "pTau217", effect == 0.25) %>%
                      mutate(reliability_of_trial_slope = round(reliability_of_trial_slope, 3))))
out(R2_trial, "R2a_trial_mean_of_platforms_bounds")

# Span needed (T16b formula) for the mean of platforms under each bound
R2_span <- bind_rows(
    T2 %>% transmute(biomarker, bound = "upper (paper model)", identifiable = TRUE,
                     ratio = (var_residual / n_arms) / var_person_slope),
    lower_full %>% transmute(biomarker, bound = "lower (shared variation not averaged)", identifiable,
                             ratio = (var_shared + var_assay / n_platforms) / var_true)) %>%
  crossing(target = c(0.7, 0.8, 0.9), k = c(2, 3, 4, 6)) %>%
  mutate(sst_needed = ratio * target / (1 - target),
         span_needed_mean_of_platforms = sqrt(sst_needed * 12 * (k - 1) / (k * (k + 1))),
         span_needed_mean_of_platforms = if_else(identifiable, span_needed_mean_of_platforms, NA_real_)) %>%
  select(biomarker, bound, identifiable, target, k, span_needed_mean_of_platforms) %>%
  arrange(factor(biomarker, biomarkers), target, k, bound)
print(as.data.frame(R2_span %>% filter(biomarker == "pTau217", k == 3) %>%
                      mutate(span_needed_mean_of_platforms = round(span_needed_mean_of_platforms, 1))))
out(R2_span, "R2b_span_needed_mean_of_platforms_bounds")



# =============================================================================================
# Section 3. Trial sizes with amyloid-positive variance components (R3)
# Original script: 08b17_paper2_revision_nov2026.R, R3
# =============================================================================================
# =============================================================================
# R3. Trial sizes with amyloid-positive variance components
# =============================================================================
say("\n== R3. Trial sizes with amyloid-positive variance ==")
amy <- long %>% filter(amyloid_bl_cl20 %in% TRUE)
r3_comp <- list(); r3_fits <- list()
for (bm in biomarkers) {
  d <- amy %>% filter(biomarker == bm) %>% mutate(arm_f = droplevels(factor(arm)))
  fr <- tryCatch(lmer(f_reduced, data = d, REML = TRUE, control = ctrl), error = function(e) NULL)
  fs <- tryCatch(lmer(f_sample,  data = d, REML = TRUE, control = ctrl), error = function(e) NULL)
  # per-platform models as for T3
  per_arm <- d %>% group_by(arm) %>% group_modify(~ {
    f <- tryCatch(lmer(as.formula(paste(YVAR, "~ plasma_time + (1 + plasma_time | RID_f)")),
                       data = .x, control = ctrl), error = function(e) NULL)
    if (is.null(f)) return(tibble(sd_residual = NA_real_, singular = NA, converged = NA))
    tibble(sd_residual = sigma(f)) %>% bind_cols(fit_info(f))
  }) %>% ungroup()
  r3_fits[[bm]] <- bind_rows(
    if (!is.null(fr)) fit_info(fr) %>% mutate(model = "reduced change model"),
    if (!is.null(fs)) fit_info(fs) %>% mutate(model = "with person-by-sample effect"),
    per_arm %>% transmute(model = paste("per-platform model:", arm), singular, converged)) %>%
    mutate(biomarker = bm, n_people = n_distinct(d$RID), .before = 1)
  r3_comp[[bm]] <- list(
    var_true = if (!is.null(fr)) vget(fr, "RID_f", "plasma_time") else NA_real_,
    var_res_joint = if (!is.null(fr)) vget(fr, "Residual") else NA_real_,
    var_true_s = if (!is.null(fs)) vget(fs, "RID_f", "plasma_time") else NA_real_,
    var_shared = if (!is.null(fs)) vget(fs, "sample_f") else NA_real_,
    var_assay  = if (!is.null(fs)) vget(fs, "Residual") else NA_real_,
    n_arms = n_distinct(d$arm), per_arm = per_arm %>% select(arm, sd_residual),
    ident = n_distinct(d$arm) >= 3 && !is.null(fs) && !isSingular(fs))
}
R3_fits <- bind_rows(r3_fits)
say("Fit diagnostics in amyloid-positive participants:")
print(as.data.frame(R3_fits))
out(R3_fits, "R3b_amyloid_positive_fit_diagnostics")

r3_rows <- list()
for (bm in names(r3_comp)) {
  cmp <- r3_comp[[bm]]
  per <- cmp$per_arm %>% mutate(biomarker = bm, var_true = cmp$var_true,
                                var_err = sd_residual^2, bound = "per platform", identifiable = TRUE)
  mean_up <- tibble(biomarker = bm, arm = "mean of platforms", var_true = cmp$var_true,
                    var_err = mean(cmp$per_arm$sd_residual^2) / cmp$n_arms, bound = "upper (paper model)",
                    identifiable = TRUE)
  mean_lo <- tibble(biomarker = bm, arm = "mean of platforms", var_true = cmp$var_true_s,
                    var_err = cmp$var_shared + cmp$var_assay / cmp$n_arms,
                    bound = "lower (shared variation not averaged)", identifiable = cmp$ident)
  r3_rows[[bm]] <- bind_rows(per, mean_up, mean_lo)
}
R3 <- bind_rows(r3_rows) %>%
  left_join(amy_slope %>% select(biomarker, arm, m_arm = mean_slope_amy_pos_z), by = c("biomarker", "arm")) %>%
  left_join(amy_mean_bm %>% rename(m_bm = mean_slope_amy_pos_z), by = "biomarker") %>%
  mutate(mean_slope_amy_pos_z = coalesce(m_arm, m_bm)) %>%
  crossing(designs) %>% crossing(effect = c(0.25, 0.50)) %>%
  mutate(var_slope_obs = var_true + var_err / sst,
         reliability_of_trial_slope = var_true / var_slope_obs,
         n_per_arm_slope_amy_pos_variance = n_slope(var_slope_obs, effect * abs(mean_slope_amy_pos_z))) %>%
  mutate(across(c(reliability_of_trial_slope, n_per_arm_slope_amy_pos_variance),
                ~ if_else(identifiable, .x, NA_real_))) %>%
  left_join(T15 %>% select(biomarker, arm, schedule, effect, n_per_arm_slope_whole_cohort_variance = n_per_arm_slope),
            by = c("biomarker", "arm", "schedule", "effect")) %>%
  mutate(n_per_arm_slope_whole_cohort_variance =
           if_else(bound == "lower (shared variation not averaged)", NA_real_, n_per_arm_slope_whole_cohort_variance)) %>%
  select(biomarker, arm, bound, identifiable, schedule, effect, reliability_of_trial_slope,
         n_per_arm_slope_amy_pos_variance, n_per_arm_slope_whole_cohort_variance) %>%
  arrange(factor(biomarker, biomarkers), schedule, effect, arm, bound)
print(as.data.frame(R3 %>% filter(biomarker == "pTau217", effect == 0.25,
                                  schedule == "24 months, 3 draws (0, 12, 24 m)") %>%
                      mutate(reliability_of_trial_slope = round(reliability_of_trial_slope, 3))))
out(R3, "R3a_trial_amyloid_positive_variance")



# =============================================================================================
# Section 4. Validity by SST fifth allowing for the precision of both slopes (R14b)
# Added November 2026. Section 1 (T17) predicts the ratio of the lowest to the highest SST fifth in
# validity against Centiloid change from plasma reliability alone (0.67 in the text). The text also gives
# the ratio when the Centiloid slope is allowed to carry measurement error too (0.45): Centiloid slope
# reliability comes from a random-slope model of Centiloid in the criterion window (plasma span +/- 180 days,
# as in 03), at each fifth's median Centiloid SST, from each participant's own scan times and, as a check,
# from evenly spaced scans at the fifth's median number of scans and span.
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR <- file.path(OUT_DIR, "paper2")
say("\n== Section 4. Validity by SST fifth, attenuation of both slopes (R14b) ==")
long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
critr  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
T2  <- read_csv(file.path(P2_DIR, "T2_change_variance_components.csv"), show_col_types = FALSE)
T17 <- read_csv(file.path(P2_DIR, "T17_validity_by_sst.csv"), show_col_types = FALSE)

# the fifths exactly as T17 (section 1, item 14), one row per person and platform
sst_person <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), span = max(plasma_time) - min(plasma_time),
            sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
bands <- slopes %>% filter(n_obs >= 2) %>% select(RID, biomarker, arm, slope) %>%
  left_join(sst_person %>% select(RID, k, span, sst), by = "RID") %>%
  left_join(critr %>% select(RID, centiloid_slope, centiloid_n, centiloid_span), by = "RID") %>%
  filter(!is.na(slope), !is.na(centiloid_slope), !is.na(sst)) %>%
  group_by(biomarker) %>% mutate(band = ntile(sst, 5)) %>% ungroup() %>%
  filter(biomarker == "pTau217")
obs_r <- bands %>% group_by(band) %>%
  summarise(r = cor(slope, centiloid_slope, method = "spearman"), .groups = "drop")

# Centiloid scans in the criterion window, as in 03 section 3
part_w <- readRDS(file.path(OUT_DIR, "participants.rds")) %>% select(RID, first_plasma, last_plasma)
cl <- readRDS(file.path(OUT_DIR, "criteria_long.rds")) %>% filter(measure == "centiloid", !is.na(value)) %>%
  inner_join(part_w, by = "RID") %>%
  filter(exam_date >= first_plasma - CRITERION_SLOPE_PAD_DAYS, exam_date <= last_plasma + CRITERION_SLOPE_PAD_DAYS) %>%
  mutate(t = as.numeric(exam_date - first_plasma) / YEAR_DAYS) %>%
  filter(RID %in% unique(bands$RID))
cl_person <- cl %>% group_by(RID) %>%
  summarise(n_scans = n(), span_c = max(t) - min(t), sst_c = sum((t - mean(t))^2), .groups = "drop")
cl_fit <- lmer(value ~ t + (1 + t | RID), data = cl %>% filter(RID %in% cl_person$RID[cl_person$n_scans >= 2]),
               REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))
vc_c <- as.data.frame(VarCorr(cl_fit))
v_slope_c <- vc_c$vcov[vc_c$grp == "RID" & vc_c$var1 == "t" & is.na(vc_c$var2)]
v_res_c   <- vc_c$vcov[vc_c$grp == "Residual"]
say("Centiloid model: true slope variance %.3f, residual variance %.3f, singular %s",
    v_slope_c, v_res_c, isSingular(cl_fit))

pT <- T2 %>% filter(biomarker == "pTau217"); stopifnot(nrow(pT) == 1)
rel_of   <- function(vs, vr, sst) vs / (vs + vr / sst)
even_sst <- function(k, s) ifelse(k >= 2, s^2 * k * (k + 1) / (12 * (k - 1)), NA_real_)
by_band <- bands %>% select(RID, band, sst) %>% left_join(cl_person, by = "RID") %>%
  group_by(band) %>%
  summarise(sst_median = median(sst), n_scans_median = median(n_scans, na.rm = TRUE),
            span_c_median = median(span_c, na.rm = TRUE), sst_c_median = median(sst_c, na.rm = TRUE), .groups = "drop") %>%
  mutate(rel_plasma = rel_of(pT$var_person_slope, pT$var_residual, sst_median),
         rel_cl_own = rel_of(v_slope_c, v_res_c, sst_c_median),
         rel_cl_even = rel_of(v_slope_c, v_res_c, even_sst(n_scans_median, span_c_median))) %>%
  left_join(obs_r, by = "band") %>%
  left_join(T17 %>% filter(biomarker == "pTau217") %>% select(band, r_T17 = r_centiloid_spearman), by = "band")
print(as.data.frame(by_band %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
lo <- by_band %>% filter(band == 1); hi <- by_band %>% filter(band == 5)
plasma_only <- sqrt(lo$rel_plasma / hi$rel_plasma)
R14b <- tibble(
  quantity = c("observed validity ratio, lowest / highest SST fifth", "predicted ratio, plasma reliability only",
               "predicted ratio, both slopes, own Centiloid scan times", "predicted ratio, both slopes, evenly spaced Centiloid scans",
               "lowest fifth, median Centiloid scans", "lowest fifth, median Centiloid span (y)",
               "highest fifth, median Centiloid scans", "highest fifth, median Centiloid span (y)",
               "Centiloid model true slope variance", "Centiloid model residual variance"),
  value = c(lo$r / hi$r, plasma_only,
            plasma_only * sqrt(lo$rel_cl_own / hi$rel_cl_own), plasma_only * sqrt(lo$rel_cl_even / hi$rel_cl_even),
            lo$n_scans_median, lo$span_c_median, hi$n_scans_median, hi$span_c_median, v_slope_c, v_res_c))
print(as.data.frame(R14b %>% mutate(value = round(value, 3))))
write_csv(R14b, file.path(P2_DIR, "R14b_validity_attenuation_both_slopes.csv"))
say("Written: R14b_validity_attenuation_both_slopes.csv in %s", P2_DIR)
