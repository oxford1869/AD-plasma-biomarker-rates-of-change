# =============================================================================================
# 07_a4learn.R
# Replication in A4/LEARN and the within-participant tests of sampling design.
#
# Sections:
#   1. A4/LEARN: random-slope change model per unit and slope reliability at the observed spread of sampling
#      times (T2); baseline to week-12 test-retest and its comparison with the long change model's residual
#      (T7, F1); within-person design test on participants with baseline, week-12 and week-240 samples (T8,
#      T8b, F8); minimal detectable change (T3, T14); validity (T5); sampling requirement (T6, T16b, F3);
#      participant bootstrap (T9; N_BOOT draws). COHORT = a4learn; writes a4learn_out/paper2/ (formerly 08b9)
#   2. the fourth (extension) sample: four slope designs from each participant's own samples, predicted and
#      refitted reliability and validity (combined/T9, T9b, F9; Figure 4C) (formerly 08b11)
#   3. the line through baseline, week 240 and the extension sample against twelve-week variance, and the
#      curvature of the mean trajectory (R8a, R8b) (formerly 08b17 R8)
# =============================================================================================


# =============================================================================================
# Section 1. A4/LEARN change model, test-retest, design test, bootstrap
# Original script: 08b9_paper2_a4learn.R
# =============================================================================================
# =============================================================================
# 08b9_paper2_a4learn.R
# Paper 2 on A4/LEARN: one platform (Lilly MSD pTau217), so the cross-platform
# decomposition of 08b does not apply. What this cohort adds instead:
#
#   1. Change decomposition on one platform, per unit (A4 placebo, LEARN, A4
#      solanezumab, placebo + LEARN, all arms) and floor rule: true-slope variance,
#      residual, and the reliability of a person's slope for the observed design
#      (var_slope / (var_slope + var_residual / SST)). With one platform the
#      platform-specific slope variance of ADNI cannot be separated; the reliability
#      here is the "one platform" quantity of Paper 2, with SST from this cohort.
#   2. Test-retest: A4 baseline and week-12 draws (about 85 days apart, before any
#      plausible true change). Within-person SD and CV on the log scale, ICC(2,1),
#      Bland-Altman limits, MDC95 between two single measurements, and a check that
#      the short-interval residual agrees with the residual of the long change model.
#   3. Design test on people's own draws: for A4 people with baseline, week-12 and
#      week-240 draws, the OLS slope is recomputed from each subset of their draws
#      (baseline + week 12; baseline + week 240; week 12 + week 240; all three). The
#      subsets differ only in SST, so the predicted reliability changes while the
#      person and the platform stay the same. Reported per subset: median SST,
#      predicted reliability, correlation with the three-draw slope, and validity
#      against concurrent Centiloid change. Also validity by observed SST bands
#      across all people with two or more draws.
#   4. Minimal detectable change per unit, on the z-log scale and as fold change.
#   5. Validity: Spearman correlation of OLS and BLUP slopes with concurrent change
#      in Centiloid, PACC, MMSE, CDR-SB, tau PET and CSF pTau/Abeta42, per unit.
#   6. Sampling requirement: reliability by number of draws and span, and the span
#      needed for reliability 0.7, 0.8, 0.9, per unit.
#   7. Person bootstrap (N_BOOT draws, default 1,000; N_BOOT=0 skips) for the slope
#      reliability per unit and the test-retest ICC.
#
# Reads   a4learn_out/frame_paper2_*.rds, participants.rds, plasma_long.rds
# Writes  a4learn_out/paper2/: T2_change_variance_components.csv,
#         T7_test_retest.csv, T8_design_subsets.csv, T8b_validity_by_sst.csv,
#         T3_minimal_detectable_change.csv, T14_mdc_fold_change.csv,
#         T5_validity_vs_criterion_change.csv, T6_sampling_requirement.csv,
#         T16b_span_needed.csv, T9_bootstrap_ci.csv, F1_test_retest.pdf,
#         F3_sampling_requirement.pdf, F8_design_subsets.pdf
# Run     by run_all.R, or with COHORT=a4learn set, source this file.
# =============================================================================

Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
if (COHORT != "a4learn") stop("08b9 is the A4/LEARN Paper 2 script; set Sys.setenv(COHORT = 'a4learn') first.")
suppressPackageStartupMessages({ library(lme4); library(ggplot2) })
N_BOOT <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))   # the paper uses 1,000 draws; 0 skips the bootstrap
set.seed(20260906)

P2_DIR <- file.path(OUT_DIR, "paper2"); dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }
YVAR <- "z_log_aligned"

part   <- readRDS(file.path(OUT_DIR, "participants.rds")) %>% select(RID, substudy, arm, has_short_interval_pair)
pl     <- readRDS(file.path(OUT_DIR, "plasma_long.rds")) %>%
  distinct(RID, visit_no, visit_code, visit_kind, on_treatment, short_interval_pair) %>% mutate(visit_code = as.character(visit_code))
long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  left_join(part, by = "RID") %>% left_join(pl, by = c("RID", "visit_no"))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>% left_join(part, by = "RID")
critr  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
define_units(long)   # MEASURES (units), m_label, unit_data(), always_floor, apply_floor(); see 00_setup.R
FLOORS <- c("as_delivered", "exclude_always_floor")

# scale of the log value per unit (z_log divided by it), to report fold changes
sd_log_of <- function(d) sd(d$log_value, na.rm = TRUE)

# ---- 1. Change decomposition per unit ----------------------------------------------------
say("\n== 1. Change decomposition (one platform) ==")
sst_of <- function(d) d %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), span = max(plasma_time) - min(plasma_time), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
fit_change1 <- function(d, label, floor) {
  d <- d %>% mutate(RID_f = factor(RID), arm_f = factor(arm))
  multi <- n_distinct(d$arm) > 1
  # arm fixed effects on level and slope in the combined units, so that the person-slope
  # variance is within-arm (the arms differ by design in amyloid burden)
  f <- as.formula(paste("y ~", if (multi) "arm_f * plasma_time" else "plasma_time", "+ (1 + plasma_time | RID_f)"))
  fit <- tryCatch(lmer(f, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
                  error = function(e) { say("  %s: change model failed: %s", label, conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  g <- function(v) { i <- which(vc$grp == "RID_f" & vc$var1 == v); if (length(i)) vc$vcov[i[1]] else NA_real_ }
  s <- sst_of(d); sst_med <- median(s$sst[s$k >= 2])
  vs <- g("plasma_time"); vr <- vc$vcov[vc$grp == "Residual"][1]
  tibble(unit = label, floor = floor, n_people = n_distinct(d$RID), n_obs = nrow(d),
         mean_slope = fixef(fit)[["plasma_time"]],
         var_person_int = g("(Intercept)"), var_person_slope = vs, var_residual = vr,
         sd_true_slope_per_year = sqrt(vs), sd_residual = sqrt(vr),
         sst_median = sst_med, median_span_3plus = median(s$span[s$k >= 3]),
         rel_slope_one_platform = vs / (vs + vr / sst_med),
         # reliability of each person's own slope at their own SST, averaged
         rel_slope_person_mean = mean(vs / (vs + vr / s$sst[s$k >= 2])),
         singular = isSingular(fit), converged = is.null(fit@optinfo$conv$lme4$messages),
         sd_log = sd_log_of(d)) %>% mutate(fit = list(fit))
}
change_rows <- list(); change_fits <- list()
for (m in MEASURES) for (fl in FLOORS) {
  d <- apply_floor(unit_data(m), fl) %>% mutate(y = .data[[YVAR]])
  r <- fit_change1(d, m_label[[m]], fl)
  if (is.null(r)) next
  if (fl == "as_delivered") change_fits[[m]] <- r$fit[[1]]
  change_rows[[paste(m, fl)]] <- r %>% select(-fit) %>% mutate(unit_key = m)
}
T2 <- bind_rows(change_rows)
print(as.data.frame(T2 %>% select(unit, floor, n_people, var_person_slope, var_residual, sst_median, rel_slope_one_platform, singular) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(T2 %>% select(-unit_key), "T2_change_variance_components")

# ---- 2. Test-retest (A4 baseline and week 12) ----------------------------------------------
say("\n== 2. Test-retest, baseline vs week 12 (A4) ==")
pairs <- long %>% filter(substudy == "A4", short_interval_pair %in% TRUE, visit_code %in% c("6", "9")) %>%
  select(RID, arm, visit_code, plasma_time, value, log_value, z = all_of(YVAR), lloq_flag) %>%
  pivot_wider(names_from = visit_code, values_from = c(plasma_time, value, log_value, z, lloq_flag)) %>%
  filter(!is.na(value_6), !is.na(value_9)) %>%
  mutate(gap_days = (plasma_time_9 - plasma_time_6) * YEAR_DAYS, d_log = log_value_9 - log_value_6, d_z = z_9 - z_6,
         any_floor = lloq_flag_6 | lloq_flag_9)
test_retest <- function(p, label) {
  if (nrow(p) < 10) return(NULL)
  dd <- p %>% select(RID, log_value_6, log_value_9) %>% pivot_longer(-RID, names_to = "occ", values_to = "lv")
  fit <- lmer(lv ~ 1 + (1 | RID) + (1 | occ), data = dd)
  v <- as.data.frame(VarCorr(fit)); vb <- v$vcov[v$grp == "RID"]; vo <- v$vcov[v$grp == "occ"]; ve <- v$vcov[v$grp == "Residual"]
  sw <- sd(p$d_log) / sqrt(2)                       # within-person SD on the log scale
  sw_z <- sd(p$d_z) / sqrt(2)                       # the same on the z-log scale used by the change models
  vt <- (var(p$log_value_6) + var(p$log_value_9)) / 2  # total variance of a single measurement
  tibble(sample = label, n_pairs = nrow(p), gap_days_median = median(p$gap_days),
         pct_any_floor = 100 * mean(p$any_floor),
         mean_diff_log = mean(p$d_log), sd_diff_log = sd(p$d_log),
         mean_diff_pct = 100 * (exp(mean(p$d_log)) - 1),
         loa_lower_pct = 100 * (exp(mean(p$d_log) - 1.96 * sd(p$d_log)) - 1),
         loa_upper_pct = 100 * (exp(mean(p$d_log) + 1.96 * sd(p$d_log)) - 1),
         sd_within_log = sw, cv_within_pct = 100 * sqrt(exp(sw^2) - 1),
         icc_2_1 = vb / (vb + vo + ve), icc_1_1_from_pairs = (vt - sw^2) / vt,
         mdc95_two_points_pct = 100 * (exp(1.96 * sqrt(2) * sw) - 1),
         pct_beyond_mdc = 100 * mean(abs(p$d_log) > 1.96 * sqrt(2) * sw),
         sd_within_z = sw_z, var_within_z = sw_z^2)
}
T7 <- bind_rows(test_retest(pairs, "A4 both arms"),
                test_retest(pairs %>% filter(arm == "Placebo"), "A4 placebo"),
                test_retest(pairs %>% filter(arm == "Solanezumab"), "A4 solanezumab"),
                test_retest(pairs %>% filter(!any_floor), "A4 both arms, no floor value"))
# compare with the residual of the long change model (z scale): residual variance from years
# of follow-up should not be smaller than the 12-week within-person variance
T7 <- T7 %>% mutate(var_residual_change_model_z = T2$var_residual[T2$unit == "A4 placebo" & T2$floor == "as_delivered"],
                    ratio_change_model_to_test_retest = var_residual_change_model_z / var_within_z)
print(as.data.frame(T7 %>% select(sample, n_pairs, gap_days_median, cv_within_pct, icc_2_1, mdc95_two_points_pct, pct_beyond_mdc, ratio_change_model_to_test_retest) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(T7, "T7_test_retest")

# ---- 3. Design test on people's own draws ---------------------------------------------------
say("\n== 3. Design test: slope from subsets of each person's draws ==")
ols <- function(t, y) if (length(t) >= 2 && var(t) > 0) cov(t, y) / var(t) else NA_real_
three <- long %>% filter(substudy == "A4", visit_code %in% c("6", "9", "66")) %>%
  select(RID, arm, visit_code, plasma_time, y = all_of(YVAR)) %>%
  group_by(RID) %>% filter(n_distinct(visit_code) == 3) %>% ungroup()
say("  A4 people with baseline, week-12 and week-240 draws: %d", n_distinct(three$RID))
subsets <- list(`baseline + week 12` = c("6", "9"), `baseline + week 240` = c("6", "66"),
                `week 12 + week 240` = c("9", "66"), `all three draws` = c("6", "9", "66"))
sub_slopes <- bind_rows(lapply(names(subsets), function(nm) three %>% filter(visit_code %in% subsets[[nm]]) %>%
  group_by(RID, arm) %>% summarise(subset = nm, k = n(), sst = sum((plasma_time - mean(plasma_time))^2),
                                   slope = ols(plasma_time, y), .groups = "drop")))
full <- sub_slopes %>% filter(subset == "all three draws") %>% select(RID, slope_full = slope)
comp <- T2 %>% filter(floor == "as_delivered") %>% select(unit, var_person_slope, var_residual)
design_row <- function(s, label, vs, vr) s %>% group_by(subset) %>%
  summarise(sample = label, n = sum(!is.na(slope)), sst_median = median(sst),
            reliability_predicted = vs / (vs + vr / median(sst)),
            r_with_three_draw_slope = cor(slope, slope_full, use = "complete.obs", method = "pearson"),
            r_centiloid_spearman = cor(slope, centiloid_slope, use = "complete.obs", method = "spearman"),
            n_centiloid = sum(complete.cases(slope, centiloid_slope)), .groups = "drop")
ss <- sub_slopes %>% left_join(full, by = "RID") %>% left_join(critr %>% select(RID, centiloid_slope), by = "RID")
T8 <- bind_rows(
  design_row(ss %>% filter(arm == "Placebo"), "A4 placebo", comp$var_person_slope[comp$unit == "A4 placebo"], comp$var_residual[comp$unit == "A4 placebo"]),
  design_row(ss %>% filter(arm == "Solanezumab"), "A4 solanezumab", comp$var_person_slope[comp$unit == "A4 solanezumab"], comp$var_residual[comp$unit == "A4 solanezumab"]),
  design_row(ss, "A4 both arms", comp$var_person_slope[comp$unit == "All three arms"], comp$var_residual[comp$unit == "All three arms"])) %>%
  mutate(subset = factor(subset, names(subsets))) %>% arrange(sample, subset) %>%
  # validity relative to the full design; if the plasma slope is the only thing that changes,
  # validity should scale with sqrt(reliability), the attenuation formula
  group_by(sample) %>% mutate(validity_ratio_to_full = r_centiloid_spearman / r_centiloid_spearman[subset == "all three draws"],
                              sqrt_reliability_ratio_to_full = sqrt(reliability_predicted / reliability_predicted[subset == "all three draws"])) %>% ungroup()
print(as.data.frame(T8 %>% select(sample, subset, n, sst_median, reliability_predicted, r_with_three_draw_slope, r_centiloid_spearman, validity_ratio_to_full, sqrt_reliability_ratio_to_full) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(T8, "T8_design_subsets")

# validity by observed SST band, all people with 2+ draws, per unit
sst_all <- sst_of(long %>% filter(!is.na(.data[[YVAR]])))
vb <- slopes %>% filter(n_obs >= 2) %>% select(RID, arm, slope) %>%
  left_join(sst_all, by = "RID") %>% left_join(critr %>% select(RID, centiloid_slope), by = "RID") %>%
  filter(!is.na(slope), !is.na(centiloid_slope))
band_of <- function(sst) cut(sst, c(-Inf, 0.1, 2, 8, 12, Inf), labels = c("< 0.1", "0.1 to 2", "2 to 8", "8 to 12", "> 12"))
T8b <- bind_rows(lapply(MEASURES, function(m) {
  ids <- unique(unit_data(m)$RID); cp <- comp %>% filter(unit == m_label[[m]])
  vb %>% filter(RID %in% ids) %>% mutate(band = band_of(sst)) %>% group_by(band) %>%
    summarise(unit = m_label[[m]], n = n(), sst_median = median(sst), k_median = median(k),
              reliability_predicted = cp$var_person_slope / (cp$var_person_slope + cp$var_residual / median(sst)),
              r_centiloid_spearman = if (n() >= 10) cor(slope, centiloid_slope, method = "spearman") else NA_real_, .groups = "drop") }))
print(as.data.frame(T8b %>% filter(unit == "Placebo + LEARN (primary)") %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(T8b, "T8b_validity_by_sst")

# ---- 4. Minimal detectable change per unit --------------------------------------------------
say("\n== 4. Minimal detectable change ==")
T3 <- T2 %>% filter(floor == "as_delivered") %>%
  mutate(sem = sd_residual, mdc95_two_points = 1.96 * sqrt(2) * sem,
         slope_se_3pts = sem / sqrt(sst_median), mdc95_slope_per_year = 1.96 * slope_se_3pts) %>%
  select(unit, unit_key, n_people, sd_residual, sd_true_slope_per_year, mean_slope, sst_median, median_span_3plus,
         mdc95_two_points, slope_se_3pts, mdc95_slope_per_year, sd_log)
rc <- bind_rows(lapply(MEASURES, function(m) slopes %>% filter(RID %in% unique(unit_data(m)$RID), n_obs >= 3) %>%
  summarise(unit_key = m, n_slopes = sum(!is.na(slope)),
            pct_reliable_increase = 100 * mean(slope > T3$mdc95_slope_per_year[T3$unit_key == m], na.rm = TRUE),
            pct_reliable_decrease = 100 * mean(slope < -T3$mdc95_slope_per_year[T3$unit_key == m], na.rm = TRUE))))
T3 <- T3 %>% left_join(rc, by = "unit_key")
out(T3 %>% select(-unit_key), "T3_minimal_detectable_change")
# on the delivered scale (fold change), as in 08b6 task 11
amy <- bind_rows(lapply(MEASURES, function(m) slopes %>% filter(RID %in% unique(unit_data(m)$RID), n_obs >= 3, amyloid_bl_cl20 %in% TRUE) %>%
  summarise(unit_key = m, mean_slope_amy_pos_z = mean(slope, na.rm = TRUE), n_amy_pos = n())))
T14 <- T3 %>% left_join(amy, by = "unit_key") %>%
  mutate(cv_within_pct = 100 * sqrt(exp((sd_residual * sd_log)^2) - 1),
         mean_pct_per_year = 100 * (exp(mean_slope * sd_log) - 1),
         mean_pct_per_year_amy_pos = 100 * (exp(mean_slope_amy_pos_z * sd_log) - 1),
         mdc95_two_points_pct = 100 * (exp(mdc95_two_points * sd_log) - 1),
         mdc95_slope_pct_per_year = 100 * (exp(mdc95_slope_per_year * sd_log) - 1),
         years_to_exceed_mdc_all = mdc95_two_points / abs(mean_slope),
         years_to_exceed_mdc_amy_pos = mdc95_two_points / abs(mean_slope_amy_pos_z)) %>%
  select(unit, n_people, sd_log, cv_within_pct, mean_pct_per_year, mean_pct_per_year_amy_pos, n_amy_pos,
         mdc95_two_points_pct, mdc95_slope_pct_per_year, years_to_exceed_mdc_all, years_to_exceed_mdc_amy_pos,
         pct_reliable_increase, pct_reliable_decrease)
print(as.data.frame(T14 %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(T14, "T14_mdc_fold_change")

# ---- 5. Validity against criterion change ---------------------------------------------------
say("\n== 5. Validity ==")
crit_cols <- intersect(c("centiloid_slope", "pacc_slope", "mmse_slope", "cdrsb_slope", "tau_pet_slope", "csf_ptau_abeta42_slope"), names(critr))
scor <- function(a, b) if (sum(complete.cases(a, b)) >= 10) cor(a, b, use = "complete.obs", method = "spearman") else NA_real_
blup <- bind_rows(lapply(names(change_fits), function(m) {
  re <- ranef(change_fits[[m]])$RID_f; fe <- fixef(change_fits[[m]])
  tibble(RID = as.integer(rownames(re)), unit_key = m, blup_slope = re[["plasma_time"]] + fe[["plasma_time"]]) }))
val <- bind_rows(lapply(MEASURES, function(m) {
  ids <- unique(unit_data(m)$RID)
  ol <- slopes %>% filter(RID %in% ids) %>% select(RID, n_obs, ols_slope = slope)
  bl <- blup %>% filter(unit_key == m) %>% select(RID, blup_slope)
  d <- ol %>% full_join(bl, by = "RID") %>% left_join(critr %>% select(RID, all_of(crit_cols)), by = "RID") %>%
    pivot_longer(all_of(crit_cols), names_to = "criterion", values_to = "crit_slope")
  bind_rows(
    d %>% filter(n_obs >= 2) %>% group_by(criterion) %>% summarise(estimate = "OLS slope, 2+ draws", n = sum(complete.cases(ols_slope, crit_slope)),
      r = scor(ols_slope, crit_slope), .groups = "drop"),
    d %>% filter(n_obs >= 3) %>% group_by(criterion) %>% summarise(estimate = "OLS slope, 3+ draws", n = sum(complete.cases(ols_slope, crit_slope)),
      r = scor(ols_slope, crit_slope), .groups = "drop"),
    d %>% group_by(criterion) %>% summarise(estimate = "BLUP slope", n = sum(complete.cases(blup_slope, crit_slope)),
      r = scor(blup_slope, crit_slope), .groups = "drop")) %>%
    mutate(unit = m_label[[m]]) }))
# same sign convention as ADNI: positive means the plasma slope tracks worsening
FLIP <- c("hcv_adj_slope", "mpacc_digit_slope", "mmse_slope")
T5 <- val %>% select(unit, criterion, estimate, n, r) %>%
  mutate(r = round(ifelse(criterion %in% FLIP, -r, r), 3))
print(as.data.frame(T5 %>% filter(criterion == "centiloid_slope")))
out(T5, "T5_validity_vs_criterion_change")

# ---- 6. Sampling requirement --------------------------------------------------------------
say("\n== 6. Sampling requirement ==")
grid <- expand_grid(k = c(2, 3, 4, 5, 6), span_years = c(2, 4, 6, 8, 10))
T6 <- T2 %>% filter(floor == "as_delivered") %>% select(unit, var_person_slope, var_residual) %>% crossing(grid) %>%
  mutate(sst = map2_dbl(k, span_years, ~ { t <- seq(0, .y, length.out = .x); sum((t - mean(t))^2) }),
         reliability_one_platform = var_person_slope / (var_person_slope + var_residual / sst)) %>%
  select(unit, k, span_years, sst, reliability_one_platform)
out(T6, "T6_sampling_requirement")
T16b <- T2 %>% filter(floor == "as_delivered") %>% select(unit, var_person_slope, var_residual) %>%
  crossing(target = c(0.7, 0.8, 0.9), k = c(2, 3, 4, 6)) %>%
  mutate(sst_needed = (var_residual / var_person_slope) * target / (1 - target),
         span_needed_years = sqrt(sst_needed * 12 * (k - 1) / (k * (k + 1)))) %>%
  select(unit, target, k, sst_needed, span_needed_years)
print(as.data.frame(T16b %>% filter(k == 3) %>% mutate(across(where(is.numeric), ~ round(.x, 1)))))
out(T16b, "T16b_span_needed")

# ---- 7. Person bootstrap -----------------------------------------------------------------------
if (N_BOOT > 0) {
  say("\n== 7. Person bootstrap (%d draws) ==", N_BOOT)
  boot_one <- function(b) {
    rows <- list()
    for (m in MEASURES) {
      d <- apply_floor(unit_data(m), "as_delivered") %>% mutate(y = .data[[YVAR]])
      ids <- unique(d$RID); smp <- sample(ids, replace = TRUE)
      db <- bind_rows(lapply(seq_along(smp), function(j) d[d$RID == smp[j], ] %>% mutate(RID = j)))
      r <- suppressMessages(fit_change1(db, m_label[[m]], "as_delivered"))
      if (!is.null(r)) rows[[m]] <- tibble(b = b, quantity = "rel_slope_one_platform", unit = m_label[[m]], value = r$rel_slope_one_platform)
    }
    pids <- unique(pairs$RID); ps <- sample(pids, replace = TRUE)
    pb <- bind_rows(lapply(seq_along(ps), function(j) pairs[pairs$RID == ps[j], ] %>% mutate(RID = j)))
    tr <- suppressMessages(test_retest(pb, "A4 both arms"))
    if (!is.null(tr)) rows[["tr"]] <- tibble(b = b, quantity = c("icc_2_1", "cv_within_pct"), unit = "A4 both arms", value = c(tr$icc_2_1, tr$cv_within_pct))
    bind_rows(rows)
  }
  draws <- bind_rows(lapply(seq_len(N_BOOT), function(b) { if (b %% 20 == 0) say("  bootstrap %d of %d", b, N_BOOT); boot_one(b) }))
  T9 <- draws %>% group_by(quantity, unit) %>%
    summarise(n_draws = n(), lo95 = quantile(value, 0.025), hi95 = quantile(value, 0.975), .groups = "drop") %>%
    left_join(bind_rows(T2 %>% filter(floor == "as_delivered") %>% transmute(quantity = "rel_slope_one_platform", unit, estimate = rel_slope_one_platform),
                        T7 %>% filter(sample == "A4 both arms") %>% transmute(quantity = "icc_2_1", unit = sample, estimate = icc_2_1),
                        T7 %>% filter(sample == "A4 both arms") %>% transmute(quantity = "cv_within_pct", unit = sample, estimate = cv_within_pct)),
              by = c("quantity", "unit")) %>% select(quantity, unit, estimate, lo95, hi95, n_draws)
  print(as.data.frame(T9 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
  out(T9, "T9_bootstrap_ci")
} else say("\nBootstrap skipped (N_BOOT = 0)")

# ---- Figures ------------------------------------------------------------------------------------
INK2 <- "#5a5751"; MUTED <- "#8d8a83"; BLUE <- "#1F6F8B"; ORANGE <- "#B0413E"
sv <- function(p, name, w, h) ggsave(file.path(P2_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = "pdf", useDingbats = FALSE, bg = "white")
# F1: Bland-Altman of the 12-week pair on the log scale
ba <- T7 %>% filter(sample == "A4 both arms")
f1 <- ggplot(pairs, aes((log_value_6 + log_value_9) / 2, d_log, colour = arm)) +
  geom_hline(yintercept = c(ba$mean_diff_log, ba$mean_diff_log + c(-1, 1) * 1.96 * ba$sd_diff_log), linetype = c(1, 2, 2), colour = INK2) +
  geom_point(size = 0.8, alpha = 0.6) +
  scale_colour_manual(values = c(Placebo = BLUE, Solanezumab = ORANGE), name = NULL) +
  labs(x = "Mean of the two log pTau217 values", y = "Week 12 minus baseline (log)",
       title = sprintf("Test-retest over %.0f days: within-person CV %.1f%%, ICC %.2f", ba$gap_days_median, ba$cv_within_pct, ba$icc_2_1))
sv(f1, "F1_test_retest", 6.5, 4.5)
# F3: sampling requirement curves per unit
f3 <- ggplot(T6, aes(span_years, reliability_one_platform, colour = factor(k))) + geom_line() + geom_point(size = 1) +
  facet_wrap(~ unit) + geom_hline(yintercept = 0.7, linetype = 2, colour = MUTED) +
  labs(x = "span (years)", y = "reliability of the slope", colour = "draws")
sv(f3, "F3_sampling_requirement", 9, 5.5)
# F8: design test, validity against SST with the attenuation prediction
f8 <- ggplot(T8, aes(sst_median, r_centiloid_spearman, colour = sample)) +
  geom_line(aes(y = r_centiloid_spearman[subset == "all three draws"][1] * sqrt_reliability_ratio_to_full), linetype = 2) +
  geom_point(size = 2.2) + geom_text(aes(label = subset), size = 2.6, vjust = -0.8, show.legend = FALSE) +
  scale_colour_manual(values = c(`A4 placebo` = BLUE, `A4 solanezumab` = ORANGE, `A4 both arms` = INK2), name = NULL) +
  labs(x = "Spread of sampling times, SST (year^2)", y = "Spearman r with Centiloid change",
       title = "Same people, same platform, different subsets of their own draws (dashed: attenuation prediction)")
sv(f8, "F8_design_subsets", 7.5, 5)

say("\nA4/LEARN Paper 2 tables and figures written to %s", P2_DIR)


# =============================================================================================
# Section 2. Fourth plasma sample: sampling span against sample count (Figure 4C)
# Original script: 08b11_paper2_fourth_sample.R
# =============================================================================================
# =============================================================================
# 08b11_paper2_fourth_sample.R
# Paper 2: is it the fourth plasma sample that helps, or the longer window?
#
# Every reliability figure in Paper 2 belongs to one design. The obvious question
# is what would improve it, and the two candidate answers, more samples and a
# longer window, normally move together. A4 separates them by accident.
#
# The design. Everyone in A4 has plasma at year 0, at about year 0.24 (week 12)
# and at about year 4.8 (week 240). A subgroup has a fourth sample from the
# extension, at a median of year 7.3. That subgroup is therefore not a different
# design: it is the same design plus one late sample, so every comparison below
# is made WITHIN the same people and selection into the extension cannot bias the
# contrast between designs. It affects only who is in the analysis.
#
# Four designs are built from each person's own samples:
#   A  first, second, third          the standard design, window 4.8 years
#   B  first, second, fourth         three samples, window 7.3 years
#   C  first, third, fourth          same window as B, middle sample moved later
#   D  all four                      window 7.3 years, four samples
# A against B isolates the window. B against D isolates the extra sample.
# B against C isolates where the middle sample sits.
#
# What is predicted, written down before the analysis runs. Reliability of a
# slope satisfies r/(1-r) = (true slope variance / residual variance) x SST, so
# the ODDS of reliability are linear in SST. Anchoring on the A4 placebo arm,
# which Paper 2 puts at 0.601 at an SST of 13.53, gives about 0.62 for design A,
# 0.75 for C, 0.79 for B and 0.81 for D. In words: the longer window is worth
# about 0.17 and the fourth sample about 0.02. The script prints the predicted
# and the observed side by side so the comparison is explicit.
#
# Caveats, which the printed output repeats:
#   - the subgroup are trial completers who entered the extension;
#   - the fourth sample falls after the treatment period, so the placebo arm is
#     the primary analysis and the arms combined are the sensitivity analysis;
#   - roughly 77 people in the placebo arm, so intervals are wide;
#   - with one assay the reliability is model-implied, not measured against a
#     second measurement of the same blood.
#
# Run by run_all.R after the A4/LEARN frames exist (it fits only small mixed
# models and takes seconds); it can also be sourced on its own.
#
# Reads   a4learn_out/frame_paper3_long.rds   the plasma series and Centiloid
#         a4learn_out/frame_paper3_person.rds the concurrent Centiloid slope
#         a4learn_out/participants.rds        the arm
# Writes  combined/T9_fourth_sample.csv       one row per design and sample
#         combined/T9b_fourth_sample_pairs.csv the three contrasts
#         combined/F9_span_against_count.pdf/.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr); library(purrr)
})
HAVE_LME4    <- requireNamespace("lme4", quietly = TRUE)
HAVE_GGPLOT2 <- requireNamespace("ggplot2", quietly = TRUE)
if (HAVE_GGPLOT2) suppressPackageStartupMessages(library(ggplot2))

Sys.setenv(COHORT = "adni"); source("00_setup.R")
out_root <- OUT_ROOT
A4_DIR   <- file.path(out_root, "a4learn_out")
COMB_DIR <- file.path(out_root, "combined")
dir.create(COMB_DIR, showWarnings = FALSE, recursive = TRUE)
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")
MEASURE <- "Lilly_plasma_ptau217"
N_BOOT  <- as.integer(Sys.getenv("P2_BOOT4", "1000"))
set.seed(as.integer(Sys.getenv("P2_SEED4", "20260908")))

say("\n============================================================")
say("08b11  Fourth plasma sample: window against sample count")
say("============================================================")

# ---- 1. the people with four samples ---------------------------------------------------------
long <- readRDS(file.path(A4_DIR, "frame_paper3_long.rds")) %>% filter(measure == MEASURE)
person <- readRDS(file.path(A4_DIR, "frame_paper3_person.rds"))
arms <- readRDS(file.path(A4_DIR, "participants.rds")) %>% select(RID, arm) %>% distinct()

series <- long %>%
  distinct(RID, plasma_time, .keep_all = TRUE) %>%
  group_by(RID) %>% arrange(plasma_time, .by_group = TRUE) %>%
  mutate(pos = row_number(), k = n()) %>% ungroup()

# the long frame already carries the arm; fall back to participants.rds if it does not
four <- series %>% filter(k == 4)
if (!"arm" %in% names(four)) four <- four %>% left_join(arms, by = "RID")
if (!nrow(four)) stop("No participant has four plasma samples in ", A4_DIR)
stopifnot("arm" %in% names(four))

say("\nParticipants with four plasma samples: %d", n_distinct(four$RID))
print(as.data.frame(four %>% distinct(RID, arm) %>% count(arm, name = "people")))
say("\nVisit times of these people, years from the first sample:")
print(as.data.frame(four %>% group_by(pos) %>%
  summarise(n = n(), median = round(median(plasma_time), 2),
            q1 = round(quantile(plasma_time, 0.25), 2),
            q3 = round(quantile(plasma_time, 0.75), 2), .groups = "drop")))
say("\nFor comparison, the people with exactly three samples (n = %d):",
    n_distinct(series$RID[series$k == 3]))
print(as.data.frame(series %>% filter(k == 3) %>% group_by(pos) %>%
  summarise(n = n(), median = round(median(plasma_time), 2), .groups = "drop")))

# ---- 2. the four designs ---------------------------------------------------------------------
DESIGNS <- list(
  `A. first, second, third (standard)`     = c(1, 2, 3),
  `B. first, second, fourth (long window)` = c(1, 2, 4),
  `C. first, third, fourth (middle moved)` = c(1, 3, 4),
  `D. all four samples`                    = c(1, 2, 3, 4))

# The slope from one subset of a person's own samples. y is the standardised log value,
# the same scale every other Paper 2 slope uses.
subset_slopes <- function(d, positions) {
  d %>% filter(pos %in% positions) %>%
    group_by(RID, arm) %>%
    filter(n() == length(positions), sd(plasma_time) > 0) %>%
    summarise(k       = n(),
              sst     = sum((plasma_time - mean(plasma_time))^2),
              span    = max(plasma_time) - min(plasma_time),
              slope   = coef(lm(z_log ~ plasma_time))[2],
              .groups = "drop")
}

slopes <- imap_dfr(DESIGNS, ~ subset_slopes(four, .x) %>% mutate(design = .y, .before = 1))
say("\nSlopes computed: %d rows over %d designs", nrow(slopes), n_distinct(slopes$design))

# ---- 3. reliability, model-implied, within each design ----------------------------------------
# Reliability of an individual slope = var(true slope) / (var(true slope) + var(residual)/SST).
# var(true slope) comes from a random-slope mixed model fitted to the samples that design uses.
reliability_of <- function(d, positions) {
  dd <- d %>% filter(pos %in% positions) %>% group_by(RID) %>% filter(n() == length(positions)) %>% ungroup()
  if (!HAVE_LME4 || n_distinct(dd$RID) < 30) return(c(var_true = NA, var_resid = NA, reliability = NA))
  m <- try(lme4::lmer(z_log ~ plasma_time + (plasma_time | RID), data = dd,
                      REML = TRUE, control = lme4::lmerControl(calc.derivs = FALSE)), silent = TRUE)
  if (inherits(m, "try-error")) return(c(var_true = NA, var_resid = NA, reliability = NA))
  vc <- as.data.frame(lme4::VarCorr(m))
  vt <- vc$vcov[vc$grp == "RID" & vc$var1 == "plasma_time" & is.na(vc$var2)]
  vr <- vc$vcov[vc$grp == "Residual"]
  sst <- dd %>% group_by(RID) %>% summarise(s = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
  c(var_true = vt, var_resid = vr, reliability = vt / (vt + vr / median(sst$s)))
}

# The prediction, from the anchor Paper 2 already reports for the A4 placebo arm.
ANCHOR_R   <- as.numeric(Sys.getenv("P2_ANCHOR_R",   "0.601"))
ANCHOR_SST <- as.numeric(Sys.getenv("P2_ANCHOR_SST", "13.533"))
odds_per_sst <- (ANCHOR_R / (1 - ANCHOR_R)) / ANCHOR_SST
predicted <- function(sst) { o <- odds_per_sst * sst; o / (1 + o) }

# ---- 4. validity, and agreement with the full four-sample slope --------------------------------
cent <- person %>% select(RID, centiloid_slope)
full <- slopes %>% filter(design == names(DESIGNS)[4]) %>% select(RID, slope_full = slope)

summarise_design <- function(s, d, label) {
  positions <- DESIGNS[[unique(s$design)]]
  rel <- reliability_of(d, positions)
  s2 <- s %>% left_join(full, by = "RID") %>% left_join(cent, by = "RID")
  ok <- !is.na(s2$centiloid_slope)
  tibble(sample = label, design = unique(s$design), n = nrow(s2),
         k = unique(s$k), sst_median = median(s2$sst), span_median = median(s2$span),
         reliability_predicted = predicted(median(s2$sst)),
         reliability_model     = unname(rel["reliability"]),
         r_with_full_slope     = suppressWarnings(cor(s2$slope, s2$slope_full,
                                                      method = "pearson", use = "complete.obs")),
         r_centiloid_spearman  = suppressWarnings(cor(s2$slope[ok], s2$centiloid_slope[ok],
                                                      method = "spearman")),
         n_centiloid           = sum(ok))
}

run_sample <- function(keep, label) {
  d <- four %>% filter(RID %in% keep)
  s <- slopes %>% filter(RID %in% keep)
  if (n_distinct(d$RID) < 20) { say("  %s: fewer than 20 people, skipped", label); return(NULL) }
  imap_dfr(DESIGNS, function(pos, nm) summarise_design(s %>% filter(design == nm), d, label))
}

placebo <- arms$RID[arms$arm == "Placebo"]
T9 <- bind_rows(
  run_sample(unique(four$RID[four$arm == "Placebo"]), "A4 placebo (primary)"),
  run_sample(unique(four$RID), "Both A4 arms (sensitivity)"))
T9 <- T9 %>%
  group_by(sample) %>%
  mutate(validity_ratio_to_D  = r_centiloid_spearman / r_centiloid_spearman[design == names(DESIGNS)[4]],
         sqrt_reliability_ratio = sqrt(reliability_predicted / reliability_predicted[design == names(DESIGNS)[4]])) %>%
  ungroup()
write_csv(T9, file.path(COMB_DIR, "T9_fourth_sample.csv"))

say("\n== Designs, predicted against observed ==")
print(as.data.frame(T9 %>%
  transmute(sample, design, n, k, sst = round(sst_median, 1), span = round(span_median, 2),
            r_predicted = round(reliability_predicted, 2),
            r_model = round(reliability_model, 2),
            r_with_D = round(r_with_full_slope, 2),
            r_centiloid = round(r_centiloid_spearman, 2),
            validity_vs_D = round(validity_ratio_to_D, 2))))

# ---- 5. the three contrasts, with a within-person bootstrap ------------------------------------
CONTRASTS <- list(
  `Window: A against B (3 samples, 4.8 against 7.3 years)`     = c(1, 2),
  `Sample count: B against D (7.3 years, 3 against 4 samples)` = c(2, 4),
  `Placement: B against C (same window, middle sample moved)`  = c(2, 3))

# Participant bootstrap of the three contrasts (added November 2026, so that the gains in
# Table S13c carry intervals like every other estimate). Draw b resamples participants with
# replacement (seed P2_SEED4 + b), gives each copy its own id so copies are separate clusters,
# recomputes within the resample the predicted reliability (from the median SST), the refitted
# reliability (the random-slope model of section 3) and the validity against Centiloid change
# for each design, and takes the three differences. Percentile 95% intervals over N_BOOT draws.
BOOT_SEED4 <- as.integer(Sys.getenv("P2_SEED4", "20260908"))
boot_contrasts <- function(keep, label) {
  d_all <- four %>% filter(RID %in% keep)
  s_all <- slopes %>% filter(RID %in% keep) %>% left_join(cent, by = "RID")
  ids <- unique(d_all$RID)
  one <- function(b) {
    set.seed(BOOT_SEED4 + b)
    key <- tibble(RID = sample(ids, length(ids), replace = TRUE)) %>% mutate(boot_id = row_number())
    d_b <- key %>% inner_join(d_all, by = "RID", relationship = "many-to-many") %>%
      mutate(RID = boot_id) %>% select(-boot_id)
    s_b <- key %>% inner_join(s_all, by = "RID", relationship = "many-to-many")
    per <- imap_dfr(DESIGNS, function(pos, nm) {
      sd <- s_b %>% filter(design == nm); ok <- !is.na(sd$centiloid_slope)
      tibble(design = nm,
             reliability_predicted = predicted(median(sd$sst)),
             reliability_model     = unname(reliability_of(d_b, pos)["reliability"]),
             r_centiloid           = suppressWarnings(cor(sd$slope[ok], sd$centiloid_slope[ok], method = "spearman")))
    })
    bind_rows(lapply(names(CONTRASTS), function(nm) {
      ij <- CONTRASTS[[nm]]
      tibble(draw = b, contrast = nm,
             predicted_gain = per$reliability_predicted[ij[2]] - per$reliability_predicted[ij[1]],
             model_gain     = per$reliability_model[ij[2]]     - per$reliability_model[ij[1]],
             validity_gain  = per$r_centiloid[ij[2]]           - per$r_centiloid[ij[1]])
    }))
  }
  draws <- bind_rows(lapply(seq_len(N_BOOT), function(b) {
    if (b %% 100 == 0) say("  %s: bootstrap draw %d of %d", label, b, N_BOOT); one(b) }))
  draws %>% group_by(contrast) %>%
    summarise(across(c(predicted_gain, model_gain, validity_gain),
                     list(lo95 = ~ quantile(.x, 0.025, na.rm = TRUE), hi95 = ~ quantile(.x, 0.975, na.rm = TRUE))),
              n_draws = sum(!is.na(model_gain)), .groups = "drop") %>%
    mutate(sample = label, .before = 1)
}

T9b <- bind_rows(lapply(unique(T9$sample), function(lab) {
  t <- T9 %>% filter(sample == lab)
  bind_rows(lapply(names(CONTRASTS), function(nm) {
    ij <- CONTRASTS[[nm]]
    tibble(sample = lab, contrast = nm,
           sst_from   = t$sst_median[ij[1]],           sst_to   = t$sst_median[ij[2]],
           predicted_from = t$reliability_predicted[ij[1]],
           predicted_to   = t$reliability_predicted[ij[2]],
           predicted_gain = t$reliability_predicted[ij[2]] - t$reliability_predicted[ij[1]],
           model_gain     = t$reliability_model[ij[2]]  - t$reliability_model[ij[1]],
           validity_gain  = t$r_centiloid_spearman[ij[2]] - t$r_centiloid_spearman[ij[1]])
  }))
}))
if (N_BOOT > 0 && HAVE_LME4) {
  say("\n== Bootstrap of the three contrasts, %d draws ==", N_BOOT)
  ci <- bind_rows(boot_contrasts(unique(four$RID[four$arm == "Placebo"]), "A4 placebo (primary)"),
                  boot_contrasts(unique(four$RID), "Both A4 arms (sensitivity)"))
  T9b <- T9b %>% left_join(ci, by = c("sample", "contrast"))
}
write_csv(T9b, file.path(COMB_DIR, "T9b_fourth_sample_pairs.csv"))

say("\n== The three contrasts ==")
print(as.data.frame(T9b %>% transmute(sample, contrast,
  sst = sprintf("%.1f to %.1f", sst_from, sst_to),
  predicted = sprintf("%.2f to %.2f (%+.2f)", predicted_from, predicted_to, predicted_gain),
  model_gain = if ("model_gain_lo95" %in% names(T9b))
    sprintf("%.2f (%.2f to %.2f)", model_gain, model_gain_lo95, model_gain_hi95) else sprintf("%.2f", model_gain),
  validity_gain = if ("validity_gain_lo95" %in% names(T9b))
    sprintf("%.2f (%.2f to %.2f)", validity_gain, validity_gain_lo95, validity_gain_hi95) else sprintf("%.2f", validity_gain))))

say("\nRead this as: the window contrast should be large and the sample-count contrast small.")
say("If it is, the recommendation is to space plasma samples further apart rather than add one.")

# ---- 6. figure ---------------------------------------------------------------------------------
if (HAVE_GGPLOT2) {
  # the shared theme is loaded by 00_setup.R
  curve <- tibble(sst = seq(5, max(T9$sst_median) * 1.15, length.out = 200)) %>%
    mutate(reliability = predicted(sst))
  # long form so that predicted and fitted are a legend rather than a sentence in the subtitle
  pts <- T9 %>%
    select(sample, design, sst_median, Predicted = reliability_predicted,
           Fitted = reliability_model) %>%
    pivot_longer(c(Predicted, Fitted), names_to = "source", values_to = "reliability") %>%
    mutate(source = factor(source, c("Predicted", "Fitted")))
  # one letter per design, placed above the predicted point of the primary sample only
  labs_df <- T9 %>%
    group_by(design) %>%
    summarise(sst_median = median(sst_median),
              reliability = max(c(reliability_predicted, reliability_model), na.rm = TRUE),
              letter = substr(first(design), 1, 1), .groups = "drop")
  F9 <- ggplot(pts, aes(sst_median, reliability)) +
    geom_line(data = curve, aes(sst, reliability), colour = "#8d8a83", linetype = 3,
              inherit.aes = FALSE) +
    geom_point(aes(colour = sample, shape = source), size = 2.6) +
    geom_text(data = labs_df, aes(label = letter), vjust = -1.4, size = 3.4,
              inherit.aes = TRUE, show.legend = FALSE) +
    scale_shape_manual(values = c(Predicted = 16, Fitted = 1), name = NULL) +
    labs(x = "Spread of sampling times (SST), years squared",
         y = "Reliability of an individual slope", colour = NULL,
         title = "A longer window, not a fourth sample, is what raises the reliability of a slope",
         subtitle = paste0("A standard design; B the same three samples over a longer window;\n",
                           "C middle sample moved later; D all four samples.\n",
                           "Dotted line: reliability implied by the spread of sampling times alone.")) +
    theme(legend.position = "bottom", legend.box = "vertical")
  ggsave(file.path(COMB_DIR, "F9_span_against_count.pdf"), F9, width = 9, height = 5,
         device = "pdf", useDingbats = FALSE, bg = "white")
  ggsave(file.path(COMB_DIR, "F9_span_against_count.png"), F9, width = 9, height = 5,
         dpi = 300, bg = "white")

  # A panel version for Figure 4C of the manuscript: the same plot with its title and
  # subtitle removed, because a panel's caption belongs in the figure legend. 08b7 reads
  # this file if it exists and appends it as panel C; if it does not exist, Figure 4 is
  # built with two panels exactly as before.
  F9_panel <- F9 + labs(title = NULL, subtitle = NULL)
  saveRDS(F9_panel, file.path(COMB_DIR, "F9_panel.rds"))
  say("  F9_panel.rds written for Figure 4C; re-run 08b7_paper2_figures.R to rebuild Figure 4")
}

say("\nWritten to %s", COMB_DIR)
say("  T9_fourth_sample.csv, T9b_fourth_sample_pairs.csv, F9_span_against_count")
say("\nCaveats to carry into the text: the subgroup are extension completers, the fourth")
say("sample falls after the treatment period, the placebo arm is small, and with one assay")
say("the reliability is model-implied rather than measured against a second measurement.")


# =============================================================================================
# Section 3. Line without week 12 against 12-week variance; curvature (R8)
# Original script: 08b17_paper2_revision_nov2026.R, R8
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
load_change_objects()
A4_OUT   <- file.path(OUT_ROOT, "a4learn_out")
N_BOOT_R8 <- as.integer(Sys.getenv("P2_BOOT_R8", "2000"))
set.seed(20261106)
# =============================================================================
# R8. A4: 12-week check without circularity, and curvature
# =============================================================================
say("\n== R8. A4: residual of a line that leaves out week 12, against 12-week variance ==")
a4 <- readRDS(file.path(A4_OUT, "frame_paper3_long.rds")) %>% filter(measure == "Lilly_plasma_ptau217")
if (!"arm" %in% names(a4)) {
  a4 <- a4 %>% left_join(readRDS(file.path(A4_OUT, "participants.rds")) %>% distinct(RID, arm), by = "RID")
}
series <- a4 %>% distinct(RID, plasma_time, .keep_all = TRUE) %>%
  group_by(RID) %>% arrange(plasma_time, .by_group = TRUE) %>%
  mutate(pos = row_number(), k = n()) %>% ungroup()
four <- series %>% filter(k == 4) %>% group_by(RID) %>%
  filter(plasma_time[pos == 2] < 0.5) %>% ungroup()          # second sample must be the week-12 draw
say("Participants with four samples and a week-12 second sample: %d", n_distinct(four$RID))
print(as.data.frame(four %>% group_by(pos) %>%
                      summarise(n = n(), median_years = round(median(plasma_time), 2), .groups = "drop")))

# Signed residual of the line through samples 1, 3, 4: the contrast orthogonal to 1 and t,
# scaled so its square is the residual sum of squares. Positive = sample 3 above the line.
person_stats <- four %>% group_by(RID, arm) %>%
  summarise(t1 = plasma_time[pos == 1], t3 = plasma_time[pos == 3], t4 = plasma_time[pos == 4],
            y1 = z_log[pos == 1], y2 = z_log[pos == 2], y3 = z_log[pos == 3], y4 = z_log[pos == 4],
            .groups = "drop") %>%
  mutate(a1 = -(t4 - t3), a3 = t4 - t1, a4 = -(t3 - t1),
         s_line = (a1 * y1 + a3 * y3 + a4 * y4) / sqrt(a1^2 + a3^2 + a4^2),
         d12 = y2 - y1)
# Both quantities centred within arm, so a mean departure (population curvature, a shift at the
# extension sample, a treatment effect) is removed from both alike and reported separately.
r8_stats <- function(ps) {
  k   <- n_distinct(ps$arm)
  s_c <- ps$s_line - ave(ps$s_line, ps$arm); d_c <- ps$d12 - ave(ps$d12, ps$arm)
  v_line <- sum(s_c^2) / (nrow(ps) - k); v_12 <- sum(d_c^2) / (nrow(ps) - k) / 2
  c(n = nrow(ps), var_residual_line_without_week12 = v_line, var_12week = v_12,
    ratio = v_line / v_12, mean_signed_residual = mean(ps$s_line))
}
r8_boot <- function(ps) {
  b <- replicate(N_BOOT_R8, r8_stats(ps[sample(nrow(ps), replace = TRUE), ])[["ratio"]])
  quantile(b, c(0.025, 0.975), na.rm = TRUE)
}
groups <- list(`all arms` = person_stats, `placebo arm` = person_stats %>% filter(arm == "Placebo"))
R8a <- imap_dfr(groups, function(ps, nm) {
  s <- r8_stats(ps); ci <- r8_boot(ps)
  tibble(group = nm, n = s[["n"]], var_residual_line_without_week12 = s[["var_residual_line_without_week12"]],
         var_12week = s[["var_12week"]], ratio = s[["ratio"]], ratio_lo95 = ci[[1]], ratio_hi95 = ci[[2]],
         mean_signed_residual = s[["mean_signed_residual"]])
})
print(as.data.frame(R8a %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
say("A ratio near 1 means the departure from a straight line over the long window is no larger than")
say("12-week variation. A ratio clearly above 1 means slower variation or curvature adds to it.")
out(R8a, "R8a_a4_line_without_week12_vs_12week")

say("\n== R8. Curvature of the mean trajectory: linear against quadratic fixed effect, all four samples ==")
say("This tests curvature in the population mean only, not curvature that differs between people.")
R8b <- imap_dfr(list(`all arms` = four, `placebo arm` = four %>% filter(arm == "Placebo")), function(d, nm) {
  lin <- if (n_distinct(d$arm) > 1) "arm * plasma_time" else "plasma_time"
  m0 <- lmer(as.formula(paste("z_log ~", lin, "+ (1 + plasma_time | RID)")), data = d, REML = FALSE, control = ctrl)
  m1 <- lmer(as.formula(paste("z_log ~", lin, "+ I(plasma_time^2) + (1 + plasma_time | RID)")), data = d,
             REML = FALSE, control = ctrl)
  a  <- anova(m0, m1)
  cf <- summary(m1)$coefficients
  tibble(group = nm, n_people = n_distinct(d$RID),
         quadratic_coef = cf["I(plasma_time^2)", "Estimate"], quadratic_se = cf["I(plasma_time^2)", "Std. Error"],
         quadratic_lo95 = quadratic_coef - 1.96 * quadratic_se, quadratic_hi95 = quadratic_coef + 1.96 * quadratic_se,
         lrt_chisq = a$Chisq[2], lrt_p = a$`Pr(>Chisq)`[2],
         singular_linear = isSingular(m0), singular_quadratic = isSingular(m1))
})
print(as.data.frame(R8b %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(R8b, "R8b_a4_curvature_four_samples")
say("The fourth sample falls after the double-blind period; check open-label treatment before")
say("interpreting the placebo-arm rows as untreated.")

