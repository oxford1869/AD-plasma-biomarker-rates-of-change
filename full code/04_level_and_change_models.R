# =============================================================================================
# 04_level_and_change_models.R
# Level and change: variance components, agreement, slope reliability and its bounds.
#
# ADNI. Sections:
#   1. level variance components and agreement; the multi-platform change model, slope reliability (one platform
#      and mean of platforms, overall and by amyloid status), minimal detectable change on the z scale, slope
#      agreement, validity, sampling requirement (T1 to T6; formerly 08b)
#   2. leave-one-platform-out pooled slopes and the split of the occasion-level variance into the part shared by
#      every platform at a blood draw and the assay-specific part (T19, T20; formerly 08b10)
#   3. the upper and lower bounds on the gain from averaging platforms (R1; formerly 08b17 R1)
#   4. reliability of each measure fitted alone (R4; formerly 08b17 R4)
#   5. pTau217 slope ICC on the three immunoassays (R6; formerly 08b17 R6)
# =============================================================================================


# =============================================================================================
# Section 1. Level and change models (T1 to T6)
# Original script: 08b_paper2_analysis.R
# =============================================================================================
# =============================================================================
# 08b_paper2_analysis.R
# Paper 2, first analysis pass: reliability and validity of individual change
# across platforms. Everything here reads frame_paper2_* from 07b.
#
# Sections
#   1. Variance decomposition of LEVEL (G-theory): person, platform,
#      person x platform, sample (person x occasion, shared by all platforms),
#      platform x occasion, batch within platform, residual. Per biomarker.
#   2. Variance decomposition of CHANGE: random slopes shared across platforms
#      (true individual change) versus platform-specific slopes (disagreement)
#      versus residual. Reliability of a person's slope on one platform and on
#      the average of all platforms. Per biomarker, then stratified by
#      baseline amyloid status.
#   3. Minimal detectable change per platform, and how many people show a
#      reliable change over their own span.
#   4. Agreement of per-person slopes across platforms (correlation matrix,
#      ICC) and of class labels (tertiles of slope) by kappa.
#   5. Validity: correlation of person slopes (single platform and pooled)
#      with concurrent change in Centiloid, CSF pTau/Ab42, MMSE, ADAS13,
#      hippocampal volume.
#   6. Sampling requirement: reliability of the slope as a function of number
#      of occasions and span, from the fitted variance components.
#
# Outputs go to pipeline_out/paper2/ as csv tables and pdf figures.
# All models use z_log_aligned (log, standardised within measure, higher =
# more abnormal). The C2N floor is handled by a sensitivity run that drops
# people always at the floor; see FLOOR_SENSITIVITY below.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages({ library(lme4); library(ggplot2) })

P2_DIR <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds"))
sw     <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes_wide.rds"))
critr  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))

YVAR <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
# Use measure (not platform) as the "platform" facet so the two C2N pTau217
# measures are separate arms; label them briefly.
long <- long %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio",
                         TRUE ~ platform),
         sample_id = paste(RID, visit_no, sep = "_"),
         RID_f = factor(RID), arm_f = factor(arm), sample_f = factor(sample_id),
         batch_f = factor(batch), visit_f = factor(visit_no),
         arm_visit = factor(paste(arm, visit_no, sep = "_")))

vc_table <- function(fit, label) {
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  vc %>% transmute(component = ifelse(is.na(var1) | var1 == "(Intercept)", grp, paste(grp, var1, sep = ":")),
                   variance = vcov, sd = sdcor) %>%
    mutate(pct = round(100 * variance / sum(variance), 1), model = label)
}

# ---- 1. Variance decomposition of level -----------------------------------------
say("\n== 1. Level decomposition ==")
level_vc <- list()
for (bm in biomarkers) {
  d <- long %>% filter(biomarker == bm)
  # person, arm, person x arm, sample (person x occasion shared across arms),
  # arm x occasion, batch (nested in arm by construction of the label)
  f <- as.formula(paste(YVAR, "~ 1 + (1|RID_f) + (1|arm_f) + (1|RID_f:arm_f) + (1|sample_f) +",
                        "(1|arm_visit) + (1|batch_f)"))
  fit <- tryCatch(lmer(f, data = d, REML = TRUE, control = lmerControl(optimizer = "bobyqa")),
                  error = function(e) { say("  %s: level model failed: %s", bm, conditionMessage(e)); NULL })
  if (is.null(fit)) next
  v <- vc_table(fit, bm) %>%
    mutate(component = recode(component,
      "RID_f" = "person", "arm_f" = "platform", "RID_f:arm_f" = "person x platform",
      "sample_f" = "sample (person x occasion)", "arm_visit" = "platform x occasion",
      "batch_f" = "batch within platform", "Residual" = "residual"))
  level_vc[[bm]] <- v
  say("  %s: person %.0f%%, person x platform %.0f%%, sample %.0f%%, batch %.0f%%, residual %.0f%%", bm,
      v$pct[v$component == "person"], v$pct[v$component == "person x platform"],
      v$pct[v$component == "sample (person x occasion)"], v$pct[v$component == "batch within platform"],
      v$pct[v$component == "residual"])
}
level_vc <- bind_rows(level_vc) %>% rename(biomarker = model)
out(level_vc, "T1_level_variance_components")

# G coefficients for a single measurement on one platform, relative and absolute
g_level <- level_vc %>%
  group_by(biomarker) %>%
  summarise(
    person = variance[component == "person"],
    rel_error = sum(variance[component %in% c("person x platform", "sample (person x occasion)", "residual")]),
    abs_error = sum(variance[component != "person"]),
    G_relative_one_platform = person / (person + rel_error),
    G_absolute_one_platform = person / (person + abs_error),
    .groups = "drop")
print(as.data.frame(g_level))
out(g_level, "T1b_G_coefficients_level")

# ---- 2. Variance decomposition of change -----------------------------------------
say("\n== 2. Change decomposition ==")
# Median over people of the spread of their sampling times. (Earlier versions
# rebound `.` inside a nested pipe and took the median over rows of `long`
# instead of over people, which gave 15.8 rather than 12.9 year^2 on the real
# data and overstated the one-platform slope reliabilities by about 0.05.)
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
  filter(k >= 2) %>% pull(sst) %>% median()
say("Median within-person spread of sampling times (SST): %.2f year^2", median_sst)

fit_change <- function(d, label) {
  # shared person intercept and slope; arm-specific person intercept and slope;
  # batch; arm fixed effects on intercept and slope
  f <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                        "(1 + plasma_time || RID_f:arm_f) + (1|batch_f)"))
  fit <- tryCatch(lmer(f, data = d, REML = TRUE,
                       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
                  error = function(e) { say("  %s: change model failed: %s", label, conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  get <- function(g, v) {
    i <- which(vc$grp == g & !is.na(vc$var1) & vc$var1 == v)
    if (length(i)) vc$vcov[i[1]] else NA_real_
  }
  # the || syntax makes lme4 rename the grouping with make.names(), so the grp names are
  # "RID_f.arm_f" and "RID_f.arm_f.1" (a dot, not a colon); match either (fixed 6 Nov 2026,
  # before which both person-by-platform variances were always recorded as 0)
  pa <- grep("^RID_f[.:]arm_f", unique(vc$grp), value = TRUE)
  pa_int <- if (length(pa)) sum(vapply(pa, function(g) get(g, "(Intercept)"), numeric(1)), na.rm = TRUE) else 0
  pa_slp <- if (length(pa)) sum(vapply(pa, function(g) get(g, "plasma_time"), numeric(1)), na.rm = TRUE) else 0
  res <- vc$vcov[vc$grp == "Residual"]
  tibble(
    model = label, n_people = n_distinct(d$RID), n_arms = n_distinct(d$arm),
    var_person_int = get("RID_f", "(Intercept)"), var_person_slope = get("RID_f", "plasma_time"),
    var_person_arm_int = pa_int, var_person_arm_slope = pa_slp,
    var_batch = { i <- which(vc$grp == "batch_f"); if (length(i)) vc$vcov[i[1]] else NA_real_ },
    var_residual = res[1],
    # Reliability of a person's observed slope. Two sources of error: a
    # systematic platform-specific slope (var_person_arm_slope) and the
    # occasion-level residual, which enters an OLS slope as residual / SST,
    # where SST is the spread of the person's sampling times. SST here is the
    # median over people in this frame (three points over about five years).
    sst_median = median_sst,
    rel_slope_one_platform_systematic_only = var_person_slope / (var_person_slope + var_person_arm_slope),
    rel_slope_one_platform = var_person_slope /
      (var_person_slope + var_person_arm_slope + var_residual / median_sst),
    # from the mean of k arms both error sources average out by k
    rel_slope_mean_of_arms = var_person_slope /
      (var_person_slope + (var_person_arm_slope + var_residual / median_sst) / n_arms),
    sd_true_slope_per_year = sqrt(var_person_slope),
    sd_platform_specific_slope = sqrt(var_person_arm_slope),
    singular = isSingular(fit),
    # a singular fit is a converged fit at a boundary; only other optimizer messages count
    converged = !any(!grepl("singular", fit@optinfo$conv$lme4$messages))
  ) %>% mutate(fit = list(fit))
}
change_all <- list(); change_fits <- list()
for (bm in biomarkers) {
  r <- fit_change(long %>% filter(biomarker == bm), bm)
  if (!is.null(r)) { change_fits[[bm]] <- r$fit[[1]]; change_all[[bm]] <- r %>% select(-fit) }
}
change_tab <- bind_rows(change_all) %>% rename(biomarker = model)
print(as.data.frame(change_tab %>% select(biomarker, n_people, n_arms, var_person_slope, var_person_arm_slope,
                                          var_residual, rel_slope_one_platform, rel_slope_mean_of_arms, singular)))
out(change_tab, "T2_change_variance_components")

# Stratified by baseline amyloid status
strat <- list()
for (bm in biomarkers) for (st in c(FALSE, TRUE)) {
  d <- long %>% filter(biomarker == bm, amyloid_bl_cl20 == st)
  r <- fit_change(d, paste0(bm, ifelse(st, " | amyloid+", " | amyloid-")))
  if (!is.null(r)) strat[[length(strat) + 1]] <- r %>% select(-fit)
}
strat_tab <- bind_rows(strat)
print(as.data.frame(strat_tab %>% select(model, n_people, var_person_slope, var_person_arm_slope,
                                         rel_slope_one_platform, rel_slope_mean_of_arms)))
out(strat_tab, "T2b_change_by_amyloid_status")

# Floor sensitivity: drop people always at the C2N floor from the pTau217 model
FLOOR_SENSITIVITY <- TRUE
if (FLOOR_SENSITIVITY) {
  always_floor <- long %>% filter(measure == "C2N_plasma_ptau217") %>%
    group_by(RID) %>% summarise(all_floor = all(lloq_flag), .groups = "drop") %>% filter(all_floor) %>% pull(RID)
  r <- fit_change(long %>% filter(biomarker == "pTau217", !RID %in% always_floor),
                  "pTau217 | excluding always-at-floor")
  if (!is.null(r)) out(r %>% select(-fit), "T2c_pTau217_floor_sensitivity")
  say("  pTau217 without %d always-at-floor people: rel one platform %.2f, mean of arms %.2f",
      length(always_floor), r$rel_slope_one_platform, r$rel_slope_mean_of_arms)
}

# ---- 3. Minimal detectable change per arm ---------------------------------------------
say("\n== 3. Minimal detectable change ==")
# Within-person, within-arm residual SD from a per-arm model with person
# random slopes (so residual excludes true change), on the z-log scale.
mdc <- long %>%
  group_by(biomarker, arm, measure) %>%
  group_modify(~ {
    fit <- tryCatch(lmer(as.formula(paste(YVAR, "~ plasma_time + (1 + plasma_time | RID_f)")),
                         data = .x, control = lmerControl(optimizer = "bobyqa")),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble(sd_residual = NA_real_))
    tibble(sd_residual = sigma(fit),
           sd_person_slope = sqrt(as.data.frame(VarCorr(fit))$vcov[2]),
           mean_slope = fixef(fit)[["plasma_time"]])
  }) %>% ungroup() %>%
  mutate(sem = sd_residual,
         mdc95_two_points = 1.96 * sqrt(2) * sem,        # difference between two measurements
         # the OLS slope SE is sem / sqrt(SST), SST = sum((t - mean t)^2); at the observed
         # median SST, the same design quantity every reliability in the pipeline uses
         # (until November 2026 this used the SST of three evenly spaced draws over the
         # median span, 2 (span/2)^2, which is 12.5 against the observed 12.855)
         median_span = median(slopes$span[slopes$n_obs == 3], na.rm = TRUE),
         slope_se_3pts = sem / sqrt(median_sst),
         mdc95_slope_per_year = 1.96 * slope_se_3pts)
# proportion of people whose own slope exceeds the MDC
rc <- slopes %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform)) %>%
  left_join(mdc %>% select(measure, mdc95_slope_per_year), by = "measure") %>%
  group_by(biomarker, arm, measure) %>%
  summarise(n = sum(!is.na(slope)),
            pct_reliable_increase = round(100 * mean(slope > mdc95_slope_per_year, na.rm = TRUE), 1),
            pct_reliable_decrease = round(100 * mean(slope < -mdc95_slope_per_year, na.rm = TRUE), 1),
            .groups = "drop")
mdc <- mdc %>% left_join(rc, by = c("biomarker", "arm", "measure"))
print(as.data.frame(mdc %>% select(biomarker, arm, sd_residual, mdc95_two_points, mdc95_slope_per_year,
                                   pct_reliable_increase, pct_reliable_decrease)))
out(mdc, "T3_minimal_detectable_change")

# ---- 4. Agreement of slopes across arms -------------------------------------------------
say("\n== 4. Agreement of slopes ==")
agree <- list()
for (bm in biomarkers) {
  m <- sw %>% filter(biomarker == bm) %>% select(starts_with("slope_")) %>% select(where(~ !all(is.na(.x))))
  if (ncol(m) < 2) next
  cm <- cor(m, use = "pairwise.complete.obs", method = "pearson")
  sm <- cor(m, use = "pairwise.complete.obs", method = "spearman")
  p <- which(upper.tri(cm), arr.ind = TRUE)
  # ICC(2,1) across arms on complete rows
  mc <- m[complete.cases(m), ]
  icc <- tryCatch({
    dd <- mc %>% mutate(id = row_number()) %>% pivot_longer(-id, names_to = "arm", values_to = "s")
    f <- lmer(s ~ 1 + (1|id) + (1|arm), data = dd)
    v <- as.data.frame(VarCorr(f)); vi <- v$vcov[v$grp == "id"]; va <- v$vcov[v$grp == "arm"]; ve <- v$vcov[v$grp == "Residual"]
    c(icc21 = vi / (vi + va + ve), icc2k = vi / (vi + (va + ve) / ncol(mc)))
  }, error = function(e) c(icc21 = NA, icc2k = NA))
  # kappa of slope tertiles between arms
  tert <- mc %>% mutate(across(everything(), ~ ntile(.x, 3)))
  kap <- sapply(seq_len(nrow(p)), function(i) {
    a <- tert[[p[i, 1]]]; b <- tert[[p[i, 2]]]
    po <- mean(a == b); pe <- sum(prop.table(table(a)) * prop.table(table(b))); (po - pe) / (1 - pe)
  })
  agree[[bm]] <- tibble(biomarker = bm,
                        pair = paste(str_remove(colnames(cm)[p[, 1]], "slope_"),
                                     str_remove(colnames(cm)[p[, 2]], "slope_"), sep = " vs "),
                        r_pearson = round(cm[p], 3), r_spearman = round(sm[p], 3),
                        kappa_tertiles = round(kap, 3), n_complete = nrow(mc),
                        icc_single_arm = round(icc[["icc21"]], 3), icc_mean_of_arms = round(icc[["icc2k"]], 3))
}
agree <- bind_rows(agree)
print(as.data.frame(agree))
out(agree, "T4_slope_agreement")

# Shared (pooled) slope per person from the multi-arm model, for validity
pooled <- list()
for (bm in names(change_fits)) {
  re <- ranef(change_fits[[bm]])$RID_f
  fe <- fixef(change_fits[[bm]])
  pooled[[bm]] <- tibble(RID = as.integer(rownames(re)), biomarker = bm,
                         pooled_slope = re[["plasma_time"]] + fe[["plasma_time"]],
                         pooled_intercept = re[["(Intercept)"]] + fe[["(Intercept)"]])
}
pooled <- bind_rows(pooled)
out(pooled, "T4b_pooled_person_slopes")

# ---- 5. Validity against criterion change ------------------------------------------------
say("\n== 5. Validity ==")
crit_cols <- c("centiloid_slope", "csf_ptau_abeta42_slope", "mmse_slope", "adas13_slope",
               "mpacc_digit_slope", "hcv_adj_slope")
crit_cols <- crit_cols[crit_cols %in% names(critr)]
val_single <- slopes %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform)) %>%
  select(RID, biomarker, arm, slope) %>%
  left_join(critr %>% select(RID, all_of(crit_cols)), by = "RID") %>%
  pivot_longer(all_of(crit_cols), names_to = "criterion", values_to = "crit_slope") %>%
  group_by(biomarker, arm, criterion) %>%
  summarise(n = sum(complete.cases(slope, crit_slope)),
            r = round(cor(slope, crit_slope, use = "complete.obs", method = "spearman"), 3), .groups = "drop") %>%
  mutate(estimate = "single measure")
val_pooled <- pooled %>%
  left_join(critr %>% select(RID, all_of(crit_cols)), by = "RID") %>%
  pivot_longer(all_of(crit_cols), names_to = "criterion", values_to = "crit_slope") %>%
  group_by(biomarker, criterion) %>%
  summarise(n = sum(complete.cases(pooled_slope, crit_slope)),
            r = round(cor(pooled_slope, crit_slope, use = "complete.obs", method = "spearman"), 3), .groups = "drop") %>%
  mutate(arm = "pooled across the five measures", estimate = "pooled")   # ADNI has no trial arms: the pooling is across the five pTau217 measures
# Orient the signs so that a positive correlation always means the plasma slope tracks
# worsening. Hippocampal volume, mPACC and MMSE fall as disease progresses, so their raw
# correlations carry the opposite sign. Figure S4 already applied this flip (08b7 line 252)
# while the table did not, so the table and its own note disagreed.
FLIP <- c("hcv_adj_slope", "mpacc_digit_slope", "mmse_slope")
validity <- bind_rows(val_single, val_pooled) %>%
  mutate(r = ifelse(criterion %in% FLIP, -r, r)) %>%
  arrange(biomarker, criterion, arm)
print(as.data.frame(validity %>% filter(criterion == "centiloid_slope")))
out(validity, "T5_validity_vs_criterion_change")

# ---- 6. Sampling requirement --------------------------------------------------------------
say("\n== 6. Sampling requirement ==")
# For a design with k equally spaced occasions over span S years, the OLS
# slope error variance on one platform is (arm-specific slope var) +
# residual / sum((t - mean t)^2). Reliability = true / (true + error).
grid <- expand_grid(k = c(2, 3, 4, 5, 6), span_years = c(2, 4, 6, 8, 10))
sampling <- change_tab %>%
  select(biomarker, var_person_slope, var_person_arm_slope, var_residual, n_arms) %>%
  crossing(grid) %>%
  mutate(sst = map2_dbl(k, span_years, ~ { t <- seq(0, .y, length.out = .x); sum((t - mean(t))^2) }),
         err_one_platform = var_person_arm_slope + var_residual / sst,
         reliability_one_platform = var_person_slope / (var_person_slope + err_one_platform),
         err_mean_of_arms = var_person_arm_slope / n_arms + var_residual / (sst * n_arms),
         reliability_mean_of_arms = var_person_slope / (var_person_slope + err_mean_of_arms)) %>%
  select(biomarker, k, span_years, reliability_one_platform, reliability_mean_of_arms)
out(sampling, "T6_sampling_requirement")
print(as.data.frame(sampling %>% filter(biomarker == "pTau217") %>%
                      mutate(across(starts_with("rel"), ~ round(.x, 2)))))

# ---- Figures ------------------------------------------------------------------------------
# F1: variance components of change
f1 <- change_tab %>%
  select(biomarker, var_person_slope, var_person_arm_slope) %>%
  pivot_longer(-biomarker, names_to = "component", values_to = "variance") %>%
  mutate(component = recode(component, var_person_slope = "shared across platforms (true change)",
                            var_person_arm_slope = "platform-specific (disagreement)")) %>%
  ggplot(aes(biomarker, variance, fill = component)) + geom_col(position = "fill") +
  labs(y = "share of slope variance", x = NULL, fill = NULL) + theme(legend.position = "bottom")
ggsave(file.path(P2_DIR, "F1_slope_variance_shares.pdf"), f1, width = 7, height = 4)
# F2: slope agreement heatmap for pTau217
if (any(agree$biomarker == "pTau217")) {
  a <- agree %>% filter(biomarker == "pTau217") %>%
    separate(pair, c("a", "b"), sep = " vs ")
  f2 <- bind_rows(a, a %>% rename(a = b, b = a)) %>%
    ggplot(aes(a, b, fill = r_pearson)) + geom_tile() +
    geom_text(aes(label = sprintf("%.2f", r_pearson)), size = 3) +
    scale_fill_gradient2(limits = c(-1, 1)) + labs(x = NULL, y = NULL, fill = "r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(P2_DIR, "F2_pTau217_slope_agreement.pdf"), f2, width = 5.5, height = 4.5)
}
# F3: sampling requirement curves
f3 <- ggplot(sampling, aes(span_years, reliability_one_platform, colour = factor(k))) +
  geom_line() + geom_point() + facet_wrap(~ biomarker) +
  geom_hline(yintercept = 0.7, linetype = 2) +
  labs(x = "span (years)", y = "reliability of one-platform slope", colour = "occasions")
ggsave(file.path(P2_DIR, "F3_sampling_requirement.pdf"), f3, width = 8, height = 5)

say("\nPaper 2 tables and figures written to %s", P2_DIR)


# =============================================================================================
# Section 2. Leave-one-platform-out and the occasion-variance split (T19, T20)
# Original script: 08b10_paper2_loo_and_variance_split.R
# =============================================================================================
# =============================================================================
# 08b10_paper2_loo_and_variance_split.R
# Two additions to Paper 2, both from material already in the frames. No model refitting beyond the
# change models, so this runs in minutes.
#
# Note on what is NOT here. An earlier version of this script also compared the plasma slope against a
# single baseline amyloid PET scan as predictors of the same outcome. That analysis was removed, not
# merely unreported: its incident-amyloid endpoint is defined by crossing a Centiloid threshold, which
# a baseline Centiloid value predicts partly by construction, so the comparison is unfair in both
# directions and is not a question this design can answer.
#
# Section 1. Leave one platform out.
#   The pooled slope is one of the paper's two headline claims, so it should not rest on any
#   single platform. The reduced change model is refitted for each biomarker with each
#   platform dropped in turn. Reported per dropped platform: the correlation of the resulting
#   pooled slope with the pooled slope from all platforms, the single-assay reliability implied by the reduced set, and its validity against concurrent
#   Centiloid change.
#
# Section 2. Is the change-model residual measurement error, or real departure from a straight line?
#   With three draws per person there is one residual degree of freedom, so within a person the two
#   cannot be told apart. Across platforms they can. Anything that moves all five platforms together
#   at the same draw, which is curvature in the person's true trajectory, that day's biology and
#   pre-analytical handling, is shared; assay noise is not. Adding a sample (person by occasion)
#   random effect to the change model therefore splits the residual into a shared part and an
#   assay-specific part. Reported per biomarker: the two variances, the shared share of the total,
#   and the slope reliability computed three ways, against the total occasion-level variance (the
#   value reported in the paper), against the assay-specific part alone (an upper bound, treating
#   everything shared as true signal) and against the shared part alone.
#
#   Read it as follows. A small shared part means the residual is assay noise, the reliability in
#   the paper is the right number, and averaging platforms helps as claimed. A large shared part
#   would mean that some of what the paper calls noise is real within-person variation, that
#   averaging platforms helps less than claimed, and that a linear slope is partly the wrong summary
#   rather than a noisy one. The A4 twelve-week re-draw answers the same question from outside and
#   the two should agree.
#
# Run   Sys.setenv(COHORT = "adni"); source("08b10_paper2_loo_and_variance_split.R")
#       ADNI only: it needs several platforms on the same people.
# Outputs (pipeline_out/paper2/)
#   T19_leave_one_platform_out.csv pooled slope with each platform dropped
#   T20_shared_vs_assay_occasion_variance.csv  split of the occasion-level variance
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
if (COHORT != "adni") stop("08b10 is ADNI only: it needs several platforms measuring the same people.\n",
                           "  Sys.setenv(COHORT = \"adni\") and run again.")
suppressPackageStartupMessages({ library(lme4); library(ggplot2) })

P2_DIR <- file.path(OUT_DIR, "paper2"); dir.create(P2_DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds"))
critr  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
YVAR <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
FOCUS_BM <- "pTau217"
CV_FOLDS <- as.integer(Sys.getenv("P2_CV_FOLDS", "10"))
set.seed(as.integer(Sys.getenv("P2_SEED", "20260907")))

long <- long %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         sample_f = factor(paste(RID, visit_no, sep = "_")),
         RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch))


# ---- section 1: leave one platform out ----------------------------------------------------------
say("\n== 1. Leave one platform out of the pooled slope ==")
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
  filter(k >= 2) %>% pull(sst) %>% median()

fit_reduced <- function(d) {
  # the reduced change model used throughout Paper 2: no platform-specific slope, since that
  # variance was estimated at zero for every biomarker
  f <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                        "(1 | RID_f:arm_f) + (1|batch_f)"))
  tryCatch(lmer(f, data = d, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
           error = function(e) NULL)
}
person_slope <- function(fit) {
  re <- ranef(fit)$RID_f
  tibble(RID = as.integer(rownames(re)), s = re[["plasma_time"]] + fixef(fit)[["plasma_time"]])
}
rel_of <- function(fit) {
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  vs <- vc$vcov[vc$grp == "RID_f" & !is.na(vc$var1) & vc$var1 == "plasma_time"][1]
  vr <- vc$vcov[vc$grp == "Residual"][1]
  vs / (vs + vr / median_sst)
}
cs <- critr %>% select(RID, centiloid_slope, centiloid_n)

loo <- list()
for (bm in biomarkers) {
  d_all <- long %>% filter(biomarker == bm)
  arms <- sort(unique(d_all$arm))
  if (length(arms) < 3) { say("  %s has %d platforms: leaving one out is not informative, skipped",
                              bm, length(arms)); next }
  f_all <- fit_reduced(d_all %>% mutate(arm_f = droplevels(factor(arm))))
  if (is.null(f_all)) { say("  %s: full model failed", bm); next }
  s_all <- person_slope(f_all)
  val <- function(s) {
    j <- s %>% inner_join(cs %>% filter(centiloid_n >= 2), by = "RID")
    if (nrow(j) < 30) return(NA_real_)
    suppressWarnings(cor(j$s, j$centiloid_slope, method = "spearman", use = "complete.obs"))
  }
  loo[[paste(bm, "none")]] <- tibble(biomarker = bm, dropped = "none (all platforms)",
    n_platforms = length(arms), r_with_full_pooled = 1,
    reliability_one_assay = rel_of(f_all), validity_vs_centiloid = val(s_all))
  for (a in arms) {
    d <- d_all %>% filter(arm != a) %>% mutate(arm_f = droplevels(factor(arm)))
    f <- fit_reduced(d)
    if (is.null(f)) { say("  %s without %s: model failed", bm, a); next }
    s <- person_slope(f)
    j <- inner_join(s_all %>% rename(s_all = s), s, by = "RID")
    loo[[paste(bm, a)]] <- tibble(biomarker = bm, dropped = a, n_platforms = length(arms) - 1,
      r_with_full_pooled = suppressWarnings(cor(j$s_all, j$s, use = "complete.obs")),
      reliability_one_assay = rel_of(f), validity_vs_centiloid = val(s))
  }
  say("  %s done (%d platforms)", bm, length(arms))
}
T19 <- bind_rows(loo) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
out(T19, "T19_leave_one_platform_out")
print(as.data.frame(T19))
say("\nThe pooled slope is robust if r_with_full_pooled stays high and the validity column")
say("does not depend on any one platform being present.")

# ---- section 2: shared against assay-specific occasion variance ---------------------------------
say("\n== 2. Is the residual assay noise or real departure from a straight line? ==")
fit_sample <- function(d) {
  # the reduced change model plus a sample effect: the person by occasion deviation shared by every
  # platform that assayed that draw
  f <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                        "(1 | RID_f:arm_f) + (1 | sample_f) + (1|batch_f)"))
  tryCatch(lmer(f, data = d, REML = TRUE,
                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))),
           error = function(e) NULL)
}
getv <- function(fit, grp, var1 = NA) {
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  i <- if (is.na(var1)) which(vc$grp == grp) else which(vc$grp == grp & !is.na(vc$var1) & vc$var1 == var1)
  if (length(i)) vc$vcov[i[1]] else NA_real_
}
share <- list()
for (bm in biomarkers) {
  d <- long %>% filter(biomarker == bm) %>% mutate(arm_f = droplevels(factor(arm)))
  n_arm <- n_distinct(d$arm)
  if (n_arm < 2) next
  f0 <- fit_reduced(d); f1 <- fit_sample(d)
  if (is.null(f0) || is.null(f1)) { say("  %s: model failed", bm); next }
  # the reported reliability comes entirely from the model the paper uses
  v_true0 <- getv(f0, "RID_f", "plasma_time"); v_res0 <- getv(f0, "Residual")
  # the decomposition comes entirely from the model with the sample effect
  v_true1 <- getv(f1, "RID_f", "plasma_time")
  v_shared <- getv(f1, "sample_f"); v_assay <- getv(f1, "Residual")
  lrt <- tryCatch({
    dv <- 2 * (as.numeric(logLik(refitML(f1))) - as.numeric(logLik(refitML(f0))))
    0.5 * pchisq(dv, df = 1, lower.tail = FALSE)
  }, error = function(e) NA_real_)
  rel <- function(vt, ve) vt / (vt + ve / median_sst)
  share[[bm]] <- tibble(
    biomarker = bm, n_platforms = n_arm,
    singular = isSingular(f1),
    identifiable = n_arm >= 3 && !isSingular(f1),
    var_true_slope_paper_model = v_true0, var_residual_paper_model = v_res0,
    var_true_slope_sample_model = v_true1,
    var_shared_at_occasion = v_shared, var_assay_specific = v_assay,
    pct_shared = 100 * v_shared / (v_shared + v_assay),
    p_lrt_shared = lrt,
    # reliability of a one-assay slope, as reported in the paper
    rel_one_assay_paper = rel(v_true0, v_res0),
    # the same quantity under the sample model, for comparison of the two models
    rel_one_assay_sample_model = rel(v_true1, v_shared + v_assay),
    # reliability of the mean of the assays, under the paper's assumption that all
    # occasion-level variance averages out, and under the assumption that the shared part does not
    rel_pooled_all_averaged = rel(v_true1, (v_shared + v_assay) / n_arm),
    rel_pooled_shared_retained = rel(v_true1, v_shared + v_assay / n_arm))
  say("  %s (%d platforms%s): shared %.0f%% of occasion variance; one-assay reliability %.2f (paper model) vs %.2f (sample model); pooled %.2f if all averages out, %.2f if the shared part does not",
      bm, n_arm, if (isSingular(f1)) ", singular fit" else "", share[[bm]]$pct_shared,
      share[[bm]]$rel_one_assay_paper, share[[bm]]$rel_one_assay_sample_model,
      share[[bm]]$rel_pooled_all_averaged, share[[bm]]$rel_pooled_shared_retained)
}
T20 <- bind_rows(share) %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
out(T20, "T20_shared_vs_assay_occasion_variance")
print(as.data.frame(T20))
say("")
say("How to read this. With three occasions per person, a person by occasion effect and a person")
say("random slope compete for the same variation, so the sample model understates the true slope")
say("variance and its shared share is an upper bound. Rows with fewer than three platforms, or with")
say("a singular fit, are not identifiable and are flagged as such; only those with identifiable TRUE")
say("carry weight. The quantity that matters for the paper is the pair of pooled reliabilities: the")
say("first assumes every occasion-level deviation is assay-specific and averages out, the second")
say("assumes the shared part does not average out. The truth lies between them, and the observed")
say("validity of the pooled slope against Centiloid change indicates which end it sits nearer.")

say("\nWritten to %s: T19, T20", P2_DIR)


# =============================================================================================
# Section 3. Averaging across platforms: upper and lower bound (R1)
# Original script: 08b17_paper2_revision_nov2026.R, R1
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
T20 <- rd("T20_shared_vs_assay_occasion_variance")
# =============================================================================
# R1. Averaging across platforms: upper and lower bound
# =============================================================================
say("\n== R1. Averaging across platforms: upper and lower bound ==")
R1 <- T2 %>% select(biomarker, n_arms, rel_one_upper = rel_slope_one_platform,
                    rel_mean_upper = rel_slope_mean_of_arms) %>%
  left_join(T20 %>% select(biomarker, identifiable, singular_sample_model = singular, pct_shared,
                           rel_one_lower = rel_one_assay_sample_model,
                           rel_mean_lower = rel_pooled_shared_retained,
                           rel_mean_lower_if_all_averaged = rel_pooled_all_averaged),
            by = "biomarker") %>%
  mutate(gain_upper = rel_mean_upper - rel_one_upper,
         gain_lower = rel_mean_lower - rel_one_lower,
         note = if_else(identifiable %in% TRUE,
                        if_else(biomarker == "pTau217",
                                "lower bound estimable; C2N pTau217 and %pTau217 share one pTau217 value, so both bounds overstate the averaging gain somewhat (see R6)",
                                "lower bound estimable"),
                        "lower bound not identifiable (fewer than three platforms or singular fit): report the upper bound with this caveat")) %>%
  mutate(biomarker = factor(biomarker, biomarkers)) %>% arrange(biomarker)
print(as.data.frame(R1 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(R1, "R1_averaging_bounds")



# =============================================================================================
# Section 4. One-platform reliability from each measure's own model (R4)
# Original script: 08b17_paper2_revision_nov2026.R, R4
# =============================================================================================
T3  <- rd("T3_minimal_detectable_change")
# =============================================================================
# R4. One-platform reliability from each measure's own single-platform model (T3 model)
# =============================================================================
say("\n== R4. One-platform reliability from single-platform models ==")
A4_SST <- as.numeric(Sys.getenv("P2_ANCHOR_SST", "13.533"))   # A4 placebo SST used in the paper
R4 <- T3 %>% select(biomarker, arm, measure, sd_person_slope, sd_residual) %>%
  mutate(rel_at_adni_median_sst = sd_person_slope^2 / (sd_person_slope^2 + sd_residual^2 / median_sst),
         rel_at_a4_placebo_sst  = sd_person_slope^2 / (sd_person_slope^2 + sd_residual^2 / A4_SST)) %>%
  left_join(T2 %>% select(biomarker, rel_joint_model_one_platform = rel_slope_one_platform), by = "biomarker") %>%
  arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(R4 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(R4, "R4_single_platform_reliability_adni")



# =============================================================================================
# Section 5. pTau217 slope ICC, immunoassays only (R6)
# Original script: 08b17_paper2_revision_nov2026.R, R6
# =============================================================================================
# =============================================================================
# R6. pTau217 slope ICC, immunoassays only
# =============================================================================
say("\n== R6. pTau217 slope ICC, immunoassays only ==")
icc_of <- function(cols) {
  m <- sw %>% filter(biomarker == "pTau217") %>% select(all_of(cols))
  mc <- m[complete.cases(m), ]
  dd <- mc %>% mutate(id = row_number()) %>% pivot_longer(-id, names_to = "arm", values_to = "s")
  f <- lmer(s ~ 1 + (1 | id) + (1 | arm), data = dd)
  v <- as.data.frame(VarCorr(f)); vi <- v$vcov[v$grp == "id"]; va <- v$vcov[v$grp == "arm"]
  ve <- v$vcov[v$grp == "Residual"]
  tibble(measures = paste(sub("slope_", "", cols), collapse = ", "), n_complete = nrow(mc),
         icc_single = vi / (vi + va + ve), icc_mean = vi / (vi + (va + ve) / length(cols)))
}
R6 <- bind_rows(icc_of(c("slope_ALZpath", "slope_C2N", "slope_C2N_ratio", "slope_Fujirebio", "slope_Janssen")),
                icc_of(c("slope_ALZpath", "slope_Fujirebio", "slope_Janssen")))
print(as.data.frame(R6 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
say("The first row should reproduce the pTau217 slope ICC in Table 2 (0.44 / 0.80).")
out(R6, "R6_ptau217_slope_icc_immunoassays")

