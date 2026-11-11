# =============================================================================================
# 08_sensitivity_analyses.R
# Sensitivity analyses.
#
# ADNI. Sections:
#   1. C2N pTau217 floor (as delivered, excluding always-at-floor participants, excluding floor observations,
#      half-floor substitution, cross-platform multiple imputation; T10, T10b), covariate adjustment (T11), scale
#      robustness (T12), level agreement per platform pair (T13) (formerly 08b5)
#   2. C2N floor counts by amyloid status, the %pTau217 comparison without floor samples, agreement of the two
#      C2N slopes (R5a to R5c; formerly 08b17 R5)
#   3. latent-process harmonisation per biomarker with lcmm::multlcmm (T8, T8b, T8c, F5; formerly 08b3)
#   4. multitrait-multimethod analysis of the Roche and Quanterix block (T9, T9b, T9c; formerly 08b4)
#   5. simulated power of the boundary-corrected likelihood-ratio test for a platform-specific slope variance
#      (T21, F13; formerly 08b12)
#   6. sample storage time: distribution, level regressed on storage years and the implied bias in the rate of
#      change, slopes regressed on storage (R9a to R9c; formerly 08b18)
# =============================================================================================


# =============================================================================================
# Section 1. Floor, covariates, scale, level agreement (T10 to T13)
# Original script: 08b5_paper2_sensitivity.R
# =============================================================================================
# =============================================================================
# 08b5_paper2_sensitivity.R
# Paper 2, tasks 7 to 10: sensitivity and robustness analyses for the change
# decomposition, plus level agreement between platforms.
#
#   7.  C2N pTau217 detection floor. The change model for pTau217 is refitted
#       under five treatments of the floor values:
#         as_delivered        the substituted 0.65 kept as a value (the main analysis)
#         exclude_always      people always at the floor removed (as in 08b)
#         exclude_floor_obs   only the floor observations removed
#         half_floor          floor values set to 0.65 / 2 on the raw scale
#         cross_platform_MI   floor values imputed from the same draw's ALZpath,
#                             Fujirebio and Janssen values, drawn from a normal
#                             truncated below log(0.65); M imputations, pooled
#   8.  Covariate-adjusted change model (age at first plasma, sex, APOE e4,
#       diagnosis at first plasma, each also interacting with time), per biomarker.
#   9.  Scale robustness: the change decomposition on z_raw (untransformed) and
#       on z_bl (baseline-standardised log) instead of z_log.
#   10. Level agreement between platforms at each visit: mean difference,
#       limits of agreement and ICC(2,1) per pair, on the log scale.
#
# Outputs T10_floor_sensitivity.csv, T11_covariate_adjusted.csv,
#         T12_scale_robustness.csv, T13_level_agreement.csv
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }
N_IMP <- as.integer(Sys.getenv("N_IMP", unset = "10"))
FLOOR_C2N <- 0.65

long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch))
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>% pull(sst) %>% median()

# Reduced change model (no platform-specific slope: confirmed zero in 08b2).
# Returns the components and reliabilities used throughout.
fit_change <- function(d, yvar, label, covariates = NULL) {
  fx <- "arm_f * plasma_time"
  if (!is.null(covariates)) fx <- paste(fx, "+", paste(paste0("(", covariates, ") * plasma_time"), collapse = " + "))
  f <- as.formula(paste(yvar, "~", fx, "+ (1 + plasma_time | RID_f) + (1 | RID_f:arm_f) + (1 | batch_f)"))
  fit <- tryCatch(lmer(f, data = d, REML = TRUE, control = ctrl),
                  error = function(e) { say("  %s failed: %s", label, conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  v <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  vs <- v$vcov[v$grp == "RID_f" & v$var1 == "plasma_time"]
  vi <- v$vcov[v$grp == "RID_f" & v$var1 == "(Intercept)"]
  re <- v$vcov[v$grp == "Residual"]; k <- n_distinct(d$arm)
  tibble(model = label, n_people = n_distinct(d$RID), n_obs = nrow(d), n_arms = k,
         var_person_int = vi, var_person_slope = vs, var_residual = re,
         rel_slope_one_platform = vs / (vs + re / median_sst),
         rel_slope_mean_of_arms = vs / (vs + re / (median_sst * k)),
         mean_slope = fixef(fit)[["plasma_time"]])
}

# ---- 7. C2N floor ------------------------------------------------------------------------
say("\n== 7. C2N pTau217 floor sensitivity ==")
p217 <- long %>% filter(biomarker == "pTau217")
c2n_floor_rows <- p217 %>% filter(measure == "C2N_plasma_ptau217", lloq_flag)
always_floor <- p217 %>% filter(measure == "C2N_plasma_ptau217") %>%
  group_by(RID) %>% summarise(all = all(lloq_flag), .groups = "drop") %>% filter(all) %>% pull(RID)
say("  floor observations: %d; people always at floor: %d", nrow(c2n_floor_rows), length(always_floor))

# z-scaling parameters of C2N pTau217 on the log scale, so modified values can
# be put back on the same z_log_aligned scale as the rest of the frame
c2n_mu <- mean(p217$log_value[p217$measure == "C2N_plasma_ptau217"], na.rm = TRUE)
c2n_sd <- sd(p217$log_value[p217$measure == "C2N_plasma_ptau217"], na.rm = TRUE)
rez <- function(logv) (logv - c2n_mu) / c2n_sd

floor_res <- list()
floor_res$a <- fit_change(p217, "z_log_aligned", "as_delivered")
floor_res$b <- fit_change(p217 %>% filter(!RID %in% always_floor), "z_log_aligned", "exclude_always_floor_people")
floor_res$c <- fit_change(p217 %>% filter(!(measure == "C2N_plasma_ptau217" & lloq_flag)), "z_log_aligned", "exclude_floor_observations")
floor_res$d <- fit_change(p217 %>% mutate(z_log_aligned = if_else(measure == "C2N_plasma_ptau217" & lloq_flag,
                                                                 rez(log(FLOOR_C2N / 2)), z_log_aligned)),
                          "z_log_aligned", "half_floor")

# Cross-platform multiple imputation. Fit log(C2N) ~ log(ALZpath) + log(Fujirebio)
# + log(Janssen) on the same draw among non-floor rows; for floor rows draw
# from N(prediction, residual SD) truncated above at log(0.65).
wide_p217 <- p217 %>% filter(arm %in% c("C2N", "ALZpath", "Fujirebio", "Janssen")) %>%
  select(RID, visit_no, arm, log_value, lloq_flag) %>%
  pivot_wider(names_from = arm, values_from = c(log_value, lloq_flag))
train <- wide_p217 %>% filter(!lloq_flag_C2N, !is.na(log_value_C2N))
imp_model <- lm(log_value_C2N ~ log_value_ALZpath + log_value_Fujirebio + log_value_Janssen, data = train)
say("  imputation model on %d non-floor draws: R2 = %.2f, residual SD = %.2f",
    nrow(train), summary(imp_model)$r.squared, sigma(imp_model))
floor_pred <- wide_p217 %>% filter(lloq_flag_C2N) %>%
  mutate(pred = predict(imp_model, newdata = .), sd = sigma(imp_model))
rtrunc_below <- function(n, mean, sd, upper) {
  # inverse-CDF draw from N(mean, sd) truncated at upper
  u <- runif(n, 0, pnorm(upper, mean, sd)); qnorm(u, mean, sd)
}
mi <- list()
set.seed(20260904)
for (m in seq_len(N_IMP)) {
  draws <- floor_pred %>%
    mutate(imp = rtrunc_below(n(), pred, sd, log(FLOOR_C2N)),
           imp = if_else(is.na(imp), log(FLOOR_C2N) - 0.5, imp)) %>%   # no predictors: below floor by a fixed margin
    select(RID, visit_no, imp)
  d_m <- p217 %>% left_join(draws, by = c("RID", "visit_no")) %>%
    mutate(z_log_aligned = if_else(measure == "C2N_plasma_ptau217" & lloq_flag & !is.na(imp),
                                   rez(imp), z_log_aligned)) %>% select(-imp)
  r <- fit_change(d_m, "z_log_aligned", paste0("MI_", m))
  if (!is.null(r)) mi[[m]] <- r
}
mi_tab <- bind_rows(mi)
floor_res$e <- mi_tab %>% summarise(model = sprintf("cross_platform_MI (averaged over %d imputations)", nrow(mi_tab)),
                                    n_people = first(n_people), n_obs = first(n_obs), n_arms = first(n_arms),
                                    across(c(var_person_int, var_person_slope, var_residual,
                                             rel_slope_one_platform, rel_slope_mean_of_arms, mean_slope), mean))
floor_tab <- bind_rows(floor_res)
print(as.data.frame(floor_tab %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(floor_tab, "T10_floor_sensitivity")
out(mi_tab, "T10b_floor_MI_draws")

# ---- 8. Covariate-adjusted change models -----------------------------------------------
say("\n== 8. Covariate-adjusted ==")
covs <- c("age_first_plasma", "PTGENDER", "APOE4_carrier", "DX_at_first_plasma")
cov_tab <- map_dfr(biomarkers, function(bm) {
  d <- long %>% filter(biomarker == bm, if_all(all_of(covs), ~ !is.na(.x))) %>%
    mutate(age_first_plasma = age_first_plasma - 72)
  bind_rows(fit_change(d, "z_log_aligned", "unadjusted (same people)") %>% mutate(biomarker = bm),
            fit_change(d, "z_log_aligned", "adjusted age, sex, APOE4, diagnosis", covs) %>% mutate(biomarker = bm))
})
print(as.data.frame(cov_tab %>% select(biomarker, model, n_people, var_person_slope, var_residual,
                                       rel_slope_one_platform, rel_slope_mean_of_arms) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(cov_tab, "T11_covariate_adjusted")

# ---- 9. Scale robustness ---------------------------------------------------------------
say("\n== 9. Scale robustness ==")
long <- long %>%
  group_by(measure) %>% mutate(z_raw_aligned = as.numeric(scale(value)) * direction) %>% ungroup()
scale_tab <- map_dfr(biomarkers, function(bm) {
  d <- long %>% filter(biomarker == bm)
  bind_rows(fit_change(d, "z_log_aligned", "z_log (main)"),
            fit_change(d, "z_raw_aligned", "z_raw (untransformed)"),
            fit_change(d, "z_bl_aligned",  "z_bl (baseline-standardised log)")) %>% mutate(biomarker = bm)
})
print(as.data.frame(scale_tab %>% select(biomarker, model, var_person_slope, var_residual,
                                         rel_slope_one_platform, rel_slope_mean_of_arms) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(scale_tab, "T12_scale_robustness")

# ---- 10. Level agreement between platforms at each visit ----------------------------------
say("\n== 10. Level agreement ==")
# Raw log scale (no direction flip) so mean_log_diff and fold_ratio read as
# "platform 1 / platform 2" in the units delivered; Pearson r is unaffected.
lvl <- long %>% select(RID, visit_no, biomarker, arm, log_value) %>%
  mutate(y = log_value)
agree <- map_dfr(biomarkers, function(bm) {
  w <- lvl %>% filter(biomarker == bm) %>% select(RID, visit_no, arm, y) %>%
    pivot_wider(names_from = arm, values_from = y)
  arms <- setdiff(names(w), c("RID", "visit_no"))
  pairs <- combn(arms, 2, simplify = FALSE)
  map_dfr(pairs, function(p) {
    x <- w[[p[1]]]; y <- w[[p[2]]]; ok <- complete.cases(x, y)
    x <- x[ok]; y <- y[ok]
    dd <- x - y
    # ICC(2,1) after removing each platform's mean (agreement on the centred log scale)
    xc <- x - mean(x); yc <- y - mean(y)
    icc <- tryCatch({
      f <- lmer(v ~ 1 + (1 | id), data = tibble(id = rep(seq_along(xc), 2), v = c(xc, yc)))
      vv <- as.data.frame(VarCorr(f)); vv$vcov[1] / sum(vv$vcov)
    }, error = function(e) NA_real_)
    tibble(biomarker = bm, pair = paste(p, collapse = " vs "), n = sum(ok),
           r_pearson = cor(x, y),
           mean_log_diff = mean(dd), sd_log_diff = sd(dd),
           loa_low = mean(dd) - 1.96 * sd(dd), loa_high = mean(dd) + 1.96 * sd(dd),
           fold_ratio_mean = exp(mean(dd)),          # platform 1 / platform 2 on the raw scale
           icc_centred = icc)
  })
})
print(as.data.frame(agree %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
out(agree, "T13_level_agreement")
say("\nSensitivity tables written to %s", P2_DIR)


# =============================================================================================
# Section 2. C2N floor by amyloid status and the two C2N measures (R5)
# Original script: 08b17_paper2_revision_nov2026.R, R5
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
load_change_objects()
# =============================================================================
# R5. C2N floor
# =============================================================================
say("\n== R5. C2N floor ==")
c2n_floor <- long %>% filter(measure == "C2N_plasma_ptau217") %>%
  distinct(RID, visit_no, lloq_flag, amyloid_bl_cl20)
always_floor <- c2n_floor %>% group_by(RID) %>% summarise(a = all(lloq_flag), .groups = "drop") %>%
  filter(a) %>% pull(RID)
ever_floor   <- c2n_floor %>% group_by(RID) %>% summarise(a = any(lloq_flag), .groups = "drop") %>%
  filter(a) %>% pull(RID)

# (a) counts by baseline amyloid status
R5a <- c2n_floor %>%
  mutate(amyloid = case_when(amyloid_bl_cl20 %in% TRUE ~ "amyloid positive",
                             amyloid_bl_cl20 %in% FALSE ~ "amyloid negative", TRUE ~ "unknown")) %>%
  group_by(amyloid) %>%
  summarise(observations = n(), observations_at_floor = sum(lloq_flag),
            pct_observations_at_floor = 100 * mean(lloq_flag),
            participants = n_distinct(RID),
            participants_ever_at_floor = n_distinct(RID[RID %in% ever_floor]),
            participants_always_at_floor = n_distinct(RID[RID %in% always_floor]), .groups = "drop") %>%
  bind_rows(summarise(., amyloid = "all", across(c(observations, observations_at_floor, participants,
                                                   participants_ever_at_floor, participants_always_at_floor), sum)) %>%
              mutate(pct_observations_at_floor = 100 * observations_at_floor / observations))
print(as.data.frame(R5a %>% mutate(pct_observations_at_floor = round(pct_observations_at_floor, 1))))
out(R5a, "R5a_c2n_floor_by_amyloid")

# (b) the two C2N measures in three samples, every column from the same code
floor_draws <- c2n_floor %>% filter(lloq_flag) %>% select(RID, visit_no) %>% mutate(floor_draw = TRUE)
c2n <- long %>% filter(biomarker == "pTau217", arm %in% c("C2N", "C2N_ratio")) %>%
  left_join(floor_draws, by = c("RID", "visit_no")) %>% mutate(floor_draw = coalesce(floor_draw, FALSE))
sd_log_tab <- long %>% filter(biomarker == "pTau217", arm %in% c("C2N", "C2N_ratio")) %>% group_by(arm) %>%
  summarise(sd_log = sd(log_value, na.rm = TRUE), direction = first(direction), .groups = "drop")
c2n_samples <- list(`all data` = c2n,
                    `excluding draws at the C2N floor` = c2n %>% filter(!floor_draw),
                    `participants never at the floor` = c2n %>% filter(!RID %in% ever_floor))
d24 <- designs %>% filter(schedule == "24 months, 3 draws (0, 12, 24 m)")
r5b <- list()
for (nm in names(c2n_samples)) for (a in c("C2N", "C2N_ratio")) {
  d <- c2n_samples[[nm]] %>% filter(arm == a)
  f <- tryCatch(lmer(as.formula(paste(YVAR, "~ plasma_time + (1 + plasma_time | RID_f)")),
                     data = d, control = ctrl), error = function(e) NULL)
  if (is.null(f)) next
  # OLS slopes in this sample, three draws, for the amyloid-positive mean rate (08b6 definition)
  ols <- d %>% group_by(RID) %>% filter(n() >= 3, sd(plasma_time) > 0) %>%
    summarise(slope = cov(plasma_time, z_log_aligned) / var(plasma_time),
              amy = first(amyloid_bl_cl20), .groups = "drop")
  s  <- sd_log_tab %>% filter(arm == a)
  sr <- sigma(f); vs <- vget(f, "RID_f", "plasma_time"); ms <- fixef(f)[["plasma_time"]]
  m_amy <- mean(ols$slope[ols$amy %in% TRUE], na.rm = TRUE)
  mdc_z <- 1.96 * sqrt(2) * sr
  v_trial <- vs + sr^2 / d24$sst
  r5b[[paste(nm, a)]] <- tibble(
    sample = nm, measure = a, n_people = n_distinct(d$RID), n_obs = nrow(d),
    singular = isSingular(f), sd_residual = sr, sd_person_slope = sqrt(vs),
    reliability_one_platform = vs / (vs + sr^2 / median_sst),
    mean_pct_per_year_all = 100 * (exp(ms * s$sd_log * s$direction) - 1),
    mean_pct_per_year_amy_pos = 100 * (exp(m_amy * s$sd_log * s$direction) - 1),
    cv_within_pct = 100 * sqrt(exp((sr * s$sd_log)^2) - 1),
    mdc95_rise_pct = 100 * (exp(mdc_z * s$sd_log) - 1),
    mdc95_fall_pct = 100 * (1 - exp(-mdc_z * s$sd_log)),
    years_to_exceed_mdc_amy_pos = mdc_z / abs(m_amy),
    n_per_arm_24m_3draws_25pct = n_slope(v_trial, 0.25 * abs(m_amy)))
}
R5b <- bind_rows(r5b)
print(as.data.frame(R5b %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(R5b, "R5b_c2n_ptau217_vs_ratio_floor_sensitivity")
say("Note: the 'all data' rows reproduce T3/T14 by construction except trial size, which here uses")
say("each measure's own single-platform slope variance so the three samples are comparable.")

# (c) agreement of the two C2N slopes
agree_c2n <- function(dd, nm) tibble(sample = nm, n = nrow(dd),
  r_pearson = cor(dd$slope_C2N, dd$slope_C2N_ratio),
  r_spearman = cor(dd$slope_C2N, dd$slope_C2N_ratio, method = "spearman"))
# OLS slopes recomputed within each sample, so the three rows match the R5b samples
c2n_slopes <- function(dd) dd %>% group_by(RID, arm) %>% filter(n() >= 2, sd(plasma_time) > 0) %>%
  summarise(slope = cov(plasma_time, z_log_aligned) / var(plasma_time), .groups = "drop") %>%
  pivot_wider(names_from = arm, values_from = slope, names_prefix = "slope_") %>%
  filter(!is.na(slope_C2N), !is.na(slope_C2N_ratio))
R5c <- imap_dfr(c2n_samples, ~ agree_c2n(c2n_slopes(.x), .y))
say("The 'all data' row should reproduce the T4 C2N pair (Pearson 0.765).")
print(as.data.frame(R5c %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(R5c, "R5c_c2n_slope_agreement_floor")



# =============================================================================================
# Section 3. Latent-process harmonisation (T8)
# Original script: 08b3_paper2_multlcmm.R
# =============================================================================================
# =============================================================================
# 08b3_paper2_multlcmm.R
# Paper 2, task 4: latent-process harmonisation with lcmm::multlcmm.
#
# For each biomarker, the platforms are treated as different measurements of
# one latent trajectory. Each platform gets its own link function onto the
# latent scale (linear first, then I-splines if LINK_SPLINES = TRUE) and its
# own measurement error. The person's latent intercept and slope are random
# effects. This gives:
#   T8   model fit per biomarker and link, with per-platform error SD on the
#        latent scale (the model-based analogue of T3)
#   T8b  each person's latent slope, compared with the lme4 pooled slope from
#        08b (T4b) and with the single-platform OLS slopes
#   T8c  reliability of change on the latent scale: latent slope variance
#        against the error a single platform contributes to a slope
#   F5   estimated link functions per platform (how each platform's log scale
#        maps onto the shared latent scale)
#
# Inputs: frame_paper2_long.rds (log_value per platform), paper2/T4b_pooled_person_slopes.csv
# Requires the lcmm package.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
if (!requireNamespace("lcmm", quietly = TRUE)) stop("Install lcmm first: install.packages('lcmm')")
suppressPackageStartupMessages({ library(lcmm); library(ggplot2) })

LINK_SPLINES <- as.logical(Sys.getenv("LINK_SPLINES", unset = "FALSE"))  # slower; run after linear works
P2_DIR <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds"))
pooled <- read_csv(file.path(P2_DIR, "T4b_pooled_person_slopes.csv"), show_col_types = FALSE)
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")

# Outcome per arm on the log scale, sign-aligned so that every arm increases
# with abnormality (Ab42/40 is negated). One row per person x visit.
long <- long %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         y = log_value * direction)

fit_one <- function(bm, link_type) {
  d <- long %>% filter(biomarker == bm)
  arms <- sort(unique(d$arm))
  wide <- d %>% select(RID, visit_no, plasma_time, arm, y) %>%
    pivot_wider(names_from = arm, values_from = y) %>%
    arrange(RID, visit_no) %>% as.data.frame()
  # centre each outcome so the latent scale is anchored near zero
  for (a in arms) wide[[a]] <- wide[[a]] - mean(wide[[a]], na.rm = TRUE)
  fml <- as.formula(paste(paste(arms, collapse = " + "), "~ plasma_time"))
  links <- if (link_type == "linear") rep("linear", length(arms)) else rep("4-equi-splines", length(arms))
  t0 <- Sys.time()
  m <- tryCatch(
    multlcmm(fml, random = ~ plasma_time, subject = "RID", data = wide,
             link = links, ng = 1, maxiter = 500, verbose = FALSE),
    error = function(e) { say("  %s (%s): failed: %s", bm, link_type, conditionMessage(e)); NULL })
  if (is.null(m)) return(NULL)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  say("  %s (%s): conv = %d, loglik = %.1f, AIC = %.1f, %.1f min", bm, link_type, m$conv, m$loglik, m$AIC, mins)
  list(bm = bm, link = link_type, model = m, arms = arms, wide = wide)
}

fits <- list()
for (bm in biomarkers) {
  say("\n== %s ==", bm)
  fits[[paste(bm, "linear")]] <- fit_one(bm, "linear")
  if (LINK_SPLINES) fits[[paste(bm, "splines")]] <- fit_one(bm, "splines")
}
fits <- Filter(Negate(is.null), fits)
if (!length(fits)) stop("No multlcmm model converged.")

# ---- T8: fit summary and per-platform error ---------------------------------------
# In multlcmm the latent process has a fixed residual variance of 1 by
# construction (identifiability); each outcome has its own measurement error
# SD on the transformed scale, reported in the summary as "std.err" per outcome.

re_variances <- function(m) {
  # Random-effect variances of the latent process. VarCovRE() returns a table
  # with one row per element (variance of intercept, covariance, variance of
  # slope) and the estimate in the first column. The intercept variance is
  # fixed to 1 in multlcmm for identifiability.
  v <- tryCatch(VarCovRE(m), error = function(e) NULL)
  if (is.null(v)) return(c(int = NA_real_, cov = NA_real_, slope = NA_real_))
  rn <- tolower(rownames(v)); est <- as.numeric(v[, 1])
  pick <- function(pat) { i <- grep(pat, rn); if (length(i)) est[i[1]] else NA_real_ }
  c(int = pick("var.*intercept|intercept.*var"),
    cov = pick("cov"),
    slope = pick("var.*plasma_time|plasma_time.*var"))
}

t8 <- map_dfr(fits, function(f) {
  m <- f$model
  est <- m$best
  err <- abs(est[grepl("std\\.err", names(est))])
  arms_named <- if (length(err) == length(m$Ynames)) m$Ynames else paste0("Y", seq_along(err))
  re <- re_variances(m)
  tibble(biomarker = f$bm, link = f$link, conv = m$conv, loglik = m$loglik, AIC = m$AIC, BIC = m$BIC,
         n_people = m$ns, n_obs = sum(!is.na(f$wide[, f$arms])),
         # the latent intercept variance is fixed at 1 in multlcmm (it defines the scale)
         var_latent_intercept = 1, cov_int_slope = re[["cov"]], var_latent_slope = re[["slope"]],
         mean_slope = est[["plasma_time"]],
         arm = arms_named, error_sd = as.numeric(err))
})
print(as.data.frame(t8 %>% select(biomarker, link, conv, AIC, var_latent_slope, arm, error_sd)))
out(t8, "T8_multlcmm_fit")

# ---- T8b: latent person slopes vs lme4 pooled and single-platform slopes ---------------
t8b <- map_dfr(fits, function(f) {
  m <- f$model
  re <- as.data.frame(m$predRE)             # columns: RID, intercept, plasma_time
  names(re) <- c("RID", "re_intercept", "re_slope")
  re %>% mutate(RID = as.integer(RID), biomarker = f$bm, link = f$link,
                latent_slope = re_slope + m$best[["plasma_time"]])
})
compare <- t8b %>%
  left_join(pooled %>% select(RID, biomarker, pooled_slope), by = c("RID", "biomarker")) %>%
  left_join(slopes %>% mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform)) %>%
              select(RID, biomarker, arm, slope) %>%
              pivot_wider(names_from = arm, values_from = slope, names_prefix = "ols_"),
            by = c("RID", "biomarker"))
out(compare, "T8b_latent_person_slopes")
safe_cor <- function(x, y) {
  ok <- complete.cases(x, y)
  if (sum(ok) < 10) return(NA_real_)
  cor(x[ok], y[ok])
}
agree <- compare %>%
  group_by(biomarker, link) %>%
  summarise(n = n(),
            r_latent_vs_lme4_pooled = safe_cor(latent_slope, pooled_slope),
            across(starts_with("ols_"), ~ safe_cor(latent_slope, .x), .names = "r_latent_vs_{.col}"),
            .groups = "drop") %>%
  select(where(~ !all(is.na(.x))))
print(as.data.frame(agree %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(agree, "T8b_latent_slope_agreement")

# ---- T8c: reliability of change on the latent scale ----------------------------------------
# On the latent scale a single platform's three-point slope has error variance
# error_sd^2 / SST (its measurement error spread over the person's sampling
# times); reliability is the latent slope variance against that. Note the
# platform error here also contains any person-specific offset between
# platforms, which the lme4 model absorbed into a separate random intercept,
# so these figures are not identical to T2.
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>% pull(sst) %>% median()
t8c <- t8 %>%
  mutate(err_slope_one_platform = error_sd^2 / median_sst,
         rel_slope_one_platform = var_latent_slope / (var_latent_slope + err_slope_one_platform)) %>%
  group_by(biomarker, link) %>%
  mutate(n_arms = n(),
         rel_slope_all_platforms = var_latent_slope /
           (var_latent_slope + mean(error_sd^2) / (median_sst * n_arms)),
         # precision-weighted version: the mean of platforms weighted by 1/error variance
         rel_slope_precision_weighted = var_latent_slope /
           (var_latent_slope + 1 / (median_sst * sum(1 / error_sd^2)))) %>%
  ungroup() %>%
  select(biomarker, link, arm, error_sd, var_latent_slope, rel_slope_one_platform,
         rel_slope_all_platforms, rel_slope_precision_weighted)
print(as.data.frame(t8c %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
out(t8c, "T8c_latent_reliability")

# ---- F5: link functions ------------------------------------------------------------------
pdf(file.path(P2_DIR, "F5_multlcmm_links.pdf"), width = 9, height = 6)
for (f in fits) {
  lk <- tryCatch(as.data.frame(f$model$estimlink), error = function(e) NULL)
  if (is.null(lk)) next
  # estimlink holds two columns per outcome, observed then transformed, in
  # outcome order; pair them by position rather than by name
  yn <- f$model$Ynames
  pieces <- map_dfr(seq_along(yn), function(i) {
    if (ncol(lk) < 2 * i) return(NULL)
    tibble(arm = yn[i], observed = lk[[2 * i - 1]], latent = lk[[2 * i]])
  })
  if (!nrow(pieces)) next
  p <- ggplot(pieces, aes(observed, latent, colour = arm)) + geom_line() +
    labs(title = paste(f$bm, ": link functions (", f$link, ")"),
         x = "observed (centred log scale)", y = "latent process") + theme_minimal()
  print(p)
}
dev.off()
saveRDS(lapply(fits, function(f) f$model), file.path(P2_DIR, "multlcmm_models.rds"))
say("\nmultlcmm tables written to %s", P2_DIR)


# =============================================================================================
# Section 4. Multitrait-multimethod analysis (T9)
# Original script: 08b4_paper2_mtmm.R
# =============================================================================================
# =============================================================================
# 08b4_paper2_mtmm.R
# Paper 2, task 5: multitrait-multimethod analysis of the Roche / Quanterix
# block (four analytes: pTau181, GFAP, NfL, Ab42/40; two vendors, same
# samples). The question: when a person's Roche and Quanterix values disagree,
# is the disagreement shared across the four analytes (a vendor or matrix
# effect for that person) or specific to each analyte?
#
# Model, on the log-z aligned scale with one row per person x visit x
# analyte x vendor:
#   y ~ analyte * vendor * plasma_time
#       + (1 | person)                  general level
#       + (1 | person:analyte)          analyte-specific level (trait)
#       + (1 | person:vendor)           vendor-specific level shared across analytes (method)
#       + (1 | person:analyte:vendor)   analyte x vendor specific level
#       + (1 | person:visit)            sample effect shared by everything
#       + (1 | person:visit:vendor)     sample x vendor: a per-draw vendor effect
#       + (1 | batch)
# and the same structure for slopes (person, person:analyte, person:vendor
# random slopes) in a second model. The size of the person:vendor terms is the
# MTMM method-factor result.
#
# Output  T9_mtmm_level_components.csv, T9b_mtmm_slope_components.csv,
#         T9c_mtmm_lrt.csv (likelihood tests of the vendor terms)
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(lme4))
P2_DIR <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }

long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
blk <- long %>%
  filter(platform %in% c("Roche", "Quanterix"), biomarker %in% c("pTau181", "GFAP", "NfL", "Ab42_40"),
         !is.na(z_log_aligned)) %>%
  transmute(RID, visit_no, plasma_time, analyte = biomarker, vendor = platform, batch,
            y = z_log_aligned,
            person = factor(RID), p_an = factor(paste(RID, biomarker)),
            p_ve = factor(paste(RID, platform)), p_an_ve = factor(paste(RID, biomarker, platform)),
            p_vi = factor(paste(RID, visit_no)), p_vi_ve = factor(paste(RID, visit_no, platform)),
            batch_f = factor(batch))
say("MTMM block: %d rows, %d people, %d analytes x %d vendors",
    nrow(blk), n_distinct(blk$RID), n_distinct(blk$analyte), n_distinct(blk$vendor))
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

vc_tab <- function(fit, label) {
  v <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2)) %>%
    transmute(term = ifelse(is.na(var1) | var1 == "(Intercept)", grp, paste(grp, var1, sep = ":")),
              variance = vcov) %>%
    mutate(pct = round(100 * variance / sum(variance), 1), model = label)
  v
}

# ---- 1. Level model -------------------------------------------------------------------
say("\n== 1. Level: is there a person x vendor effect shared across analytes? ==")
f_lvl <- y ~ analyte * vendor + (1 | person) + (1 | p_an) + (1 | p_ve) + (1 | p_an_ve) +
  (1 | p_vi) + (1 | p_vi_ve) + (1 | batch_f)
m_lvl <- lmer(f_lvl, data = blk, REML = TRUE, control = ctrl)
lvl <- vc_tab(m_lvl, "level") %>%
  mutate(term = recode(term, person = "person (general level)", p_an = "person x analyte (trait)",
                       p_ve = "person x vendor (method, shared across analytes)",
                       p_an_ve = "person x analyte x vendor", p_vi = "person x visit (sample)",
                       p_vi_ve = "person x visit x vendor (per-draw vendor effect)",
                       batch_f = "batch", Residual = "residual"))
print(as.data.frame(lvl))
out(lvl, "T9_mtmm_level_components")

# Likelihood tests of the two vendor terms (ML fits, boundary-corrected p)
lrt_row <- function(full, reduced, label) {
  st <- 2 * (as.numeric(logLik(full)) - as.numeric(logLik(reduced)))
  tibble(test = label, LR_stat = st, p_boundary = 0.5 * pchisq(st, 1, lower.tail = FALSE))
}
m_lvl_ml  <- lmer(f_lvl, data = blk, REML = FALSE, control = ctrl)
m_no_pve  <- update(m_lvl_ml, . ~ . - (1 | p_ve))
m_no_pvv  <- update(m_lvl_ml, . ~ . - (1 | p_vi_ve))
lrt <- bind_rows(lrt_row(m_lvl_ml, m_no_pve, "person x vendor (level)"),
                 lrt_row(m_lvl_ml, m_no_pvv, "person x visit x vendor (level)"))

# ---- 2. Slope model -------------------------------------------------------------------
say("\n== 2. Change: is there a person x vendor slope shared across analytes? ==")
f_slp <- y ~ analyte * vendor * plasma_time + (1 + plasma_time || person) +
  (1 + plasma_time || p_an) + (1 + plasma_time || p_ve) + (1 | p_an_ve) + (1 | p_vi) + (1 | batch_f)
m_slp <- lmer(f_slp, data = blk, REML = TRUE, control = ctrl)
slp <- vc_tab(m_slp, "slope") %>%
  # with the || syntax lme4 labels the slope part "<grp>:plasma_time" and the
  # intercept part "<grp>.1"
  mutate(term = str_replace_all(term, c("^person:plasma_time$" = "person slope (general change)",
                                        "^person\\.1$" = "person level",
                                        "^p_an:plasma_time$" = "person x analyte slope (trait change)",
                                        "^p_an\\.1$" = "person x analyte level",
                                        "^p_ve:plasma_time$" = "person x vendor slope (method change)",
                                        "^p_ve\\.1$" = "person x vendor level",
                                        "^p_an_ve$" = "person x analyte x vendor level",
                                        "^p_vi$" = "person x visit (sample)",
                                        "^batch_f$" = "batch", "^Residual$" = "residual")))
print(as.data.frame(slp))
out(slp, "T9b_mtmm_slope_components")
m_slp_ml <- lmer(f_slp, data = blk, REML = FALSE, control = ctrl)
m_slp_no <- lmer(update(f_slp, . ~ . - (1 + plasma_time || p_ve) + (1 | p_ve)), data = blk, REML = FALSE, control = ctrl)
lrt <- bind_rows(lrt, lrt_row(m_slp_ml, m_slp_no, "person x vendor slope (change)"))
print(as.data.frame(lrt))
out(lrt, "T9c_mtmm_lrt")

# ---- 3. Plain-language summary -----------------------------------------------------------
get <- function(df, pat) { x <- df$variance[grepl(pat, df$term)]; if (length(x)) x[1] else NA_real_ }
say("\nLevel: person x vendor shared across analytes = %.1f%% of variance; analyte-specific vendor = %.1f%%; per-draw vendor = %.1f%%",
    lvl$pct[grepl("method", lvl$term)], lvl$pct[grepl("analyte x vendor", lvl$term)],
    lvl$pct[grepl("per-draw", lvl$term)])
say("Change: person x vendor slope = %.2e against general person slope %.2e and analyte-specific slope %.2e",
    get(slp, "method change"), get(slp, "general change"), get(slp, "trait change"))
say("\nMTMM tables written to %s", P2_DIR)


# =============================================================================================
# Section 5. Power for a platform-specific slope variance (T21)
# Original script: 08b12_paper2_power_platform_slope.R
# =============================================================================================
# =============================================================================
# 08b12_paper2_power_platform_slope.R
#
# What magnitude of platform-specific slope variance could this design have
# detected?
#
# The paper's headline rests on accepting a null at a variance boundary: the
# variance of platform-specific slopes was estimated at zero for every
# biomarker, and the paper reads that as "the disagreement between platforms
# about change is random rather than systematic". A null at a boundary is only
# informative if the design had the power to reject it. This script measures
# that power.
#
# Method. Data are simulated from the reduced model's own estimated components
# for a chosen biomarker, using the real people, the real sampling times and the
# real platform structure, so the design is exactly the observed one. A
# platform-specific slope variance is then injected at a range of magnitudes,
# expressed as a fraction of the estimated true between-person slope variance.
# For each magnitude the full and reduced models are refitted and the
# boundary-corrected likelihood-ratio test is applied, exactly as in 08b2
# section 1. Power is the proportion of replications rejecting at p < 0.05.
#
# Output: T21_power_platform_slope.csv and F13_power_platform_slope.pdf/.png,
# plus a one-line summary giving the smallest ratio detectable with 80% power.
#
# Usage, from the pipeline folder:
#   Sys.setenv(POWER_BM = "pTau217", POWER_NREP = "200")
#   source("08b12_paper2_power_platform_slope.R")
# A 200-replication run over 10 magnitudes is roughly 2000 model pairs; start
# with POWER_NREP=25 to time it on your machine before committing to the full run.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages({ library(lme4); library(dplyr); library(ggplot2) })
if (!exists("BLUE"))   BLUE <- "#1F6F8B"
if (!exists("ORANGE")) ORANGE <- "#B0413E"
if (!exists("GRID"))   GRID <- "#e6e3dc"

SCRIPT_VERSION <- "0.1.0"
P2_DIR <- file.path(OUT_DIR, "paper2"); dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { readr::write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }
say <- function(fmt, ...) cat(sprintf(fmt, ...), "\n")

BM     <- Sys.getenv("POWER_BM", "pTau217")
NREP   <- as.integer(Sys.getenv("POWER_NREP", "200"))
SEED   <- as.integer(Sys.getenv("POWER_SEED", "20260908"))
# platform-specific slope SD as a fraction of the true between-person slope SD
RATIOS <- as.numeric(strsplit(Sys.getenv("POWER_RATIOS", "0,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8"), ",")[[1]])   # the grid of Table S17

YVAR <- "z_log_aligned"
long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch))
d <- long %>% filter(biomarker == BM)
stopifnot(nrow(d) > 0)

ctrl   <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
f_full <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                           "(1 + plasma_time || RID_f:arm_f) + (1|batch_f)"))
f_red  <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                           "(1 | RID_f:arm_f) + (1|batch_f)"))

# ---- 1. the observed reduced fit supplies the truth we simulate from -------------------
say("Fitting the reduced model to the observed %s data (%d rows, %d people, %d platforms)",
    BM, nrow(d), n_distinct(d$RID_f), n_distinct(d$arm_f))
m0 <- lmer(f_red, data = d, REML = TRUE, control = ctrl)
vc <- as.data.frame(VarCorr(m0))
gv <- function(grp, v1, v2 = NA) {
  r <- vc[vc$grp == grp & (is.na(vc$var1) | vc$var1 == v1) &
          (if (is.na(v2)) is.na(vc$var2) else (!is.na(vc$var2) & vc$var2 == v2)), "vcov"]
  if (!length(r)) 0 else sum(r)
}
V_int    <- gv("RID_f", "(Intercept)")
V_slope  <- gv("RID_f", "plasma_time")
COV_is   <- gv("RID_f", "(Intercept)", "plasma_time")
V_pa_int <- gv("RID_f:arm_f", "(Intercept)")
V_batch  <- gv("batch_f", "(Intercept)")
V_resid  <- attr(VarCorr(m0), "sc")^2
FIX      <- fixef(m0); X <- model.matrix(m0)
say("  true slope variance %.5f (SD %.4f); residual %.4f", V_slope, sqrt(V_slope), V_resid)

# ---- 2. simulate, inject a platform-specific slope, refit, test -------------------------
set.seed(SEED)
key_p  <- as.integer(d$RID_f)
key_pa <- as.integer(interaction(d$RID_f, d$arm_f, drop = TRUE))
key_b  <- as.integer(d$batch_f)
np <- max(key_p); npa <- max(key_pa); nb <- max(key_b)
tt <- d$plasma_time
mu <- as.vector(X %*% FIX)

L <- chol(matrix(c(V_int, COV_is, COV_is, V_slope), 2, 2) + diag(1e-10, 2))

one_rep <- function(sd_pa_slope) {
  b  <- matrix(rnorm(2 * np), np, 2) %*% L                 # person intercept and slope
  u  <- rnorm(npa, 0, sqrt(V_pa_int))                      # person by platform intercept
  s  <- rnorm(npa, 0, sd_pa_slope)                         # person by platform SLOPE, the truth
  bt <- rnorm(nb, 0, sqrt(V_batch))
  y  <- mu + b[key_p, 1] + b[key_p, 2] * tt + u[key_pa] + s[key_pa] * tt +
        bt[key_b] + rnorm(length(tt), 0, sqrt(V_resid))
  dd <- d; dd[[YVAR]] <- y
  f1 <- try(lmer(f_full, data = dd, REML = FALSE, control = ctrl), silent = TRUE)
  f2 <- try(lmer(f_red,  data = dd, REML = FALSE, control = ctrl), silent = TRUE)
  if (inherits(f1, "try-error") || inherits(f2, "try-error")) return(NA_real_)
  stat <- 2 * (as.numeric(logLik(f1)) - as.numeric(logLik(f2)))
  if (!is.finite(stat) || stat < 0) stat <- 0
  0.5 * pchisq(stat, df = 1, lower.tail = FALSE)           # 50:50 boundary mixture, as 08b2
}

res <- list()
for (r in RATIOS) {
  sd_pa <- r * sqrt(V_slope)
  p <- vapply(seq_len(NREP), function(i) one_rep(sd_pa), numeric(1))
  ok <- sum(!is.na(p))
  res[[as.character(r)]] <- tibble(
    biomarker = BM, ratio = r, sd_platform_slope = sd_pa,
    var_platform_slope = sd_pa^2, var_true_slope = V_slope,
    n_rep = NREP, n_converged = ok,
    power = mean(p[!is.na(p)] < 0.05))
  say("  ratio %.2f (platform slope SD %.4f): power %.3f from %d fits",
      r, sd_pa, res[[as.character(r)]]$power, ok)
}
T21 <- bind_rows(res)
out(T21, "T21_power_platform_slope")

# smallest ratio reaching 80% power, by linear interpolation between the bracketing points
hit <- NA_real_
for (i in seq_len(nrow(T21) - 1)) {
  a <- T21[i, ]; b <- T21[i + 1, ]
  if (a$power < 0.8 && b$power >= 0.8) {
    hit <- a$ratio + (0.8 - a$power) * (b$ratio - a$ratio) / (b$power - a$power); break
  }
}
say("\nType I error at ratio 0 (should be near 0.05): %.3f", T21$power[T21$ratio == 0])
if (is.na(hit)) {
  say("80%% power not reached within the range tested; widen POWER_RATIOS.")
} else {
  say("Smallest detectable platform-specific slope SD with 80%% power: %.2f of the true", hit)
  say("  between-person slope SD, i.e. %.4f standardised units per year.", hit * sqrt(V_slope))
}

p13 <- ggplot(T21, aes(ratio, power)) +
  geom_hline(yintercept = c(0.05, 0.8), linetype = "22", colour = GRID) +
  geom_line(colour = BLUE, linewidth = 0.7) + geom_point(colour = BLUE, size = 2) +
  { if (!is.na(hit)) geom_vline(xintercept = hit, colour = ORANGE, linetype = 3) } +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Platform-specific slope SD, as a fraction of the true between-person slope SD",
       y = "Power to reject a zero platform-specific slope variance",
       title = sprintf("What size of platform-specific slope could this design have detected? (%s)", BM),
       subtitle = sprintf("Boundary-corrected likelihood-ratio test at p < 0.05, %d replications per point", NREP))
f <- file.path(P2_DIR, "F13_power_platform_slope")
ggsave(paste0(f, ".pdf"), p13, width = 8, height = 4.5, device = "pdf", useDingbats = FALSE, bg = "white")
ggsave(paste0(f, ".png"), p13, width = 8, height = 4.5, dpi = 300, bg = "white")
say("\nWritten to %s", P2_DIR)
say("  T21_power_platform_slope.csv, F13_power_platform_slope.pdf and .png")


# =============================================================================================
# Section 6. Storage time (R9)
# Original script: 08b18_paper2_storage_check.R
# =============================================================================================
# =============================================================================
# 08b18_paper2_storage_check.R
# Paper 2, November 2026 revision, analysis A5: does storage time bias the rates of change?
#
# All ADNI plasma was assayed years after collection, so within every participant the first
# sample was stored longer than the last. If an analyte degrades at a steady rate d per year
# of storage, each measured value is lowered by d x (years stored), and because years stored
# falls by exactly one per year of follow-up, every participant's measured slope is raised by
# d. That bias is the same for everyone, so it cannot be seen by relating slopes to storage
# time; it can be seen between participants as a relation between a sample's level and how
# long it was stored. Participants were enrolled over many years, so storage at the first
# sample varies widely between them.
#
# Three outputs:
#   R9a  storage time at the first and last sample: distribution per platform.
#   R9b  level against storage, between participants: for each measure, the standardised log
#        value of the FIRST sample regressed on its storage years, adjusted for age at that
#        sample, baseline amyloid status and diagnosis. The coefficient is -d in SD units;
#        the implied slope bias is +(-coef) per year, given as a per cent of the mean annual
#        rate (all participants and amyloid positive) so its size can be judged directly.
#        Reported with the 95% CI of the coefficient and of the implied bias. Storage years is
#        the calendar date of collection (every sample was assayed in January to March 2024, in
#        collection order), so it is indistinguishable from enrolment era (ADNI1/GO/2/3) and
#        from assay batch; adding either as a covariate would absorb the storage term. The
#        coefficient cannot be attributed to degradation alone; a pattern that differs between
#        analytes in the direction known from stability studies is the only evidence that can
#        be weighed, and era and degradation of opposite sign could cancel.
#   R9c  slope against storage, between participants: each measure's three-sample OLS slope
#        and the pooled slope regressed on storage years at the first sample, with the same
#        covariates. This detects only non-linear degradation; a null here does not rule out
#        the steady-rate bias that R9b estimates. It also picks up era differences in the true
#        rate of change (later cohorts were followed for shorter spans), and pooled slopes are
#        shrunken towards the mean.
#
# Run in R:
# Fits no mixed models; runs in seconds. Changes no existing file.
# Reads pipeline_out/frame_paper2_long.rds, frame_paper2_slopes.rds, paper2/T4b, T14.
# Writes pipeline_out/paper2/R9a, R9b, R9c.
# =============================================================================

Sys.setenv(COHORT = "adni")
source("00_setup.R")
if (COHORT != "adni") stop("COHORT must be 'adni' for this script.")

P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)

rd  <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
out <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
YEAR <- 365.25

say("08b18 started %s", format(Sys.time()))

long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         storage_years = as.numeric(as.Date(RUNDATE) - as.Date(visit_date)) / YEAR,
         age_at_sample = age_first_plasma + plasma_time)
stopifnot(all(!is.na(long$storage_years)), all(long$storage_years > 0))

# ---- R9a. storage time at first and last sample ------------------------------------------------
say("\n== R9a. Storage time (years from collection to assay run) ==")
R9a <- long %>% group_by(RID, biomarker, arm) %>%
  summarise(first = storage_years[which.min(plasma_time)], last = storage_years[which.max(plasma_time)],
            .groups = "drop") %>%
  group_by(biomarker, arm) %>%
  summarise(n_people = n(),
            first_median = median(first), first_q1 = quantile(first, 0.25), first_q3 = quantile(first, 0.75),
            first_min = min(first), first_max = max(first), first_sd = sd(first),
            last_median = median(last), last_min = min(last), last_max = max(last), .groups = "drop") %>%
  arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(R9a %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
out(R9a, "R9a_storage_time_distribution")

# ---- R9b. level of the first sample against its storage time -----------------------------------
say("\n== R9b. First-sample level against storage years, between participants ==")
scale_tab <- long %>% group_by(biomarker, arm) %>%
  summarise(sd_log = sd(log_value, na.rm = TRUE), direction = first(direction), .groups = "drop")
T14 <- rd("T14_mdc_fold_change") %>% select(biomarker, arm, mean_pct_per_year, mean_pct_per_year_amy_pos)
first_s <- long %>% group_by(RID, biomarker, arm) %>% filter(plasma_time == min(plasma_time)) %>% ungroup() %>%
  filter(!is.na(amyloid_bl_cl20), !is.na(age_at_sample), !is.na(DX_at_first_plasma))
r9b <- list()
for (bm in biomarkers) for (a in sort(unique(first_s$arm[first_s$biomarker == bm]))) {
  d <- first_s %>% filter(biomarker == bm, arm == a)
  m <- lm(z_log_aligned ~ storage_years + age_at_sample + amyloid_bl_cl20 + factor(DX_at_first_plasma), data = d)
  cf <- summary(m)$coefficients["storage_years", ]
  ci <- confint(m)["storage_years", ]
  s  <- scale_tab %>% filter(biomarker == bm, arm == a)
  r9b[[paste(bm, a)]] <- tibble(
    biomarker = bm, arm = a, n = nrow(d), storage_range_years = diff(range(d$storage_years)),
    coef_per_storage_year_sd = cf[["Estimate"]], se = cf[["Std. Error"]],
    lo95 = ci[[1]], hi95 = ci[[2]], p = cf[["Pr(>|t|)"]],
    # coef is on the aligned z scale; measured slopes are raised by -coef SD/year regardless of
    # direction. On the delivered scale the level falls with storage when coef * direction < 0.
    level_pct_per_storage_year = 100 * (exp(cf[["Estimate"]] * s$sd_log * s$direction) - 1),
    implied_slope_bias_sd_per_year = -cf[["Estimate"]],
    implied_slope_bias_pct_per_year = 100 * (exp(-cf[["Estimate"]] * s$sd_log * s$direction) - 1))
}
# Bias as a share of the mean rate, computed on the z scale (inverting 08b6's per-cent formula)
# so the ratio is exact, with its 95% interval. Not interpretable for Ab42/40, whose mean rate
# is near zero.
R9b <- bind_rows(r9b) %>% left_join(T14, by = c("biomarker", "arm")) %>%
  left_join(scale_tab, by = c("biomarker", "arm")) %>%
  mutate(mean_slope_z = log(1 + mean_pct_per_year / 100) / (sd_log * direction),
         mean_slope_z_amy_pos = log(1 + mean_pct_per_year_amy_pos / 100) / (sd_log * direction),
         bias_as_pct_of_mean_rate_all = 100 * implied_slope_bias_sd_per_year / mean_slope_z,
         bias_as_pct_of_mean_rate_amy_pos = 100 * implied_slope_bias_sd_per_year / mean_slope_z_amy_pos,
         bias_lo95_pct_of_mean_rate_amy_pos = 100 * (-hi95) / mean_slope_z_amy_pos,
         bias_hi95_pct_of_mean_rate_amy_pos = 100 * (-lo95) / mean_slope_z_amy_pos) %>%
  arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(R9b %>% select(biomarker, arm, n, coef_per_storage_year_sd, lo95, hi95, p,
                                   level_pct_per_storage_year, implied_slope_bias_pct_per_year,
                                   mean_pct_per_year_amy_pos, bias_as_pct_of_mean_rate_amy_pos,
                                   bias_lo95_pct_of_mean_rate_amy_pos, bias_hi95_pct_of_mean_rate_amy_pos) %>%
                      mutate(across(where(is.numeric), ~ signif(.x, 3)))))
say("Read: the bias_lo95 and bias_hi95 columns bound the steady-rate bias as a share of the mean")
say("amyloid-positive rate, if enrolment era and degradation do not cancel. Storage years equals the")
say("collection date, so era, batch and degradation cannot be separated; compare the pattern across")
say("analytes. The ratio is not interpretable for Ab42_40, whose mean rate is near zero.")
out(R9b, "R9b_level_vs_storage")

# ---- R9c. slopes against storage at the first sample --------------------------------------------
say("\n== R9c. Slopes against storage years at the first sample, between participants ==")
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform)) %>%
  filter(n_obs >= 3)
pooled <- rd("T4b_pooled_person_slopes") %>% mutate(arm = "pooled") %>%
  select(RID, biomarker, arm, slope = pooled_slope)
cov1 <- first_s %>% distinct(RID, biomarker, arm, storage_years, age_at_sample, amyloid_bl_cl20, DX_at_first_plasma)
cov_pool <- cov1 %>% group_by(RID, biomarker) %>%
  summarise(storage_years = mean(storage_years), age_at_sample = first(age_at_sample),
            amyloid_bl_cl20 = first(amyloid_bl_cl20), DX_at_first_plasma = first(DX_at_first_plasma),
            .groups = "drop") %>% mutate(arm = "pooled")
sl <- bind_rows(slopes %>% select(RID, biomarker, arm, slope) %>% inner_join(cov1, by = c("RID", "biomarker", "arm")),
                pooled %>% inner_join(cov_pool, by = c("RID", "biomarker", "arm")))
r9c <- list()
for (bm in biomarkers) for (a in sort(unique(sl$arm[sl$biomarker == bm]))) {
  d <- sl %>% filter(biomarker == bm, arm == a)
  m <- lm(slope ~ storage_years + age_at_sample + amyloid_bl_cl20 + factor(DX_at_first_plasma), data = d)
  cf <- summary(m)$coefficients["storage_years", ]; ci <- confint(m)["storage_years", ]
  r9c[[paste(bm, a)]] <- tibble(biomarker = bm, arm = a, n = nrow(d),
    coef_slope_per_storage_year = cf[["Estimate"]], lo95 = ci[[1]], hi95 = ci[[2]], p = cf[["Pr(>|t|)"]])
}
T3 <- rd("T3_minimal_detectable_change") %>% select(biomarker, arm, sd_person_slope)
R9c <- bind_rows(r9c) %>% left_join(T3, by = c("biomarker", "arm")) %>% arrange(factor(biomarker, biomarkers), arm)
print(as.data.frame(R9c %>% mutate(across(where(is.numeric), ~ signif(.x, 3)))))
say("Slopes are in SD units per year; sd_person_slope (from T3) gives the scale of true between-person")
say("differences in rate for each measure.")
out(R9c, "R9c_slope_vs_storage")

say("\n08b18 finished %s. Files in %s", format(Sys.time()), REV_DIR)
