# =============================================================================================
# 09_bootstrap_intervals.R
# The participant bootstrap and every interval and paired comparison derived from it.
#
# ADNI (the A4/LEARN bootstrap is in 07, section 1). 1,000 participant resamples (seed 20260904 + draw), saved
# draw by draw and resumed if interrupted. N_BOOT = 0 skips the resampling in sections 1 and 11; sections 3 to 6
# and 8 to 10 need the saved draws, so a test run without a bootstrap should stop at STEPS 08. Sections:
#   1. likelihood-ratio test of the platform-specific slope variance (T2d), residual diagnostics (T2e, T2f, F4), and
#      the bootstrap of the G coefficients, reliabilities, slope ICCs and validity correlations (T7) (formerly 08b2)
#   2. intervals for the minimal detectable change, years to exceed it and trial sample sizes (T3b, T15b; formerly 08b16)
#   3. intervals for the averaging lower bound, each measure fitted alone, the amyloid-positive trial sizes and the
#      immunoassay ICC (R10; formerly 08b19)
#   4. intervals for the sampling span needed to reach a target reliability (R11; formerly 08b20)
#   5. paired differences and ratios within resamples for the head-to-head comparisons in the text, and the
#      boundary components of the singular amyloid-positive fits (R12, R12b; formerly 08b21)
#   6. the amyloid-positive minus amyloid-negative difference in pTau217 validity, within draws (printed only;
#      formerly p2_amyloid_interval)
#   7. observed-data validity of the pooled slope by amyloid group and Pearson versions of the design-test
#      correlations, both cohorts (T7c, T17p; A4 T8p, T8bp; formerly p2_observed_checks)
#   8. the paired intervals quoted in the text and in Table S26, computed within the saved draws: averaging gain,
#      between-biomarker and amyloid-group differences (R13a), platform differences in CV, MDC and years to exceed
#      (R13b), and the shares of singular sample-model fits (R13c)
#   9. the reduction in trial size from five rather than three samples over 24 months, from the section 3 draws (R15a)
#  10. the pooled minus single-platform pTau217 validity against Centiloid change, paired within the section 1
#      resamples (R15b)
#  11. intervals for the whole-cohort trial sizes of the pTau217 mean of platforms under the lower bound, from the
#      same resamples refitted (R17)
# =============================================================================================


# =============================================================================================
# Section 1. Platform-specific slope test, diagnostics, the participant bootstrap (T2d to T7)
# Original script: 08b2_paper2_checks_bootstrap.R
# =============================================================================================
# =============================================================================
# 08b2_paper2_checks_bootstrap.R
# Paper 2, tasks 1 to 3 of the task list.
#
#   1. Likelihood-ratio test of the platform-specific slope variance: the
#      change model with and without the (plasma_time | person:platform)
#      term, per biomarker, fitted by ML. Reports the test and which model is
#      used for the bootstrap.
#   2. Residual diagnostics for the change model: residual SD by platform, by
#      fitted level (tertiles), and by floor status; QQ summary. Writes one
#      pdf of residual plots per biomarker.
#   3. Cluster bootstrap by person for the quantities the paper reports:
#      G coefficients (T1b), slope reliabilities (T2, T2b), slope ICCs (T4),
#      pooled-slope validity correlations (T5). Percentile 95% intervals.
#
# The bootstrap is slow (five mixed models per draw). It saves every draw to
# pipeline_out/paper2/bootstrap_draws.rds as it goes and resumes from there,
# so it can be stopped and restarted. Run it from a terminal so RStudio stays
# free:
# With N_BOOT=0 only sections 1 and 2 run (a few minutes).
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages({ library(lme4); library(ggplot2); library(tidyr) })
# the shared palette and theme, so Supplementary Figure S5 matches the rest of the set
if (!exists("GRID")) GRID <- "#e6e3dc"

N_BOOT   <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))   # the paper uses 1,000 draws; 0 runs the checks only
BOOT_SEED <- 20260904
P2_DIR <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
critr  <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
YVAR <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
long <- long %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch))
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
  pull(sst) %>% median()

ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
f_full <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                           "(1 + plasma_time || RID_f:arm_f) + (1|batch_f)"))
f_red  <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                           "(1 | RID_f:arm_f) + (1|batch_f)"))

# ---- 1. Likelihood-ratio test for the platform-specific slope ---------------------
say("\n== 1. Is there any platform-specific slope variance? ==")
lrt <- list()
for (bm in biomarkers) {
  d <- long %>% filter(biomarker == bm)
  m1 <- lmer(f_full, data = d, REML = FALSE, control = ctrl)
  m0 <- lmer(f_red,  data = d, REML = FALSE, control = ctrl)
  ll1 <- as.numeric(logLik(m1)); ll0 <- as.numeric(logLik(m0))
  stat <- 2 * (ll1 - ll0)
  # variance on the boundary: p from a 50:50 mixture of chi-square(0) and chi-square(1)
  p <- 0.5 * pchisq(stat, df = 1, lower.tail = FALSE)
  vc <- as.data.frame(VarCorr(m1))
  # with || lme4 names the grouping "RID_f.arm_f" (dot), so match either (fixed 6 Nov 2026)
  pa_slope <- vc$vcov[grepl("^RID_f[.:]arm_f", vc$grp) & !is.na(vc$var1) & vc$var1 == "plasma_time" & is.na(vc$var2)]
  lrt[[bm]] <- tibble(biomarker = bm, var_person_arm_slope_full = sum(pa_slope),
                      logLik_full = ll1, logLik_reduced = ll0, LR_stat = stat, p_boundary = p,
                      AIC_full = AIC(m1), AIC_reduced = AIC(m0))
  say("  %s: LR = %.3f, p = %.3f (boundary-corrected); platform-specific slope var = %.2e",
      bm, stat, p, sum(pa_slope))
}
lrt <- bind_rows(lrt)
out(lrt, "T2d_platform_slope_LRT")
USE_REDUCED <- all(lrt$p_boundary > 0.05)
say("Bootstrap will use the %s model", ifelse(USE_REDUCED, "reduced (no platform-specific slope)", "full"))

# ---- 2. Residual diagnostics ----------------------------------------------------------
say("\n== 2. Residual diagnostics ==")
diag_rows <- list()
res_rows  <- list()          # per-observation residuals, kept so the figure can be redrawn
pdf(file.path(P2_DIR, "F4_residual_diagnostics.pdf"), width = 10, height = 7)
for (bm in biomarkers) {
  d <- long %>% filter(biomarker == bm)
  fit <- lmer(if (USE_REDUCED) f_red else f_full, data = d, REML = TRUE, control = ctrl)
  d$res <- resid(fit); d$fit <- fitted(fit)
  d <- d %>% mutate(fit_tertile = ntile(fit, 3),
                    floor = if_else(lloq_flag, "at floor", "measured"))
  by_arm   <- d %>% group_by(arm) %>% summarise(sd_res = sd(res), n = n(), .groups = "drop") %>% mutate(facet = "platform", level = arm)
  by_fit   <- d %>% group_by(fit_tertile) %>% summarise(sd_res = sd(res), n = n(), .groups = "drop") %>%
    mutate(facet = "fitted tertile", level = paste0("T", fit_tertile)) %>% select(-fit_tertile)
  by_floor <- d %>% group_by(floor) %>% summarise(sd_res = sd(res), n = n(), .groups = "drop") %>% mutate(facet = "floor", level = floor) %>% select(-floor)
  qq <- tibble(facet = "shape", level = c("skewness", "excess kurtosis"),
               sd_res = c(mean((d$res - mean(d$res))^3) / sd(d$res)^3,
                          mean((d$res - mean(d$res))^4) / sd(d$res)^4 - 3), n = nrow(d))
  diag_rows[[bm]] <- bind_rows(by_arm, by_fit, by_floor, qq) %>% mutate(biomarker = bm)
  res_rows[[bm]]  <- d %>% transmute(biomarker = bm, arm, plasma_time, fit, res)
  p1 <- ggplot(d, aes(fit, res)) + geom_point(alpha = .15, size = .6) + geom_smooth(se = FALSE, colour = "red") +
    facet_wrap(~ arm) + labs(title = paste(bm, ": residual vs fitted by platform"), x = "fitted", y = "residual") + theme_minimal()
  p2 <- ggplot(d, aes(sample = res)) + stat_qq(alpha = .3, size = .6) + stat_qq_line() +
    facet_wrap(~ arm) + labs(title = paste(bm, ": normal QQ by platform")) + theme_minimal()
  print(p1); print(p2)
}
dev.off()
res_all <- bind_rows(res_rows)
out(res_all, "T2f_residuals")

# The ten-page PDF above stays as an internal check. Supplementary Figure S5 is this single
# styled page: residuals against the fitted value and against time, one column per biomarker.
{
  BM  <- c(pTau217 = "pTau217", pTau181 = "pTau181", GFAP = "GFAP", NfL = "NfL", Ab42_40 = "Ab42/40")
  ARM <- function(a) sub("^C2N_ratio$", "C2N %pTau217", as.character(a))
  rp <- res_all %>%
    mutate(bm = factor(BM[biomarker], BM), platform = factor(ARM(arm))) %>%
    tidyr::pivot_longer(c(fit, plasma_time), names_to = "x_kind", values_to = "x") %>%
    mutate(x_lab = if_else(x_kind == "fit", "fitted value", "years from the first sample"),
           # one panel per biomarker and per kind of x, because the two share no x scale
           panel = factor(paste0(bm, "\n", x_lab),
                          paste0(rep(levels(bm), each = 2),
                                 c("\nfitted value", "\nyears from the first sample"))))
  # a handful of residuals reach 8; the axis is held at the range that holds 99.8 per cent of them
  CLIP <- 2.5
  n_out <- sum(abs(res_all$res) > CLIP)
  pS3 <- ggplot(rp, aes(x, res, colour = platform)) +
    geom_hline(yintercept = 0, colour = GRID, linewidth = 0.4) +
    geom_point(alpha = 0.10, size = 0.3, stroke = 0) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.6) +
    facet_wrap(~ panel, nrow = 2, dir = "v", scales = "free_x") +
    coord_cartesian(ylim = c(-CLIP, CLIP)) +
    labs(x = NULL, y = "Residual, standardised units", colour = NULL,
         title = "Residual diagnostics of the multi-platform change model",
         subtitle = sprintf(paste("Points are single measurements; lines are loess fits within a platform.",
                                  "%d of %d residuals fall outside the axis."), n_out, nrow(res_all))) +
    theme(legend.position = "bottom", strip.text = element_text(lineheight = 1.1))
  fS3 <- file.path(P2_DIR, "figures", "supplementary", "FigureS5_residual_diagnostics")
  dir.create(dirname(fS3), showWarnings = FALSE, recursive = TRUE)
  ggsave(paste0(fS3, ".png"), pS3, width = 11, height = 6, dpi = 300, bg = "white")
  ggsave(paste0(fS3, ".pdf"), pS3, width = 11, height = 6, device = "pdf", useDingbats = FALSE, bg = "white")
  say("  Supplementary Figure S5 written to %s.png", fS3)
}

diag_tab <- bind_rows(diag_rows) %>% select(biomarker, facet, level, sd_res, n)
print(as.data.frame(diag_tab %>% filter(facet != "shape") %>% mutate(sd_res = round(sd_res, 3))))
out(diag_tab, "T2e_residual_diagnostics")

# ---- 3. Cluster bootstrap ------------------------------------------------------------------
if (N_BOOT > 0) {
  say("\n== 3. Cluster bootstrap by person, %d draws ==", N_BOOT)
  crit_cols <- c("centiloid_slope", "csf_ptau_abeta42_slope", "mmse_slope", "adas13_slope", "mpacc_digit_slope", "hcv_adj_slope")
  crit_cols <- crit_cols[crit_cols %in% names(critr)]
  # Signs oriented as in T5 (04, section 1): hippocampal volume, mPACC and MMSE fall as disease
  # progresses, so their correlations are reversed and a positive value always means the plasma
  # slope tracks worsening. Added November 2026 with mPACC, so that every criterion in Figure S6
  # has an interval and the intervals carry the same sign as the point estimates.
  FLIP <- c("hcv_adj_slope", "mpacc_digit_slope", "mmse_slope")
  f_boot <- if (USE_REDUCED) f_red else f_full
  f_level <- as.formula(paste(YVAR, "~ 1 + (1|RID_f) + (1|arm_f) + (1|RID_f:arm_f) + (1|sample_f) + (1|arm_visit) + (1|batch_f)"))
  long <- long %>% mutate(sample_f = factor(paste(RID, visit_no)), arm_visit = factor(paste(arm, visit_no)))
  rids <- unique(long$RID)

  one_draw <- function(b) {
    set.seed(BOOT_SEED + b)
    picked <- sample(rids, length(rids), replace = TRUE)
    # give each resampled copy of a person a new id so copies are separate clusters
    boot_ids <- tibble(RID = picked, boot_id = seq_along(picked))
    res <- list()
    for (bm in biomarkers) {
      d <- long %>% filter(biomarker == bm) %>%
        inner_join(boot_ids, by = "RID", relationship = "many-to-many") %>%
        mutate(RID_f = factor(boot_id), sample_f = factor(paste(boot_id, visit_no)))
      # level G
      g <- tryCatch({
        fl <- lmer(f_level, data = d, REML = TRUE, control = ctrl)
        v <- as.data.frame(VarCorr(fl)); getv <- function(g) { x <- v$vcov[v$grp == g]; if (length(x)) x[1] else 0 }
        person <- getv("RID_f")
        person / (person + getv("RID_f:arm_f") + getv("sample_f") + getv("Residual"))
      }, error = function(e) NA_real_)
      # change model, all people and by amyloid status
      ch <- function(dd, lab) tryCatch({
        fc <- lmer(f_boot, data = dd, REML = TRUE, control = ctrl)
        v <- as.data.frame(VarCorr(fc)) %>% filter(is.na(var2))
        vs <- v$vcov[v$grp == "RID_f" & v$var1 == "plasma_time"]
        pa <- sum(v$vcov[grepl("^RID_f[.:]arm_f", v$grp) & !is.na(v$var1) & v$var1 == "plasma_time"])
        re <- v$vcov[v$grp == "Residual"]; k <- n_distinct(dd$arm)
        re_bl <- ranef(fc)$RID_f
        pooled <- tibble(boot_id = as.integer(rownames(re_bl)), pooled_slope = re_bl[["plasma_time"]])
        val <- pooled %>% left_join(dd %>% distinct(boot_id, RID), by = "boot_id") %>%
          left_join(critr %>% select(RID, all_of(crit_cols)), by = "RID")
        rs <- sapply(crit_cols, function(cc) cor(val$pooled_slope, val[[cc]], use = "complete.obs", method = "spearman"))
        rs <- ifelse(crit_cols %in% FLIP, -rs, rs)
        tibble(quantity = c("rel_slope_one_platform", "rel_slope_mean_of_arms", paste0("validity_", crit_cols)),
               value = c(vs / (vs + pa + re / median_sst), vs / (vs + (pa + re / median_sst) / k), rs),
               subset = lab)
      }, error = function(e) NULL)
      # slope ICC across arms from per-person OLS slopes
      icc <- tryCatch({
        sl <- slope_by(d, YVAR, "plasma_time", boot_id, arm) %>% select(boot_id, arm, slope) %>%
          pivot_wider(names_from = arm, values_from = slope)
        mc <- sl[complete.cases(sl), -1]
        dd <- mc %>% mutate(id = row_number()) %>% pivot_longer(-id, names_to = "arm", values_to = "s")
        fi <- lmer(s ~ 1 + (1|id) + (1|arm), data = dd, control = ctrl)
        v <- as.data.frame(VarCorr(fi)); vi <- v$vcov[v$grp == "id"]; va <- v$vcov[v$grp == "arm"]; ve <- v$vcov[v$grp == "Residual"]
        c(vi / (vi + va + ve), vi / (vi + (va + ve) / ncol(mc)))
      }, error = function(e) c(NA, NA))
      res[[bm]] <- bind_rows(
        tibble(quantity = "G_relative_one_platform", value = g, subset = "all"),
        ch(d, "all"),
        ch(d %>% filter(amyloid_bl_cl20 %in% FALSE), "amyloid-"),
        ch(d %>% filter(amyloid_bl_cl20 %in% TRUE),  "amyloid+"),
        tibble(quantity = c("icc_single_arm", "icc_mean_of_arms"), value = icc, subset = "all")
      ) %>% mutate(biomarker = bm)
    }
    bind_rows(res) %>% mutate(draw = b)
  }

  draws_file <- file.path(P2_DIR, "bootstrap_draws.rds")
  draws <- if (file.exists(draws_file)) readRDS(draws_file) else tibble()
  # Saved draws from before November 2026 lack the mPACC validity and carry the unoriented signs,
  # so they are not resumed: the file is renamed and every draw is recomputed with the same seeds.
  if (nrow(draws) && !any(draws$quantity == "validity_mpacc_digit_slope")) {
    old <- sub("\\.rds$", paste0("_superseded_", format(Sys.time(), "%Y%m%d"), ".rds"), draws_file)
    file.rename(draws_file, old)
    say("  saved draws predate the November 2026 change to the validity quantities; moved to %s and recomputing", basename(old))
    draws <- tibble()
  }
  done <- if (nrow(draws)) unique(draws$draw) else integer(0)
  todo <- setdiff(seq_len(N_BOOT), done)
  say("  %d draws already saved, %d to do", length(done), length(todo))
  t0 <- Sys.time()
  for (b in todo) {
    draws <- bind_rows(draws, one_draw(b))
    saveRDS(draws, draws_file)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    say("  draw %d done (%.1f min elapsed, ~%.0f min remaining)", b, el,
        el / (which(todo == b)) * (length(todo) - which(todo == b)))
  }

  # Summarise: percentile intervals
  ci <- draws %>%
    group_by(biomarker, subset, quantity) %>%
    summarise(n_draws = sum(!is.na(value)),
              boot_mean = mean(value, na.rm = TRUE),
              ci_low = quantile(value, .025, na.rm = TRUE),
              ci_high = quantile(value, .975, na.rm = TRUE), .groups = "drop")
  out(ci, "T7_bootstrap_intervals")
  print(as.data.frame(ci %>% filter(quantity %in% c("rel_slope_one_platform", "rel_slope_mean_of_arms", "icc_single_arm")) %>%
                        mutate(across(c(boot_mean, ci_low, ci_high), ~ round(.x, 3)))))
  say("\nBootstrap intervals written to %s", file.path(P2_DIR, "T7_bootstrap_intervals.csv"))
} else {
  say("\nN_BOOT = 0: bootstrap skipped.")
}


# =============================================================================================
# Section 2. Intervals for the minimal detectable change and trial sizes (T3b, T15b)
# Original script: 08b16_paper2_mdc_bootstrap.R
# =============================================================================================
# =============================================================================
# 08b16_paper2_mdc_bootstrap.R                                        v0.1.0
#
# Confidence intervals for the three quantities the paper reports as point
# estimates: the minimal detectable change, the years an average change needs to
# exceed it, and the trial sample size. Everything else already has an interval
# from 08b2.
#
# The estimators are those of 08b (section 3) and 08b6 (sections 11 and 12), so
# the bootstrap mean should sit on the published point estimate. Within a draw:
#
#   1. participants are resampled with replacement, each copy given a new id so
#      that copies are separate clusters, using the same seeds as 08b2, so a
#      draw here is the same resample as the draw of that number there;
#   2. the multi-platform change model is refitted per biomarker (the reduced
#      model, as in the 08b2 bootstrap, since the platform-specific slope
#      variance is zero for every biomarker) for var(true slope) and the
#      residual, which the trial size needs;
#   3. a per-arm model with participant intercepts and slopes is refitted for
#      the within-person residual SD and the mean rate, which the minimal
#      detectable change needs;
#   4. per-participant least-squares slopes, the SD of the log value and the
#      median span are recomputed on the resample, so every input to the
#      delivered-scale conversion is resampled rather than held fixed.
#
# Writes  pipeline_out/paper2/T3b_mdc_bootstrap.csv   (per biomarker and platform)
#         pipeline_out/paper2/T15b_trial_bootstrap.csv (per schedule)
#         pipeline_out/paper2/mdc_bootstrap_draws.rds  (the cache)
#
# Every draw is saved as it is computed and the script resumes from the cache,
# so it can be stopped with ctrl-C and restarted. Run it from a terminal:
# P2_MDC_BOOT=0 computes the point estimates only and writes nothing.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
if (COHORT != "adni") stop("08b16 is the ADNI script; set Sys.setenv(COHORT = 'adni').")
suppressPackageStartupMessages(library(lme4))

N_BOOT    <- as.integer(Sys.getenv("P2_MDC_BOOT", unset = "1000"))
BOOT_SEED <- 20260904                      # the seed 08b2 uses, so draws match
P2_DIR    <- file.path(OUT_DIR, "paper2")
dir.create(P2_DIR, showWarnings = FALSE)
out <- function(df, name) { write_csv(df, file.path(P2_DIR, paste0(name, ".csv"))); invisible(df) }
YVAR <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
zsum2 <- (qnorm(0.975) + qnorm(0.80))^2
DESIGNS <- tibble(
  schedule = c("18 months, 4 draws (0, 6, 12, 18 m)",
               "24 months, 3 draws (0, 12, 24 m)",
               "24 months, 5 draws (0, 6, 12, 18, 24 m)"),
  span = c(1.5, 2, 2),
  sst  = c(sum((c(0, .5, 1, 1.5) - mean(c(0, .5, 1, 1.5)))^2),
           sum((c(0, 1, 2)      - mean(c(0, 1, 2)))^2),
           sum((c(0, .5, 1, 1.5, 2) - mean(c(0, .5, 1, 1.5, 2)))^2)))
EFFECT <- 0.25

long0 <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))

# ---- one dataset in, one row per quantity out ---------------------------------------
# `d` must carry: id (the clustering unit), plasma_time, z_log_aligned, log_value,
# value, direction, biomarker, arm, measure, batch, amyloid_bl_cl20.
compute_one <- function(d) {
  d <- d %>% mutate(RID_f = factor(id), arm_f = factor(arm), batch_f = factor(batch))

  # Per-participant least-squares slopes, one row per participant and measure, with
  # the same inclusion rule as slope_by() in 00_setup.R: a slope needs two or more
  # observations and a non-zero span.
  sl <- d %>%
    filter(!is.na(.data[[YVAR]]), !is.na(plasma_time)) %>%
    group_by(biomarker, arm, measure, id) %>%
    summarise(n_obs = n(),
              span  = max(plasma_time) - min(plasma_time),
              slope = if (n() >= 2 && (max(plasma_time) - min(plasma_time)) > 0)
                        cov(plasma_time, .data[[YVAR]]) / var(plasma_time) else NA_real_,
              amy   = first(amyloid_bl_cl20), .groups = "drop")

  # The median span is taken over every participant and measure with exactly three
  # observations, pooled across biomarkers, as in 08b section 3; the median SST over
  # participants with two or more draws, recomputed on the resample, is what the slope
  # MDC uses (as 04 section 3 since November 2026).
  median_span <- median(sl$span[sl$n_obs == 3], na.rm = TRUE)
  sst_med <- d %>% distinct(id, plasma_time) %>% group_by(id) %>%
    summarise(k = n(), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
    filter(k >= 2) %>% pull(sst) %>% median()

  # The mean amyloid-positive rate uses participants with three or more observations,
  # as in 08b6 section 11.
  amy_slope <- sl %>% filter(n_obs >= 3, amy %in% TRUE) %>% group_by(biomarker, arm) %>%
    summarise(mean_slope_amy_pos_z = mean(slope, na.rm = TRUE), .groups = "drop")

  # the scale factor, resampled with everything else
  scale_tab <- d %>% group_by(biomarker, arm, measure, direction) %>%
    summarise(sd_log = sd(log_value, na.rm = TRUE), .groups = "drop")

  # per-arm within-person residual SD and mean rate
  per_arm <- d %>% group_by(biomarker, arm, measure) %>%
    group_modify(~ {
      fit <- tryCatch(lmer(as.formula(paste(YVAR, "~ plasma_time + (1 + plasma_time | RID_f)")),
                           data = .x, control = ctrl), error = function(e) NULL)
      if (is.null(fit)) return(tibble(sd_residual = NA_real_, mean_slope = NA_real_))
      tibble(sd_residual = sigma(fit), mean_slope = fixef(fit)[["plasma_time"]])
    }) %>% ungroup()

  # the multi-platform change model, for the trial size
  chg <- lapply(biomarkers, function(bm) {
    dd <- d %>% filter(biomarker == bm)
    f <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                          "(1 | RID_f:arm_f) + (1|batch_f)"))
    fit <- tryCatch(lmer(f, data = dd, REML = TRUE, control = ctrl), error = function(e) NULL)
    if (is.null(fit)) return(tibble(biomarker = bm, var_person_slope = NA_real_,
                                    var_residual = NA_real_, n_arms = n_distinct(dd$arm)))
    vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
    vs <- vc$vcov[vc$grp == "RID_f" & !is.na(vc$var1) & vc$var1 == "plasma_time"]
    tibble(biomarker = bm, var_person_slope = if (length(vs)) vs[1] else NA_real_,
           var_residual = vc$vcov[vc$grp == "Residual"][1], n_arms = n_distinct(dd$arm))
  }) %>% bind_rows()

  # ---- the reported quantities, exactly as 08b6 derives them
  mdc <- per_arm %>%
    left_join(scale_tab, by = c("biomarker", "arm", "measure")) %>%
    left_join(amy_slope, by = c("biomarker", "arm")) %>%
    mutate(slope_se_3pts        = sd_residual / sqrt(sst_med),
           mdc95_two_points     = 1.96 * sqrt(2) * sd_residual,
           mdc95_slope_per_year = 1.96 * slope_se_3pts,
           mdc95_two_points_pct = 100 * (exp(mdc95_two_points * sd_log) - 1),
           mdc95_slope_pct_per_year = 100 * (exp(mdc95_slope_per_year * sd_log) - 1),
           mean_pct_per_year_amy_pos = 100 * (exp(mean_slope_amy_pos_z * sd_log * direction) - 1),
           years_to_exceed_mdc_amy_pos = mdc95_two_points / abs(mean_slope_amy_pos_z),
           cv_within_pct = 100 * sqrt(exp((sd_residual * sd_log)^2) - 1))
  rc <- sl %>% left_join(mdc %>% select(measure, mdc95_slope_per_year), by = "measure") %>%
    group_by(biomarker, arm) %>%
    summarise(pct_reliable_increase = 100 * mean(slope >  mdc95_slope_per_year, na.rm = TRUE),
              pct_reliable_decrease = 100 * mean(slope < -mdc95_slope_per_year, na.rm = TRUE),
              .groups = "drop")   # na.rm drops participants with no estimable slope
  mdc <- mdc %>% left_join(rc, by = c("biomarker", "arm"))

  # trial size, per arm and for the mean of platforms
  arm_comp <- mdc %>% select(biomarker, arm, sd_residual, mean_slope_amy_pos_z) %>%
    left_join(chg %>% select(biomarker, var_person_slope, n_arms), by = "biomarker")
  mean_comp <- arm_comp %>% group_by(biomarker) %>%
    summarise(arm = "mean of platforms",
              sd_residual = sqrt(mean(sd_residual^2, na.rm = TRUE) / first(n_arms)),
              var_person_slope = first(var_person_slope), n_arms = first(n_arms),
              mean_slope_amy_pos_z = mean(mean_slope_amy_pos_z, na.rm = TRUE), .groups = "drop")
  trial <- bind_rows(arm_comp, mean_comp) %>%
    tidyr::crossing(DESIGNS) %>%
    mutate(var_slope_obs = var_person_slope + sd_residual^2 / sst,
           delta_slope   = EFFECT * abs(mean_slope_amy_pos_z),
           n_per_arm_slope = ceiling(2 * zsum2 * var_slope_obs / delta_slope^2),
           reliability_of_trial_slope = var_person_slope / var_slope_obs)

  bind_rows(
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "mdc95_two_points_pct", value = mdc95_two_points_pct),
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "mdc95_slope_pct_per_year", value = mdc95_slope_pct_per_year),
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "years_to_exceed_mdc_amy_pos", value = years_to_exceed_mdc_amy_pos),
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "cv_within_pct", value = cv_within_pct),
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "mean_pct_per_year_amy_pos", value = mean_pct_per_year_amy_pos),
    mdc %>% transmute(biomarker, arm, schedule = NA_character_,
                      quantity = "pct_reliable_increase", value = pct_reliable_increase),
    trial %>% transmute(biomarker, arm, schedule,
                        quantity = "n_per_arm_slope", value = as.numeric(n_per_arm_slope)),
    trial %>% transmute(biomarker, arm, schedule,
                        quantity = "reliability_of_trial_slope", value = reliability_of_trial_slope))
}

# ---- the point estimate ---------------------------------------------------------------
say("\n== point estimates on the observed data ==")
point <- compute_one(long0 %>% mutate(id = RID)) %>% rename(estimate = value)
print(as.data.frame(point %>% filter(biomarker == "pTau217", is.na(schedule)) %>%
                      tidyr::pivot_wider(names_from = quantity, values_from = estimate) %>%
                      mutate(across(where(is.numeric), ~ round(.x, 2)))))

RUN_BOOT <- N_BOOT > 0
# Do the recomputed point estimates reproduce the published table? If they do not,
# the bootstrap would be an interval around the wrong centre, so check this first.
pub <- file.path(P2_DIR, "T14_mdc_fold_change.csv")
if (file.exists(pub)) {
  T14 <- read_csv(pub, show_col_types = FALSE)
  chk <- point %>% filter(quantity %in% c("mdc95_two_points_pct", "mdc95_slope_pct_per_year",
                                          "years_to_exceed_mdc_amy_pos", "cv_within_pct",
                                          "mean_pct_per_year_amy_pos", "pct_reliable_increase")) %>%
    tidyr::pivot_wider(names_from = quantity, values_from = estimate) %>%
    inner_join(T14, by = c("biomarker", "arm"), suffix = c("_new", "_published"))
  # 08b rounds the reliable-change percentages to one decimal before writing T3, so the
  # published value carries that rounding and can differ from the recomputation by up
  # to half of the last digit. Every other quantity is stored at full precision.
  ROUNDED <- c("pct_reliable_increase", "pct_reliable_decrease")
  cmp <- lapply(c("mdc95_two_points_pct", "mdc95_slope_pct_per_year", "years_to_exceed_mdc_amy_pos",
                  "cv_within_pct", "mean_pct_per_year_amy_pos", "pct_reliable_increase"), function(q) {
    a <- chk[[paste0(q, "_new")]]; b <- chk[[paste0(q, "_published")]]
    tol <- if (q %in% ROUNDED) 0.05 else 1e-6
    d <- max(abs(a - b), na.rm = TRUE)
    tibble(quantity = q, max_abs_difference = d, tolerance = tol, passes = d <= tol)
  }) %>% bind_rows()
  say("\n== recomputed point estimates against T14_mdc_fold_change.csv ==")
  print(as.data.frame(cmp %>% mutate(max_abs_difference = signif(max_abs_difference, 3))))
  if (!all(cmp$passes))
    say("  WARNING: the recomputation does not reproduce the published values; do not run the bootstrap until this is understood.")
  else say("  The recomputation reproduces the published values, the reliable-change percentages to the one decimal at which they are stored.")
} else say("\nT14_mdc_fold_change.csv not found, so the point estimates could not be checked.")

if (!RUN_BOOT) say("\nP2_MDC_BOOT = 0: point estimates only, nothing written.")

if (RUN_BOOT) {

  # ---- the bootstrap --------------------------------------------------------------------
  say("\n== person bootstrap, %d draws ==", N_BOOT)
  rids <- unique(long0$RID)
  one_draw <- function(b) {
    set.seed(BOOT_SEED + b)
    picked <- sample(rids, length(rids), replace = TRUE)
    ids <- tibble(RID = picked, id = seq_along(picked))
    d <- long0 %>% inner_join(ids, by = "RID", relationship = "many-to-many")
    compute_one(d) %>% mutate(draw = b)
  }
  cache <- file.path(P2_DIR, "mdc_bootstrap_draws.rds")
  draws <- if (file.exists(cache)) readRDS(cache) else tibble()
  done  <- if (nrow(draws)) unique(draws$draw) else integer(0)
  todo  <- setdiff(seq_len(N_BOOT), done)
  say("  %d draws already saved, %d to do", length(done), length(todo))
  t0 <- Sys.time()
  for (i in seq_along(todo)) {
    b <- todo[i]
    draws <- bind_rows(draws, suppressMessages(one_draw(b)))
    saveRDS(draws, cache)
    el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    if (i %% 10 == 0 || i == 1)
      say("  draw %d done (%.1f min elapsed, ~%.0f min remaining)", b, el,
          el / i * (length(todo) - i))
  }

  # ---- percentile intervals -------------------------------------------------------------
  # One row per draw and quantity, in case a draw was written more than once (two
  # copies of the script running at the same time would do that; the draws themselves
  # are identical, since each seeds itself from its own draw number).
  ci <- draws %>% filter(draw <= N_BOOT) %>%
    distinct(draw, biomarker, arm, schedule, quantity, .keep_all = TRUE) %>%
    group_by(biomarker, arm, schedule, quantity) %>%
    summarise(n_draws = sum(is.finite(value)),
              boot_mean = mean(value[is.finite(value)], na.rm = TRUE),
              ci_low  = quantile(value[is.finite(value)], .025, na.rm = TRUE),
              ci_high = quantile(value[is.finite(value)], .975, na.rm = TRUE), .groups = "drop") %>%
    left_join(point, by = c("biomarker", "arm", "schedule", "quantity")) %>%
    select(biomarker, arm, schedule, quantity, estimate, n_draws, boot_mean, ci_low, ci_high)

  out(ci %>% filter(is.na(schedule)) %>% select(-schedule), "T3b_mdc_bootstrap")
  out(ci %>% filter(!is.na(schedule)), "T15b_trial_bootstrap")
  say("\nWritten: T3b_mdc_bootstrap.csv and T15b_trial_bootstrap.csv in %s", P2_DIR)
  print(as.data.frame(ci %>% filter(biomarker == "pTau217", is.na(schedule),
                                    quantity %in% c("mdc95_two_points_pct", "years_to_exceed_mdc_amy_pos")) %>%
                        mutate(across(where(is.numeric), ~ round(.x, 2)))))

}   # end of the bootstrap block


# =============================================================================================
# Section 3. Intervals for the bounds, fitted-alone reliabilities, amyloid-positive trial sizes, immunoassay ICC (R10)
# Original script: 08b19_paper2_revision_bootstrap.R
# =============================================================================================
# =============================================================================
# 08b19_paper2_revision_bootstrap.R
# Paper 2, November 2026 revision: percentile bootstrap intervals for the quantities that
# 08b17 reported as point estimates, so that every reliability and trial size in the paper
# has a 95% interval.
#
# Quantities, one set per bootstrap draw (participants resampled with replacement, the same
# seeds as 08b2 so draw b here uses the same participants as draw b of the main bootstrap):
#   B1  pTau217 lower bound of averaging (model with a person-by-sample effect, whole cohort):
#       reliability of one platform and of the mean of platforms with the shared part kept,
#       and the sampling span needed for the mean of platforms to reach 0.7, 0.8 and 0.9
#       (Table S16b, lower-bound column). pTau217 only: for the other biomarkers the sample
#       model is not identifiable (08b17 R1).
#   B2  Reliability of each of the 15 measures fitted alone (08b17 R4), at the ADNI median
#       SST and at the A4 placebo SST.
#   B3  Trial sizes per arm on amyloid-positive variance components (08b17 R3): every
#       platform, and the mean of platforms under the upper and (pTau217) lower bound, for
#       the three schedules and both effect sizes.
#   B4  pTau217 slope ICC, immunoassays only (08b17 R6).
# Draw 0 is the full data without resampling; its values must reproduce R1, R3a, R4 and R6
# (except that the amyloid-positive lower bound is reported here even where R3a left it blank).
# Draws in which the sample model is singular are kept (a boundary estimate is a legitimate
# draw); the proportion singular is reported.
#
# Requires 08b17 to have been run (its R1, R3a, R4 and R6 CSVs are read for the draw-0 check).
# Run in R (the console stays busy; expect one to two minutes per draw, so 1,000 draws take
# more than one night; stop and resume as often as needed, or set N_BOOT to 500 first, since the
# first 500 seeds are the same and a later run to 1,000 only adds draws):
# The first three draws are timed and the projected total is printed; stop with Esc whenever
# needed and rerun later: every draw is saved to R10_bootstrap_draws.rds and the script
# resumes from there. To change the number of draws: Sys.setenv(N_BOOT = "500") first.
# Writes: pipeline_out/paper2/R10_bootstrap_draws.rds,
#         R10_bootstrap_intervals.csv (estimate, lo95, hi95 per quantity).
# Changes no existing file.
# =============================================================================

Sys.setenv(COHORT = "adni")
source("00_setup.R")
if (COHORT != "adni") stop("COHORT must be 'adni' for this script.")
suppressPackageStartupMessages(library(lme4))

P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)

N_BOOT    <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
BOOT_SEED <- 20260904                       # as 08b2, so draws pair with the main bootstrap
YVAR       <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
zsum2      <- (qnorm(0.975) + qnorm(0.80))^2
A4_SST     <- as.numeric(Sys.getenv("P2_ANCHOR_SST", "13.533"))
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
out  <- function(df, name) { write_csv(df, file.path(REV_DIR, paste0(name, ".csv"))); invisible(df) }

say("08b19 started %s. N_BOOT = %d. Outputs: %s", format(Sys.time()), N_BOOT, REV_DIR)

long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         batch_f = factor(batch))
rids <- unique(long$RID)

# Median SST as in 08b, 08b17: over people with at least two draws, full data (held fixed)
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
  filter(k >= 2) %>% pull(sst) %>% median()
say("Median SST %.3f (paper: 12.855)", median_sst)

designs <- tribble(~schedule, ~times,
                   "18 months, 4 draws (0, 6, 12, 18 m)",     c(0, 0.5, 1, 1.5),
                   "24 months, 3 draws (0, 12, 24 m)",        c(0, 1, 2),
                   "24 months, 5 draws (0, 6, 12, 18, 24 m)", c(0, 0.5, 1, 1.5, 2)) %>%
  mutate(k = map_int(times, length), sst = map_dbl(times, ~ sum((.x - mean(.x))^2))) %>% select(-times)
effects <- c(0.25, 0.50)

f_reduced <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                              "(1 | RID_f:arm_f) + (1 | batch_f)"))
f_sample  <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                              "(1 | RID_f:arm_f) + (1 | sample_f) + (1 | batch_f)"))
f_single  <- as.formula(paste(YVAR, "~ plasma_time + (1 + plasma_time | RID_f)"))

n_slope <- function(v, delta) ceiling(2 * zsum2 * v / delta^2)
vget <- function(fit, grp, var1 = NA) {
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  i <- if (is.na(var1)) which(vc$grp == grp) else which(vc$grp == grp & !is.na(vc$var1) & vc$var1 == var1)
  if (length(i)) vc$vcov[i[1]] else NA_real_
}
pick1 <- function(x) if (length(x)) x[1] else NA_real_
safe_lmer <- function(f, d) tryCatch(suppressMessages(lmer(f, data = d, REML = TRUE, control = ctrl)),
                                     error = function(e) NULL)
# min_n = 3 for the amyloid-positive mean rate (08b6 definition); min_n = 2 for the slope ICC,
# matching frame_paper2_slopes_wide (slope_by in 00_setup.R).
ols_slopes <- function(d, min_n) d %>% group_by(boot_id, biomarker, arm) %>%
  filter(n() >= min_n, sd(plasma_time) > 0) %>%
  summarise(slope = cov(plasma_time, z_log_aligned) / var(plasma_time),
            amy = first(amyloid_bl_cl20), .groups = "drop")

# ---- one draw -----------------------------------------------------------------------------------
one_draw <- function(b) {
  if (b == 0) {
    boot_ids <- tibble(RID = rids, boot_id = seq_along(rids))
  } else {
    set.seed(BOOT_SEED + b)
    picked <- sample(rids, length(rids), replace = TRUE)
    boot_ids <- tibble(RID = picked, boot_id = seq_along(picked))
  }
  d_all <- long %>% inner_join(boot_ids, by = "RID", relationship = "many-to-many") %>%
    mutate(RID_f = factor(boot_id), sample_f = factor(paste(boot_id, visit_no)))
  res <- list()
  row <- function(quantity, value, biomarker = NA, arm = NA, schedule = NA, effect = NA)
    tibble(quantity = quantity, biomarker = biomarker, arm = arm, schedule = schedule,
           effect = effect, value = value)

  # B1. pTau217 sample model, whole cohort
  d <- d_all %>% filter(biomarker == "pTau217") %>% mutate(arm_f = droplevels(factor(arm)))
  k <- n_distinct(d$arm)
  fs <- safe_lmer(f_sample, d)
  if (!is.null(fs)) {
    vt <- vget(fs, "RID_f", "plasma_time"); vsh <- vget(fs, "sample_f"); va <- vget(fs, "Residual")
    res$b1 <- bind_rows(
      row("B1_sample_model_singular", as.numeric(isSingular(fs)), "pTau217"),
      row("B1_pct_shared_of_occasion", 100 * vsh / (vsh + va), "pTau217"),
      row("B1_rel_one_platform_lower", vt / (vt + (vsh + va) / median_sst), "pTau217"),
      row("B1_rel_mean_of_platforms_lower", vt / (vt + (vsh + va / k) / median_sst), "pTau217"),
      crossing(target = c(0.7, 0.8, 0.9), kk = c(2, 3, 4, 6)) %>%
        mutate(ratio = (vsh + va / k) / vt, sst_needed = ratio * target / (1 - target),
               value = sqrt(sst_needed * 12 * (kk - 1) / (kk * (kk + 1))),
               quantity = sprintf("B1_span_needed_lower_target%.1f_k%d", target, kk),
               biomarker = "pTau217") %>% select(quantity, biomarker, value))
  }

  # B2. each measure fitted alone, whole cohort
  res$b2 <- d_all %>% group_by(biomarker, arm) %>% group_modify(~ {
    f <- safe_lmer(f_single, .x)
    if (is.null(f)) return(tibble(quantity = character(), value = numeric()))
    vs <- vget(f, "RID_f", "plasma_time"); sr <- sigma(f)
    tibble(quantity = c("B2_rel_single_at_adni_sst", "B2_rel_single_at_a4_sst", "B2_single_singular"),
           value = c(vs / (vs + sr^2 / median_sst), vs / (vs + sr^2 / A4_SST), as.numeric(isSingular(f))))
  }) %>% ungroup()

  # B3. trial sizes on amyloid-positive components
  amy <- d_all %>% filter(amyloid_bl_cl20 %in% TRUE)
  sl  <- ols_slopes(d_all, 3) %>% filter(amy %in% TRUE)
  m_arm <- sl %>% group_by(biomarker, arm) %>% summarise(m = mean(slope, na.rm = TRUE), .groups = "drop")
  m_bm  <- m_arm %>% group_by(biomarker) %>% summarise(m = mean(m, na.rm = TRUE), .groups = "drop")
  b3 <- list()
  for (bm in biomarkers) {
    d <- amy %>% filter(biomarker == bm) %>% mutate(arm_f = droplevels(factor(arm)))
    k <- n_distinct(d$arm)
    fr <- safe_lmer(f_reduced, d)
    if (is.null(fr)) next
    vt <- vget(fr, "RID_f", "plasma_time")
    per <- d %>% group_by(arm) %>% group_modify(~ {
      f <- safe_lmer(f_single, .x)
      tibble(var_err = if (is.null(f)) NA_real_ else sigma(f)^2)
    }) %>% ungroup() %>% left_join(m_arm %>% filter(biomarker == bm) %>% select(arm, m), by = "arm") %>%
      mutate(var_true = vt, bound = "per platform")
    comps <- bind_rows(per,
                       tibble(arm = "mean of platforms", var_err = mean(per$var_err) / k,
                              m = pick1(m_bm$m[m_bm$biomarker == bm]), var_true = vt, bound = "upper"))
    if (bm == "pTau217") {
      fsa <- safe_lmer(f_sample, d)
      if (!is.null(fsa)) {
        comps <- bind_rows(comps, tibble(arm = "mean of platforms",
                                         var_err = vget(fsa, "sample_f") + vget(fsa, "Residual") / k,
                                         m = pick1(m_bm$m[m_bm$biomarker == bm]),
                                         var_true = vget(fsa, "RID_f", "plasma_time"), bound = "lower"))
        b3[[paste(bm, "sing")]] <- row("B3_sample_model_amy_pos_singular", as.numeric(isSingular(fsa)), bm)
      }
    }
    b3[[bm]] <- comps %>% crossing(designs) %>% crossing(effect = effects) %>%
      mutate(v = var_true + var_err / sst,
             rel = var_true / v,
             n = n_slope(v, effect * abs(m)),
             arm = if_else(bound == "per platform", arm, paste(arm, bound)), biomarker = bm) %>%
      select(biomarker, arm, schedule, effect, rel, n) %>%
      pivot_longer(c(rel, n), names_to = "quantity", values_to = "value") %>%
      mutate(quantity = if_else(quantity == "n", "B3_n_per_arm_amy_pos", "B3_reliability_trial_slope"))
  }
  res$b3 <- bind_rows(b3)

  # B4. pTau217 slope ICC, immunoassays only (and all five as a check)
  sw <- ols_slopes(d_all, 2) %>% filter(biomarker == "pTau217") %>% select(boot_id, arm, slope) %>%
    pivot_wider(names_from = arm, values_from = slope)
  icc_of <- function(cols, lab) {
    if (!all(cols %in% names(sw))) return(tibble())
    mc <- sw %>% select(all_of(cols)); mc <- mc[complete.cases(mc), ]
    dd <- mc %>% mutate(id = row_number()) %>% pivot_longer(-id, names_to = "arm", values_to = "s")
    f <- tryCatch(lmer(s ~ 1 + (1 | id) + (1 | arm), data = dd, control = ctrl), error = function(e) NULL)
    if (is.null(f)) return(tibble())
    v <- as.data.frame(VarCorr(f)); vi <- v$vcov[v$grp == "id"]; va <- v$vcov[v$grp == "arm"]
    ve <- v$vcov[v$grp == "Residual"]
    row(c("B4_icc_single", "B4_icc_mean"), c(vi / (vi + va + ve), vi / (vi + (va + ve) / length(cols))),
        "pTau217", lab)
  }
  res$b4 <- bind_rows(icc_of(c("ALZpath", "Fujirebio", "Janssen"), "immunoassays"),
                      icc_of(c("ALZpath", "C2N", "C2N_ratio", "Fujirebio", "Janssen"), "all five"))

  bind_rows(res) %>% mutate(draw = b)
}

# ---- draw 0: point estimates, checked against 08b17 ----------------------------------------------
say("\n== Draw 0 (full data) ==")
t0 <- Sys.time()
est <- one_draw(0)
say("Draw 0 took %.1f s", as.numeric(difftime(Sys.time(), t0, units = "secs")))
rd_rev <- function(name) read_csv(file.path(REV_DIR, paste0(name, ".csv")), show_col_types = FALSE)
chk <- list(
  est %>% filter(quantity == "B1_rel_mean_of_platforms_lower") %>% pull(value),
  rd_rev("R1_averaging_bounds") %>% filter(biomarker == "pTau217") %>% pull(rel_mean_lower),
  est %>% filter(quantity == "B2_rel_single_at_adni_sst", biomarker == "pTau217", arm == "Fujirebio") %>% pull(value),
  rd_rev("R4_single_platform_reliability_adni") %>% filter(biomarker == "pTau217", arm == "Fujirebio") %>% pull(rel_at_adni_median_sst),
  est %>% filter(quantity == "B3_n_per_arm_amy_pos", biomarker == "pTau217", arm == "mean of platforms upper",
                 schedule == "24 months, 3 draws (0, 12, 24 m)", effect == 0.25) %>% pull(value),
  rd_rev("R3a_trial_amyloid_positive_variance") %>%
    filter(biomarker == "pTau217", arm == "mean of platforms", bound == "upper (paper model)",
           schedule == "24 months, 3 draws (0, 12, 24 m)", effect == 0.25) %>% pull(n_per_arm_slope_amy_pos_variance),
  est %>% filter(quantity == "B4_icc_single", arm == "immunoassays") %>% pull(value),
  rd_rev("R6_ptau217_slope_icc_immunoassays") %>% slice(2) %>% pull(icc_single))
chk <- lapply(chk, pick1)
say("Checks against 08b17 (this script / 08b17): lower bound %.4f / %.4f; Fujirebio single %.4f / %.4f;",
    chk[[1]], chk[[2]], chk[[3]], chk[[4]])
say("  mean-of-platforms upper n (24 m, 3 draws, 25%%) %d / %d; immunoassay ICC %.4f / %.4f",
    chk[[5]], chk[[6]], chk[[7]], chk[[8]])
differs <- function(a, b, tol = 1e-3) !(length(a) == 1 && length(b) == 1 && isTRUE(abs(a - b) <= tol))
if (differs(chk[[1]], chk[[2]]) || differs(chk[[3]], chk[[4]]) || differs(chk[[5]], chk[[6]], 1) ||
    differs(chk[[7]], chk[[8]]))
  warning("Draw 0 does not reproduce 08b17; check before using the intervals.")

# ---- bootstrap draws, resumable -------------------------------------------------------------------
draws_file <- file.path(REV_DIR, "R10_bootstrap_draws.rds")
draws <- if (file.exists(draws_file)) readRDS(draws_file) else tibble()
done  <- if (nrow(draws)) unique(draws$draw) else integer(0)
todo  <- setdiff(seq_len(N_BOOT), done)
say("\n== Bootstrap: %d draws already saved, %d to do ==", length(done), length(todo))
t0 <- Sys.time()
for (i in seq_along(todo)) {
  b <- todo[i]
  draws <- bind_rows(draws, one_draw(b))
  tmp <- paste0(draws_file, ".tmp"); saveRDS(draws, tmp); file.rename(tmp, draws_file)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (i <= 3 || i %% 25 == 0 || i == length(todo))
    say("  draw %d done (%.1f min elapsed, %.0f min projected remaining)", b, el, el / i * (length(todo) - i))
  if (i == 3) say("  Projected total for %d draws: %.0f min. Press Esc to stop; rerun to resume.",
                  length(todo), el / 3 * length(todo))
}

# ---- summary: estimate from draw 0 with percentile interval ------------------------------------------
ci <- draws %>% filter(draw > 0) %>%
  group_by(quantity, biomarker, arm, schedule, effect) %>%
  summarise(n_draws = sum(!is.na(value)), lo95 = quantile(value, 0.025, na.rm = TRUE),
            hi95 = quantile(value, 0.975, na.rm = TRUE), boot_median = median(value, na.rm = TRUE),
            .groups = "drop")
R10 <- est %>% select(quantity, biomarker, arm, schedule, effect, estimate = value) %>%
  left_join(ci, by = c("quantity", "biomarker", "arm", "schedule", "effect")) %>%
  arrange(quantity, factor(biomarker, biomarkers), arm, schedule, effect)
out(R10, "R10_bootstrap_intervals")
say("\nSingular fits among draws: sample model whole cohort %.0f%%, amyloid-positive %.0f%%",
    100 * mean(draws$value[draws$quantity == "B1_sample_model_singular" & draws$draw > 0]),
    100 * mean(draws$value[draws$quantity == "B3_sample_model_amy_pos_singular" & draws$draw > 0]))
show <- R10 %>% filter(quantity %in% c("B1_rel_one_platform_lower", "B1_rel_mean_of_platforms_lower",
                                       "B2_rel_single_at_adni_sst", "B4_icc_single", "B4_icc_mean") |
                         (quantity == "B3_n_per_arm_amy_pos" & schedule == "24 months, 3 draws (0, 12, 24 m)" &
                            effect == 0.25)) %>%
  mutate(across(c(estimate, lo95, hi95), ~ signif(.x, 3)))
print(as.data.frame(show %>% select(quantity, biomarker, arm, estimate, lo95, hi95, n_draws)))
say("\n08b19 finished %s. Intervals in %s", format(Sys.time()), file.path(REV_DIR, "R10_bootstrap_intervals.csv"))


# =============================================================================================
# Section 4. Intervals for the span needed (R11)
# Original script: 08b20_paper2_span_intervals.R
# =============================================================================================
# =============================================================================
# 08b20_paper2_span_intervals.R
# Paper 2, November 2026 revision: bootstrap intervals for the sampling span needed to reach
# a target reliability (Table S16b, one-platform and upper-bound columns), derived from the
# saved draws of the main bootstrap (08b2, pipeline_out/paper2/bootstrap_draws.rds).
#
# No refitting. The span formula of 08b6 depends on the variance components only through the
# ratio var(residual) / var(true slope), and that ratio is a function of the reliability at the
# observed median SST that 08b2 saved for every draw:
#   one platform:      ratio = (1 / rel_one - 1) * median_sst
#   mean of platforms: ratio / n_arms = (1 / rel_mean - 1) * median_sst
#   sst_needed = ratio * target / (1 - target);  span = sqrt(sst_needed * 12 (k - 1) / (k (k + 1)))
# The point estimates recomputed from T2 must equal T16b exactly, which the script checks.
# The lower-bound spans already have intervals from 08b19 (R10).
#
# Run in R (seconds):
# Writes pipeline_out/paper2/R11_span_intervals.csv. Changes no
# existing file.
# =============================================================================

Sys.setenv(COHORT = "adni")
source("00_setup.R")
if (COHORT != "adni") stop("COHORT must be 'adni' for this script.")

P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
dir.create(REV_DIR, showWarnings = FALSE, recursive = TRUE)
say("08b20 started %s", format(Sys.time()))

biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
# median SST exactly as 08b2 computes it for the saved draws
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>% pull(sst) %>% median()
say("Median SST %.3f", median_sst)

T2   <- read_csv(file.path(P2_DIR, "T2_change_variance_components.csv"), show_col_types = FALSE)
T16b <- read_csv(file.path(P2_DIR, "T16b_span_needed.csv"), show_col_types = FALSE) %>%
  mutate(target = round(target, 2), k = as.numeric(k))
draws <- readRDS(file.path(P2_DIR, "bootstrap_draws.rds")) %>%
  filter(subset == "all", quantity %in% c("rel_slope_one_platform", "rel_slope_mean_of_arms"))
say("Draws available: %d", n_distinct(draws$draw))

span_from_rel <- function(rel, target, k) {
  ratio <- (1 / rel - 1) * median_sst
  sqrt(ratio * target / (1 - target) * 12 * (k - 1) / (k * (k + 1)))
}
grid <- crossing(target = c(0.7, 0.8, 0.9), k = c(2, 3, 4, 6))

# Point estimates from T2, checked against T16b
pt <- T2 %>% select(biomarker, rel_one = rel_slope_one_platform, rel_mean = rel_slope_mean_of_arms) %>%
  crossing(grid) %>%
  mutate(span_one = span_from_rel(rel_one, target, k), span_mean = span_from_rel(rel_mean, target, k)) %>%
  left_join(T16b, by = c("biomarker", "target", "k"))
if (anyNA(pt$span_needed_one_platform)) stop("T16b join failed for some biomarker/target/k rows.")
mx <- max(abs(pt$span_one - pt$span_needed_one_platform), abs(pt$span_mean - pt$span_needed_mean_of_platforms), na.rm = TRUE)
say("Largest difference between recomputed and T16b spans: %.2e years (must be ~0)", mx)
if (mx > 1e-6) warning("Recomputed spans do not reproduce T16b; check median_sst and T2 columns.")

# Intervals from the draws
bd <- draws %>% select(draw, biomarker, quantity, value) %>%
  pivot_wider(names_from = quantity, values_from = value) %>%
  crossing(grid) %>%
  mutate(span_one = span_from_rel(rel_slope_one_platform, target, k),
         span_mean = span_from_rel(rel_slope_mean_of_arms, target, k)) %>%
  group_by(biomarker, target, k) %>%
  summarise(n_draws = sum(!is.na(span_one)),
            span_one_lo95 = quantile(span_one, 0.025, na.rm = TRUE), span_one_hi95 = quantile(span_one, 0.975, na.rm = TRUE),
            span_mean_lo95 = quantile(span_mean, 0.025, na.rm = TRUE), span_mean_hi95 = quantile(span_mean, 0.975, na.rm = TRUE),
            .groups = "drop")
R11 <- pt %>% select(biomarker, target, k, span_one_platform = span_needed_one_platform,
                     span_mean_of_platforms_upper = span_needed_mean_of_platforms) %>%
  left_join(bd, by = c("biomarker", "target", "k")) %>%
  arrange(factor(biomarker, biomarkers), target, k)
print(as.data.frame(R11 %>% filter(k == 3) %>% mutate(across(where(is.numeric), ~ round(.x, 1)))))
write_csv(R11, file.path(REV_DIR, "R11_span_intervals.csv"))
say("\n08b20 finished %s. Written: %s", format(Sys.time()), file.path(REV_DIR, "R11_span_intervals.csv"))


# =============================================================================================
# Section 5. Paired differences within resamples (R12, R12b)
# Original script: 08b21_paper2_paired_differences.R
# =============================================================================================
# =============================================================================
# 08b21_paper2_paired_differences.R
# Paper 2, November 2026 revision: paired bootstrap differences for the head-to-head
# comparisons the text makes, and the boundary components of the singular amyloid-positive fits.
#
# Every quantity below was estimated in the same 1,000 participant resamples (08b2 for the
# joint model, 08b19 for the new quantities; both use set.seed(20260904 + b) and the same
# participant list), so a difference or ratio computed within each draw has a valid
# percentile interval, and it is the right test of "fewer than" or "four times as many",
# rather than comparing two overlapping intervals by eye.
#
#   D1  Trial sizes per arm, 24 months, 3 draws, 25% slowing, amyloid-positive components:
#       mean of platforms (upper bound) minus Fujirebio; mean (lower bound) minus Fujirebio;
#       Roche minus Quanterix pTau181; ratio C2N %pTau217 / Fujirebio; ratio Roche / Quanterix
#       pTau181; and the same for GFAP and NfL platform pairs.
#   D2  Reliability: each pTau217 measure fitted alone minus the joint-model one-platform value;
#       the one-platform lower bound (sample model) minus the joint-model value; the pooled lower
#       bound minus the pooled upper bound; GFAP and NfL single fits minus their joint values.
#   D3  pTau217 slope ICC, immunoassays only minus all five measures (single and mean).
#   D4  Diagnostics: variance components of the amyloid-positive reduced change model for GFAP,
#       NfL and Ab42/40 (singular in 08b17 R3b), to state which component sat at the boundary.
#
# Run in R (D4 refits three small models; under a minute):
# Reads pipeline_out/paper2/bootstrap_draws.rds and T2; pipeline_out/paper2/R10_bootstrap_draws.rds.
# Writes pipeline_out/paper2/R12_paired_differences.csv, R12b_amy_pos_singular_components.csv, log.
# =============================================================================

Sys.setenv(COHORT = "adni")
source("00_setup.R")
if (COHORT != "adni") stop("COHORT must be 'adni' for this script.")
suppressPackageStartupMessages(library(lme4))

P2_DIR  <- file.path(OUT_DIR, "paper2")
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
say("08b21 started %s", format(Sys.time()))
YVAR <- "z_log_aligned"
S24  <- "24 months, 3 draws (0, 12, 24 m)"

# ---- draws ---------------------------------------------------------------------------------------
j <- readRDS(file.path(P2_DIR, "bootstrap_draws.rds")) %>% filter(subset == "all", draw > 0) %>%
  select(draw, biomarker, quantity, value)
r <- readRDS(file.path(REV_DIR, "R10_bootstrap_draws.rds")) %>% filter(draw > 0)
# draw 0 (full data) was not saved in the draws file; take it from the intervals CSV
r0 <- read_csv(file.path(REV_DIR, "R10_bootstrap_intervals.csv"), show_col_types = FALSE) %>%
  transmute(quantity, biomarker, arm, schedule, effect, value = estimate, draw = 0L)
r <- bind_rows(r0 %>% mutate(across(c(biomarker, arm, schedule), as.character)),
               r %>% mutate(across(c(biomarker, arm, schedule), as.character)))
say("Draws: joint model %d, revision %d (draw 0 of the revision set is the full data)",
    n_distinct(j$draw), n_distinct(r$draw[r$draw > 0]))
T2 <- read_csv(file.path(P2_DIR, "T2_change_variance_components.csv"), show_col_types = FALSE)

# Point estimates: full data (draw 0 of R10; T2 for the joint model)
pick <- function(df, q, bm = NA, arm = NA, sched = NA, eff = NA) {
  d <- df %>% filter(quantity == q)
  if (!is.na(bm))    { b_ <- bm; d <- d %>% filter(.data$biomarker == b_) }
  if (!is.na(arm))   { a_ <- arm; d <- d %>% filter(.data$arm == a_) }
  if (!is.na(sched)) { s_ <- sched; d <- d %>% filter(.data$schedule == s_) }
  if (!is.na(eff))   { e_ <- eff; d <- d %>% filter(abs(.data$effect - e_) < 1e-9) }
  d %>% select(draw, value)
}
pair <- function(label, a, b, op = c("minus", "ratio")) {
  op <- match.arg(op)
  x <- inner_join(a, b, by = "draw", suffix = c("_a", "_b")) %>%
    mutate(d = if (op == "minus") value_a - value_b else value_a / value_b)
  est <- x$d[x$draw == 0]; bd <- x$d[x$draw > 0]
  tibble(comparison = label, operation = op, estimate = if (length(est)) est else NA_real_,
         n_draws = sum(!is.na(bd)), lo95 = quantile(bd, 0.025, na.rm = TRUE), hi95 = quantile(bd, 0.975, na.rm = TRUE),
         share_of_draws_a_below_b = mean(x$value_a[x$draw > 0] < x$value_b[x$draw > 0], na.rm = TRUE))
}
# joint-model draws need a draw-0 row from T2
joint <- function(bm, q) bind_rows(
  tibble(draw = 0L, value = { v <- T2[[q]][T2$biomarker == bm]; if (length(v)) v[1] else NA_real_ }),
  j %>% filter(biomarker == bm, quantity == q) %>% select(draw, value))
n_of <- function(bm, arm) pick(r, "B3_n_per_arm_amy_pos", bm, arm, S24, 0.25)
rel_single <- function(bm, arm) pick(r, "B2_rel_single_at_adni_sst", bm, arm)

say("\n== D1. Trial sizes per arm (24 months, 3 draws, 25%% slowing), paired within draws ==")
D1 <- bind_rows(
  pair("pTau217: mean of platforms (upper bound) minus Fujirebio", n_of("pTau217", "mean of platforms upper"), n_of("pTau217", "Fujirebio")),
  pair("pTau217: mean of platforms (lower bound) minus Fujirebio", n_of("pTau217", "mean of platforms lower"), n_of("pTau217", "Fujirebio")),
  pair("pTau217: mean of platforms (upper bound) / Fujirebio", n_of("pTau217", "mean of platforms upper"), n_of("pTau217", "Fujirebio"), "ratio"),
  pair("pTau217: C2N %pTau217 / Fujirebio", n_of("pTau217", "C2N_ratio"), n_of("pTau217", "Fujirebio"), "ratio"),
  pair("pTau217: ALZpath minus Fujirebio", n_of("pTau217", "ALZpath"), n_of("pTau217", "Fujirebio")),
  pair("pTau181: Roche minus Quanterix", n_of("pTau181", "Roche"), n_of("pTau181", "Quanterix")),
  pair("pTau181: Roche / Quanterix", n_of("pTau181", "Roche"), n_of("pTau181", "Quanterix"), "ratio"),
  pair("GFAP: Roche minus Quanterix", n_of("GFAP", "Roche"), n_of("GFAP", "Quanterix")),
  pair("NfL: Roche minus Quanterix", n_of("NfL", "Roche"), n_of("NfL", "Quanterix")))
print(as.data.frame(D1 %>% mutate(across(c(estimate, lo95, hi95), ~ signif(.x, 4)))))

say("\n== D2. Reliability, paired within draws ==")
D2 <- bind_rows(
  bind_rows(lapply(c("ALZpath", "C2N", "C2N_ratio", "Fujirebio", "Janssen"), function(a)
    pair(paste0("pTau217 ", a, " fitted alone minus joint-model one platform"), rel_single("pTau217", a), joint("pTau217", "rel_slope_one_platform")))),
  pair("pTau217 one-platform lower bound (sample model) minus joint-model one platform", pick(r, "B1_rel_one_platform_lower", "pTau217"), joint("pTau217", "rel_slope_one_platform")),
  pair("pTau217 pooled lower bound minus pooled upper bound", pick(r, "B1_rel_mean_of_platforms_lower", "pTau217"), joint("pTau217", "rel_slope_mean_of_arms")),
  pair("pTau217 pooled lower bound minus joint-model one platform", pick(r, "B1_rel_mean_of_platforms_lower", "pTau217"), joint("pTau217", "rel_slope_one_platform")),
  bind_rows(lapply(list(c("GFAP", "Roche"), c("GFAP", "Quanterix"), c("NfL", "Roche"), c("NfL", "Quanterix"), c("pTau181", "Roche"), c("pTau181", "Quanterix")),
                   function(v) pair(paste0(v[1], " ", v[2], " fitted alone minus joint-model one platform"), rel_single(v[1], v[2]), joint(v[1], "rel_slope_one_platform")))))
print(as.data.frame(D2 %>% mutate(across(c(estimate, lo95, hi95), ~ round(.x, 3)))))

say("\n== D3. pTau217 slope ICC, immunoassays only minus all five ==")
D3 <- bind_rows(
  pair("ICC single: immunoassays minus all five", pick(r, "B4_icc_single", "pTau217", "immunoassays"), pick(r, "B4_icc_single", "pTau217", "all five")),
  pair("ICC mean: immunoassays minus all five", pick(r, "B4_icc_mean", "pTau217", "immunoassays"), pick(r, "B4_icc_mean", "pTau217", "all five")))
print(as.data.frame(D3 %>% mutate(across(c(estimate, lo95, hi95), ~ round(.x, 3)))))
R12 <- bind_rows(D1 %>% mutate(block = "D1 trial size"), D2 %>% mutate(block = "D2 reliability"), D3 %>% mutate(block = "D3 ICC"))
write_csv(R12, file.path(REV_DIR, "R12_paired_differences.csv"))
say("Read: 'share_of_draws_a_below_b' is the proportion of resamples in which the first quantity was smaller.")

say("\n== D4. Boundary components of the singular amyloid-positive reduced fits ==")
long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         RID_f = factor(RID), batch_f = factor(batch))  # unused levels dropped by lme4 in the subset fits
f_reduced <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) + (1 | RID_f:arm_f) + (1 | batch_f)"))
ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
D4 <- bind_rows(lapply(c("GFAP", "NfL", "Ab42_40"), function(bm) {
  d <- long %>% filter(amyloid_bl_cl20 %in% TRUE, biomarker == bm) %>% mutate(arm_f = droplevels(factor(arm)))
  f <- lmer(f_reduced, data = d, REML = TRUE, control = ctrl)
  as.data.frame(VarCorr(f)) %>% transmute(biomarker = bm, grp, var1, var2, vcov, sdcor, singular = isSingular(f))
}))
print(as.data.frame(D4 %>% mutate(across(c(vcov, sdcor), ~ signif(.x, 3)))))
say("A component with vcov 0, or a correlation of 1 or -1, is the one at the boundary. If it is the")
say("participant slope variance, the GFAP or NfL trial sizes rest on a boundary estimate; if it is batch")
say("or the participant-by-platform intercept, the trial sizes are unaffected.")
write_csv(D4, file.path(REV_DIR, "R12b_amy_pos_singular_components.csv"))
say("\n08b21 finished %s", format(Sys.time()))


# =============================================================================================
# Section 6. Amyloid-positive minus amyloid-negative validity difference, within draws
# Original script: p2_amyloid_interval.R
# =============================================================================================
# p2_amyloid_interval.R                                                   Paper 2
# Prints the 95% percentile interval, to three decimals, of the amyloid-positive minus
# amyloid-negative difference in pTau217 validity against Centiloid change, computed within
# each of the saved bootstrap draws (the same construction as Supplementary Table S26).
# Reads pipeline_out/paper2/bootstrap_draws.rds; writes nothing.
suppressPackageStartupMessages({ library(dplyr); library(tidyr) })
Sys.setenv(COHORT = "adni"); source("00_setup.R")
f <- file.path(OUT_ROOT, "pipeline_out", "paper2", "bootstrap_draws.rds")
d <- readRDS(f) %>%
  filter(biomarker == "pTau217", quantity == "validity_centiloid_slope", subset %in% c("amyloid+", "amyloid-")) %>%
  select(draw, subset, value) %>% pivot_wider(names_from = subset, values_from = value) %>%
  mutate(diff = `amyloid+` - `amyloid-`)
q <- quantile(d$diff, c(0.025, 0.975), na.rm = TRUE)
cat(sprintf("Draws used: %d\n", sum(!is.na(d$diff))))
cat(sprintf("95%% interval, three decimals: %.3f to %.3f\n", q[1], q[2]))
cat(sprintf("Check, two decimals (Table S26 shows -0.42 to -0.03): %.2f to %.2f\n", q[1], q[2]))


# =============================================================================================
# Section 8. Paired intervals from the saved draws (R13a, R13b, R13c)
# Added November 2026
# =============================================================================================
# The comparisons quoted in the text and in Supplementary Table S26 were computed within the saved
# draws and until now only printed. This section writes them, so that every interval in the paper
# has a file: R13a from the main bootstrap (section 1), R13b from the MDC bootstrap (section 2),
# R13c the shares of singular sample-model fits in the R10 bootstrap (section 3). Draw 0, where
# present, is the full data and is excluded from the intervals.
Sys.setenv(COHORT = "adni"); source("00_setup.R")
P2_DIR <- file.path(OUT_DIR, "paper2")
pair_summary <- function(df, label) {
  x <- df$diff[df$draw > 0]
  tibble(comparison = label, n_draws = sum(!is.na(x)), boot_mean = mean(x, na.rm = TRUE),
         lo95 = unname(quantile(x, 0.025, na.rm = TRUE)), hi95 = unname(quantile(x, 0.975, na.rm = TRUE)))
}
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")

# ---- R13a: from the main bootstrap draws ----------------------------------------------------
bd <- readRDS(file.path(P2_DIR, "bootstrap_draws.rds")) %>% filter(draw > 0)
wide_q <- function(q, sub) bd %>% filter(quantity == q, subset == sub) %>% select(draw, biomarker, value) %>%
  pivot_wider(names_from = biomarker, values_from = value)
r13a <- list()
# A. gain from averaging: mean of platforms minus one platform, per biomarker (all participants)
g <- bd %>% filter(subset == "all", quantity %in% c("rel_slope_one_platform", "rel_slope_mean_of_arms")) %>%
  select(draw, biomarker, quantity, value) %>% pivot_wider(names_from = quantity, values_from = value) %>%
  mutate(diff = rel_slope_mean_of_arms - rel_slope_one_platform)
for (bm in biomarkers)
  r13a[[length(r13a) + 1]] <- pair_summary(g %>% filter(biomarker == bm),
                                           sprintf("%s: mean of platforms minus one platform", bm)) %>% mutate(block = "A", biomarker = bm)
# B. between biomarkers (all participants): each minus pTau217, and Ab42_40 minus NfL
w  <- wide_q("rel_slope_one_platform", "all")
wm <- wide_q("rel_slope_mean_of_arms", "all")
for (bm in setdiff(biomarkers, "pTau217")) {
  r13a[[length(r13a) + 1]] <- pair_summary(w %>% transmute(draw, diff = .data[[bm]] - pTau217),
                                           sprintf("%s minus pTau217: one-platform reliability", bm)) %>% mutate(block = "B", biomarker = bm)
  r13a[[length(r13a) + 1]] <- pair_summary(wm %>% transmute(draw, diff = .data[[bm]] - pTau217),
                                           sprintf("%s minus pTau217: mean-of-platforms reliability", bm)) %>% mutate(block = "B", biomarker = bm)
}
r13a[[length(r13a) + 1]] <- pair_summary(w %>% transmute(draw, diff = Ab42_40 - NfL),
                                         "Ab42_40 minus NfL: one-platform reliability") %>% mutate(block = "B", biomarker = "Ab42_40")
r13a[[length(r13a) + 1]] <- pair_summary(w %>% transmute(draw, diff = NfL - GFAP),
                                         "NfL minus GFAP: one-platform reliability") %>% mutate(block = "B", biomarker = "NfL")
# C. amyloid-positive minus amyloid-negative, per biomarker: one-platform and mean reliability, validity against Centiloid
for (q in c("rel_slope_one_platform", "rel_slope_mean_of_arms", "validity_centiloid_slope")) {
  a <- bd %>% filter(quantity == q, subset %in% c("amyloid+", "amyloid-")) %>% select(draw, biomarker, subset, value) %>%
    pivot_wider(names_from = subset, values_from = value) %>% mutate(diff = `amyloid+` - `amyloid-`)
  for (bm in biomarkers)
    r13a[[length(r13a) + 1]] <- pair_summary(a %>% filter(biomarker == bm),
                                             sprintf("%s: amyloid-positive minus amyloid-negative, %s", bm, q)) %>% mutate(block = "C", biomarker = bm)
}
R13a <- bind_rows(r13a) %>% select(block, biomarker, comparison, n_draws, boot_mean, lo95, hi95)
write_csv(R13a, file.path(P2_DIR, "R13a_paired_reliability_from_draws.csv"))
say("\n== R13a: paired intervals from the main bootstrap ==")
print(as.data.frame(R13a %>% mutate(across(c(boot_mean, lo95, hi95), ~ round(.x, 3)))))

# ---- R13b: from the MDC bootstrap draws ---------------------------------------------------
md <- readRDS(file.path(P2_DIR, "mdc_bootstrap_draws.rds")) %>% filter(draw > 0, is.na(schedule)) %>%
  distinct(draw, biomarker, arm, quantity, .keep_all = TRUE)
pairs <- list(c("pTau217", "ALZpath", "Fujirebio"), c("pTau217", "ALZpath", "Janssen"), c("pTau217", "ALZpath", "C2N"),
              c("pTau217", "ALZpath", "C2N_ratio"), c("pTau217", "Janssen", "Fujirebio"), c("pTau217", "Janssen", "C2N"),
              c("pTau217", "Janssen", "C2N_ratio"), c("pTau217", "Fujirebio", "C2N"), c("pTau217", "Fujirebio", "C2N_ratio"),
              c("pTau217", "C2N", "C2N_ratio"), c("pTau181", "Roche", "Quanterix"), c("GFAP", "Roche", "Quanterix"),
              c("NfL", "Roche", "Quanterix"))
r13b <- list()
for (q in c("cv_within_pct", "mdc95_two_points_pct", "mdc95_slope_pct_per_year", "years_to_exceed_mdc_amy_pos")) {
  w <- md %>% filter(quantity == q) %>% select(draw, biomarker, arm, value)
  for (p in pairs) {
    x <- w %>% filter(biomarker == p[1], arm %in% p[2:3]) %>% pivot_wider(names_from = arm, values_from = value)
    if (!all(p[2:3] %in% names(x))) next
    x <- x %>% mutate(diff = .data[[p[3]]] - .data[[p[2]]])
    r13b[[length(r13b) + 1]] <- pair_summary(x, sprintf("%s: %s minus %s, %s", p[1], p[3], p[2], q)) %>%
      mutate(biomarker = p[1], quantity = q, arm_a = p[3], arm_b = p[2])
  }
}
R13b <- bind_rows(r13b) %>% select(biomarker, quantity, arm_a, arm_b, comparison, n_draws, boot_mean, lo95, hi95)
write_csv(R13b, file.path(P2_DIR, "R13b_paired_mdc_from_draws.csv"))
say("\n== R13b: paired intervals from the MDC bootstrap (pTau217 rows) ==")
print(as.data.frame(R13b %>% filter(biomarker == "pTau217") %>%
                      mutate(across(c(boot_mean, lo95, hi95), ~ round(.x, 2))) %>% select(comparison, boot_mean, lo95, hi95)))

# ---- R13c: shares of singular sample-model fits in the R10 bootstrap -------------------------
rd10 <- readRDS(file.path(P2_DIR, "R10_bootstrap_draws.rds")) %>% filter(draw > 0)
R13c <- rd10 %>% filter(quantity %in% c("B1_sample_model_singular", "B3_sample_model_amy_pos_singular")) %>%
  group_by(quantity, biomarker) %>% summarise(n_draws = n(), pct_singular = 100 * mean(value), .groups = "drop")
write_csv(R13c, file.path(P2_DIR, "R13c_singular_shares.csv"))
say("\n== R13c: singular sample-model fits, per cent of draws ==")
print(as.data.frame(R13c))


# =============================================================================================
# Section 7. Observed-data validity by amyloid group and Pearson design-test correlations (T7c, T17p, T8p, T8bp)
# Original script: p2_observed_checks.R
# =============================================================================================
# =============================================================================
# p2_observed_checks.R                                                 Paper 2
#
# Observed-data values that the manuscript needs but the pipeline produced only inside the
# bootstrap, and Pearson versions of the Spearman design-test correlations (28 Sep 2026).
# Reads the existing data frames and output tables. Writes four NEW files and changes no
# existing file.
#
#   1. ADNI, validity of the pooled slope by amyloid group on the observed data, fitted exactly
#      as in the bootstrap of 08b2 (change model refitted within each group, Spearman
#      correlation of the pooled slope with each reference measure; Pearson added).
#      Writes pipeline_out/paper2/T7c_observed_by_amyloid.csv
#   2. ADNI, Table S24 (T17, validity by SST quintile) with Pearson beside Spearman.
#      Writes pipeline_out/paper2/T17p_validity_by_sst_pearson.csv
#   3. A4/LEARN, Tables S23b (T8) and S23c (T8b) with Pearson beside Spearman.
#      Writes a4learn_out/paper2/T8p_design_subsets_pearson.csv and
#             a4learn_out/paper2/T8bp_validity_by_sst_pearson.csv
#
# Time    a few minutes (ten mixed-model fits)
# =============================================================================
suppressPackageStartupMessages(library(lme4))
wr <- function(df, dir, name) { f <- file.path(dir, paste0(name, ".csv")); readr::write_csv(df, f); cat("Wrote", f, "\n") }

# ---- 1 and 2: ADNI ------------------------------------------------------------------------
Sys.setenv(COHORT = "adni"); source("00_setup.R")
P2 <- file.path(OUT_DIR, "paper2")
YVAR <- "z_log_aligned"
biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
long  <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch))
critr <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>% pull(sst) %>% median()
ctrl  <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
f_full <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) + (1 + plasma_time || RID_f:arm_f) + (1|batch_f)"))
f_red  <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) + (1 | RID_f:arm_f) + (1|batch_f)"))
lrt <- read_csv(file.path(P2, "T2d_platform_slope_LRT.csv"), show_col_types = FALSE)
f_use <- if (all(lrt$p_boundary > 0.05)) f_red else f_full
cat("Model:", if (identical(f_use, f_red)) "reduced (as in the bootstrap)" else "full", "\n")
crit_cols <- intersect(c("centiloid_slope", "csf_ptau_abeta42_slope", "mmse_slope", "adas13_slope", "mpacc_digit_slope", "hcv_adj_slope"), names(critr))
FLIP <- c("hcv_adj_slope", "mpacc_digit_slope", "mmse_slope")   # oriented as in T5 and the bootstrap (section 1)
ch_obs <- function(dd, bm, lab) {
  fc <- lmer(f_use, data = dd, REML = TRUE, control = ctrl)
  v  <- as.data.frame(VarCorr(fc)) %>% filter(is.na(var2))
  vs <- v$vcov[v$grp == "RID_f" & v$var1 == "plasma_time"]
  pa <- sum(v$vcov[grepl("^RID_f:arm_f", v$grp) & !is.na(v$var1) & v$var1 == "plasma_time"])
  re <- v$vcov[v$grp == "Residual"]; k <- n_distinct(dd$arm)
  rb <- ranef(fc)$RID_f
  val <- tibble(RID = as.integer(rownames(rb)), pooled_slope = rb[["plasma_time"]]) %>%
    left_join(critr %>% select(RID, all_of(crit_cols)), by = "RID")
  bind_rows(
    tibble(biomarker = bm, subset = lab, n = nrow(val), quantity = c("rel_slope_one_platform", "rel_slope_mean_of_arms"),
           method = "-", value = c(vs / (vs + pa + re / median_sst), vs / (vs + (pa + re / median_sst) / k))),
    bind_rows(lapply(crit_cols, function(cc) tibble(biomarker = bm, subset = lab, n = sum(complete.cases(val$pooled_slope, val[[cc]])),
      quantity = paste0("validity_", cc), method = c("spearman", "pearson"),
      value = (if (cc %in% FLIP) -1 else 1) *
              c(cor(val$pooled_slope, val[[cc]], use = "complete.obs", method = "spearman"),
                cor(val$pooled_slope, val[[cc]], use = "complete.obs", method = "pearson"))))))
}
obs <- bind_rows(lapply(biomarkers, function(bm) {
  cat("  ", bm, "\n"); d <- long %>% filter(biomarker == bm)
  bind_rows(ch_obs(d, bm, "all"),
            ch_obs(d %>% filter(amyloid_bl_cl20 %in% FALSE), bm, "amyloid-"),
            ch_obs(d %>% filter(amyloid_bl_cl20 %in% TRUE),  bm, "amyloid+")) }))
wr(obs, P2, "T7c_observed_by_amyloid")
cat("\nCheck against existing outputs: pTau217 one-platform reliability should be 0.40 (all), 0.35 (amyloid-), 0.49 (amyloid+);",
    "validity against Centiloid (all) should be 0.39.\n")
print(as.data.frame(obs %>% filter(biomarker == "pTau217", quantity %in% c("rel_slope_one_platform", "validity_centiloid_slope")) %>%
                      mutate(value = round(value, 3))))

slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
T2 <- read_csv(file.path(P2, "T2_change_variance_components.csv"), show_col_types = FALSE)
sst_person <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), span = max(plasma_time) - min(plasma_time), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
T17p <- slopes %>% filter(n_obs >= 2) %>% select(RID, biomarker, arm, slope) %>%
  left_join(sst_person %>% select(RID, k, span, sst), by = "RID") %>%
  left_join(critr %>% select(RID, centiloid_slope), by = "RID") %>%
  filter(!is.na(slope), !is.na(centiloid_slope), !is.na(sst)) %>%
  group_by(biomarker) %>% mutate(band = ntile(sst, 5)) %>% ungroup() %>%
  group_by(biomarker, band) %>%
  summarise(n_person_platform = n(), sst_median = median(sst),
            r_centiloid_spearman = cor(slope, centiloid_slope, method = "spearman"),
            r_centiloid_pearson  = cor(slope, centiloid_slope, method = "pearson"), .groups = "drop") %>%
  left_join(T2 %>% select(biomarker, var_person_slope, var_residual), by = "biomarker") %>%
  mutate(reliability_predicted = var_person_slope / (var_person_slope + var_residual / sst_median)) %>%
  group_by(biomarker) %>%
  mutate(ratio_spearman = r_centiloid_spearman / r_centiloid_spearman[band == max(band)],
         ratio_pearson  = r_centiloid_pearson  / r_centiloid_pearson[band == max(band)],
         sqrt_reliability_ratio = sqrt(reliability_predicted / reliability_predicted[band == max(band)])) %>% ungroup()
wr(T17p, P2, "T17p_validity_by_sst_pearson")

# ---- 3: A4/LEARN --------------------------------------------------------------------------
Sys.setenv(COHORT = "a4learn"); source("00_setup.R")
P2 <- file.path(OUT_DIR, "paper2")
part  <- readRDS(file.path(OUT_DIR, "participants.rds")) %>% select(RID, substudy, arm, has_short_interval_pair)
pl    <- readRDS(file.path(OUT_DIR, "plasma_long.rds")) %>%
  distinct(RID, visit_no, visit_code, visit_kind, on_treatment, short_interval_pair) %>% mutate(visit_code = as.character(visit_code))
long  <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>% left_join(part, by = "RID") %>% left_join(pl, by = c("RID", "visit_no"))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>% left_join(part, by = "RID")
critr <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds"))
define_units(long)
T2 <- read_csv(file.path(P2, "T2_change_variance_components.csv"), show_col_types = FALSE)
comp <- T2 %>% filter(floor == "as_delivered") %>% select(unit, var_person_slope, var_residual)
ols <- function(t, y) if (length(t) >= 2 && var(t) > 0) cov(t, y) / var(t) else NA_real_
three <- long %>% filter(substudy == "A4", visit_code %in% c("6", "9", "66")) %>%
  select(RID, arm, visit_code, plasma_time, y = all_of(YVAR)) %>%
  group_by(RID) %>% filter(n_distinct(visit_code) == 3) %>% ungroup()
subsets <- list(`baseline + week 12` = c("6", "9"), `baseline + week 240` = c("6", "66"),
                `week 12 + week 240` = c("9", "66"), `all three draws` = c("6", "9", "66"))
sub_slopes <- bind_rows(lapply(names(subsets), function(nm) three %>% filter(visit_code %in% subsets[[nm]]) %>%
  group_by(RID, arm) %>% summarise(subset = nm, k = n(), sst = sum((plasma_time - mean(plasma_time))^2),
                                   slope = ols(plasma_time, y), .groups = "drop")))
full <- sub_slopes %>% filter(subset == "all three draws") %>% select(RID, slope_full = slope)
ss <- sub_slopes %>% left_join(full, by = "RID") %>% left_join(critr %>% select(RID, centiloid_slope), by = "RID")
design_row <- function(s, label, vs, vr) s %>% group_by(subset) %>%
  summarise(sample = label, n = sum(!is.na(slope)), sst_median = median(sst),
            reliability_predicted = vs / (vs + vr / median(sst)),
            r_centiloid_spearman = cor(slope, centiloid_slope, use = "complete.obs", method = "spearman"),
            r_centiloid_pearson  = cor(slope, centiloid_slope, use = "complete.obs", method = "pearson"), .groups = "drop")
T8p <- bind_rows(
  design_row(ss %>% filter(arm == "Placebo"), "A4 placebo", comp$var_person_slope[comp$unit == "A4 placebo"], comp$var_residual[comp$unit == "A4 placebo"]),
  design_row(ss %>% filter(arm == "Solanezumab"), "A4 solanezumab", comp$var_person_slope[comp$unit == "A4 solanezumab"], comp$var_residual[comp$unit == "A4 solanezumab"]),
  design_row(ss, "A4 both arms", comp$var_person_slope[comp$unit == "All three arms"], comp$var_residual[comp$unit == "All three arms"])) %>%
  group_by(sample) %>%
  mutate(ratio_spearman = r_centiloid_spearman / r_centiloid_spearman[subset == "all three draws"],
         ratio_pearson  = r_centiloid_pearson  / r_centiloid_pearson[subset == "all three draws"],
         sqrt_reliability_ratio = sqrt(reliability_predicted / reliability_predicted[subset == "all three draws"])) %>% ungroup()
wr(T8p, P2, "T8p_design_subsets_pearson")
sst_of <- function(d) d %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
  summarise(k = n(), span = max(plasma_time) - min(plasma_time), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop")
sst_all <- sst_of(long %>% filter(!is.na(.data[[YVAR]])))
vb <- slopes %>% filter(n_obs >= 2) %>% select(RID, arm, slope) %>%
  left_join(sst_all, by = "RID") %>% left_join(critr %>% select(RID, centiloid_slope), by = "RID") %>%
  filter(!is.na(slope), !is.na(centiloid_slope))
band_of <- function(sst) cut(sst, c(-Inf, 0.1, 2, 8, 12, Inf), labels = c("< 0.1", "0.1 to 2", "2 to 8", "8 to 12", "> 12"))
T8bp <- bind_rows(lapply(MEASURES, function(m) {
  ids <- unique(unit_data(m)$RID); cp <- comp %>% filter(unit == m_label[[m]])
  vb %>% filter(RID %in% ids) %>% mutate(band = band_of(sst)) %>% group_by(band) %>%
    summarise(unit = m_label[[m]], n = n(), sst_median = median(sst),
              reliability_predicted = cp$var_person_slope / (cp$var_person_slope + cp$var_residual / median(sst)),
              r_centiloid_spearman = if (n() >= 10) cor(slope, centiloid_slope, method = "spearman") else NA_real_,
              r_centiloid_pearson  = if (n() >= 10) cor(slope, centiloid_slope, method = "pearson")  else NA_real_, .groups = "drop") }))
wr(T8bp, P2, "T8bp_validity_by_sst_pearson")
cat("\nA4 design test, Spearman against Pearson ratios:\n")
print(as.data.frame(T8p %>% select(sample, subset, ratio_spearman, ratio_pearson, sqrt_reliability_ratio) %>% mutate(across(where(is.numeric), ~ round(.x, 2)))))
cat("\nDone. Four files written; no existing file was changed.\n")
Sys.setenv(COHORT = "adni")


# =============================================================================================
# Section 9. Five against three samples over 24 months, amyloid-positive components (R15a)
# Added November 2026. The percentage reduction in participants per arm when a 24-month trial takes five
# samples instead of three (25% slowing), per pTau217 measure and for the mean of platforms under both
# bounds, with percentile intervals from the section 3 draws (R10_bootstrap_draws.rds). Point estimates are
# from R3a (06, section 3). Nothing is refitted.
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages({ library(lme4); library(tidyr) })
P2_DIR <- file.path(OUT_DIR, "paper2")
need <- file.path(P2_DIR, c("R10_bootstrap_draws.rds", "bootstrap_draws.rds"))
if (!all(file.exists(need))) stop("Sections 9 and 10 need the saved bootstrap draws (",
  paste(basename(need[!file.exists(need)]), collapse = ", "), "); run sections 1 and 3 with N_BOOT > 0 first.")
ci95 <- function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE, names = FALSE)
say("\n== Section 9. Five against three samples over 24 months (R15a) ==")
S3 <- "24 months, 3 draws (0, 12, 24 m)"
S5 <- "24 months, 5 draws (0, 6, 12, 18, 24 m)"
r10 <- readRDS(file.path(P2_DIR, "R10_bootstrap_draws.rds")) %>%
  filter(quantity == "B3_n_per_arm_amy_pos", biomarker == "pTau217",
         abs(as.numeric(effect) - 0.25) < 1e-9, schedule %in% c(S3, S5)) %>%
  select(draw, arm, schedule, value) %>%
  distinct(draw, arm, schedule, .keep_all = TRUE) %>%
  pivot_wider(names_from = schedule, values_from = value) %>%
  mutate(reduction_pct = 100 * (1 - .data[[S5]] / .data[[S3]]))
r3a <- read_csv(file.path(P2_DIR, "R3a_trial_amyloid_positive_variance.csv"), show_col_types = FALSE) %>%
  filter(biomarker == "pTau217", abs(as.numeric(effect) - 0.25) < 1e-9, schedule %in% c(S3, S5)) %>%
  mutate(arm = case_when(bound == "per platform" ~ arm,
                         grepl("^upper", bound) ~ "mean of platforms upper",
                         grepl("^lower", bound) ~ "mean of platforms lower")) %>%
  select(arm, schedule, n = n_per_arm_slope_amy_pos_variance) %>%
  pivot_wider(names_from = schedule, values_from = n) %>%
  mutate(estimate = 100 * (1 - .data[[S5]] / .data[[S3]]))
R15a <- r10 %>% filter(draw > 0) %>% group_by(arm) %>%
  summarise(n_draws = sum(!is.na(reduction_pct)), lo95 = ci95(reduction_pct)[1], hi95 = ci95(reduction_pct)[2],
            .groups = "drop") %>%
  left_join(r3a %>% select(arm, n_three = all_of(S3), n_five = all_of(S5), estimate), by = "arm") %>%
  select(arm, n_three, n_five, estimate, lo95, hi95, n_draws)
print(as.data.frame(R15a %>% mutate(across(where(is.numeric), ~ round(.x, 1)))))
write_csv(R15a, file.path(P2_DIR, "R15a_five_against_three_samples.csv"))
say("Written: R15a_five_against_three_samples.csv in %s", P2_DIR)


# =============================================================================================
# Section 10. Pooled minus single-platform pTau217 validity against Centiloid change (R15b)
# Added November 2026. The pooled value of each draw is read from the section 1 cache (bootstrap_draws.rds);
# the single-platform value is recomputed on the same resample, rebuilt from the same seed (20260904 + draw)
# and the same participant list. Draws 1 to 3 are first checked against the cached slope ICC, and the section
# stops if they do not match, since the intervals would then not be paired.
# =============================================================================================
say("\n== Section 10. Pooled minus single pTau217 platform, validity against Centiloid (R15b) ==")
BOOT_SEED <- 20260904
long_b <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
rids_b <- unique(long_b$RID)                  # the participant list section 1 samples from
slopes_b <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform)) %>%
  filter(biomarker == "pTau217") %>% select(RID, arm, slope)
critr_b <- readRDS(file.path(OUT_DIR, "frame_paper2_criteria.rds")) %>% select(RID, centiloid_slope)
T5_b <- read_csv(file.path(P2_DIR, "T5_validity_vs_criterion_change.csv"), show_col_types = FALSE) %>%
  filter(biomarker == "pTau217", criterion == "centiloid_slope")
cache_b <- readRDS(file.path(P2_DIR, "bootstrap_draws.rds"))
pooled_b <- cache_b %>% filter(biomarker == "pTau217", subset == "all", quantity == "validity_centiloid_slope") %>%
  distinct(draw, .keep_all = TRUE) %>% select(draw, pooled = value)
single_validity <- function(picked) {
  tibble(RID = picked) %>%
    inner_join(slopes_b, by = "RID", relationship = "many-to-many") %>%
    left_join(critr_b, by = "RID") %>%
    group_by(arm) %>%
    summarise(r = cor(slope, centiloid_slope, use = "complete.obs", method = "spearman"), .groups = "drop")
}
obs_single <- single_validity(unique(slopes_b$RID)) %>%
  left_join(T5_b %>% filter(estimate == "single measure") %>% select(arm, r_T5 = r), by = "arm")
print(as.data.frame(obs_single %>% mutate(r = round(r, 3))))

ctrl_b <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5), calc.derivs = FALSE)
icc_cache <- cache_b %>% filter(biomarker == "pTau217", quantity == "icc_single_arm") %>%
  distinct(draw, .keep_all = TRUE) %>% select(draw, icc_cache = value)
pair_chk <- lapply(1:3, function(b) {
  set.seed(BOOT_SEED + b)
  picked <- sample(rids_b, length(rids_b), replace = TRUE)
  boot_ids <- tibble(RID = picked, boot_id = seq_along(picked))
  d <- long_b %>% filter(biomarker == "pTau217") %>%
    inner_join(boot_ids, by = "RID", relationship = "many-to-many")
  sl <- slope_by(d, "z_log_aligned", "plasma_time", boot_id, arm) %>% select(boot_id, arm, slope) %>%
    pivot_wider(names_from = arm, values_from = slope)
  mc <- sl[complete.cases(sl), -1]
  dd <- mc %>% mutate(id = row_number()) %>% pivot_longer(-id, names_to = "arm", values_to = "s")
  fi <- lmer(s ~ 1 + (1|id) + (1|arm), data = dd, control = ctrl_b)
  v <- as.data.frame(VarCorr(fi))
  tibble(draw = b, icc_rebuilt = v$vcov[v$grp == "id"] / sum(v$vcov))
}) %>% bind_rows() %>% left_join(icc_cache, by = "draw")
print(as.data.frame(pair_chk %>% mutate(across(starts_with("icc"), ~ round(.x, 4)))))
if (any(abs(pair_chk$icc_rebuilt - pair_chk$icc_cache) > 1e-4, na.rm = TRUE) || anyNA(pair_chk$icc_cache))
  stop("Section 10: the rebuilt resamples do not match the cached draws, so the intervals would not be paired.")

draws_b <- lapply(pooled_b$draw[pooled_b$draw > 0], function(b) {
  set.seed(BOOT_SEED + b)
  picked <- sample(rids_b, length(rids_b), replace = TRUE)
  single_validity(picked) %>% mutate(draw = b)
}) %>% bind_rows() %>% left_join(pooled_b, by = "draw") %>% mutate(diff = pooled - r)
obs_pooled <- T5_b %>% filter(estimate == "pooled") %>% pull(r)
R15b <- draws_b %>% group_by(arm) %>%
  summarise(n_draws = sum(!is.na(diff)), lo95 = ci95(diff)[1], hi95 = ci95(diff)[2], .groups = "drop") %>%
  left_join(obs_single %>% select(arm, single = r_T5), by = "arm") %>%
  mutate(pooled = obs_pooled, estimate = pooled - single) %>%
  select(arm, pooled, single, estimate, lo95, hi95, n_draws)
print(as.data.frame(R15b %>% mutate(across(where(is.numeric), ~ round(.x, 3)))))
write_csv(R15b, file.path(P2_DIR, "R15b_pooled_minus_single_validity.csv"))
say("Written: R15b_pooled_minus_single_validity.csv in %s", P2_DIR)


# =============================================================================================
# Section 11. Trial sizes on whole-cohort components, mean of platforms, lower bound (R17)
# Added November 2026. The whole-cohort sensitivity column of the trial-size table gives the pTau217 mean of
# platforms under the lower bound (R2a, 06 section 2); this section adds its intervals. Per draw (the same
# participant resamples and seeds as section 3): the pTau217 sample model on the whole cohort (true slope,
# shared occasion and assay-specific variance), the mean amyloid-positive rate (mean over the five measures
# of the mean least-squares slope in amyloid-positive participants with three or more samples), and
# n per arm = ceiling(2 (z0.975 + z0.80)^2 v / (effect x |rate|)^2), v = var(true slope) +
# (var(shared) + var(assay-specific)/5) / SST of the schedule. Draw 0 must reproduce R2a. Saved every 10 draws
# and resumed if interrupted. Writes R17_trial_lower_bound_whole_cohort.csv.
# =============================================================================================
Sys.setenv(COHORT = "adni"); source("00_setup.R")
if (COHORT != "adni") stop("COHORT must be 'adni'.")
suppressPackageStartupMessages(library(lme4))

REV_DIR    <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
draws_file <- file.path(REV_DIR, "R17_trial_lower_bound_whole_cohort_draws.rds")
N_BOOT     <- as.integer(Sys.getenv("N_BOOT", unset = "1000"))
BOOT_SEED  <- 20260904
YVAR       <- "z_log_aligned"
zsum2      <- (qnorm(0.975) + qnorm(0.80))^2
ctrl       <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))

long_all <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds"))
rids <- unique(long_all$RID)          # same participant list and order as 09, so draws pair
long <- long_all %>%
  mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
         batch_f = factor(batch)) %>%
  filter(biomarker == "pTau217")
k_platforms <- n_distinct(long$arm)
stopifnot(k_platforms == 5)

designs <- tibble(schedule = c("18 months, 4 draws (0, 6, 12, 18 m)",
                               "24 months, 3 draws (0, 12, 24 m)",
                               "24 months, 5 draws (0, 6, 12, 18, 24 m)"),
                  sst = c(sum((c(0, 0.5, 1, 1.5) - 0.75)^2),
                          sum((c(0, 1, 2) - 1)^2),
                          sum((c(0, 0.5, 1, 1.5, 2) - 1)^2)))
effects <- c(0.25, 0.50)

f_sample <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                             "(1 | RID_f:arm_f) + (1 | sample_f) + (1 | batch_f)"))
vget <- function(fit, grp, var1 = NA) {
  vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
  i <- if (is.na(var1)) which(vc$grp == grp) else which(vc$grp == grp & !is.na(vc$var1) & vc$var1 == var1)
  if (length(i)) vc$vcov[i[1]] else NA_real_
}

one_draw <- function(b) {
  boot_ids <- if (b == 0) tibble(RID = rids, boot_id = seq_along(rids)) else {
    set.seed(BOOT_SEED + b)
    picked <- sample(rids, length(rids), replace = TRUE)
    tibble(RID = picked, boot_id = seq_along(picked))
  }
  d <- long %>% inner_join(boot_ids, by = "RID", relationship = "many-to-many") %>%
    mutate(RID_f = factor(boot_id), sample_f = factor(paste(boot_id, visit_no)),
           arm_f = droplevels(factor(arm)))
  fs <- tryCatch(suppressMessages(lmer(f_sample, data = d, REML = TRUE, control = ctrl)),
                 error = function(e) NULL)
  if (is.null(fs)) return(NULL)
  vt <- vget(fs, "RID_f", "plasma_time"); vsh <- vget(fs, "sample_f"); va <- vget(fs, "Residual")
  # mean amyloid-positive rate: per-measure mean OLS slope (>= 3 samples), averaged over measures
  m <- d %>% group_by(boot_id, arm) %>%
    filter(n() >= 3, sd(plasma_time) > 0) %>%
    summarise(slope = cov(plasma_time, z_log_aligned) / var(plasma_time),
              amy = first(amyloid_bl_cl20), .groups = "drop") %>%
    filter(amy %in% TRUE) %>%
    group_by(arm) %>% summarise(m = mean(slope, na.rm = TRUE), .groups = "drop") %>%
    summarise(m = mean(m, na.rm = TRUE)) %>% pull(m)
  designs %>% crossing(effect = effects) %>%
    mutate(v = vt + (vsh + va / k_platforms) / sst,
           reliability = vt / v,
           n_per_arm = ceiling(2 * zsum2 * v / (effect * abs(m))^2),
           singular = isSingular(fs), draw = b)
}

if (N_BOOT > 0) {
# ---- draw 0: must reproduce R2a ----------------------------------------------------------
est <- one_draw(0)
R2a <- read_csv(file.path(REV_DIR, "R2a_trial_mean_of_platforms_bounds.csv"), show_col_types = FALSE) %>%
  filter(biomarker == "pTau217", grepl("^lower", bound)) %>%
  select(schedule, effect, n_R2a = n_per_arm_slope)
chk <- est %>% left_join(R2a, by = c("schedule", "effect"))
print(as.data.frame(chk %>% select(schedule, effect, n_per_arm, n_R2a)))
stopifnot(all(abs(chk$n_per_arm - chk$n_R2a) <= 1))
say("Draw 0 reproduces R2a. Starting %d bootstrap draws.", N_BOOT)

# ---- bootstrap, resumable -------------------------------------------------------------------
draws <- if (file.exists(draws_file)) readRDS(draws_file) else est
done  <- unique(draws$draw)
todo  <- setdiff(seq_len(N_BOOT), done)
t0 <- Sys.time()
for (i in seq_along(todo)) {
  b <- todo[i]
  r <- one_draw(b)
  if (!is.null(r)) draws <- bind_rows(draws, r)
  if (i == 3) say("Projected time for %d draws: %.0f minutes", length(todo),
                  as.numeric(difftime(Sys.time(), t0, units = "mins")) / 3 * length(todo))
  if (i %% 10 == 0 || i == length(todo)) { tmp <- paste0(draws_file, ".tmp"); saveRDS(draws, tmp); file.rename(tmp, draws_file); say("Saved %d of %d draws", i, length(todo)) }
}

# ---- intervals ----------------------------------------------------------------------------------
boot <- draws %>% filter(draw > 0)
ci <- boot %>% group_by(schedule, effect) %>%
  summarise(n_draws = n(), pct_singular = round(100 * mean(singular), 1),
            n_lo95 = quantile(n_per_arm, 0.025), n_hi95 = quantile(n_per_arm, 0.975),
            rel_lo95 = quantile(reliability, 0.025), rel_hi95 = quantile(reliability, 0.975),
            .groups = "drop")
res <- est %>% select(schedule, effect, reliability, n_per_arm) %>%
  left_join(ci, by = c("schedule", "effect"))
write_csv(res, file.path(REV_DIR, "R17_trial_lower_bound_whole_cohort.csv"))
say("\nWhole-cohort components, pTau217 mean of platforms lower bound (25%% slowing):")
print(as.data.frame(res %>% filter(effect == 0.25) %>%
                      transmute(schedule, n = sprintf("%s (%s to %s)", format(n_per_arm, big.mark = ",", trim = TRUE),
                                                      format(round(n_lo95), big.mark = ",", trim = TRUE),
                                                      format(round(n_hi95), big.mark = ",", trim = TRUE)),
                                pct_singular)))
say("Written: %s", file.path(REV_DIR, "R17_trial_lower_bound_whole_cohort.csv"))
} else say("N_BOOT = 0: section 11 skipped.")
