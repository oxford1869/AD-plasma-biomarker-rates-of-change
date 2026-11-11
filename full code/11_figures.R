# =============================================================================================
# 11_figures.R
# All figures.
#
# ADNI, reading the finished tables of both cohorts. Working figures (pipeline_out/paper2/figures/), journal
# versions with legends (figures/final/), the manuscript composites Figure 1 to Figure 4 (figures/main/) and
# Supplementary Figures S1 to S5 (figures/supplementary/), with manuscript_figure_legends.md. Formerly 08b7.
# No model fits.
# =============================================================================================

# =============================================================================
# 08b7_paper2_figures.R
# Paper 2, task 14: figures. Reads the frames from 07b and the tables written
# by 08b, 08b2, 08b5 and 08b6. No model fits.
#
#   F7  spaghetti plots: every person's series on every platform, per biomarker
#       (one file per biomarker; panels are platforms; colour is baseline amyloid)
#   F8  per-person three-draw slopes by platform with the minimal detectable
#       change band (the slope that noise alone can produce)
#   F9a slope agreement between platforms (Pearson r, T4)
#   F9b level agreement between platforms (centred ICC, T13)
#   F10a level variance components per biomarker (T1)
#   F10b reliability of an individual slope, one platform and mean of platforms,
#        with bootstrap intervals (T2, T2b, T7), all people and by amyloid status
#   F11 validity: correlation of the pooled slope with concurrent criterion
#       change, with bootstrap intervals (T5, T7)
#   F12 years between two draws before the mean amyloid-positive change
#       exceeds the two-point mdc95, per platform (T14)
#
# Design rules used throughout (kept deliberately plain for a journal):
#   colour carries one meaning per figure and is never used for identity that
#   position already gives; biomarkers and platforms are panels or positions,
#   not hues. Two series colours only (blue and orange, checked for colour
#   vision deficiency), a neutral grey for "no information", a single-hue blue
#   ramp for magnitude. Text is always ink, never a series colour. Hairline
#   horizontal gridlines, no vertical ones, no panel boxes, no shaded strips.
#   Files are written as PDF (vector) and PNG (300 dpi) with the same name.
# =============================================================================

Sys.setenv(COHORT = "adni"); source("00_setup.R")
suppressPackageStartupMessages(library(ggplot2))
P2_DIR  <- file.path(OUT_DIR, "paper2")
FIG_DIR <- file.path(P2_DIR, "figures"); dir.create(FIG_DIR, showWarnings = FALSE)
rd <- function(name) read_csv(file.path(P2_DIR, paste0(name, ".csv")), show_col_types = FALSE)
# Every figure is written twice: a working version with title and subtitle on
# the figure (figures/), and a journal version (figures/final/) with the title
# and subtitle removed and the full legend printed under the plot, wrapped to
# the figure width. The same legends are collected in figures/final/figure_legends.md.
FINAL_DIR <- file.path(FIG_DIR, "final"); dir.create(FINAL_DIR, showWarnings = FALSE)
LEG <- list()   # filled below, before the figures are drawn

# capabilities("cairo") can report TRUE on a machine where the cairo library then fails to
# load, which produced one "failed to load cairo DLL" warning per PDF. Test it once, for real.
# The failure is signalled as a warning rather than an error, so the test has to catch both.
cairo_works <- function() {
  f <- tempfile(fileext = ".pdf")
  on.exit({ if (grDevices::dev.cur() > 1) try(grDevices::dev.off(), silent = TRUE); unlink(f) }, add = TRUE)
  warned <- FALSE
  withCallingHandlers({ grDevices::cairo_pdf(f); grid::grid.rect() },
                      warning = function(w) { warned <<- TRUE; invokeRestart("muffleWarning") })
  !warned
}
CAIRO_OK <- isTRUE(capabilities("cairo")) && isTRUE(tryCatch(cairo_works(), error = function(e) FALSE))
if (!CAIRO_OK) say("cairo is unavailable on this machine; PDFs use the base device")
sv <- function(p, name, w, h, legend_key = name) {
  dev <- if (CAIRO_OK) grDevices::cairo_pdf else "pdf"
  ggsave(file.path(FIG_DIR, paste0(name, ".pdf")), p, width = w, height = h, device = dev, bg = "white")
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
  cap <- LEG[[legend_key]]
  if (is.null(cap)) cap <- ""
  wrapped <- paste(strwrap(cap, width = round(w * 17)), collapse = "\n")   # about 17 characters per inch at 7.5 pt
  n_lines <- length(strsplit(wrapped, "\n")[[1]])
  q <- p + labs(title = NULL, subtitle = NULL, caption = wrapped) +
    theme(plot.margin = margin(8, 12, 6, 8),
          plot.caption = element_text(size = 7.5, colour = INK2, hjust = 0, lineheight = 1.15,
                                      margin = margin(t = 10)))
  h_final <- h - 0.45 + 0.13 * n_lines
  ggsave(file.path(FINAL_DIR, paste0(name, ".pdf")), q, width = w, height = h_final, device = dev, bg = "white")
  ggsave(file.path(FINAL_DIR, paste0(name, ".png")), q, width = w, height = h_final, dpi = 300, bg = "white")
}

# ---- palette and theme --------------------------------------------------------------------
# Defined once in 00_setup.R so every script in the pipeline draws in the same style.

biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
bm_label   <- c(pTau217 = "pTau217", pTau181 = "pTau181", GFAP = "GFAP", NfL = "NfL", Ab42_40 = "Ab42/40")
arm_levels <- c("ALZpath", "C2N", "C2N %pTau217", "Fujirebio", "Janssen", "Roche", "Quanterix")
arm_of  <- function(measure, platform) case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N %pTau217", TRUE ~ platform)
fix_arm <- function(x) factor(str_replace(x, "^C2N_ratio$", "C2N %pTau217"), arm_levels)
amy_lab <- function(x) factor(case_when(x %in% TRUE ~ "amyloid positive", x %in% FALSE ~ "amyloid negative",
                                        TRUE ~ "no PET"), names(pal_amy))
lab_bm  <- as_labeller(bm_label)

long   <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
  mutate(arm = fix_arm(arm_of(measure, platform)), amyloid = amy_lab(amyloid_bl_cl20),
         biomarker = factor(biomarker, biomarkers))
slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
  mutate(arm = fix_arm(arm_of(measure, platform)), amyloid = amy_lab(amyloid_bl_cl20),
         biomarker = factor(biomarker, biomarkers))

# ---- figure legends (printed under the journal versions, and collected in figure_legends.md) ----
T7 <- rd("T7_bootstrap_intervals")
n_people <- n_distinct(long$RID)
n_core   <- long %>% count(RID, biomarker, arm, name = "k") %>% filter(k >= 3) %>% distinct(RID) %>% nrow()
n_boot   <- suppressWarnings(max(T7$n_draws, na.rm = TRUE))
for (bm in biomarkers) LEG[[paste0("F7_spaghetti_", bm)]] <- sprintf(
  "Figure. Individual plasma %s trajectories on every platform. Each thin line is one participant's series of measurements on the platform named above the panel (N = %d participants, up to three samples over a median of five years). Values are log-transformed and standardised within each measure, with the sign set so that higher values indicate greater abnormality%s. Heavy lines are ordinary least-squares trends by amyloid PET status at the first plasma sample (above 20 Centiloids). Because every platform assayed the same aliquots, panels differ only in how each platform measured the same samples. The vertical axis is truncated at ±4 standardised units; a small number of values lie outside it.",
  bm_label[[bm]], n_people, if (bm == "Ab42_40") " (the ratio is reversed)" else "")
LEG[["F8_slopes_with_mdc_band"]] <- sprintf(
  "Figure. Individual rates of change by platform, against the change that noise alone can produce. Points are per-participant slopes from three plasma samples (standardised log units per year; N = %d participants with three samples), coloured by amyloid PET status at the first sample; boxes show the median and interquartile range. The shaded band is the 95%% minimal detectable change for a three-sample slope on that platform (\u00b11.96 \u00d7 the standard error of a slope estimated from within-person, within-platform residual variation over the median sampling span), that is, the range of slopes expected in 95%% of people whose true rate of change is zero. Slopes inside the band cannot be distinguished from no change.", n_core)
LEG[["F9a_slope_agreement_heatmaps"]] <-
  "Figure. Agreement between platforms about an individual's rate of change. Pearson correlation between per-participant three-sample slopes estimated separately on each pair of platforms measuring the same biomarker on the same samples. Diagonal cells are blank. Higher values mean that the two platforms rank people's rates of change similarly."
LEG[["F9b_level_agreement_heatmaps"]] <-
  "Figure. Agreement between platforms about an individual's level. Intraclass correlation between platforms measuring the same biomarker at the same visit, computed on the log scale after removing each platform's mean, so that fixed calibration differences between platforms do not count as disagreement. Diagonal cells are blank."
LEG[["F10a_level_variance_components"]] <-
  "Figure. Sources of variance in a single plasma measurement. Variance components from a crossed random-effects model per biomarker fitted to all platforms and all visits (person; platform; person \u00d7 platform; sample, meaning the person \u00d7 visit effect shared by all platforms assaying that aliquot; platform \u00d7 visit; batch within platform; residual), expressed as a percentage of the total. The person component is the between-person variance in true level; every other component is error for a single measurement on one platform. Segments below 5% are unlabelled."
LEG[["F10b_slope_reliability"]] <- sprintf(
  "Figure. Reliability of an individual's rate of change. Reliability is the share of the variance of an observed three-sample slope that is true between-person variance in rate of change, from the multi-platform change model: var(true slope) / [var(true slope) + residual variance / SST], where SST is the median within-person spread of sampling times. Shown for a slope from one platform and for the mean of all platforms measuring that biomarker, in all participants and by amyloid PET status at the first sample. Filled symbols are the joint-model values, which treat variation shared across platforms at a blood draw as true change and are upper bounds; open symbols (pTau217, all participants) are the lower bounds from the model with a participant-by-sample effect; grey bands span the reliability of each measure fitted alone, from the lowest to the highest platform. Vertical lines are 95%% percentile intervals from %d bootstrap resamples of participants.", n_boot)
LEG[["F11_validity_pooled_slope"]] <- sprintf(
  "Figure. Validity of the pooled plasma rate of change against concurrent change in four criteria. Spearman correlation between each participant's pooled plasma slope (the shared person slope from the multi-platform change model) and the participant's concurrent rate of change in amyloid PET (Centiloid), CSF pTau181/Ab42, adjusted hippocampal volume, ADAS-Cog13, mPACC and MMSE. Signs are set so that a positive correlation always means the plasma slope tracks worsening on the criterion. Vertical lines are 95%% percentile intervals from %d bootstrap resamples of participants.", n_boot)
LEG[["F12_years_to_exceed_mdc"]] <-
  "Figure. Time needed before a change exceeds measurement noise. For each platform, the number of years between two plasma samples that a participant changing at the mean rate observed in amyloid-positive participants would need before the difference between the two measurements exceeds the 95% minimal detectable change for that platform. Bars are truncated at 30 years (pale bars) where the mean rate of change is too small for the difference to exceed noise within a realistic interval."
ABBREV <- "Abbreviations: ADAS-Cog13, Alzheimer's Disease Assessment Scale, cognitive subscale, 13 items; CSF, cerebrospinal fluid; GFAP, glial fibrillary acidic protein; mPACC, modified Preclinical Alzheimer Cognitive Composite; MMSE, Mini-Mental State Examination; NfL, neurofilament light chain; PET, positron emission tomography; pTau, phosphorylated tau; SST, sum of squared centred sampling times."

# ---- F7 spaghetti ---------------------------------------------------------------------------
say("F7 spaghetti plots")
p7 <- list()
for (bm in biomarkers) {
  d <- long %>% filter(biomarker == bm)
  n_arm <- n_distinct(d$arm)
  p <- ggplot(d, aes(plasma_time, z_log_aligned, group = RID, colour = amyloid)) +
    geom_line(alpha = 0.22, linewidth = 0.25) +
    geom_smooth(aes(group = amyloid), method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.1) +
    facet_wrap(~ arm, nrow = 1) +
    scale_colour_manual(values = pal_amy) +
    scale_x_continuous(breaks = seq(0, 10, 2)) +
    coord_cartesian(ylim = c(-4, 4)) +   # a handful of extreme values would otherwise set the axis
    labs(x = "Years from first plasma sample", y = "Standardised log value",
         title = paste0(bm_label[[bm]], ": every person's series on every platform"),
         subtitle = "Standardised so that higher is more abnormal. Thin lines are individuals; heavy lines are the group linear trend by baseline amyloid PET status") +
    guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1.2)))
  sv(p, paste0("F7_spaghetti_", bm), w = 1.9 * n_arm + 1.4, h = 3.9)
  p7[[bm]] <- p
}

# ---- F8 slopes with the MDC band -------------------------------------------------------------
say("F8 slope distributions with MDC bands")
T3 <- rd("T3_minimal_detectable_change") %>%
  mutate(arm = fix_arm(arm), biomarker = factor(biomarker, biomarkers))
s3 <- slopes %>% filter(n_obs >= 3, !is.na(slope))
p8 <- ggplot() +
  geom_tile(data = T3, aes(x = arm, y = 0, height = 2 * mdc95_slope_per_year), width = 0.9, fill = BAND) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.3) +
  geom_jitter(data = s3, aes(arm, slope, colour = amyloid), width = 0.24, height = 0, size = 0.45, alpha = 0.5,
              shape = 16) +
  geom_boxplot(data = s3, aes(arm, slope), width = 0.32, outlier.shape = NA, fill = NA,
               colour = INK, linewidth = 0.4) +
  facet_wrap(~ biomarker, scales = "free_x", nrow = 1, labeller = lab_bm) +
  scale_colour_manual(values = pal_amy) +
  coord_cartesian(ylim = c(-0.6, 0.8)) +
  labs(x = NULL, y = "Three-draw slope (standardised log units per year)",
       title = "Individual rates of change by platform",
       subtitle = "Shaded band: the range of slopes that measurement noise alone produces in 95% of people with no true change (mdc95)") +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
sv(p8, "F8_slopes_with_mdc_band", w = 11, h = 4.2)

# ---- F9 agreement heatmaps ------------------------------------------------------------------
say("F9 agreement heatmaps")
sym_pairs <- function(df, value) {
  a <- df %>% separate(pair, c("a", "b"), sep = " vs ") %>%
    mutate(a = fix_arm(a), b = fix_arm(b)) %>% select(biomarker, a, b, v = all_of(value))
  arms <- bind_rows(a %>% distinct(biomarker, arm = a), a %>% distinct(biomarker, arm = b)) %>% distinct()
  diag <- arms %>% transmute(biomarker, a = arm, b = arm, v = NA_real_)
  bind_rows(a, a %>% rename(a = b, b = a), diag) %>%
    mutate(biomarker = factor(biomarker, biomarkers), label = if_else(is.na(v), "", sprintf("%.2f", v)),
           ink = if_else(!is.na(v) & v > 0.6, "white", INK))
}
heat <- function(df, title, subtitle, fill_lab) {
  ggplot(df, aes(a, b, fill = v)) +
    geom_tile(colour = "white", linewidth = 1.2) +
    geom_text(aes(label = label, colour = ink), size = 2.7) +
    scale_colour_identity() +
    facet_wrap(~ biomarker, scales = "free", nrow = 1, labeller = lab_bm) +
    scale_fill_gradientn(colours = SEQ, limits = c(0, 1), breaks = c(0, 0.5, 1), name = fill_lab,
                         na.value = BAND, oob = scales::squish) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme(panel.grid = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
          axis.text.x.top = element_text(angle = 35, hjust = 0, vjust = 0),
          axis.text = element_text(colour = INK2),
          legend.position = "right", legend.title = element_text(size = 8, colour = INK2),
          legend.key.height = unit(18, "pt"), legend.key.width = unit(8, "pt"))
}
T4  <- rd("T4_slope_agreement")
p9a <- heat(sym_pairs(T4, "r_pearson"), "Agreement about change between platforms",
            "Pearson correlation between per-person three-draw slopes", "r")
sv(p9a, "F9a_slope_agreement_heatmaps", w = 12, h = 3.4)
if (file.exists(file.path(P2_DIR, "T13_level_agreement.csv"))) {
  T13 <- rd("T13_level_agreement")
  p9b <- heat(sym_pairs(T13, "icc_centred"), "Agreement about level between platforms",
              "Intraclass correlation between platforms measuring the same sample, log scale, mean offsets removed", "ICC")
  sv(p9b, "F9b_level_agreement_heatmaps", w = 12, h = 3.4)
}

# ---- F10 variance components and reliability ------------------------------------------------
say("F10 variance components and reliability")
comp_levels <- c("person", "person x platform", "platform", "sample (person x occasion)",
                 "platform x occasion", "batch within platform", "residual")
comp_labels <- c("Person (true level)", "Person x platform", "Platform", "Sample (person x occasion)",
                 "Platform x occasion", "Batch within platform", "Residual")
# the signal in blue; every error component in a graded neutral so the eye reads signal against noise
comp_fill <- c(BLUE, "#8f8b82", "#d8d4cb", "#6b675f", "#e8e5dd", "#aeaaa1", "#c9c5bc")
names(comp_fill) <- comp_labels
dark_fill <- comp_labels[c(1, 2, 4)]   # segments whose labels need white ink
T1 <- rd("T1_level_variance_components") %>%
  mutate(biomarker = factor(biomarker, biomarkers),
         component = factor(comp_labels[match(component, comp_levels)], rev(comp_labels))) %>%
  group_by(biomarker) %>% arrange(desc(component), .by_group = TRUE) %>%   # stacking order: person at the bottom
  mutate(pos = cumsum(pct) - pct / 2, lab = if_else(pct >= 5, sprintf("%.0f", pct), "")) %>% ungroup()
p10a <- ggplot(T1, aes(biomarker, pct, fill = component)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.6) +
  geom_text(aes(y = pos, label = lab, colour = if_else(component %in% dark_fill, "white", INK)), size = 2.7) +
  scale_colour_identity() +
  scale_fill_manual(values = comp_fill, breaks = comp_labels) +
  scale_x_discrete(labels = bm_label) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)), labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Share of variance in a single measurement",
       title = "Where the variance in one plasma measurement comes from",
       subtitle = "Blue is between-person variance in true level; every other component is error for a single measurement on one platform") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank())
sv(p10a, "F10a_level_variance_components", w = 7.5, h = 4.6)

T7  <- rd("T7_bootstrap_intervals")
T2  <- rd("T2_change_variance_components")
T2b <- rd("T2b_change_by_amyloid_status") %>%
  separate(model, c("biomarker", "subset"), sep = " \\| ") %>%
  select(biomarker, subset, rel_slope_one_platform, rel_slope_mean_of_arms)
rel_pt <- bind_rows(T2 %>% mutate(subset = "all") %>% select(biomarker, subset, rel_slope_one_platform, rel_slope_mean_of_arms),
                    T2b) %>%
  pivot_longer(starts_with("rel_"), names_to = "quantity", values_to = "estimate")
rel_ci <- T7 %>% filter(quantity %in% c("rel_slope_one_platform", "rel_slope_mean_of_arms")) %>%
  select(biomarker, subset, quantity, ci_low, ci_high)
rel <- rel_pt %>% left_join(rel_ci, by = c("biomarker", "subset", "quantity")) %>%
  mutate(biomarker = factor(biomarker, biomarkers),
         quantity = factor(recode(quantity, rel_slope_one_platform = "one platform",
                                  rel_slope_mean_of_arms = "mean of all platforms"),
                           c("one platform", "mean of all platforms")),
         subset = factor(recode(subset, all = "All people", `amyloid-` = "Amyloid negative", `amyloid+` = "Amyloid positive"),
                         c("All people", "Amyloid negative", "Amyloid positive")))
# Revision, November 2026: the joint-model reliabilities are upper bounds. R10 (08b19) adds, for
# pTau217 in all participants, the lower bounds from the model with a participant-by-sample
# effect (open symbols), and for every biomarker the range of each measure fitted alone (grey
# band), which is the like-for-like comparison with the second cohort. Without R10 the panel is
# drawn as before.
REV_DIR <- file.path(OUT_ROOT, "pipeline_out", "paper2")   # the R tables sit with the T tables
r10_file <- file.path(REV_DIR, "R10_bootstrap_intervals.csv")
rel <- rel %>% mutate(bound = "joint model (upper bound)")
rel_alone <- NULL
if (file.exists(r10_file)) {
  R10 <- read_csv(r10_file, show_col_types = FALSE)
  rel_lower <- R10 %>%
    filter(quantity %in% c("B1_rel_one_platform_lower", "B1_rel_mean_of_platforms_lower")) %>%
    transmute(biomarker = factor(biomarker, biomarkers), subset = factor("All people", levels(rel$subset)),
              quantity = factor(if_else(quantity == "B1_rel_one_platform_lower", "one platform", "mean of all platforms"),
                                levels(rel$quantity)),
              estimate, ci_low = lo95, ci_high = hi95, bound = "lower bound")
  rel <- bind_rows(rel, rel_lower)
  rel_alone <- R10 %>% filter(quantity == "B2_rel_single_at_adni_sst") %>%
    group_by(biomarker) %>% summarise(lo = min(estimate), hi = max(estimate), .groups = "drop") %>%
    mutate(biomarker = factor(biomarker, biomarkers), subset = factor("All people", levels(rel$subset)))
}
p10b <- ggplot(rel, aes(biomarker, estimate, colour = quantity, shape = bound)) +
  geom_hline(yintercept = c(0.5, 0.7, 0.9), colour = GRID, linewidth = 0.3) +
  { if (!is.null(rel_alone)) geom_linerange(data = rel_alone, aes(x = biomarker, ymin = lo, ymax = hi),
                                            inherit.aes = FALSE, colour = "#b8b5ad", linewidth = 4, alpha = 0.55) } +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.55), linewidth = 0.6) +
  geom_point(position = position_dodge(width = 0.55), size = 2.2) +
  scale_shape_manual(values = c("joint model (upper bound)" = 16, "lower bound" = 1), name = NULL) +
  facet_wrap(~ subset, nrow = 1) +
  scale_x_discrete(labels = bm_label) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_colour_manual(values = c("one platform" = ORANGE, "mean of all platforms" = BLUE)) +
  labs(x = NULL, y = "Reliability of a three-draw slope",
       title = "How much of an individual's observed rate of change is real",
       subtitle = paste0("Share of slope variance that is true between-person change; ",
                         "95% percentile intervals from ", format(n_boot, big.mark = ","),
                         " person-level bootstrap resamples; open symbols are lower bounds, grey bands each measure fitted alone")) +
  theme(panel.grid.major.y = element_blank())
sv(p10b, "F10b_slope_reliability", w = 9.5, h = 3.8)

# ---- F11 validity ----------------------------------------------------------------------------
say("F11 validity")
crit_lab <- c(centiloid_slope = "Amyloid\nPET", csf_ptau_abeta42_slope = "CSF\npTau/Ab42",
              hcv_adj_slope = "Hippo-\ncampal\nvolume", adas13_slope = "ADAS-\nCog13", mpacc_digit_slope = "mPACC",
              mmse_slope = "MMSE")
# signs are oriented in 08b when T5 is written, so nothing to flip here
T5 <- rd("T5_validity_vs_criterion_change")
val <- T5 %>% filter(estimate == "pooled") %>% select(biomarker, criterion, n, r) %>%
  left_join(T7 %>% filter(subset == "all", str_starts(quantity, "validity_")) %>%
              mutate(criterion = str_remove(quantity, "^validity_")) %>% select(biomarker, criterion, ci_low, ci_high),
            by = c("biomarker", "criterion")) %>%
  filter(criterion %in% names(crit_lab)) %>%
  mutate(lo = pmin(ci_low, ci_high), hi = pmax(ci_low, ci_high),
         biomarker = factor(biomarker, biomarkers), criterion = factor(crit_lab[criterion], crit_lab))
p11 <- ggplot(val, aes(criterion, r)) +
  geom_hline(yintercept = 0, colour = AXIS, linewidth = 0.4) +
  geom_linerange(aes(ymin = lo, ymax = hi), colour = BLUE, linewidth = 0.6) +
  geom_point(colour = BLUE, size = 2.2, shape = 16) +
  geom_text(aes(y = hi, label = sprintf("%.2f", r)), vjust = -0.6, size = 2.5, colour = INK2) +
  facet_wrap(~ biomarker, nrow = 1, labeller = lab_bm) +
  scale_y_continuous(limits = c(-0.3, 0.65), breaks = seq(-0.2, 0.6, 0.2)) +
  labs(x = NULL, y = "Spearman r with concurrent criterion change",
       title = "Validity of the pooled plasma rate of change against four criteria",
       subtitle = "Positive means the plasma slope tracks worsening on the criterion; 95% bootstrap intervals") +
  theme(axis.text.x = element_text(size = 6.5, lineheight = 0.85), panel.spacing.x = unit(18, "pt"))
sv(p11, "F11_validity_pooled_slope", w = 12, h = 4)

# ---- F12 years needed --------------------------------------------------------------------------
p12 <- NULL
if (file.exists(file.path(P2_DIR, "T14_mdc_fold_change.csv"))) {
  say("F12 years to exceed MDC")
  CAP <- 30
  T14 <- rd("T14_mdc_fold_change") %>%
    mutate(arm = fix_arm(arm), biomarker = factor(biomarker, biomarkers),
           yrs = pmin(years_to_exceed_mdc_amy_pos, CAP), capped = years_to_exceed_mdc_amy_pos > CAP,
           lab = if_else(capped, paste0("> ", CAP), sprintf("%.0f", years_to_exceed_mdc_amy_pos)))
  p12 <- ggplot(T14, aes(arm, yrs)) +
    geom_col(width = 0.6, fill = BLUE, alpha = if_else(T14$capped, 0.35, 1)) +
    geom_text(aes(label = lab), vjust = -0.5, size = 2.7, colour = INK2) +
    facet_wrap(~ biomarker, scales = "free_x", nrow = 1, labeller = lab_bm) +
    scale_y_continuous(limits = c(0, CAP * 1.12), breaks = seq(0, CAP, 10), expand = expansion(mult = c(0, 0))) +
    labs(x = NULL, y = "Years between two draws",
         title = "How long before a change is larger than noise",
         subtitle = "Years of average amyloid-positive change needed before two measurements on the same platform differ by more than the mdc95; pale bars exceed 30 years") +
    theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.line.x = element_blank(), axis.ticks.x = element_blank())
  sv(p12, "F12_years_to_exceed_mdc", w = 11, h = 3.8)
}

# =============================================================================
# Composite figures for the manuscript (main: figures/main/, supplementary:
# figures/supplementary/). Panels are the working plots above with titles
# removed; the legend text goes in the manuscript, not on the figure.
# =============================================================================
say("Composite figures")
MAIN_DIR <- file.path(FIG_DIR, "main"); dir.create(MAIN_DIR, showWarnings = FALSE)
SUPP_DIR <- file.path(FIG_DIR, "supplementary"); dir.create(SUPP_DIR, showWarnings = FALSE)
suppressPackageStartupMessages(library(patchwork))
strip_title <- function(p) p + labs(title = NULL, subtitle = NULL)
save_fig <- function(p, dir, name, w, h) {
  dev <- if (CAIRO_OK) grDevices::cairo_pdf else "pdf"
  ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = dev, bg = "white")
  ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white")
}
tag_theme <- theme(plot.tag = element_text(size = 11, face = "bold", colour = INK),
                   plot.tag.position = c(0, 1))

# --- Figure 1A: design schematic (which platform measured which biomarker, and on how many people)
design <- long %>% group_by(biomarker, arm) %>%
  summarise(n_people = n_distinct(RID), n_obs = n(), .groups = "drop") %>%
  mutate(lab = sprintf("%d participants\n%d samples", n_people, n_obs),
         biomarker = factor(biomarker, biomarkers), arm = factor(arm, rev(arm_levels)))
p1a <- ggplot(design, aes(biomarker, arm)) +
  geom_tile(fill = SEQ[2], colour = "white", linewidth = 1.5, width = 0.95, height = 0.9) +
  geom_text(aes(label = lab), size = 2.4, colour = INK, lineheight = 0.9) +
  scale_x_discrete(labels = bm_label, position = "top") +
  labs(x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text = element_text(colour = INK, size = 8.5), axis.text.x.top = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8, colour = INK2))
fig1 <- (p1a / strip_title(p7[["pTau217"]])) + plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = "A") & tag_theme
say("  Figure 1"); save_fig(fig1, MAIN_DIR, "Figure1_design_and_trajectories", w = 10, h = 7.4)

# --- Figure 2: level components and slope reliability
fig2 <- (strip_title(p10a) + guides(fill = guide_legend(ncol = 1)) +
           theme(legend.position = "bottom", legend.key.height = unit(8, "pt"), legend.text = element_text(size = 7.5)) |
         strip_title(p10b) + theme(axis.text.x = element_text(angle = 35, hjust = 1))) +
  plot_layout(widths = c(1, 1.9)) + plot_annotation(tag_levels = "A") & tag_theme
say("  Figure 2"); save_fig(fig2, MAIN_DIR, "Figure2_level_vs_change", w = 13, h = 5.2)

# --- Figure 3: slopes with the MDC band, and years to exceed it
if (!is.null(p12)) {
  fig3 <- (strip_title(p8) / strip_title(p12)) + plot_layout(heights = c(1.25, 1)) +
    plot_annotation(tag_levels = "A") & tag_theme
  say("  Figure 3"); save_fig(fig3, MAIN_DIR, "Figure3_individual_change_vs_noise", w = 11, h = 8)
}

# --- Figure 4: sampling design and trial sample size
if (file.exists(file.path(P2_DIR, "T16_sampling_schedules.csv")) && file.exists(file.path(P2_DIR, "T15_trial_sample_size.csv"))) {
  T2c <- rd("T2_change_variance_components"); T16 <- rd("T16_sampling_schedules"); T15 <- rd("T15_trial_sample_size")
  curve_df <- T2c %>% select(biomarker, var_person_slope, var_residual) %>%
    crossing(sst = seq(0.5, 40, by = 0.25)) %>%
    mutate(rel = var_person_slope / (var_person_slope + var_residual / sst), biomarker = factor(biomarker, biomarkers))
  curve_lab <- curve_df %>% filter(sst == max(sst))
  sched_pts <- T16 %>% filter(biomarker == "pTau217",
                              schedule %in% c("2 draws, 2 y apart", "3 draws evenly over 4 y", "5 draws annually over 4 y",
                                              "ADNI observed (median SST)", "4 draws: replicate at each end, 4 y", "2 draws, 6 y apart")) %>%
    mutate(short = recode(schedule, "2 draws, 2 y apart" = "2 samples, 2 y", "3 draws evenly over 4 y" = "2 or 3 samples, 4 y",
                          "5 draws annually over 4 y" = "5 annual samples, 4 y", "ADNI observed (median SST)" = "ADNI design",
                          "4 draws: replicate at each end, 4 y" = "2 + 2 replicate samples, 4 y", "2 draws, 6 y apart" = "2 samples, 6 y"))
  p4a <- ggplot(curve_df, aes(sst, rel, group = biomarker)) +
    geom_hline(yintercept = c(0.7, 0.9), colour = GRID, linewidth = 0.4, linetype = "22") +
    geom_line(aes(colour = biomarker == "pTau217", linewidth = biomarker == "pTau217")) +
    ggrepel::geom_text_repel(data = curve_lab, aes(label = bm_label[as.character(biomarker)]),
                             hjust = 0, nudge_x = 1.2, direction = "y", size = 2.6, colour = INK2,
                             min.segment.length = Inf, box.padding = 0.08, seed = 1) +
    geom_point(data = sched_pts, aes(sst, rel_one_platform), colour = ORANGE, size = 2.2, inherit.aes = FALSE) +
    # the shortest schedule sits near the axis, so its label is nudged up rather than down (a label
    # nudged below zero was dropped by the axis limit and the point was left unlabelled)
    ggrepel::geom_text_repel(data = sched_pts, aes(sst, rel_one_platform, label = short), colour = ORANGE, size = 2.4,
                             inherit.aes = FALSE, nudge_x = 9,
                             nudge_y = ifelse(sched_pts$short == "2 samples, 2 y", 0.14, -0.12), direction = "y", hjust = 0,
                             segment.colour = "#d9c7c2", segment.size = 0.3, min.segment.length = 0, seed = 1,
                             box.padding = 0.25) +
    scale_colour_manual(values = c(`TRUE` = BLUE, `FALSE` = "#b8b5ad"), guide = "none") +
    scale_linewidth_manual(values = c(`TRUE` = 1, `FALSE` = 0.5), guide = "none") +
    scale_x_continuous(limits = c(0, 46), breaks = seq(0, 40, 10)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Spread of sampling times, SST (sum of squared centred years)", y = "Reliability of a one-platform slope")
  # Revision, November 2026: trial sizes on amyloid-positive variance components with bootstrap
  # intervals (R10, 08b19), the mean of platforms under both bounds. Without R10, as before.
  if (file.exists(r10_file)) {
    trial <- read_csv(r10_file, show_col_types = FALSE) %>%
      filter(quantity == "B3_n_per_arm_amy_pos", schedule == "24 months, 3 draws (0, 12, 24 m)",
             effect == 0.25, biomarker != "Ab42_40") %>%
      mutate(arm_lab = case_when(arm == "mean of platforms upper" ~ "Mean, upper bound",
                                 arm == "mean of platforms lower" ~ "Mean, lower bound",
                                 TRUE ~ as.character(fix_arm(arm))),
             is_mean = str_starts(arm, "mean"), biomarker = factor(biomarker, biomarkers),
             arm_lab = factor(arm_lab, c(arm_levels, "Mean, upper bound", "Mean, lower bound")),
             n_per_arm_slope = estimate)
  } else {
    trial <- T15 %>% filter(schedule == "24 months, 3 draws (0, 12, 24 m)", effect == 0.25, biomarker != "Ab42_40") %>%
      mutate(arm_lab = if_else(arm == "mean of platforms", "Mean of platforms", as.character(fix_arm(arm))),
             is_mean = arm == "mean of platforms", biomarker = factor(biomarker, biomarkers),
             arm_lab = factor(arm_lab, c(arm_levels, "Mean of platforms")), lo95 = NA_real_, hi95 = NA_real_)
  }
  p4b <- ggplot(trial, aes(arm_lab, n_per_arm_slope, fill = is_mean)) +
    geom_col(width = 0.65) +
    geom_linerange(aes(ymin = lo95, ymax = hi95), colour = INK2, linewidth = 0.45, na.rm = TRUE) +
    ggrepel::geom_text_repel(aes(label = format(round(n_per_arm_slope), big.mark = ",")),
                             direction = "y", nudge_y = 180, min.segment.length = Inf,
                             box.padding = 0.06, point.padding = 0, size = 2.2, colour = INK2, seed = 1) +
    facet_grid(~ biomarker, scales = "free_x", space = "free_x", labeller = lab_bm) +
    scale_fill_manual(values = c(`FALSE` = BLUE, `TRUE` = ORANGE), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.14)), labels = function(x) format(x, big.mark = ",")) +
    labs(x = NULL, y = "Participants per arm (95% interval)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6.5), axis.line.x = element_blank(), axis.ticks.x = element_blank())
  # Panel C, the within-participant test of the same axes as panel A. It is built here from the
  # table 08b11_paper2_fourth_sample.R writes, not from a saved plot object: a ggplot built in
  # another session does not survive being combined with a theme in this one, which is what the
  # "annotation$theme is not a valid theme" failure was. Optional: without the table, Figure 4
  # is the two panels it has always been.
  t9_file <- file.path(OUT_ROOT, "combined", "T9_fourth_sample.csv")
  t9 <- if (file.exists(t9_file)) tryCatch(read_csv(t9_file, show_col_types = FALSE),
                                           error = function(e) NULL) else NULL
  if (is.null(t9) || !nrow(t9)) {
    fig4 <- (p4a | p4b) + plot_layout(widths = c(1, 1.5)) +
      plot_annotation(tag_levels = "A") & tag_theme
    say("  Figure 4, two panels"); save_fig(fig4, MAIN_DIR, "Figure4_design_and_trial", w = 13, h = 4.6)
  } else {
    # Panel C must be readable against panel A, so it takes A's axes exactly. It also must
    # not be mistaken for A: the curve in C is the A4 assay's relation, which sits above the
    # ADNI pTau217 curve in A because that assay is more reliable at the same design, so the
    # axis label and an in-panel annotation both say which cohort the panel describes.
    # The predicted curve, anchored on the A4 placebo arm's reliability at its own design.
    a_r <- 0.601; a_sst <- 13.533
    odds_per_sst <- (a_r / (1 - a_r)) / a_sst
    c_curve <- tibble(sst = seq(3, 46, length.out = 300)) %>%
      mutate(rel = (odds_per_sst * sst) / (1 + odds_per_sst * sst))
    c_pts <- t9 %>%
      select(sample, design, sst_median, Predicted = reliability_predicted,
             Fitted = reliability_model) %>%
      pivot_longer(c(Predicted, Fitted), names_to = "source", values_to = "rel") %>%
      mutate(source = factor(source, c("Predicted", "Fitted")))
    c_labs <- t9 %>% group_by(design) %>%
      summarise(sst_median = median(sst_median),
                rel = max(c(reliability_predicted, reliability_model), na.rm = TRUE),
                letter = substr(first(design), 1, 1), .groups = "drop")
    p4c <- ggplot(c_pts, aes(sst_median, rel)) +
      geom_line(data = c_curve, aes(sst, rel), colour = MUTED, linetype = 3, inherit.aes = FALSE) +
      geom_point(aes(colour = sample, shape = source), size = 2.2) +
      geom_text(data = c_labs, aes(label = letter), vjust = -1.2, size = 3, colour = INK,
                inherit.aes = TRUE, show.legend = FALSE) +
      scale_colour_manual(values = c(BLUE, ORANGE)) +
      scale_shape_manual(values = c(Predicted = 16, Fitted = 1)) +
      scale_x_continuous(limits = c(0, 46), breaks = seq(0, 40, 10)) +
      scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
      annotate("text", x = 1, y = 0.98, hjust = 0, vjust = 1, size = 2.4, colour = INK2,
               label = "A4, one assay, which is more\nreliable than ADNI pTau217\nat the same design") +
      labs(x = "Spread of sampling times, SST\n(sum of squared centred years)",
           y = "Reliability of a slope, A4") +
      # The two guides sit side by side rather than stacked, so the panel keeps its height in
      # a single row of three. Nothing here sets a key size: under ggplot2 4.0 a legend key
      # smaller than the glyph it holds draws into a zero-dimension viewport, which is what
      # made this figure fail. Titles are removed with legend.title rather than name = NULL.
      theme(legend.position = "bottom", legend.box = "horizontal",
            legend.title = element_blank(), legend.text = element_text(size = 6.5))
    # One row: A the prediction, B its consequence for trials, C the prediction tested
    # within participants. A and C are narrower than B, which carries four facets.
    fig4 <- (p4a | p4b | p4c) + plot_layout(widths = c(1.2, 1, 1.2)) +
      plot_annotation(tag_levels = "A") & tag_theme
    say("  Figure 4, three panels")
    # patchwork's print method catches every error raised while a figure is drawn and reports
    # it as "the Plots window may be too small", which hides the real cause. Building the
    # gtable first shows the actual message, and saving the gtable rather than the patchwork
    # object bypasses that print method altogether.
    fig4g <- tryCatch(patchwork::patchworkGrob(fig4), error = function(e) {
      say("  ! Figure 4 could not be assembled: %s", conditionMessage(e)); NULL })
    if (!is.null(fig4g)) {
      ok <- tryCatch({ save_fig(fig4g, MAIN_DIR, "Figure4_design_and_trial", w = 15, h = 5); TRUE },
                     error = function(e) { say("  ! Figure 4 could not be drawn: %s",
                                               conditionMessage(e)); FALSE })
    } else ok <- FALSE
    if (!ok) {
      # Fall back to the two panels, so the rest of the script finishes. Under a different
      # name: a three-panel Figure4_design_and_trial assembled elsewhere must not be
      # silently overwritten by the fallback.
      say("  falling back to the two-panel Figure 4, saved as Figure4_two_panel")
      fig4b <- (p4a | p4b) + plot_layout(widths = c(1, 1.5)) +
        plot_annotation(tag_levels = "A") & tag_theme
      save_fig(fig4b, MAIN_DIR, "Figure4_two_panel", w = 13, h = 4.6)
      save_fig(p4c, MAIN_DIR, "Figure4C_fourth_sample", w = 6, h = 4.6)
    }
  }
}

# --- Supplementary figures
supp_bm <- setdiff(biomarkers, "pTau217")
# Numbered by first mention in the manuscript: S1 level agreement, S2 slope
# agreement, S3 residual diagnostics (from 08b2), S4 validity, S5 trajectories
# of the other four biomarkers (cited in the Figure 1 legend).
if (exists("p9b")) save_fig(strip_title(p9b), SUPP_DIR, "FigureS1_level_agreement", w = 12, h = 3.2)
save_fig(strip_title(p9a), SUPP_DIR, "FigureS2_slope_agreement", w = 12, h = 3.2)
f4 <- file.path(P2_DIR, "F4_residual_diagnostics.pdf")
if (file.exists(f4)) file.copy(f4, file.path(SUPP_DIR, "FigureS3_residual_diagnostics.pdf"), overwrite = TRUE)
save_fig(strip_title(p11), SUPP_DIR, "FigureS4_validity", w = 12, h = 3.8)
figS5 <- wrap_plots(lapply(supp_bm, function(bm) strip_title(p7[[bm]]) + theme(legend.position = "none")), ncol = 1) +
  plot_annotation(tag_levels = "A") & tag_theme
save_fig(figS5, SUPP_DIR, "FigureS5_trajectories_other_biomarkers", w = 8, h = 12)

LEG_MAIN <- c(
  sprintf("Figure 1. Study design and individual pTau217 trajectories. (A) Platforms and biomarkers: each cell gives the number of participants and plasma samples measured on that platform for that biomarker. All %d participants provided three samples over a median of five years, and every platform assayed an aliquot of the same blood draw. C2N reports pTau217 both as a concentration and as the percentage of phosphorylated to non-phosphorylated tau (%%pTau217). (B) Every participant's pTau217 series on each platform. Values are log-transformed and standardised within measure so that higher is more abnormal; thin lines are individual participants, heavy lines the ordinary least-squares trend by amyloid PET status at the first sample (above 20 Centiloids). The vertical axis is truncated at ±4 standardised units. The corresponding panels for the other four biomarkers are in Supplementary Figure S5.", n_people),
  paste0("Figure 2. Platforms agree about level but not about change. (A) Sources of variance in a single measurement, from a crossed random-effects model per biomarker: the participant component (blue) is between-participant variance in true level; every other component is error for one measurement on one platform (sample is the participant × visit effect shared by all platforms assaying that aliquot). Segments under 5% are unlabelled. (B) Reliability of an individual's three-sample rate of change, defined as the share of the variance of an observed slope that is true between-participant variance in rate of change, for a slope from one platform and for the mean of all platforms measuring that biomarker, in all participants and by amyloid PET status at the first sample. Vertical lines are 95% percentile intervals from ", format(n_boot, big.mark = ","),
         " bootstrap resamples of participants. Filled symbols are the joint-model values, which treat variation shared across platforms at a blood draw as true change and are upper bounds; open symbols (pTau217, all participants) are the lower bounds from the model with a participant-by-sample effect; grey bands span the reliability of each measure fitted alone, from the lowest to the highest platform. Because every platform measured the same blood sample, averaging removes only the assay-specific part of the occasion-level variation."),
  sprintf("Figure 3. Individual change against measurement noise. (A) Per-participant three-sample slopes by platform (standardised log units per year; N = %d participants), coloured by amyloid PET status at the first sample, with boxes showing the median and interquartile range. The shaded band is the 95%% minimal detectable change for a three-sample slope on that platform, the range of slopes expected in 95%% of people whose true rate of change is zero; slopes inside the band cannot be distinguished from no change. (B) Years between two plasma samples that a participant changing at the mean rate observed in amyloid-positive participants would need before the difference between the two measurements exceeds the 95%% minimal detectable change for that platform; pale bars exceed 30 years.", n_core),
  "Figure 4. Consequences for sampling design and for trials. (A) Reliability of a one-platform slope as a function of the spread of sampling times (SST, the sum of squared deviations of a participant's sampling times from their mean), which is the only feature of a schedule that reliability depends on. Curves are the five biomarkers (pTau217 in blue); points mark selected schedules for pTau217, including the observed ADNI design. Dotted lines mark reliabilities of 0.7 and 0.9. (B) Participants per arm required by a two-arm trial using the plasma slope as its outcome to detect a 25% slowing of the mean amyloid-positive rate of change over 24 months with three samples (two-sided α 0.05, 80% power), for each platform and for the mean of all platforms measuring that biomarker (orange) under the upper bound, in which averaging removes all occasion-level variation, and, for pTau217, the lower bound, in which the part shared across platforms at a blood draw remains. Vertical lines are 95% percentile intervals from 1,000 bootstrap resamples of participants. The calculation uses the variance components and mean rate of amyloid-positive ADNI participants and describes an ADNI-like population; Aβ42/40 is omitted because no platform gives a usable slope outcome. (C) The same axes tested within participants. In the A4 subgroup who gave a fourth plasma sample, each participant's rate of change was recomputed from four subsets of their own samples: A the standard three samples over 4.8 years, B the same three samples over 7.3 years, C the same interval with the middle sample placed later, D all four samples. Filled points are the reliability the spread of sampling times predicts, open points the reliability refitted within these participants. The dotted line is the predicted relation for the A4 assay, which lies above the ADNI pTau217 curve in panel A because that assay is more reliable at the same design; the axes are otherwise identical to panel A.",
  "Supplementary Figure S1. Agreement between platforms about an individual's level: intraclass correlation between platforms measuring the same biomarker at the same visit, on the log scale after removing each platform's mean. Diagonal cells are blank.",
  "Supplementary Figure S2. Agreement between platforms about an individual's rate of change: Pearson correlation between per-participant three-sample slopes for each pair of platforms measuring the same biomarker. Diagonal cells are blank.",
  "Supplementary Figure S3. Residual diagnostics of the multi-platform change model (09_bootstrap_intervals.R, section 1, F4).",
  "Supplementary Figure S5. Individual trajectories on every platform for pTau181, GFAP, NfL and Aβ42/40, drawn as in Figure 1B.",
  "Supplementary Figure S4. Validity of the pooled rate of change: Spearman correlation between each participant's pooled plasma slope and the concurrent rate of change in amyloid PET (Centiloid), CSF pTau181/Aβ42, adjusted hippocampal volume, ADAS-Cog13, mPACC and MMSE, with 95% bootstrap intervals. Signs are set so that a positive correlation means the plasma slope tracks worsening on the reference measure.")
writeLines(c("# Manuscript figure legends", "", LEG_MAIN, "", ABBREV), file.path(MAIN_DIR, "manuscript_figure_legends.md"))

writeLines(c("# Figure legends (as printed under the journal versions in figures/final/)", "",
             sprintf("Generated by 08b7_paper2_figures.R, pipeline version %s.", PIPELINE_VERSION), "",
             unlist(lapply(names(LEG), function(k) c(paste0("## ", k), LEG[[k]], ""))), ABBREV),
           file.path(FINAL_DIR, "figure_legends.md"))

say("\nWorking figures: %s; journal versions: %s; manuscript composites: %s and %s", FIG_DIR, FINAL_DIR, MAIN_DIR, SUPP_DIR)
