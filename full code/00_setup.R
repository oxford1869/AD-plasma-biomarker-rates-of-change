# =============================================================================================
# 00_setup.R
# Shared settings for the whole pipeline: paths, constants, helpers, the figure theme, the analysis units and the driver helper.
#
# Sourced at the start of every section of every script. Four parts:
#   1. settings (the former 00_config.R): paths, cohort switch, raw file names, constants, helper functions
#   2. figure theme (the former 00_theme.R)
#   3. define_units(long) (the former 08c0_paper3_units.R) and load_change_objects() (the shared setup of
#      the former 08b17), both called by the sections that need them
#   4. run_steps(), used by run_all.R
# Paths: BASE_DIR (raw ADNI files), OUT_ROOT (where pipeline_out/, a4learn_out/ and combined/ are written) and
# A4_DIR (the A4/LEARN release); each has a default and an environment variable (ADNI_JM_DIR, ADNI_OUT_DIR, A4_DIR).
# =============================================================================================

# ---- 1. settings ----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(lubridate)
  library(purrr)
  library(tibble)
})

# ---- Paths -------------------------------------------------------------------
# BASE_DIR : where the raw ADNI CSVs live. Read only; nothing is written here.
# OUT_ROOT : where pipeline_out/ is created. Everything the pipeline writes
#            goes under OUT_ROOT/pipeline_out/.
# Both can be overridden with environment variables ADNI_JM_DIR and ADNI_OUT_DIR.
BASE_DIR <- Sys.getenv(
  "ADNI_JM_DIR",
  unset = "/Users/yiyang/Library/CloudStorage/OneDrive-Nexus365/Desktop/ADNI/JM data"
)
OUT_ROOT <- Sys.getenv(
  "ADNI_OUT_DIR",
  unset = "/Users/yiyang/Desktop/3 fast"
)
# COHORT selects which data the downstream scripts (07b to 08c) read:
#   "adni"    (default) pipeline_out/ built by 01 and 03 from the ADNI FNIH files
#   "a4learn" a4learn_out/ built by 02 and 03 from the A4/LEARN release
# Set with Sys.setenv(COHORT = "a4learn") before sourcing a script.
COHORT <- tolower(Sys.getenv("COHORT", unset = "adni"))
if (!COHORT %in% c("adni", "a4learn")) stop("COHORT must be 'adni' or 'a4learn'")
A4_DIR <- Sys.getenv("A4_DIR", unset = "/Users/yiyang/Desktop/3 fast/A4&learn")   # read only
OUT_DIR    <- file.path(OUT_ROOT, if (COHORT == "adni") "pipeline_out" else "a4learn_out")
CHECK_DIR  <- file.path(OUT_DIR, "checks")
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(CHECK_DIR, showWarnings = FALSE, recursive = TRUE)

RAW_FILES <- list(
  plasma   = file.path(BASE_DIR, "FNIHBC_BLOOD_BIOMARKER_TRAJECTORIES.csv"),
  merge    = file.path(BASE_DIR, "ADNIMERGE_DATA.csv"),
  amyloid  = file.path(BASE_DIR, "UCBERKELEY_AMY_6MM.CSV"),
  csf      = file.path(BASE_DIR, "UPENNBIOMK_ROCHE_ELECSYS.csv"),
  lab      = file.path(BASE_DIR, "LABDATA.CSV"),
  demog    = file.path(BASE_DIR, "PTDEMOG.csv"),
  mmse     = file.path(BASE_DIR, "MMSE.CSV"),
  vitals   = file.path(BASE_DIR, "VITALS.csv")     # optional: height and weight for BMI
)
A4_FILES <- list(
  ptau217  = file.path(A4_DIR, "External Data", "biomarker_pTau217.csv"),
  roche    = file.path(A4_DIR, "External Data", "biomarker_Plasma_Roche_Results.csv"),
  amyloid  = file.path(A4_DIR, "External Data", "imaging_SUVR_amyloid.csv"),
  tau      = file.path(A4_DIR, "External Data", "imaging_SUVR_tau.csv"),
  lab      = file.path(A4_DIR, "External Data", "clrm_lab.csv"),
  subjinfo = file.path(A4_DIR, "Derived Data", "SUBJINFO.csv"),
  sv       = file.path(A4_DIR, "Derived Data", "SV.csv"),
  ds       = file.path(A4_DIR, "Derived Data", "DS.csv"),
  pacc     = file.path(A4_DIR, "Derived Data", "PACC.csv"),
  ptdemog  = file.path(A4_DIR, "Raw Data", "ptdemog.csv"),
  vitals   = file.path(A4_DIR, "Raw Data", "vitals.csv"),
  cdr      = file.path(A4_DIR, "Raw Data", "cdr.csv"),
  mmse     = file.path(A4_DIR, "Raw Data", "mmse.csv"),
  visits   = file.path(A4_DIR, "Documents", "Data Dictionaries", "visits_datadic.csv")
)
# A4/LEARN give every date as days since T0 (A4: randomisation; LEARN: baseline
# registry exam). The pipeline needs Date columns, so days are mapped onto a
# nominal calendar anchored at A4_EPOCH; only differences between dates are used.
A4_EPOCH <- as.Date("2000-01-01")
a4_date  <- function(days_t0) A4_EPOCH + as.integer(round(as.numeric(days_t0)))
# Centiloids for the Avid florbetapir pipeline (whole cerebellum reference):
# Navitsky et al. 2018, Alzheimer's & Dementia 14:1565. CL = 183 x SUVR - 177.
a4_centiloid <- function(suvr) 183 * suvr - 177
# An integer id for the analysis scripts (they expect numeric RID): the digits of the BID
a4_rid <- function(bid) as.integer(stringr::str_remove(bid, "^B"))

# ---- Constants ---------------------------------------------------------------
PIPELINE_VERSION <- "0.34.0"
YEAR_DAYS        <- 365.25
# Plasma dates within this many days of each other, for the same person, are
# treated as one visit. On the real file every person's distinct dates match
# their distinct SAMPLEID prefixes, and one person has two genuine samples 27
# days apart, so the window is kept narrow: it only absorbs a lab entering a
# date a few days off, never two separate draws.
VISIT_GAP_DAYS   <- 7
CENTILOID_CUTS   <- c(cl20 = 20, cl37 = 37)
# Windows used by 06_master.R when attaching the nearest assessment to a
# plasma visit. Each paper can override these when it builds its own frame.
MATCH_WINDOW_DAYS <- c(criteria = 180, modifiers = 365, dx = 180)

# ---- Helpers -----------------------------------------------------------------

# Dates in the ADNI exports are day/month/year. Try that first, then fall back.
parse_date_dmy <- function(x) {
  as.Date(lubridate::parse_date_time(x, orders = c("dmy", "ymd", "mdy"), quiet = TRUE))
}
# ADNIMERGE and the UPenn/Berkeley files are usually year-month-day. Try that
# first; the plasma file is the only one that is day/month/year.
parse_date_ymd <- function(x) {
  as.Date(lubridate::parse_date_time(x, orders = c("ymd", "dmy", "mdy"), quiet = TRUE))
}

# Screening and follow-up codes collapse to baseline for joins to ADNIMERGE.
recode_viscode <- function(v) {
  v <- as.character(v)
  dplyr::case_when(
    v %in% c("sc", "f", "scmri") ~ "bl",
    TRUE ~ v
  )
}

# Platform label from the PLASMA_BIOMARKER prefix.
platform_from_prefix <- function(pb) {
  dplyr::case_when(
    str_starts(pb, "C2N_")     ~ "C2N",
    str_starts(pb, "Fuji_")    ~ "Fujirebio",
    str_starts(pb, "Roche_")   ~ "Roche",
    str_starts(pb, "AlzPath_") ~ "ALZpath",
    str_starts(pb, "Janssen_") ~ "Janssen",
    str_starts(pb, "QX_")      ~ "Quanterix",
    str_starts(pb, "Lilly_")   ~ "Lilly MSD",
    TRUE ~ NA_character_
  )
}

# ---- Measure lookup ----------------------------------------------------------
# One row per PLASMA_BIOMARKER value in the raw file.
#   biomarker : the construct this measure belongs to
#   role      : "measure" enters the analyses; "component" is an ingredient
#               of a ratio and is kept for deriving ratios and for checks
#   component : which ingredient, for ratio derivation
#   direction : +1 if higher means more abnormal, -1 if lower does
MEASURE_LOOKUP <- tibble::tribble(
  ~PLASMA_BIOMARKER,            ~biomarker, ~role,       ~component, ~direction,
  "C2N_plasma_ptau217",         "pTau217",  "measure",   NA,         1,
  "C2N_plasma_ptau217_ratio",   "pTau217",  "measure",   NA,         1,
  "C2N_plasma_nptau217",        "pTau217",  "component", "npTau217", 1,
  "Fuji_plasma_ptau217",        "pTau217",  "measure",   NA,         1,
  "AlzPath_plasma_ptau217",     "pTau217",  "measure",   NA,         1,
  "Janssen_plasma_ptau217",     "pTau217",  "measure",   NA,         1,
  "Lilly_plasma_ptau217",       "pTau217",  "measure",   NA,         1,   # A4/LEARN only
  "Roche_plasma_ptau181",       "pTau181",  "measure",   NA,         1,
  "QX_plasma_ptau181",          "pTau181",  "measure",   NA,         1,
  "Roche_plasma_GFAP",          "GFAP",     "measure",   NA,         1,
  "QX_plasma_GFAP",             "GFAP",     "measure",   NA,         1,
  "Roche_plasma_NfL",           "NfL",      "measure",   NA,         1,
  "QX_plasma_NfL",              "NfL",      "measure",   NA,         1,
  "C2N_plasma_Abeta42_Abeta40", "Ab42_40",  "measure",   NA,        -1,
  "C2N_plasma_Abeta40",         "Ab42_40",  "component", "Ab40",     1,
  "C2N_plasma_Abeta42",         "Ab42_40",  "component", "Ab42",     1,
  "Fuji_plasma_Ab40",           "Ab42_40",  "component", "Ab40",     1,
  "Fuji_plasma_Ab42",           "Ab42_40",  "component", "Ab42",     1,
  "Roche_plasma_Ab40",          "Ab42_40",  "component", "Ab40",     1,
  "Roche_plasma_Ab42",          "Ab42_40",  "component", "Ab42",     1,
  "QX_plasma_Ab40",             "Ab42_40",  "component", "Ab40",     1,
  "QX_plasma_Ab42",             "Ab42_40",  "component", "Ab42",     1,
  # Ratios derived in 01_plasma_long.R (not in the raw file); listed so that
  # every measure in the master has a biomarker and a direction.
  "Fuji_plasma_Ab42_Ab40",      "Ab42_40",  "measure",   NA,        -1,
  "Roche_plasma_Ab42_Ab40",     "Ab42_40",  "measure",   NA,        -1,
  "QX_plasma_Ab42_Ab40",        "Ab42_40",  "measure",   NA,        -1
)

# Platforms whose Ab42/40 ratio is not supplied by the lab and must be derived.
DERIVE_RATIO_FOR <- c("Fujirebio", "Roche", "Quanterix")

# Read any ADNI CSV with every column as character, so nothing is silently
# coerced. Column names are upper-cased on the way in because ADNI exports
# are inconsistent (VISCODE2 vs viscode2, SCANDATE vs scandate).
read_adni <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  x <- readr::read_csv(path, col_types = readr::cols(.default = readr::col_character()),
                       show_col_types = FALSE, progress = FALSE)
  names(x) <- toupper(names(x))
  x
}

# Return the first of several candidate column names that exists in df,
# or stop with a message listing what is there. Used because ADNI files
# name the same thing differently (EXAMDATE / VISDATE / SCANDATE).
pick_col <- function(df, candidates, required = TRUE, label = candidates[1]) {
  hit <- candidates[toupper(candidates) %in% names(df)]
  if (length(hit) == 0) {
    if (required) stop(sprintf("None of [%s] found for %s. Columns are: %s",
                               paste(candidates, collapse = ", "), label,
                               paste(names(df), collapse = ", ")))
    return(NA_character_)
  }
  toupper(hit[1])
}

# Numeric parse that keeps a comparator ("<8", ">1700") as a flag. Returns a
# tibble with value, comparator, and the raw string.
parse_censored <- function(x) {
  x <- as.character(x)
  comp <- stringr::str_trim(stringr::str_extract(x, "^\\s*[<>]"))
  val  <- ifelse(stringr::str_detect(dplyr::coalesce(x, ""), "[0-9]"),
                 suppressWarnings(readr::parse_number(x)), NA_real_)
  tibble::tibble(value = val, comparator = comp)
}

# For each row of `target` (RID + a date), find the plasma visit of the same
# person closest in time. Returns target with visit_no, visit_date, gap_days
# (target date minus plasma date; negative = before the plasma sample) and
# plasma_time of the target relative to the first plasma sample.
match_to_plasma <- function(target, date_col, plasma_visits) {
  stopifnot(all(c("RID", "visit_no", "visit_date", "first_plasma") %in% names(plasma_visits)))
  target$.row <- seq_len(nrow(target))
  target$.d   <- target[[date_col]]
  matched <- target %>%
    dplyr::select(.row, RID, .d) %>%
    dplyr::left_join(plasma_visits, by = "RID", relationship = "many-to-many") %>%
    dplyr::mutate(gap_days = as.integer(.d - visit_date)) %>%
    dplyr::group_by(.row) %>%
    dplyr::slice_min(abs(gap_days), n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(plasma_time = as.numeric(.d - first_plasma) / YEAR_DAYS) %>%
    dplyr::select(.row, visit_no, visit_date, gap_days, plasma_time)
  target %>%
    dplyr::left_join(matched, by = ".row") %>%
    dplyr::select(-.row, -.d)
}

# eGFR from serum creatinine (mg/dL), age (years) and sex ("Female"/"Male").
# CKD-EPI 2021 (race-free) is the primary equation; CKD-EPI 2009 without the
# race coefficient is kept for comparison with earlier work.
egfr_ckdepi_2021 <- function(scr, age, sex) {
  female <- sex == "Female"
  k      <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.241, -0.302)
  142 * pmin(scr / k, 1)^alpha * pmax(scr / k, 1)^-1.200 * 0.9938^age * ifelse(female, 1.012, 1)
}
egfr_ckdepi_2009 <- function(scr, age, sex) {
  female <- sex == "Female"
  k      <- ifelse(female, 0.7, 0.9)
  alpha  <- ifelse(female, -0.329, -0.411)
  141 * pmin(scr / k, 1)^alpha * pmax(scr / k, 1)^-1.209 * 0.993^age * ifelse(female, 1.018, 1)
}

# Paper-specific settings used by the 07 frame scripts
EGFR_SENS_GAP_YEARS <- 1        # sensitivity analysis: creatinine within this of visit 1
CRITERION_SLOPE_PAD_DAYS <- 180 # assessments this far outside the plasma span still count
                                # towards a person's criterion rate of change

# Per-person least-squares slope of y on t within groups. Returns one row per
# group with n, intercept, slope, slope_se, first/last t and the span. Groups
# with fewer than 2 points get NA slope.
slope_by <- function(df, y, t, ...) {
  df %>%
    dplyr::filter(!is.na(.data[[y]]), !is.na(.data[[t]])) %>%
    dplyr::group_by(...) %>%
    dplyr::summarise(
      n_obs   = dplyr::n(),
      t_first = min(.data[[t]]), t_last = max(.data[[t]]),
      span    = t_last - t_first,
      y_first = .data[[y]][which.min(.data[[t]])],
      y_last  = .data[[y]][which.max(.data[[t]])],
      slope   = if (dplyr::n() >= 2 && span > 0) stats::cov(.data[[t]], .data[[y]]) / stats::var(.data[[t]]) else NA_real_,
      intercept = if (!is.na(slope)) mean(.data[[y]]) - slope * mean(.data[[t]]) else NA_real_,
      slope_se  = if (dplyr::n() >= 3 && span > 0) {
        r <- .data[[y]] - (intercept + slope * .data[[t]])
        sqrt(sum(r^2) / (dplyr::n() - 2) / sum((.data[[t]] - mean(.data[[t]]))^2))
      } else NA_real_,
      .groups = "drop"
    )
}

# Log transform for positive concentrations and ratios; NA otherwise.
safe_log <- function(x) ifelse(!is.na(x) & x > 0, log(x), NA_real_)

# ---- Small utilities ---------------------------------------------------------
save_rds <- function(obj, name) {
  attr(obj, "pipeline_version") <- PIPELINE_VERSION
  attr(obj, "built_at") <- Sys.time()
  saveRDS(obj, file.path(OUT_DIR, paste0(name, ".rds")))
  invisible(obj)
}

say <- function(...) cat(sprintf(...), "\n")


# ---- 2. figure theme ------------------------------------------------------------------------

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  message("00_theme.R: ggplot2 is not installed, figure styling skipped")
} else {
  suppressPackageStartupMessages(library(ggplot2))

  INK    <- "#1c1c1c"; INK2 <- "#5a5751"; MUTED <- "#8d8a83"
  GRID   <- "#e6e3dc"; AXIS <- "#c9c6be"; BAND <- "#efece5"
  BLUE   <- "#1F6F8B"; ORANGE <- "#B0413E"; GREY <- "#a39f97"
  SEQ    <- c("#f4f1ea", "#d5dde2", "#a9bcc8", "#7d9bad", "#527b92", "#2f5b74", "#1b3b50")

  # categorical palette: the two series colours first, then four that stay
  # distinguishable from them and from each other in greyscale
  PAL_CAT <- c("#1F6F8B", "#B0413E", "#6E8B3D", "#8C6A9E", "#C98A2E", "#4F4F4F", "#7FA8B8", "#D08C8A")

  pal_amy    <- c("amyloid negative" = BLUE, "amyloid positive" = ORANGE, "no PET" = GREY)
  # cohorts and units appear side by side in the combined figures
  pal_cohort <- c("ADNI" = BLUE, "adni" = BLUE, "A4/LEARN" = ORANGE, "a4learn" = ORANGE)

  theme_paper <- function(base_size = 9) {
    theme_minimal(base_size = base_size) +
      theme(text = element_text(colour = INK),
            plot.title = element_text(size = base_size + 1, face = "plain", colour = INK, hjust = 0,
                                      margin = margin(b = 8)),
            plot.subtitle = element_text(size = base_size - 0.5, colour = INK2, hjust = 0, margin = margin(b = 8)),
            plot.title.position = "plot",
            plot.caption = element_text(size = base_size - 2, colour = MUTED, hjust = 0),
            plot.caption.position = "plot",
            axis.title = element_text(size = base_size - 0.5, colour = INK2),
            axis.title.y = element_text(margin = margin(r = 6)),
            axis.title.x = element_text(margin = margin(t = 6)),
            axis.text = element_text(size = base_size - 1, colour = MUTED),
            axis.line.x = element_line(colour = AXIS, linewidth = 0.3),
            axis.ticks.x = element_line(colour = AXIS, linewidth = 0.3), axis.ticks.length.x = unit(2, "pt"),
            panel.grid.major.y = element_line(colour = GRID, linewidth = 0.3),
            panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
            panel.spacing.x = unit(14, "pt"), panel.spacing.y = unit(10, "pt"),
            strip.text = element_text(size = base_size, colour = INK, face = "bold", hjust = 0,
                                      margin = margin(b = 4)),
            strip.background = element_blank(), strip.placement = "outside",
            legend.position = "bottom", legend.title = element_blank(),
            legend.text = element_text(size = base_size - 1, colour = INK2),
            legend.key.height = unit(8, "pt"), legend.key.width = unit(14, "pt"),
            legend.margin = margin(t = 2),
            plot.margin = margin(10, 12, 8, 10),
            plot.background = element_rect(fill = "white", colour = NA))
  }

  scale_colour_paper <- function(...) ggplot2::scale_colour_manual(values = unname(PAL_CAT), ...)
  scale_color_paper  <- scale_colour_paper
  scale_fill_paper   <- function(...) ggplot2::scale_fill_manual(values = unname(PAL_CAT), ...)

  theme_set(theme_paper())
  options(ggplot2.discrete.colour   = unname(PAL_CAT),
          ggplot2.discrete.fill     = unname(PAL_CAT),
          ggplot2.continuous.colour = function(...) ggplot2::scale_colour_gradientn(colours = SEQ, ...),
          ggplot2.continuous.fill   = function(...) ggplot2::scale_fill_gradientn(colours = SEQ, ...))
}

# =============================================================================
# Figure and table labels (added 0.31.5)
#
# Every script until now printed the internal name of a summary on the axis of a
# figure, so a reader saw hlme_3class or threshold_q80_3group. Those names are
# fine in a CSV, where they are keys, and wrong in a figure, where they are the
# only text a reader has. The maps below turn each internal name into the wording
# used in the manuscript, and nothing else changes: the data keep their internal
# names throughout, and only the label shown is translated.
#
# Objects: pretty_summary(), pretty_world(), pretty_family(), SUMMARY_LEVELS
# =============================================================================

SUMMARY_LABELS <- c(
  level_bl                    = "Level at the first sample",
  level_mean                  = "Mean level",
  level_last                  = "Level at the last sample",
  ols_slope                   = "OLS rate of change",
  blup_slope                  = "BLUP rate of change",
  gamm_slope                  = "GAMM rate of change",
  pooled_slope                = "Pooled rate of change, five measures",
  lcmm_slope                  = "Latent class model rate of change",
  hlme_2class                 = "Latent class assignment, 2 classes",
  hlme_3class                 = "Latent class assignment, 3 classes",
  hlme_4class                 = "Latent class assignment, 4 classes",
  hlme_bic                    = "Latent class assignment, classes chosen by BIC",
  hlme_2class_slope           = "Latent class model rate of change, 2 classes",
  hlme_3class_slope           = "Latent class model rate of change, 3 classes",
  hlme_4class_slope           = "Latent class model rate of change, 4 classes",
  hlme_2class_post            = "Posterior class probabilities, 2 classes",
  hlme_3class_post            = "Posterior class probabilities, 3 classes",
  hlme_4class_post            = "Posterior class probabilities, 4 classes",
  kmeans_2class               = "k-means cluster assignment, 2 clusters",
  kmeans_3class               = "k-means cluster assignment, 3 clusters",
  kmeans_4class               = "k-means cluster assignment, 4 clusters",
  threshold_q60_3group        = "Threshold group, 60th percentile",
  threshold_q67_3group        = "Threshold group, 67th percentile",
  threshold_q75_3group        = "Threshold group, 75th percentile",
  threshold_q80_3group        = "Threshold group, 80th percentile",
  threshold_quantile_3group   = "Threshold group, sample quantile",
  threshold_quantile_baseline = "Above the threshold at the first sample",
  threshold_mixture_3group    = "Threshold group, mixture crossing point",
  threshold_published_3group  = "Threshold group, published cut-off",
  threshold_3group            = "Threshold group, Youden cut-off",
  threshold_converter         = "Threshold crossing, Youden cut-off"
)

# the order the manuscript uses: level, then continuous rate, then the three
# categorical families. reorder() inside a figure still overrides this.
SUMMARY_LEVELS <- names(SUMMARY_LABELS)

# the threshold summaries are a rule crossed with a form, so they are built rather than
# listed: any new percentile or new form is labelled correctly without editing the map
THRESH_RULE <- c(q60 = "60th percentile", q67 = "67th percentile", q75 = "75th percentile",
                 q80 = "80th percentile", quantile = "sample quantile",
                 mixture = "mixture crossing point", published = "published cut-off",
                 youden = "Youden cut-off")

pretty_summary <- function(x) {
  x <- as.character(x)
  lab <- unname(SUMMARY_LABELS[x])        # the explicit map wins wherever it has an entry
  miss <- is.na(lab)
  if (any(miss)) {
    y <- x[miss]
    m <- regmatches(y, regexec("^threshold_([a-z0-9]+)_(3group|baseline|converter)$", y))
    built <- vapply(seq_along(y), function(i) {
      g <- m[[i]]
      if (length(g) != 3) return(NA_character_)
      rule <- unname(THRESH_RULE[g[2]])
      if (is.na(rule)) return(NA_character_)
      switch(g[3],
             `3group`  = sprintf("Threshold group, %s", rule),
             baseline  = sprintf("Above the %s at the first sample", rule),
             converter = sprintf("Threshold crossing, %s", rule))
    }, character(1))
    # anything still unmatched is made readable rather than left as an internal key
    plain <- gsub("_", " ", sub("_post$", " posterior probabilities",
                                sub("_slope$", " rate of change", y)))
    y <- ifelse(is.na(built), plain, built)
    lab[miss] <- y
  }
  lab
}

WORLD_LABELS <- c(
  `2class`      = "Two classes in rate of change",
  `3class`      = "Three classes in rate of change",
  continuous    = "No classes, Gaussian rate",
  skewed        = "No classes, skewed rate",
  censored      = "No classes, detection floor",
  skewcens      = "No classes, skewed rate and floor",
  level2        = "No classes, two-component level",
  level2cens    = "No classes, two-component level and floor",
  levelskew     = "No classes, skewed level",
  levelskewcens = "No classes, skewed level and floor"
)
WORLD_LEVELS <- names(WORLD_LABELS)

pretty_world <- function(x) {
  x <- as.character(x)
  lab <- unname(WORLD_LABELS[x]); lab[is.na(lab)] <- x[is.na(lab)]; lab
}

pretty_family <- function(x) {
  x <- as.character(x)
  m <- c(level = "Level", `continuous summary` = "Continuous rate",
         `latent class` = "Latent class", `distance-based cluster` = "Distance-based cluster",
         `threshold crossing` = "Threshold group", other = "Other")
  lab <- unname(m[x]); lab[is.na(lab)] <- x[is.na(lab)]; lab
}

# ---- the other specification-grid dimensions (added 0.31.6) ------------------
# The specification curve prints one row of levels per dimension, and until now only the
# summary row was translated, so a reader still saw floor_na, baseline_negative_cl37 and
# age_sex_apoe_educ. pretty_level() translates a value given the dimension it belongs to.

LEVEL_LABELS <- list(
  floor = c(as_delivered         = "Retained as delivered",
            exclude_always_floor = "Participants always at the floor excluded",
            floor_na             = "Floor values set to missing"),
  outcome = c(centiloid_slope  = "Concurrent Centiloid rate of change",
              centiloid_change = "Centiloid change, last minus first",
              centiloid_last   = "Last Centiloid, adjusted for the first"),
  sample = c(all_core               = "All participants",
             baseline_negative_cl20 = "Amyloid negative at 20 Centiloid",
             baseline_negative_cl37 = "Amyloid negative at 37 Centiloid",
             baseline_cn            = "Cognitively normal at the first sample"),
  covariates = c(none               = "None",
                 age_sex            = "Age and sex",
                 age_sex_apoe       = "Age, sex, APOE e4",
                 age_sex_apoe_educ  = "Age, sex, APOE e4, education",
                 age_sex_amyloid_bl = "Age, sex, baseline Centiloid"),
  class_order = c(`-`       = "Not applicable",
                  intercept = "Ordered by level",
                  slope     = "Ordered by rate of change"),
  measure_label = c(C2N                    = "C2N pTau217",
                    `Pooled (5 platforms)` = "Pooled across the five measures",
                    `pooled (all arms)`    = "Pooled across the five measures")
)

pretty_level <- function(dimension, value) {
  value <- as.character(value); dimension <- as.character(dimension)
  out <- value
  for (d in unique(dimension)) {
    i <- dimension == d
    if (identical(d, "summary")) { out[i] <- pretty_summary(value[i]); next }
    m <- LEVEL_LABELS[[d]]
    if (is.null(m)) next
    lab <- unname(m[value[i]])
    out[i] <- ifelse(is.na(lab), value[i], lab)
  }
  out
}

# class_order appears on its own as a shape or colour legend, not only inside the grid
pretty_class_order <- function(x) pretty_level("class_order", x)


# ---- 3. analysis units and shared change-model objects --------------------------------------
# define_units(): the analysis units by cohort and the floor rule. ADNI: the five pTau217 platforms on the
# same people. A4/LEARN: one platform, so the units are subsamples (A4 placebo, LEARN, A4 solanezumab,
# placebo plus LEARN, all three arms). Requires the frame `long` (frame_paper3_long) already read.
define_units <- function(long) {
  # Objects defined: MEASURES, m_label, unit_data(), ref_data(), always_floor, apply_floor(), UNITS_SHARE_PEOPLE,
  #   CONS_UNITS, UNIT_NOUN, UNITS_TITLE, FLOOR_LABEL, REF_LABEL, PUBLISHED_CUTS, hlme_person_slope(); all placed in
  #   the global environment, as when this code was a script sourced after the frame was read.
  if (COHORT == "adni") {
    MEASURES <- c("C2N_plasma_ptau217", "C2N_plasma_ptau217_ratio", "Fuji_plasma_ptau217",
                  "AlzPath_plasma_ptau217", "Janssen_plasma_ptau217")
    m_label <- c(C2N_plasma_ptau217 = "C2N", C2N_plasma_ptau217_ratio = "C2N %pTau217", Fuji_plasma_ptau217 = "Fujirebio",
                 AlzPath_plasma_ptau217 = "ALZpath", Janssen_plasma_ptau217 = "Janssen", pooled = "Pooled (5 platforms)")
    p217 <- long %>% filter(measure %in% MEASURES)
    unit_data <- function(m) p217 %>% filter(measure == m)
    always_floor <- p217 %>% filter(measure == "C2N_plasma_ptau217") %>% group_by(RID) %>%
      summarise(af = all(lloq_flag), .groups = "drop") %>% filter(af) %>% pull(RID)
    UNITS_SHARE_PEOPLE <- TRUE
    # Sensitivity: C2N contributes twice, as an absolute concentration and as %pTau217 from the
    # same assay on the same blood, so the five measures are four assays. P3_DROP_MEASURES drops
    # named measures, e.g. P3_DROP_MEASURES=C2N_plasma_ptau217_ratio for the four-assay run.
    drop_m <- trimws(strsplit(Sys.getenv("P3_DROP_MEASURES", ""), ",")[[1]])
    drop_m <- drop_m[nzchar(drop_m)]
    if (length(drop_m)) {
      unknown <- setdiff(drop_m, MEASURES)
      if (length(unknown)) stop("P3_DROP_MEASURES names measures that are not units: ",
                                paste(unknown, collapse = ", "))
      MEASURES <- setdiff(MEASURES, drop_m)
      p217 <- p217 %>% filter(measure %in% MEASURES)
      say("P3_DROP_UNITS: dropped %s; %d units remain", paste(m_label[drop_m], collapse = ", "),
          length(MEASURES))
    }
    CONS_UNITS  <- MEASURES
    UNIT_NOUN   <- "platform"
    UNITS_TITLE <- sprintf("%s platforms",
                           c("one", "two", "three", "four", "five")[length(MEASURES)])
    FLOOR_LABEL <- "C2N floor"
  } else {
    MEASURES <- c("A4_placebo", "LEARN", "A4_solanezumab", "placebo_LEARN", "all_arms")
    m_label <- c(A4_placebo = "A4 placebo", LEARN = "LEARN", A4_solanezumab = "A4 solanezumab",
                 placebo_LEARN = "Placebo + LEARN (primary)", all_arms = "All three arms", pooled = "Pooled")
    p217 <- long %>% filter(measure == "Lilly_plasma_ptau217")
    unit_data <- function(m) switch(m,
      A4_placebo     = p217 %>% filter(arm == "Placebo"),
      LEARN          = p217 %>% filter(arm == "LEARN"),
      A4_solanezumab = p217 %>% filter(arm == "Solanezumab"),
      placebo_LEARN  = p217 %>% filter(arm %in% c("Placebo", "LEARN")),
      all_arms       = p217)
    always_floor <- p217 %>% group_by(RID) %>%
      summarise(af = all(lloq_flag), .groups = "drop") %>% filter(af) %>% pull(RID)
    UNITS_SHARE_PEOPLE <- FALSE
    # P3_DROP_UNITS does for the subsamples what P3_DROP_MEASURES does for the ADNI assays, so a
    # quick check run can be limited to a few of them. It has its own name because the two cohorts
    # have different unit names and the ADNI branch rejects a name it does not recognise.
    # Inert unless the variable is set.
    drop_m <- trimws(strsplit(Sys.getenv("P3_DROP_UNITS", ""), ",")[[1]])
    drop_m <- drop_m[nzchar(drop_m)]
    unknown_u <- setdiff(drop_m, MEASURES)
    if (length(unknown_u)) stop("P3_DROP_UNITS names subsamples that are not units: ",
                                paste(unknown_u, collapse = ", "))
    if (length(drop_m)) {
      MEASURES <- setdiff(MEASURES, drop_m)
      say("P3_DROP_UNITS: dropped %s; %d units remain",
          paste(m_label[drop_m], collapse = ", "), length(MEASURES))
    }
    CONS_UNITS  <- intersect(c("A4_placebo", "LEARN", "A4_solanezumab"), MEASURES)
    UNIT_NOUN   <- "subsample"
    UNITS_TITLE <- "A4 placebo, LEARN and A4 solanezumab (one platform)"
    FLOOR_LABEL <- "Lilly MSD floor"
  }
  # Key of the unit that a "one platform on all people" figure should show in the cohort
  UNIT_FOCUS <- if (COHORT == "adni") "C2N_plasma_ptau217" else "placebo_LEARN"

  # Reference sample a threshold cut-off is derived on, before being applied unchanged to
  # every unit. A laboratory fixes a cut-off once and then applies that same value to
  # whoever walks in; it does not recompute it inside each new group of patients. Deriving
  # it per unit is what made the A4/LEARN quantile cut-off land on a different absolute
  # concentration in every subsample, which is a property of the procedure and not of the
  # biomarker.
  #   ADNI: the units are five assays on different scales, so a cut-off cannot be shared
  #     between them. The reference is the whole ADNI sample on that assay, which is what
  #     unit_data() already returns, so nothing changes here.
  #   A4/LEARN: the units are subsamples of one assay on one scale, so they can and should
  #     share a cut-off. The reference is the whole cohort on that assay.
  ref_data <- if (COHORT == "adni") function(m) unit_data(m) else function(m) p217
  REF_LABEL <- if (COHORT == "adni") "whole ADNI sample, per assay" else "whole A4/LEARN sample, one assay"

  # Published cut-offs, in the assay's own units, from cohorts outside this analysis.
  # Included as a descriptive check ("does a value published elsewhere land near the one we
  # would have chosen?"), not as a row in the cross-unit comparison, because they exist for
  # only some of the assays and a comparison across units needs all of them.
  #   Deliberately excluded: the ADNI-specific cut-offs in Schindler et al. 2024 (4.06% for
  #   C2N %pTau217, 0.158 pg/mL for Fujirebio). They were derived in this cohort against
  #   this outcome, so using them would be the Youden cut-off with a citation attached.
  #   Also excluded: C2N absolute concentration (no cut-off calibrated against amyloid
  #   exists), Janssen (published values come from a different commercial implementation of
  #   the antibody pair, transferability unverified), and Lilly MSD in A4/LEARN (reported in
  #   U/mL, so the published pg/mL values do not apply).
  PUBLISHED_CUTS <- if (COHORT == "adni") list(
    AlzPath_plasma_ptau217   = list(value = 0.42, units = "pg/mL",
                                    source = "Ashton et al., JAMA Neurol 2024;81(3):255-263"),
    Fuji_plasma_ptau217      = list(value = 0.27, units = "pg/mL",
                                    source = "Palmqvist et al., Nat Med 2025;31(6):2036-2043"),
    C2N_plasma_ptau217_ratio = list(value = 3.26, units = "%",
                                    source = "Palmqvist et al., JAMA 2024;332(15):1245-1257")
  ) else list()

  apply_floor <- function(d, floor) {
    switch(floor,
      as_delivered         = d %>% mutate(y = z_log),
      exclude_always_floor = d %>% filter(!RID %in% always_floor) %>% mutate(y = z_log),
      floor_na             = d %>% mutate(y = if_else(lloq_flag, NA_real_, z_log)) %>% filter(!is.na(y)))
  }
  # Person-level predicted slope from a fitted hlme object: the posterior-weighted class slope
  # plus that person's random slope. This is the latent class model used as a continuous model,
  # so the classes play no part in the prediction. Returns tibble(RID, slope, has_random_effect),
  # or NULL if the fit does not expose what is needed.
  hlme_person_slope <- function(m, timevar = "plasma_time") {
    tryCatch({
      pp <- m$pprob
      if (is.null(pp) || !nrow(pp)) return(NULL)
      ids <- as.integer(pp[[1]])
      pcols <- grep("^prob", names(pp), value = TRUE)
      if (!length(pcols)) return(NULL)
      P <- as.matrix(pp[pcols]); P <- P / pmax(rowSums(P), 1e-12)
      b <- m$best
      if (is.null(names(b))) return(NULL)
      sl <- b[grepl(paste0("^", timevar), names(b))]
      if (!length(sl)) return(NULL)
      fixed <- if (length(sl) == ncol(P)) as.numeric(P %*% as.numeric(sl)) else rep(as.numeric(sl[[1]]), nrow(P))
      re <- m$predRE; rs <- rep(0, length(ids)); has_re <- FALSE
      if (!is.null(re) && timevar %in% names(re)) {
        v <- re[[timevar]][match(ids, as.integer(re[[1]]))]
        rs <- ifelse(is.na(v), 0, v); has_re <- TRUE
      }
      tibble(RID = ids, slope = fixed + rs, has_random_effect = has_re)
    }, error = function(e) NULL)
  }

  say("Paper 3 units (%s): %s; people always at the floor: %d", COHORT, paste(m_label[MEASURES], collapse = ", "), length(always_floor))
  list2env(mget(setdiff(ls(), c("long"))), envir = globalenv())
  invisible(NULL)
}

load_change_objects <- function() {
  # The objects shared by the change-model sections of 04 to 08 (the frames with factor columns,
  # the median spread of sampling times, the trial designs, the amyloid-positive mean rates, the
  # two model formulas and three helpers). Placed in the global environment. Reads nothing but
  # the frames; the tables each section needs are read in that section.
  YVAR       <- "z_log_aligned"
  biomarkers <- c("pTau217", "pTau181", "GFAP", "NfL", "Ab42_40")
  zsum2      <- (qnorm(0.975) + qnorm(0.80))^2
  ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))


  long <- readRDS(file.path(OUT_DIR, "frame_paper2_long.rds")) %>%
    mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform),
           RID_f = factor(RID), arm_f = factor(arm), batch_f = factor(batch),
           sample_f = factor(paste(RID, visit_no, sep = "_")))
  slopes <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes.rds")) %>%
    mutate(arm = case_when(measure == "C2N_plasma_ptau217_ratio" ~ "C2N_ratio", TRUE ~ platform))
  sw  <- readRDS(file.path(OUT_DIR, "frame_paper2_slopes_wide.rds"))

  # Median SST over people with at least two draws: the same definition as 08b, 08b6, 08b10.
  median_sst <- long %>% distinct(RID, plasma_time) %>% group_by(RID) %>%
    summarise(k = n(), sst = sum((plasma_time - mean(plasma_time))^2), .groups = "drop") %>%
    filter(k >= 2) %>% pull(sst) %>% median()
  say("Median SST %.3f (paper: 12.855)", median_sst)

  # Same trial designs as 08b6
  designs <- tribble(~schedule, ~times,
                     "18 months, 4 draws (0, 6, 12, 18 m)",     c(0, 0.5, 1, 1.5),
                     "24 months, 3 draws (0, 12, 24 m)",        c(0, 1, 2),
                     "24 months, 5 draws (0, 6, 12, 18, 24 m)", c(0, 0.5, 1, 1.5, 2)) %>%
    mutate(k = map_int(times, length), span = map_dbl(times, max),
           sst = map_dbl(times, ~ sum((.x - mean(.x))^2))) %>% select(-times)

  # Mean OLS slope in amyloid-positive people, as in 08b6
  amy_slope <- slopes %>% filter(n_obs >= 3, amyloid_bl_cl20 %in% TRUE) %>%
    group_by(biomarker, arm) %>%
    summarise(mean_slope_amy_pos_z = mean(slope, na.rm = TRUE), n_amy_pos = n(), .groups = "drop")
  amy_mean_bm <- amy_slope %>% group_by(biomarker) %>%
    summarise(mean_slope_amy_pos_z = mean(mean_slope_amy_pos_z, na.rm = TRUE), .groups = "drop")

  # Model formulas, as in 08b10
  f_reduced <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                                "(1 | RID_f:arm_f) + (1 | batch_f)"))
  f_sample  <- as.formula(paste(YVAR, "~ arm_f * plasma_time + (1 + plasma_time | RID_f) +",
                                "(1 | RID_f:arm_f) + (1 | sample_f) + (1 | batch_f)"))

  n_slope  <- function(v, delta) ceiling(2 * zsum2 * v / delta^2)
  fit_info <- function(fit) {
    msg <- fit@optinfo$conv$lme4$messages
    tibble(singular = isSingular(fit), converged = !any(!grepl("singular", msg)),
           messages = paste(msg, collapse = "; "))
  }
  vget <- function(fit, grp, var1 = NA) {
    vc <- as.data.frame(VarCorr(fit)) %>% filter(is.na(var2))
    i <- if (is.na(var1)) which(vc$grp == grp) else which(vc$grp == grp & !is.na(vc$var1) & vc$var1 == var1)
    if (length(i)) vc$vcov[i[1]] else NA_real_
  }
  list2env(mget(ls()), envir = globalenv())
  invisible(NULL)
}

# ---- the driver helper ------------------------------------------------------------------------
# Runs a list of steps, each list(cohort =, file =), announcing each with the time, and leaves the
# session on the ADNI default at the end. Used by run_all.R.
run_steps <- function(steps) {
  for (s in steps) {
    cat("\n\n############################################################\n## ", s$file,
        "  [COHORT = ", s$cohort, "]   ", format(Sys.time()),
        "\n############################################################\n", sep = "")
    Sys.setenv(COHORT = s$cohort)
    source(s$file, echo = FALSE)
  }
  Sys.setenv(COHORT = "adni")
}
