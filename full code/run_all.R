# =============================================================================================
# run_all.R
# The driver. Open in RStudio and press Source, or source() it by its full path. Runs the whole
# pipeline in order, for both cohorts, writing one console log to <OUT_ROOT>/pipeline_out/.
#
# STEPS   which files to run, by their two-digit prefix. The default is everything; to run a subset,
#         set STEPS in the console before sourcing (STEPS <- c("09", "10", "11")), or edit the line below.
#         Data build 01 to 03 takes a few minutes; 04 to 11 several hours, most of it the A4/LEARN
#         bootstrap in 07 and the power simulation in 08. The ADNI bootstraps in 09 save every draw
#         and resume if interrupted; the A4/LEARN one in 07 does not.
# N_BOOT  bootstrap draws for 07 (A4/LEARN) and 09 (ADNI); the paper uses 1,000, the default.
#         N_BOOT = 0 skips the resampling in 07 and in section 1 of 09; the later sections of 09,
#         and 10 and 11, need the saved draws, so a run without a bootstrap should set STEPS to
#         "01" to "08".
# Paths   set ADNI_JM_DIR, ADNI_OUT_DIR and A4_DIR with Sys.setenv() before sourcing, or edit the
#         defaults in 00_setup.R.
# =============================================================================================
if (!exists("STEPS")) STEPS <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11")   # set STEPS before sourcing to run a subset

# ---- move to the folder this file lives in, whatever the working directory is ----
this_file <- tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) NA_character_)
if (is.na(this_file)) { a <- grep("^--file=", commandArgs(), value = TRUE); if (length(a)) this_file <- normalizePath(sub("^--file=", "", a[1])) }
if (is.na(this_file) && requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
  this_file <- tryCatch(normalizePath(rstudioapi::getSourceEditorContext()$path), error = function(e) NA_character_)
if (!is.na(this_file) && nzchar(this_file)) setwd(dirname(this_file))
if (!file.exists("00_setup.R")) stop("Cannot find 00_setup.R. setwd() to the code folder and source run_all.R again.")
Sys.setenv(COHORT = "adni"); source("00_setup.R")

log_dir <- file.path(OUT_ROOT, "pipeline_out"); dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(log_dir, paste0("run_all_log_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt"))
while (sink.number() > 0) sink()
sink(log_file, split = TRUE)

files <- list.files(pattern = "^[0-9]{2}_.*\\.R$")
plan <- list(
  list(cohort = "adni",    file = "01_data_adni.R"),
  list(cohort = "a4learn", file = "02_data_a4learn.R"),
  list(cohort = "adni",    file = "03_analysis_frames.R"),
  list(cohort = "a4learn", file = "03_analysis_frames.R"),
  list(cohort = "adni",    file = "04_level_and_change_models.R"),
  list(cohort = "adni",    file = "05_minimal_detectable_change.R"),
  list(cohort = "adni",    file = "06_sampling_design_and_trials.R"),
  list(cohort = "adni",    file = "07_a4learn.R"),             # sets the cohort itself, section by section
  list(cohort = "adni",    file = "08_sensitivity_analyses.R"),
  list(cohort = "adni",    file = "09_bootstrap_intervals.R"),
  list(cohort = "adni",    file = "10_tables.R"),
  list(cohort = "adni",    file = "11_figures.R"))
plan <- Filter(function(s) substr(s$file, 1, 2) %in% STEPS, plan)
missing <- setdiff(unique(vapply(plan, `[[`, "", "file")), files)
if (length(missing)) stop("Not in this folder: ", paste(missing, collapse = ", "))

cat("Started at", format(Sys.time()), "\nLog:", log_file, "\n")
run_steps(plan)
cat("\nFinished at", format(Sys.time()), "\n",
    "  ADNI      pipeline_out/ and pipeline_out/paper2/\n",
    "  A4/LEARN  a4learn_out/ and a4learn_out/paper2/\n",
    "  Combined  combined/\n", sep = "")
while (sink.number() > 0) sink()
