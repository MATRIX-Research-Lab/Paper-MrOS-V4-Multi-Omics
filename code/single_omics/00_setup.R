set.seed(123)
options(stringsAsFactors = FALSE)

# Resolve the repository from the working directory or the location of the
# entrypoint. Do not depend on the clone directory name or on an ignored
# .Rproj file; both are common sources of clean-clone failures.
find_repository_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  ofile <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  starts <- c(getwd(), if (length(file_arg)) dirname(file_arg[[1]]) else NULL,
              if (!is.null(ofile)) dirname(ofile) else NULL)
  starts <- unique(starts[!is.na(starts) & nzchar(starts)])
  is_repo <- function(d) {
    file.exists(file.path(d, "code", "single_omics", "00_setup.R")) &&
      file.exists(file.path(d, "code", "single_omics", "helpers.R"))
  }
  for (s in starts) {
    d <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(20L)) {
      if (is_repo(d)) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Could not locate the repository root. Start inside the checkout or run an entrypoint by its path.")
}

project_root <- find_repository_root()
setwd(project_root)

# There are two sets of layer-specific model outputs and they are not
# interchangeable.
#
#   model_outputs        The primary analysis. These are what 01_-05_ read, and
#                        what every published primary figure and table is
#                        computed from.
#   model_outputs_refit  A later re-fit of the same three layers. 06_ reads
#                        these, and only these, to assemble the matched cohort:
#                        they carry the provenance manifests that the matched
#                        stage checks, which the primary objects predate.
#
# The two agree on the analytic samples -- same participants, same feature
# blocks -- so the matched cohort is the same either way. Keeping them apart is
# what stops a render from silently mixing them.
paths <- list(
  data = file.path(project_root, "data"),
  raw = file.path(project_root, "data", "raw_data"),
  processed = file.path(project_root, "data", "processed_data"),
  model_outputs = file.path(project_root, "data", "model_outputs", "single_omics"),
  model_outputs_refit = file.path(project_root, "data", "model_outputs", "single_omics_rerun"),
  results = file.path(project_root, "results"),
  manifests = file.path(project_root, "results", "manifests"),
  figures = file.path(project_root, "results", "figures"),
  tables = file.path(project_root, "results", "tables")
)

dir.create(paths$results, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$manifests, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)
dir.create(paths$tables, recursive = TRUE, showWarnings = FALSE)

# Frozen analysis contract shared by preprocessing, primary-model fitting,
# preflight validation, and the matched-cohort refit. Keeping these definitions
# in one place prevents a silent change in feature boundaries or expected sample
# counts in only one entry point.
analysis_covariates <- c(
  "race", "edu", "ol_health", "bmi", "mstat", "smoke", "diab", "hbp",
  "cancer", "tmm_score", "gds", "pase", "total_meds", "age_v4"
)
analysis_non_feature_columns <- c(
  "ID", "status", "ATTEND_V5", "D3CR_MB_COHORT", "batch", "site", "grip",
  "chair", "gait400m", "gait6m", "cr_cmm", "abx", "b4thd", "b4fnd",
  "b4lsd", "faprev4", "hqdrfefl", "hqdtfefl", "hqptfefl", "fr_chs4",
  "age_v5", "type"
)
analysis_contract <- list(
  microbiome = c(n = 765L, deaths = 439L, features = 147L),
  metabolomics = c(n = 447L, deaths = 216L, features = 1574L),
  proteomics = c(n = 446L, deaths = 215L, features = 7596L)
)
analysis_matched_metadata_columns <- unique(c(
  analysis_non_feature_columns, analysis_covariates
))

# `tune_grid()` can parallelize either over resamples or over the full
# resample-by-hyperparameter grid. The ordinary/default mode is intentionally
# conservative; the 64-CPU cluster driver opts into the finer-grained mode.
analysis_tune_parallel_over <- function(default = "resamples") {
  value <- tolower(Sys.getenv("MROS_TUNE_PARALLEL_OVER", unset = default))
  if (!value %in% c("resamples", "everything")) {
    stop("MROS_TUNE_PARALLEL_OVER must be 'resamples' or 'everything'")
  }
  value
}

analysis_tune_control <- function(save_pred = FALSE, verbose = FALSE) {
  tune::control_grid(
    save_pred = save_pred,
    verbose = verbose,
    parallel_over = analysis_tune_parallel_over()
  )
}

# Establish the parallel backend for a fitting stage, and return its label.
#
# This must be forked workers (future::multicore), never future::multisession.
# The limma screen is applied through step_limma_select(), a custom recipes step
# whose S3 methods (prep./bake.) are defined in the global environment by
# feature_selection.R. multisession workers are separate R processes: future
# exports the globals it can detect statically, but S3 methods are resolved by
# dispatch at run time and never reach the worker, so every tuning call dies
# with
#
#   no applicable method for 'prep' applied to an object of class
#   "step_limma_select"
#
# and tune::collect_metrics() then aborts with "All models failed". Forked
# workers inherit the parent's global environment, so dispatch works.
#
# Failing loudly here is deliberate: falling back to multisession would waste a
# multi-hour allocation before the error surfaced, and falling back to
# sequential would silently turn a 20-worker job into a single-core one.
analysis_parallel_plan <- function(n_workers) {
  n_workers <- suppressWarnings(as.integer(n_workers))
  if (is.na(n_workers) || n_workers < 1L) {
    stop("n_workers must be a positive integer")
  }
  if (n_workers == 1L) {
    future::plan(future::sequential)
    return("sequential")
  }
  supported <- requireNamespace("parallelly", quietly = TRUE) &&
    isTRUE(parallelly::supportsMulticore())
  if (!supported) {
    stop(
      "Parallel fitting requires forked workers (future::multicore), which this ",
      "platform does not support. The limma recipe step cannot dispatch inside ",
      "future::multisession workers. Re-run with n_workers = 1, or run on Linux ",
      "outside RStudio."
    )
  }
  future::plan(future::multicore, workers = n_workers)
  "multicore"
}

# Attach, in the parent, the packages tune would otherwise attach inside every
# worker.
#
# Before evaluating anything, tune calls required_pkgs() on the workflow and
# runs library() for each result in each worker: limma for the screening step,
# and the engine package for the model. A forked worker inherits whatever is
# already attached in the parent, so those library() calls then return without
# touching the filesystem at all.
#
# That matters because the conda environment lives on a parallel filesystem.
# Sixty forked R processes resolving packages out of a library of roughly a
# quarter of a million small files loads the metadata server heavily, and a
# single failed read surfaces as
#
#   there is no package called 'limma'
#
# which reads like a broken environment but is a transient I/O failure: the
# same code, on the same environment, succeeds in the other layers running
# concurrently. Preloading takes the filesystem out of that path.
analysis_preload_worker_packages <- function(
    packages = c("limma", "glmnet", "ranger", "xgboost")) {
  attached <- character()
  for (p in packages) {
    if (requireNamespace(p, quietly = TRUE)) {
      suppressPackageStartupMessages(
        library(p, character.only = TRUE, quietly = TRUE)
      )
      attached <- c(attached, p)
    }
  }
  invisible(attached)
}

# Upper bound for xgboost's `min_n` on a given analysis set.
#
# parsnip maps min_n to xgboost's min_child_weight, which is a minimum sum of
# *Hessian weights* in a child node, not a count of observations. For binary
# logistic the Hessian per observation is p(1-p) <= 0.25, so a node holding k
# observations carries at most 0.25k. A split needs both children to clear the
# threshold, and row subsampling shrinks what is available, so the largest
# usable value is about
#
#   n_analysis * min(sample_size) * 0.25 / 2
#
# Above that, no split can satisfy the constraint: every tree collapses to a
# single leaf, the model predicts one constant probability for everybody, and
# its AUC is exactly 0.500.
#
# The fixed range c(10, 120) this replaces sits above the bound for analysis
# sets of a few hundred men. It produced exactly that failure -- constant
# XGBoost models with min_n 48 (metabolomics) and 107 (proteomics), zero split
# nodes, AUC 0.500 -- while the larger microbiome layer escaped at min_n 33.
analysis_xgb_min_n_range <- function(n_analysis, min_sample_prop = 0.5) {
  n_analysis <- suppressWarnings(as.integer(n_analysis))
  if (is.na(n_analysis) || n_analysis < 1L) {
    stop("n_analysis must be a positive integer")
  }
  upper <- as.integer(floor(n_analysis * min_sample_prop * 0.25 / 2))
  c(2L, max(4L, upper))
}

assert_modeling_api <- function() {
  if (!"parallel_over" %in% names(formals(tune::control_grid))) {
    stop("The installed tune package is too old: control_grid() lacks parallel_over")
  }
  required_functions <- list(
    "tune::collect_predictions" = tune::collect_predictions,
    "dials::grid_space_filling" = dials::grid_space_filling,
    "recipes::add_step" = recipes::add_step
  )
  if (any(!vapply(required_functions, is.function, logical(1)))) {
    stop("The installed tidymodels components lack required workflow APIs")
  }
  invisible(TRUE)
}

# Figures produced by code/cross_omics/02_make_figures.R.
# run_all.R never generates these; it only confirms they are still present, so a
# render cannot quietly drop them from results/ the way it once did.
# Figures this workflow does not produce. Each entry point owns its own outputs;
# run_all.R only confirms these survived the render, so results/ cannot quietly
# lose them the way it once did.
cross_omics_artifacts <- c(
  "Fig5_micro-metab_correlation.pdf.pdf",
  "Supplementary_Fig3_micro_prot.pdf",
  "Supplementary_Fig4_meta_prot.pdf"
)

latex_artifacts <- c(
  "Supplementary_Fig5_primary_modeling_workflow.pdf",
  "Supplementary_Fig6_matched_cohort_workflow.pdf"
)

# Producer label per external artifact, used when recording their status.
external_artifact_source <- c(
  setNames(rep("code/cross_omics/02_make_figures.R", length(cross_omics_artifacts)), cross_omics_artifacts),
  setNames(rep("code/latex_workflow/", length(latex_artifacts)), latex_artifacts)
)

expected_artifacts <- c(
  "Fig1_flowchart.pdf",
  "Fig2_roc.pdf",
  "Fig3_bubble_omics.pdf",
  "Fig4_bubble_ensamble.pdf",
  "Fig5_micro-metab_correlation.pdf.pdf",
  "Supplementary_Fig1_omics_importance.pdf",
  "Supplementary_Fig2_ensamble_importance.pdf",
  "Supplementary_Fig3_micro_prot.pdf",
  "Supplementary_Fig4_meta_prot.pdf",
  "Table1.pdf",
  "Table1.docx",
  "Table2_AUC_CI.pdf",
  "SuppTable_model_comparisons.pdf",
  "SuppTable_calibration.pdf",
  "SuppTable_age_proxy.pdf",
  "SuppTable_complete_case_vs_rest.pdf",
  "SuppTable_characteristics_microbiome.pdf",
  "SuppTable_characteristics_metabolomics.pdf",
  "SuppTable_characteristics_proteomics.pdf",
  "SuppTable_en_xgb_overlap.pdf",
  "SuppTable_ensemble_overlap.pdf",
  "Supplementary_Fig5_primary_modeling_workflow.pdf",
  "Supplementary_Fig6_matched_cohort_workflow.pdf",
  # matched-cohort refit, generated by 07_
  "SuppTable_matched_auc.pdf",
  "SuppTable_matched_comparisons.pdf"
)

# `tidyverse` and `tidymodels` are the meta-packages that omics_pred_pipeline.R
# and run_omics_pred.R attach with library(). They must be declared, not just
# their components: assert_packages() checks exactly this list, so a missing
# meta-package otherwise slips past every gate and fails at library() once the
# job is already running.
required_packages <- c(
  "dplyr", "tidyr", "readr", "purrr", "stringr", "tibble",
  "haven", "readxl", "openxlsx", "flextable", "officer", "broom",
  "ggplot2", "scales", "ggbeeswarm", "ggpubr", "cowplot", "ggtext",
  "phyloseq", "microbiome", "ANCOMBC", "SomaDataIO",
  "tidyverse", "tidymodels", "glmnet", "ranger", "xgboost", "limma",
  "kernelshap", "shapviz", "pROC", "future"
)

primary_packages <- c(
  "dplyr", "tibble", "purrr", "tidyverse", "tidymodels", "glmnet", "ranger",
  "xgboost", "limma", "kernelshap", "shapviz", "broom", "future"
)

matched_packages <- c("dplyr", "tibble", "tidymodels", "glmnet", "ranger",
                      "xgboost", "limma", "future", "pROC")

render_packages <- setdiff(
  required_packages,
  c("haven", "readxl", "phyloseq", "microbiome", "ANCOMBC", "SomaDataIO")
)

issues <- data.frame(
  step = character(),
  category = character(),
  severity = character(),
  detail = character(),
  stringsAsFactors = FALSE
)

artifact_status <- data.frame(
  artifact = expected_artifacts,
  generated = FALSE,
  status = "not_checked",
  detail = "",
  stringsAsFactors = FALSE
)

record_issue <- function(step, category, detail, severity = "warning") {
  issue <- data.frame(
    step = step,
    category = category,
    severity = severity,
    detail = detail,
    stringsAsFactors = FALSE
  )
  assign("issues", rbind(get("issues", envir = globalenv()), issue), envir = globalenv())
}

record_artifact <- function(artifact, generated, status, detail = "") {
  tab <- get("artifact_status", envir = globalenv())
  idx <- match(artifact, tab$artifact)
  if (is.na(idx)) {
    tab <- rbind(
      tab,
      data.frame(
        artifact = artifact,
        generated = generated,
        status = status,
        detail = detail,
        stringsAsFactors = FALSE
      )
    )
  } else {
    tab$generated[idx] <- generated
    tab$status[idx] <- status
    tab$detail[idx] <- detail
  }
  assign("artifact_status", tab, envir = globalenv())
}

artifact_output_file <- function(artifact) {
  if (grepl("^Table", artifact)) {
    return(file.path(paths$tables, artifact))
  }
  file.path(paths$figures, artifact)
}

source(file.path(project_root, "code", "single_omics", "helpers.R"))
source(file.path(project_root, "code", "single_omics", "analysis_provenance.R"))
source(file.path(project_root, "code", "single_omics", "feature_selection.R"))
