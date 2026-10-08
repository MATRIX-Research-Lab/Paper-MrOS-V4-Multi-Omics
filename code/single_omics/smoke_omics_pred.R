# Smoke test for the primary fitting pipeline: 3 omics x all fs_method
# (none | limma | mrmr | sis), on a subsample with capped CV/grid and SHAP
# mocked. It exercises the code paths; it does not reproduce manuscript AUCs.
#
# Optional filters (character vectors):
#   SMOKE_LAYERS <- c("metabolomics")          # default: all three
#   SMOKE_FS     <- c("none", "limma")          # default: none/limma/mrmr/sis
#
# Terminal:
#   Rscript code/single_omics/smoke_omics_pred.R

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidymodels)
})

find_smoke_paths <- function() {
  pipeline_rel <- file.path("code", "single_omics", "omics_pred_pipeline.R")
  cands <- character()

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    ap <- tryCatch(rstudioapi::getActiveProject(), error = function(e) NULL)
    if (!is.null(ap) && nzchar(ap)) cands <- c(cands, ap)
  }

  this_file <- NULL
  ofile <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
  if (!is.null(ofile) && nzchar(ofile)) {
    this_file <- normalizePath(ofile, mustWork = FALSE)
  }
  cmd <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  if (length(fa) && nzchar(fa[[1]])) {
    this_file <- normalizePath(fa[[1]], mustWork = FALSE)
  }
  if (!is.null(this_file) && file.exists(this_file)) {
    cands <- c(cands, dirname(this_file), dirname(dirname(this_file)),
               dirname(dirname(dirname(this_file))),
               dirname(dirname(dirname(dirname(this_file)))))
  }

  d <- normalizePath(getwd(), mustWork = FALSE)
  for (i in 1:10) {
    cands <- c(cands, d)
    d2 <- dirname(d)
    if (identical(d2, d)) break
    d <- d2
  }

  cands <- unique(cands[!is.na(cands) & nzchar(cands)])

  for (root in cands) {
    if (file.exists(file.path(root, "omics_pred_pipeline.R"))) {
      return(list(
        root = normalizePath(file.path(root, "..", "..", "..")),
        pipeline = normalizePath(file.path(root, "omics_pred_pipeline.R"))
      ))
    }
    pipe <- file.path(root, pipeline_rel)
    if (file.exists(pipe)) {
      return(list(root = normalizePath(root), pipeline = normalizePath(pipe)))
    }
    hit <- file.path(root, "Workstation-MrOS-Microbiome-dev 3", pipeline_rel)
    if (file.exists(hit)) {
      pr <- normalizePath(file.path(root, "Workstation-MrOS-Microbiome-dev 3"))
      return(list(root = pr, pipeline = normalizePath(file.path(pr, pipeline_rel))))
    }
  }

  stop(
    "Cannot find project root / omics_pred_pipeline.R.\n",
    "Fix: Session → Set Working Directory → To Project Directory\n",
    "Project folder should be: Workstation-MrOS-Microbiome-dev 3\n",
    "Current getwd(): ", getwd()
  )
}

resolve_rds <- function(fname, data_dirs) {
  for (d in data_dirs) {
    p <- file.path(d, fname)
    if (file.exists(p)) return(normalizePath(p))
  }
  NA_character_
}

paths <- find_smoke_paths()
root <- paths$root
setwd(root)
message("[SMOKE] project root = ", root)
message("[SMOKE] pipeline     = ", paths$pipeline)

source(paths$pipeline, local = FALSE)

data_dirs <- file.path(root, "data", "processed_data")
data_dirs <- data_dirs[dir.exists(data_dirs)]
if (!length(data_dirs)) {
  stop("No data/processed_data/ found under the repository root: ", root)
}
message("[SMOKE] data dirs    = ", paste(data_dirs, collapse = " | "))

cov_cols_default <- c(
  "race", "edu", "ol_health", "bmi", "mstat", "smoke", "diab",
  "hbp", "cancer", "tmm_score", "gds", "pase", "total_meds", "age_v4"
)
non_feature_cols <- c(
  "ID", "status", "batch", "site", "grip", "chair", "gait400m", "gait6m", "cr_cmm",
  "abx", "b4thd", "b4fnd", "b4lsd", "faprev4",
  "hqdrfefl", "hqdtfefl", "hqptfefl", "fr_chs4", "age_v5", "type"
)

# Manuscript-facing tables + all FS engine methods
layer_defs <- list(
  microbiome = list(
    label = "Microbiome-ANCOMBC",
    files = c("df_status_micro_ancombc.rds", "df_status_micro.rds")
  ),
  metabolomics = list(
    label = "Metabolomics",
    files = c("df_status_metab.rds")
  ),
  proteomics = list(
    label = "Proteomics",
    files = c("df_status_prot.rds", "df_status_prot_clean.rds")
  )
)

fs_cases <- list(
  list(fs_method = "none",  limma_fdr = 0.05, mrmr_top = 10L),
  list(fs_method = "limma", limma_fdr = 0.05, mrmr_top = 10L),
  list(fs_method = "mrmr",  limma_fdr = 0.05, mrmr_top = 10L),
  list(fs_method = "sis",   limma_fdr = 0.05, mrmr_top = 10L)
)

want_layers <- if (exists("SMOKE_LAYERS", inherits = TRUE)) {
  intersect(names(layer_defs), as.character(get("SMOKE_LAYERS", inherits = TRUE)))
} else {
  names(layer_defs)
}
want_fs <- if (exists("SMOKE_FS", inherits = TRUE)) {
  intersect(c("none", "limma", "mrmr", "sis"), as.character(get("SMOKE_FS", inherits = TRUE)))
} else {
  c("none", "limma", "mrmr", "sis")
}
fs_cases <- Filter(function(x) x$fs_method %in% want_fs, fs_cases)

if (!length(want_layers)) stop("SMOKE_LAYERS empty / unknown. Use: microbiome, metabolomics, proteomics")
if (!length(fs_cases)) stop("SMOKE_FS empty / unknown. Use: none, limma, mrmr, sis")

# speed patches for smoke only
suppressWarnings(try(untrace(dials::grid_space_filling), silent = TRUE))
suppressWarnings(try(untrace(rsample::vfold_cv), silent = TRUE))
trace(dials::grid_space_filling, tracer = quote({
  if (!missing(size) && !is.null(size) && is.finite(size[1])) size <- min(as.integer(size[1]), 2L)
}), print = FALSE)
trace(rsample::vfold_cv, tracer = quote({
  if (!missing(v) && !is.null(v) && is.finite(v[1]) && as.integer(v[1]) > 2L) v <- 2L
}), print = FALSE)

if (requireNamespace("kernelshap", quietly = TRUE)) {
  ns <- asNamespace("kernelshap")
  unlockBinding("kernelshap", ns)
  assign("kernelshap", function(object, X, bg_X = NULL, ...) {
    X <- as.data.frame(X)
    structure(list(S = matrix(0, nrow(X), ncol(X))), class = "ks_smoke")
  }, envir = ns)
  lockBinding("kernelshap", ns)
}
if (requireNamespace("shapviz", quietly = TRUE)) {
  ns <- asNamespace("shapviz")
  unlockBinding("shapviz", ns); unlockBinding("sv_importance", ns)
  assign("shapviz", function(...) list(smoke = TRUE), envir = ns)
  assign("sv_importance", function(...) ggplot2::ggplot() + ggplot2::theme_void(), envir = ns)
  lockBinding("shapviz", ns); lockBinding("sv_importance", ns)
}

prepare_layer <- function(rds_path, seed = 2026L, n_per_class = 30L, p_feat = 40L) {
  df0 <- readRDS(rds_path) %>%
    dplyr::filter(status %in% c("Active", "Deceased")) %>%
    dplyr::mutate(status = factor(status, levels = c("Deceased", "Active")))
  cov_cols <- intersect(cov_cols_default, names(df0))
  fcols_all <- setdiff(names(df0), unique(c(cov_cols, non_feature_cols)))
  set.seed(seed)
  df_s <- df0 %>%
    dplyr::group_by(status) %>%
    dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(nrow(.x), n_per_class))) %>%
    dplyr::ungroup()
  avail <- intersect(fcols_all, names(df_s))
  if (!length(avail)) stop("No omics feature columns in ", rds_path)
  fcols <- sample(avail, min(p_feat, length(avail)))
  df_s <- df_s[, unique(c("ID", "status", cov_cols, fcols)), drop = FALSE]
  list(df = df_s, cov_cols = cov_cols, feature_cols = fcols)
}

pass <- character()
fail <- list()
skipped <- character()

for (ly in want_layers) {
  def <- layer_defs[[ly]]
  rds <- NA_character_
  for (fn in def$files) {
    hit <- resolve_rds(fn, data_dirs)
    if (!is.na(hit)) { rds <- hit; break }
  }
  if (is.na(rds)) {
    skipped <- c(skipped, ly)
    message("[SMOKE SKIP] layer=", ly, " | missing RDS: ", paste(def$files, collapse = ", "))
    next
  }
  message("[SMOKE] layer=", ly, " | data=", rds)
  prep <- tryCatch(prepare_layer(rds), error = function(e) e)
  if (inherits(prep, "error")) {
    key <- paste0(ly, "/prep")
    fail[[key]] <- conditionMessage(prep)
    message("[SMOKE FAIL] ", key, " | ", conditionMessage(prep))
    next
  }

  event_class <- levels(prep$df$status)[2]
  for (case in fs_cases) {
    fs <- case$fs_method
    key <- paste0(ly, "/", fs)
    message(
      "[SMOKE] ", key,
      " | n=", nrow(prep$df), " p=", length(prep$feature_cols),
      " | limma_fdr=", case$limma_fdr, " mrmr_top=", case$mrmr_top
    )
    ok <- tryCatch({
      res <- run_omics_ensemble_pipeline(
        df = prep$df,
        outcome = "status",
        id_col = "ID",
        event_level = "second",
        event_class = event_class,
        event_prob_col = paste0(".pred_", event_class),
        n_workers = 2,
        seed = 2026,
        cov_cols = prep$cov_cols,
        feature_cols = prep$feature_cols,
        omics_label = paste0("smoke_", ly, "_", fs),
        fs_method = fs,
        limma_fdr = case$limma_fdr,
        mrmr_top = case$mrmr_top,
        sis_cap = 20L
      )
      stopifnot(!is.null(res$auc), length(res$auc) == 5)
      ests <- vapply(res$auc, function(x) {
        if (is.data.frame(x) && ".estimate" %in% names(x)) as.numeric(x$.estimate[[1]]) else NA_real_
      }, numeric(1))
      message("[SMOKE PASS] ", key, " | ", paste(sprintf("%s=%.3f", names(ests), ests), collapse = " | "))
      TRUE
    }, error = function(e) {
      fail[[key]] <<- conditionMessage(e)
      message("[SMOKE FAIL] ", key, " | ", conditionMessage(e))
      FALSE
    })
    if (isTRUE(ok)) pass <- c(pass, key)
  }
}

suppressWarnings(try(untrace(dials::grid_space_filling), silent = TRUE))
suppressWarnings(try(untrace(rsample::vfold_cv), silent = TRUE))

message(
  "[SMOKE SUMMARY] passed=", if (length(pass)) paste(pass, collapse = ",") else "none",
  " | failed=", if (length(fail)) paste(names(fail), collapse = ",") else "none",
  " | skipped=", if (length(skipped)) paste(skipped, collapse = ",") else "none"
)

if (!length(pass) && !length(fail)) {
  stop("Nothing ran. No omics RDS found for requested layers.")
}
if (length(fail)) {
  stop(
    "Smoke failed:\n",
    paste(sprintf("  - %s: %s", names(fail), unlist(fail)), collapse = "\n")
  )
}
# Require all three layers when no filter was set
if (!exists("SMOKE_LAYERS", inherits = TRUE) && length(skipped)) {
  stop("Missing omics layers (needed for full smoke): ", paste(skipped, collapse = ", "))
}
invisible(pass)
