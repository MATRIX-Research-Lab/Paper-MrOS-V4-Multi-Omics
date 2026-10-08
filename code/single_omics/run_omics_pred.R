#!/usr/bin/env Rscript
# =============================================================================
# Run single-omics FS × model combinations via omics_pred_pipeline.R
#
# From repo root:
#   Rscript code/single_omics/run_omics_pred.R --plan-only
#   Rscript code/single_omics/run_omics_pred.R --quick
#   Rscript code/single_omics/run_omics_pred.R
#   Rscript code/single_omics/run_omics_pred.R --layer proteomics
#
# Needs data/processed_data/*.rds (gitignored).
# =============================================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---- specs (merged from former all_fs_specs.R) ----
cov_cols_default <- c(
  "race", "edu", "ol_health", "bmi", "mstat", "smoke", "diab",
  "hbp", "cancer", "tmm_score", "gds", "pase", "total_meds", "age_v4"
)
name_map_default <- c(
  race = "Race", edu = "Education", ol_health = "Self-Rated Overall Health",
  bmi = "BMI", mstat = "Marital Status", smoke = "Smoking Status",
  diab = "Diabetes", hbp = "High Blood Pressure", cancer = "Cancer",
  tmm_score = "Teng 3MS Score", gds = "Geriatric Depression Scale",
  pase = "PASE Score", total_meds = "Total Medications", age_v4 = "Age"
)
cat_prefixes_default <- c("race", "edu", "mstat", "smoke", "diab", "hbp", "cancer")
non_feature_cols_default <- c(
  "ID", "status", "batch", "site", "grip", "chair", "gait400m", "gait6m", "cr_cmm",
  "abx", "b4thd", "b4fnd", "b4lsd", "faprev4",
  "hqdrfefl", "hqdtfefl", "hqptfefl", "fr_chs4", "age_v5", "type"
)

fs_specs_full <- list(
  noFS      = list(label = "No FS",          fs_method = "none"),
  limma_005 = list(label = "limma FDR 0.05", fs_method = "limma", limma_fdr = 0.05),
  limma_010 = list(label = "limma FDR 0.10", fs_method = "limma", limma_fdr = 0.10),
  mrmr_50   = list(label = "mRMR top 50",    fs_method = "mrmr",  mrmr_top = 50),
  mrmr_100  = list(label = "mRMR top 100",   fs_method = "mrmr",  mrmr_top = 100),
  mrmr_150  = list(label = "mRMR top 150",   fs_method = "mrmr",  mrmr_top = 150),
  mrmr_200  = list(label = "mRMR top 200",   fs_method = "mrmr",  mrmr_top = 200),
  sis_def   = list(label = "SIS default",    fs_method = "sis")
)
fs_specs_quick <- fs_specs_full[c("noFS", "limma_005", "mrmr_50", "sis_def")]

# Primary microbiome table uses the same FS engine as metab/prot.
# Other preprocess tables stay No-FS (abundance already filtered upstream).
micro_jobs_full <- list(
  list(prefix = "Microbiome-ANCOMBC", path = "df_status_micro_ancombc.rds",
       specs = fs_specs_full),
  list(prefix = "Microbiome-DESeq2", path = "df_status_micro_deseq.rds",
       specs = list(nofs = list(label = "DESeq2 (No FS)", fs_method = "none"))),
  list(prefix = "Microbiome-CLR", path = "df_status_micro_clr_all.rds",
       specs = list(nofs = list(label = "CLR (No FS)", fs_method = "none"))),
  list(prefix = "Microbiome-ANCOMBC-raw005", path = "df_status_micro_ancombc_raw005.rds",
       specs = list(nofs = list(label = "ANCOM-BC rawP<0.05 (No FS)", fs_method = "none"))),
  list(prefix = "Microbiome-ANCOMBC-fdr005", path = "df_status_micro_ancombc_fdr005.rds",
       specs = list(nofs = list(label = "ANCOM-BC FDR<0.05 (No FS)", fs_method = "none"))),
  list(prefix = "Microbiome-ANCOMBC-fdr010", path = "df_status_micro_ancombc_fdr010.rds",
       specs = list(nofs = list(label = "ANCOM-BC FDR<0.10 (No FS)", fs_method = "none")))
)
micro_jobs_quick <- list(
  list(prefix = "Microbiome-ANCOMBC", path = "df_status_micro_ancombc.rds",
       specs = fs_specs_quick)
)

# ---- tiny helpers ----
find_project_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = TRUE)
  for (i in seq_len(16)) {
    if (file.exists(file.path(d, "code", "single_omics", "omics_pred_pipeline.R")) &&
        file.exists(file.path(d, "proj.Rproj"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Run from Workstation-MrOS-Microbiome-dev root.")
}

prepare_status_outcome <- function(df, outcome = "status") {
  df <- df %>% dplyr::filter(.data[[outcome]] %in% c("Active", "Deceased"))
  df[[outcome]] <- factor(df[[outcome]], levels = c("Deceased", "Active"))
  df
}

omics_feature_cols <- function(df, cov_cols) {
  setdiff(names(df), unique(c(cov_cols, non_feature_cols_default)))
}

run_one_setting <- function(spec, df, outcome, id_col, event_level, event_class,
                            event_prob_col, cov_cols, feature_cols, omics_prefix,
                            seed, n_workers, out_rds, overwrite = FALSE) {
  dir.create(dirname(out_rds), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(out_rds) && !overwrite) {
    message("[CACHE HIT] ", out_rds)
    return(readRDS(out_rds))
  }
  message("[RUN] ", omics_prefix, " | ", spec$label)
  res <- run_omics_ensemble_pipeline(
    df = df, outcome = outcome, id_col = id_col,
    event_level = event_level, event_class = event_class,
    event_prob_col = event_prob_col, n_workers = n_workers, seed = seed,
    cov_cols = cov_cols, feature_cols = feature_cols,
    omics_label = paste0(omics_prefix, " (", spec$label, ")"),
    fs_method = spec$fs_method %||% "none",
    limma_fdr = spec$limma_fdr %||% 0.05,
    mrmr_top  = spec$mrmr_top  %||% 50,
    sis_top   = spec$sis_top   %||% NULL,
    sis_cap   = spec$sis_cap   %||% 200,
    sis_method = spec$sis_method %||% "spearman"
  )
  saveRDS(res, out_rds)
  message("[SAVED] ", out_rds)
  res
}

save_plots_from_res <- function(res, fig_dir, prefix, height = 5, width = 6.25) {
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  plots <- res$plots
  if (is.null(plots)) return(invisible(FALSE))
  for (k in names(plots)) {
    p <- plots[[k]]
    if (is.null(p)) next
    tryCatch(
      ggplot2::ggsave(
        file.path(fig_dir, paste0(prefix, "_", sub("^p_", "", k), ".pdf")),
        p, height = height, width = width
      ),
      error = function(e) message("[PLOT SKIP] ", k, ": ", conditionMessage(e))
    )
  }
  invisible(TRUE)
}

run_omics_specs <- function(df, outcome, id_col, event_level, event_class,
                            event_prob_col, cov_cols, feature_cols, omics_prefix,
                            tag, specs, out_dir, fig_dir, seed = 2026,
                            n_workers = 5, overwrite = FALSE) {
  results <- list()
  for (nm in names(specs)) {
    spec <- specs[[nm]]
    prefix <- paste0(tolower(omics_prefix), "_", nm, "_", tag)
    out_rds <- file.path(out_dir, paste0(prefix, ".rds"))
    fig_subdir <- file.path(fig_dir, tolower(omics_prefix), tag)
    res <- run_one_setting(
      spec, df, outcome, id_col, event_level, event_class, event_prob_col,
      cov_cols, feature_cols, omics_prefix, seed, n_workers, out_rds, overwrite
    )
    save_plots_from_res(res, fig_subdir, prefix)
    results[[nm]] <- res
    rm(res); gc()
  }
  results
}

# ---- CLI ----
parse_args <- function(args) {
  out <- list(layer = "all", quick = FALSE, plan_only = FALSE,
              overwrite = FALSE, seed = 2026L, n_workers = 5L)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a %in% c("--layer", "-l") && i < length(args)) {
      out$layer <- tolower(args[[i + 1L]]); i <- i + 2L
    } else if (a %in% c("--quick", "-q")) {
      out$quick <- TRUE; i <- i + 1L
    } else if (a %in% c("--plan-only", "-n")) {
      out$plan_only <- TRUE; i <- i + 1L
    } else if (a %in% c("--overwrite", "-f")) {
      out$overwrite <- TRUE; i <- i + 1L
    } else if (a == "--seed" && i < length(args)) {
      out$seed <- as.integer(args[[i + 1L]]); i <- i + 2L
    } else if (a == "--n-workers" && i < length(args)) {
      out$n_workers <- as.integer(args[[i + 1L]]); i <- i + 2L
    } else if (a %in% c("--help", "-h")) {
      cat("Usage: Rscript code/single_omics/run_omics_pred.R [--layer ...] [--quick] [--plan-only] [--overwrite]\n")
      quit(save = "no", status = 0)
    } else stop("Unknown argument: ", a, call. = FALSE)
  }
  out
}

# Only run CLI when executed as a script (not when sourced by smoke)
is_main <- !interactive() && any(grepl("run_omics_pred\\.R$", commandArgs(trailingOnly = FALSE)))
# Rscript always non-interactive; detect via --file=
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
is_main <- length(file_arg) > 0 && grepl("run_omics_pred\\.R$", file_arg[[1]])

if (is_main) {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  suppressPackageStartupMessages({ library(tidyverse); library(tidymodels) })

  this_file <- normalizePath(file_arg[[1]])
  mod_dir <- dirname(this_file)
  root <- find_project_root(mod_dir)
  setwd(root)
  source(file.path(mod_dir, "omics_pred_pipeline.R"))

  fs_specs <- if (args$quick) fs_specs_quick else fs_specs_full
  micro_jobs <- if (args$quick) micro_jobs_quick else micro_jobs_full
  out_dir <- file.path(root, "results", "omics_pred", "outputs")
  fig_dir <- file.path(root, "results", "omics_pred", "figures")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  processed <- file.path(root, "data", "processed_data")
  if (!dir.exists(processed)) {
    stop("Missing data/processed_data/ under the repository root: ", processed)
  }
  message("[DATA] ", processed)

  plan_rows <- list()
  for (job in micro_jobs) {
    for (sn in names(job$specs)) {
      plan_rows[[length(plan_rows) + 1]] <- tibble::tibble(
        layer = job$prefix, file = job$path, spec = sn,
        fs = job$specs[[sn]]$fs_method
      )
    }
  }
  for (layer in c("Metabolomics", "Proteomics")) {
    for (sn in names(fs_specs)) {
      plan_rows[[length(plan_rows) + 1]] <- tibble::tibble(
        layer = layer,
        file = if (layer == "Metabolomics") "df_status_metab.rds" else "df_status_prot.rds",
        spec = sn, fs = fs_specs[[sn]]$fs_method
      )
    }
  }
  plan <- dplyr::bind_rows(plan_rows)
  if (!identical(args$layer, "all")) {
    keep <- switch(
      args$layer,
      microbiome = grepl("^Microbiome", plan$layer),
      metabolomics = plan$layer == "Metabolomics",
      proteomics = plan$layer == "Proteomics",
      stop("bad --layer")
    )
    plan <- plan[keep, , drop = FALSE]
  }
  message("[ROOT] ", root)
  message("[JOBS] ", nrow(plan), " | RF|EN|XGB|RF+EN|RF+XGB each")
  print(plan)
  if (args$plan_only) quit(save = "no", status = 0)

  id_col <- "ID"; outcome <- "status"; event_level <- "second"
  cov_cols <- cov_cols_default

  run_layer <- function(path, prefix, specs) {
    if (!file.exists(path)) { message("[SKIP] ", path); return(invisible(NULL)) }
    df <- prepare_status_outcome(readRDS(path))
    event_class <- levels(df[[outcome]])[2]
    run_omics_specs(
      df = df, outcome = outcome, id_col = id_col, event_level = event_level,
      event_class = event_class, event_prob_col = paste0(".pred_", event_class),
      cov_cols = intersect(cov_cols, names(df)),
      feature_cols = omics_feature_cols(df, cov_cols),
      omics_prefix = prefix, tag = "status", specs = specs,
      out_dir = out_dir, fig_dir = fig_dir, seed = args$seed,
      n_workers = args$n_workers, overwrite = args$overwrite
    )
  }

  if (args$layer %in% c("all", "microbiome")) {
    for (job in micro_jobs) {
      run_layer(file.path(processed, job$path), job$prefix, job$specs)
    }
  }
  if (args$layer %in% c("all", "metabolomics")) {
    run_layer(file.path(processed, "df_status_metab.rds"), "Metabolomics", fs_specs)
  }
  if (args$layer %in% c("all", "proteomics")) {
    pp <- file.path(processed, "df_status_prot.rds")
    if (!file.exists(pp)) pp <- file.path(processed, "df_status_prot_clean.rds")
    run_layer(pp, "Proteomics", fs_specs)
  }
  message("[DONE] ", out_dir)
}
