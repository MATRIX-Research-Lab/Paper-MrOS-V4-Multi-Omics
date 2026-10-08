#!/usr/bin/env Rscript
# Read-only checks for the three expensive/release stages.
#
# Examples:
#   Rscript code/single_omics/preflight.R --stage primary --layer microbiome
#   Rscript code/single_omics/preflight.R --stage matched
#   Rscript code/single_omics/preflight.R --stage render

find_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  starts <- c(getwd(), if (length(fa)) dirname(fa[[1]]) else NULL)
  for (s in unique(starts)) {
    d <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(20L)) {
      if (file.exists(file.path(d, "code", "single_omics", "00_setup.R")) &&
          file.exists(file.path(d, "code", "single_omics", "helpers.R"))) {
        return(d)
      }
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Could not locate the repository root")
}

parse_args <- function(args) {
  out <- list(stage = NULL, layer = "all", allow_count_change = FALSE)
  i <- 1L
  while (i <= length(args)) {
    a <- args[[i]]
    if (a == "--stage" && i < length(args)) {
      out$stage <- tolower(args[[i + 1L]]); i <- i + 2L
    } else if (a == "--layer" && i < length(args)) {
      out$layer <- tolower(args[[i + 1L]]); i <- i + 2L
    } else if (a == "--allow-count-change") {
      out$allow_count_change <- TRUE; i <- i + 1L
    } else if (a %in% c("--help", "-h")) {
      cat("Usage: Rscript code/single_omics/preflight.R --stage primary|matched|render [--layer ...]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", a, call. = FALSE)
    }
  }
  if (is.null(out$stage) || !out$stage %in% c("primary", "matched", "render")) {
    stop("--stage must be primary, matched, or render")
  }
  if (!out$layer %in% c("all", "microbiome", "metabolomics", "proteomics")) {
    stop("--layer must be all, microbiome, metabolomics, or proteomics")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- find_root()
setwd(project_root)
source(file.path(project_root, "code", "single_omics", "00_setup.R"))

primary_inputs <- function() {
  c(
    microbiome = file.path(paths$processed, "df_status_micro_ancombc.rds"),
    metabolomics = file.path(paths$processed, "df_status_metab.rds"),
    proteomics = file.path(paths$processed, "df_status_prot.rds")
  )
}

check_primary_input <- function(path, layer) {
  if (!file.exists(path)) stop(layer, ": missing modeling frame: ", path)
  df <- readRDS(path)
  if (!is.data.frame(df)) stop(layer, ": modeling input is not a data.frame: ", path)
  required <- c("ID", "status", analysis_covariates)
  missing_required <- setdiff(required, names(df))
  if (length(missing_required)) {
    stop(layer, ": modeling input is missing required columns: ",
         paste(missing_required, collapse = ", "))
  }
  df <- df[df$status %in% c("Active", "Deceased"), , drop = FALSE]
  df$status <- factor(df$status, levels = c("Deceased", "Active"))
  feature_cols <- setdiff(names(df), unique(c(
    analysis_covariates, analysis_non_feature_columns
  )))
  validate_analytic_frame(df, "ID", "status", analysis_covariates, feature_cols,
                          paste0("preflight / ", layer))
  n <- nrow(df)
  deaths <- sum(df$status == "Deceased")
  expected <- analysis_contract[[layer]]
  if (!args$allow_count_change && n != expected[["n"]]) {
    stop(layer, ": expected n = ", expected[["n"]], " but found ", n,
         ". Pass --allow-count-change only after explicitly updating the analysis specification.")
  }
  if (!args$allow_count_change && deaths != expected[["deaths"]]) {
    stop(layer, ": expected deaths = ", expected[["deaths"]], " but found ", deaths,
         ". Pass --allow-count-change only after explicitly updating the analysis specification.")
  }
  if (length(feature_cols) != expected[["features"]]) {
    stop(layer, ": expected ", expected[["features"]], " features but found ", length(feature_cols), ".")
  }
  data.frame(layer = layer, path = path, n = n, deaths = deaths,
             features = length(feature_cols), stringsAsFactors = FALSE)
}

if (args$stage == "primary") {
  assert_packages(primary_packages, "primary preflight")
  assert_modeling_api()
  suppressPackageStartupMessages(library(dplyr))
  layers <- if (args$layer == "all") names(analysis_contract) else args$layer
  inp <- primary_inputs()
  report <- do.call(rbind, lapply(layers, function(layer) check_primary_input(inp[[layer]], layer)))
  print(report, row.names = FALSE)
  message("Primary preflight passed. The runner will use limma FDR 0.05 for these layers.")
}

if (args$stage %in% c("matched", "render")) {
  # The two stages read different objects, for the reason set out in 00_setup.R:
  # 06_ builds the matched cohort from the re-fitted layers, which carry the
  # provenance manifests this checks; the render reads the primary analysis,
  # which predates them. Both must describe the same 332 men.
  is_matched <- args$stage == "matched"
  if (is_matched) {
    assert_packages(matched_packages, "matched preflight")
    assert_modeling_api()
  } else {
    assert_packages(render_packages, "render preflight")
  }
  suppressPackageStartupMessages(library(dplyr))
  primary <- load_prediction_model_outputs(
    paths,
    dir = if (is_matched) paths$model_outputs_refit else paths$model_outputs,
    require_provenance = is_matched,
    expect_score = if (is_matched) "deceased" else "auto"
  )
  source_files <- attr(primary, "source_files")
  if (!all(basename(source_files) == c("microbiome_model_outputs.rds",
                                       "metabolomics_model_outputs.rds",
                                       "proteomics_model_outputs.rds"))) {
    stop("Matched/render stages require the three canonical primary outputs.")
  }
  if (is_matched) {
    for (i in seq_along(primary)) {
      cfg <- primary[[i]]$provenance$signature$config
      if (!identical(cfg$analysis_mode, "primary") ||
          !identical(cfg$specification$fs_method, "limma") ||
          !identical(cfg$event_level, "first") ||
          !identical(cfg$event_class, "Deceased")) {
        stop(names(primary)[[i]],
             ": canonical output is not from the prespecified death-risk limma analysis")
      }
    }
  }
  full <- lapply(primary, function(x) {
    dplyr::bind_rows(rsample::training(x$split), rsample::testing(x$split))
  })
  ids <- lapply(full, function(d) as.character(d$ID))
  if (any(vapply(ids, function(x) anyDuplicated(x) > 0L, logical(1)))) {
    stop("Primary outputs contain duplicate IDs")
  }
  overlap <- Reduce(intersect, ids)
  if (length(overlap) != 332L) stop("Expected 332 matched participants, found ", length(overlap))
  status <- lapply(full, function(d) as.character(d$status[match(overlap, d$ID)]))
  if (!all(vapply(status[-1], identical, logical(1), status[[1]]))) {
    stop("Vital status differs across primary layers in the matched cohort")
  }
  if (sum(status[[1]] == "Deceased") != 169L) {
    stop("Expected 169 matched-cohort deaths, found ", sum(status[[1]] == "Deceased"))
  }
  print(data.frame(n_matched = length(overlap), deaths = sum(status[[1]] == "Deceased"),
                   stringsAsFactors = FALSE), row.names = FALSE)
  if (is_matched) {
    message("Matched preflight passed. Run 06_fit_matched_cohort_models.R next.")
  } else {
    matched_file <- file.path(paths$data, "model_outputs", "matched_cohort",
                              "matched_cohort_model_outputs.rds")
    if (!file.exists(matched_file)) stop("Render preflight: missing matched output: ", matched_file)
    mc <- load_matched_cohort_outputs(paths, require_provenance = TRUE)
    if (mc$cohort$n != 332L || mc$cohort$n_deaths != 169L) {
      stop("Matched output has unexpected cohort dimensions")
    }
    message("Render preflight passed for canonical primary and matched outputs.")
  }
}
