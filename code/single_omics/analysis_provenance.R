# Shared provenance, validation, and atomic-output helpers.
#
# These functions are deliberately dependency-light so that every entrypoint
# can use them before loading the model-specific packages.

repo_git_state <- function(root = project_root) {
  commit <- tryCatch(
    system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE),
    error = function(e) character()
  )
  dirty <- tryCatch(
    length(system2("git", c("-C", shQuote(root), "status", "--porcelain"), stdout = TRUE, stderr = FALSE)) > 0L,
    error = function(e) NA
  )
  list(
    commit = if (length(commit)) trimws(commit[[1]]) else NA_character_,
    dirty = dirty
  )
}

file_manifest <- function(files, root = project_root) {
  files <- unique(as.character(files))
  files <- files[nzchar(files)]
  missing <- files[!file.exists(files)]
  if (length(missing)) {
    stop("Cannot build a provenance manifest; missing files: ",
         paste(missing, collapse = ", "))
  }
  files <- unique(normalizePath(files, winslash = "/", mustWork = TRUE))
  if (!length(files)) {
    return(data.frame(path = character(), md5 = character(), bytes = numeric(),
                      stringsAsFactors = FALSE))
  }
  rel <- function(x) {
    prefix <- paste0(normalizePath(root, winslash = "/", mustWork = TRUE), "/")
    sub(prefix, "", x, fixed = TRUE)
  }
  data.frame(
    path = vapply(files, rel, character(1)),
    md5 = unname(as.character(tools::md5sum(files))),
    bytes = as.numeric(file.info(files)$size),
    stringsAsFactors = FALSE
  )
}

package_versions <- function(packages) {
  packages <- unique(as.character(packages))
  out <- vapply(packages, function(p) {
    if (!requireNamespace(p, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(p))
  }, character(1))
  unname(out)
}

assert_packages <- function(packages, stage = "workflow") {
  packages <- unique(as.character(packages))
  ok <- vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  if (any(!ok)) {
    stop(stage, " requires missing R packages: ", paste(packages[!ok], collapse = ", "),
         ". Install the analysis environment before running this stage.")
  }
  invisible(TRUE)
}

make_provenance <- function(stage, inputs = character(), code_files = character(),
                            config = list(), seed = NULL, packages = character()) {
  git <- repo_git_state()
  r_version <- paste(R.version$major, R.version$minor, sep = ".")
  input_manifest <- file_manifest(inputs)
  code_manifest <- file_manifest(code_files)
  if (!length(packages)) packages <- character()
  signature <- list(
    stage = stage,
    seed = seed,
    config = config,
    r_version = r_version,
    packages = package_versions(packages),
    package_names = as.character(packages),
    input_manifest = input_manifest,
    code_manifest = code_manifest,
    git_commit = git$commit
  )
  list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    git = git,
    signature = signature
  )
}

provenance_matches <- function(object, expected) {
  !is.null(object$provenance) &&
    !is.null(object$provenance$signature) &&
    identical(object$provenance$signature, expected$signature)
}

atomic_save_rds <- function(object, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(object, tmp)
  if (!file.rename(tmp, file)) stop("Could not atomically replace ", file)
  invisible(file)
}

assert_unique_ids <- function(df, id_col = "ID", label = "data") {
  if (!id_col %in% names(df)) stop(label, ": missing ID column ", id_col)
  ids <- as.character(df[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids))) stop(label, ": missing participant IDs")
  if (anyDuplicated(ids)) stop(label, ": duplicated participant IDs")
  invisible(TRUE)
}

validate_analytic_frame <- function(df, id_col = "ID", outcome = "status",
                                    cov_cols = character(), feature_cols = character(),
                                    label = "analytic data") {
  if (!is.data.frame(df)) stop(label, ": expected a data.frame")
  assert_unique_ids(df, id_col, label)
  required <- unique(c(id_col, outcome, cov_cols, feature_cols))
  missing <- setdiff(required, names(df))
  if (length(missing)) stop(label, ": missing columns: ", paste(missing, collapse = ", "))
  if (anyNA(df[[outcome]])) stop(label, ": outcome contains missing values")
  if (length(feature_cols)) {
    non_numeric <- feature_cols[!vapply(df[feature_cols], is.numeric, logical(1))]
    if (length(non_numeric)) stop(label, ": non-numeric omics features: ", paste(head(non_numeric, 10), collapse = ", "))
    infinite <- feature_cols[vapply(df[feature_cols], function(x) any(is.infinite(x)), logical(1))]
    if (length(infinite)) stop(label, ": infinite omics values: ", paste(head(infinite, 10), collapse = ", "))
  }
  numeric_covariates <- cov_cols[vapply(df[cov_cols], is.numeric, logical(1))]
  if (length(numeric_covariates)) {
    infinite <- numeric_covariates[vapply(df[numeric_covariates], function(x) any(is.infinite(x)), logical(1))]
    if (length(infinite)) stop(label, ": infinite clinical values: ", paste(infinite, collapse = ", "))
  }
  invisible(TRUE)
}

# `expect_score` selects how the stored p1..p5 / p1..p3 score column relates to
# the predicted death risk. The two saved generations of primary outputs differ
# here and both are internally consistent:
#
#   "deceased"  p == .pred_Deceased. The re-fitted outputs under
#               data/model_outputs/single_omics_rerun/, which 06_ reads.
#   "active"    p == 1 - .pred_Deceased, i.e. the survival-oriented score. The
#               primary outputs under data/model_outputs/single_omics/. Their
#               readers flip it explicitly -- see out_of_sample_predictions() --
#               and every AUC is computed from .pred_Deceased regardless, so the
#               published numbers are unaffected by which convention is stored.
#   "auto"      accept either, and record which was found on the returned value.
#
# "auto" is the default so the primary objects get every other check here --
# test-set ID order, truth alignment, probability ranges, out-of-fold row keys --
# rather than failing at the score column and being skipped entirely.
validate_prediction_output <- function(x, label = "model output",
                                       expect_score = c("auto", "deceased", "active")) {
  expect_score <- match.arg(expect_score)
  if (!is.list(x)) stop(label, ": expected a list")
  required <- c("split", "feature_selection", "auc", "fits", "predictions")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, ": missing components: ", paste(missing, collapse = ", "))
  if (is.null(x$split$data) || is.null(x$split$in_id)) stop(label, ": incomplete rsample split")
  if (!all(c("ID", "status") %in% names(x$split$data))) {
    stop(label, ": split data must contain ID and status")
  }
  assert_unique_ids(x$split$data, "ID", label)
  in_id <- as.integer(x$split$in_id)
  if (length(in_id) < 2L || anyNA(in_id) || any(in_id < 1L) || any(in_id > nrow(x$split$data)) || anyDuplicated(in_id)) {
    stop(label, ": invalid training-row indices in rsample split")
  }
  test_id <- as.character(x$split$data$ID[setdiff(seq_len(nrow(x$split$data)), in_id)])
  if (!length(test_id)) stop(label, ": empty held-out test set")
  # Resolved from the first prediction frame when expect_score is "auto", then
  # required to hold for every remaining one: a object that mixes conventions
  # across models is broken in a way neither convention would explain.
  score_convention <- if (identical(expect_score, "auto")) NA_character_ else expect_score
  check_score <- function(p, p_deceased, where) {
    is_deceased <- isTRUE(all.equal(as.numeric(p), as.numeric(p_deceased),
                                    tolerance = 1e-12, check.attributes = FALSE))
    is_active <- isTRUE(all.equal(as.numeric(p), 1 - as.numeric(p_deceased),
                                  tolerance = 1e-12, check.attributes = FALSE))
    found <- if (is_deceased) "deceased" else if (is_active) "active" else NA_character_
    if (is.na(found)) {
      stop(where, ": stored model score is neither the predicted death risk nor its complement")
    }
    if (is.na(score_convention)) {
      score_convention <<- found
    } else if (!identical(found, score_convention)) {
      stop(where, ": stored model score is '", found, "' but '", score_convention,
           "' was expected")
    }
    invisible(found)
  }

  pred_names <- c("prob_test_cov", "prob_test_micro_en", "prob_test_micro",
                  "prob_test_4", "prob_test_5")
  if (!all(pred_names %in% names(x$predictions))) {
    stop(label, ": missing primary prediction components")
  }
  expected_truth <- as.character(x$split$data$status[setdiff(seq_len(nrow(x$split$data)), in_id)])
  for (nm in pred_names) {
    d <- x$predictions[[nm]]
    if (!all(c("truth", ".pred_Deceased") %in% names(d))) {
      stop(label, ": ", nm, " lacks truth/.pred_Deceased")
    }
    # The stacked frames in the primary outputs carry no ID column; the
    # re-fitted ones do. Where IDs exist they must match the held-out test set
    # exactly and in order. Where they do not, row order against the split is
    # the only key available, and it is what the readers rely on -- so check
    # the row count and let the truth-label check below carry the alignment.
    if ("ID" %in% names(d)) {
      assert_unique_ids(d, "ID", paste0(label, "/", nm))
      if (!identical(as.character(d$ID), test_id)) {
        stop(label, "/", nm, ": prediction IDs do not exactly match the held-out test set in order")
      }
    } else if (nrow(d) != length(test_id)) {
      stop(label, "/", nm, ": has no ID column and does not have one row per held-out participant")
    }
    if (!identical(as.character(d$truth), expected_truth)) {
      stop(label, "/", nm, ": prediction truth labels do not match the held-out test set")
    }
    if (anyNA(d$.pred_Deceased) || any(!is.finite(d$.pred_Deceased)) ||
        any(d$.pred_Deceased < 0 | d$.pred_Deceased > 1)) {
      stop(label, "/", nm, ": invalid predicted probabilities")
    }
    p_cols <- grep("^p[1-5]$", names(d), value = TRUE)
    if (length(p_cols) != 1L) {
      stop(label, "/", nm, ": expected exactly one stored model score column")
    }
    check_score(d[[p_cols]], d$.pred_Deceased, paste0(label, "/", nm))
  }
  train_id <- as.character(x$split$data$ID[in_id])
  train_truth <- as.character(x$split$data$status[in_id])
  oof_names <- c("pred_oof_cov", "pred_oof_micro_en", "pred_oof_micro_xgb")
  if (!all(oof_names %in% names(x$predictions))) {
    stop(label, ": missing nested OOF prediction components")
  }
  for (nm in oof_names) {
    d <- x$predictions[[nm]]
    if (!".train_row" %in% names(d)) {
      stop(label, "/", nm, ": OOF predictions lack .train_row")
    }
    # .train_row keys each out-of-fold prediction back to a training row. Check
    # the mapping, not the row order: these frames are written fold by fold in
    # the primary outputs and in training order in the re-fitted ones, and every
    # reader indexes through .train_row rather than assuming a layout.
    row_id <- as.integer(d$.train_row)
    if (!identical(sort(row_id), seq_along(train_id))) {
      stop(label, "/", nm, ": OOF row keys are not a permutation of the training rows")
    }
    # Only the re-fitted outputs store truth alongside the omics OOF scores; in
    # the primary outputs it is carried by .train_row alone. Verify it where it
    # is there rather than requiring the newer layout.
    if ("truth" %in% names(d) &&
        !identical(as.character(d$truth), train_truth[row_id])) {
      stop(label, "/", nm, ": OOF truth labels are not aligned to the training rows they key to")
    }
    p_cols <- grep("^p[1-3]$", names(d), value = TRUE)
    if (length(p_cols) != 1L || anyNA(d[[p_cols]]) ||
        any(!is.finite(d[[p_cols]])) || any(d[[p_cols]] < 0 | d[[p_cols]] > 1)) {
      stop(label, "/", nm, ": invalid OOF predicted probabilities")
    }
    # The OOF frames carry no .pred_Deceased to compare against, so the
    # convention resolved from the test-set frames above is what applies.
  }
  invisible(structure(TRUE, score_convention = score_convention))
}

write_manifest_csv <- function(provenance, file) {
  rows <- list()
  add <- function(kind, d) {
    if (!is.null(d) && nrow(d)) {
      rows[[length(rows) + 1L]] <<- data.frame(
        kind = kind, d, check.names = FALSE, stringsAsFactors = FALSE
      )
    }
  }
  add("input", provenance$signature$input_manifest)
  add("code", provenance$signature$code_manifest)
  if (!length(rows)) return(invisible(FALSE))
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(out, file, row.names = FALSE)
  invisible(TRUE)
}
