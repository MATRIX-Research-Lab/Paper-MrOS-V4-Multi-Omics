# Matched-cohort refit: every omics layer fitted on the same 332 men, on
# identical folds, against a single clinical model.
#
# Why this exists. The primary analysis fits each layer inside its own analytic
# sample (765 / 447 / 446 men), each with its own train/test split. That makes
# the layer-specific AUCs incomparable, because the case mix differs, and it
# means there is no single clinical benchmark: the clinical model is refitted in
# each subset, so it has three different AUCs. Restricting the *saved*
# predictions to the overlapping men mitigates the first problem but not the
# second, and mixes predictions of different provenance. Refitting everything on
# the matched cohort removes both problems at once.
#
# What is different from the primary analysis:
#   * one cohort (the 332 men with all three layers), one set of folds, shared
#     by the clinical model and by all three layers;
#   * genuinely nested cross-validation. Hyperparameters are tuned on an inner
#     3-fold split of each outer fold's analysis set, and the outer assessment
#     fold is used only for evaluation. Feature screening likewise happens
#     inside the analysis set, so the assessment fold never influences it.
#
# Nothing here touches the primary analysis: the layer-specific model outputs in
# data/model_outputs/single_omics are read only.
#
# Outputs: data/model_outputs/matched_cohort/*.rds, consumed by
# 07_make_matched_cohort_tables.R.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidymodels)
})

env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out) || out < 1L) stop(name, " must be a positive integer")
  out
}

env_bool <- function(name, default) {
  value <- tolower(Sys.getenv(name, unset = ""))
  if (!nzchar(value)) return(isTRUE(default))
  if (value %in% c("1", "true", "yes", "y")) return(TRUE)
  if (value %in% c("0", "false", "no", "n")) return(FALSE)
  stop(name, " must be true/false")
}

MATCHED <- list(
  seed       = env_int("MROS_MATCHED_SEED", 2026L),
  outer_v    = 5L,
  inner_v    = 3L,
  n_repeats  = env_int("MROS_MATCHED_REPEATS", 5L),
  n_workers  = env_int("MROS_MATCHED_WORKERS", Sys.getenv("SLURM_CPUS_PER_TASK", "8")),
  tune_parallel_over = analysis_tune_parallel_over(),
  limma_fdr  = 0.05,
  limma_fallback_n = 50L,
  # The matched release specification screens only the high-dimensional
  # proteome; microbiome and metabolome use their complete frozen feature
  # blocks. Keeping this in the design object makes the choice explicit.
  screened_layers = "proteomics",
  # Reuse per-job checkpoints from an earlier run instead of refitting. Set to
  # FALSE to force a clean run.
  resume = env_bool("MROS_MATCHED_RESUME", TRUE)
)

matched_dir <- file.path(paths$data, "model_outputs", "matched_cohort")
dir.create(matched_dir, recursive = TRUE, showWarnings = FALSE)

# Columns that are phenotype/metadata rather than omics features. Identical in
# all three saved analytic samples; the omics block is whatever remains.
MATCHED_META_COLS <- analysis_matched_metadata_columns

# The 14 covariates of the clinical model, as in the primary analysis.
MATCHED_COV_COLS <- analysis_covariates

# Expected feature-block sizes, from Methods. Asserted, not assumed.
MATCHED_EXPECTED_FEATURES <- vapply(
  analysis_contract, function(x) x[["features"]], integer(1)
)

# --- cohort assembly -------------------------------------------------------
#
# The saved split objects carry the complete analytic sample including every
# feature column, so the matched cohort is assembled from them directly. Nothing
# is recomputed: the microbiome block is already ANCOM-BC2 bias-corrected and the
# metabolomics and proteomics blocks are already in their modelling form, so the
# refit sees byte-identical inputs to the primary analysis.
build_matched_cohort <- function(model_outputs) {
  full <- lapply(model_outputs, function(x) {
    dplyr::bind_rows(rsample::training(x$split), rsample::testing(x$split))
  })
  for (om in names(full)) assert_unique_ids(full[[om]], "ID", paste0("matched / ", om))
  ids <- Reduce(intersect, lapply(full, function(d) as.character(d$ID)))
  if (length(ids) != 332L) {
    stop("Expected 332 complete-case participants, found ", length(ids),
         ". Check the primary input objects and participant-ID harmonisation.")
  }

  feature_cols <- lapply(full, function(d) setdiff(names(d), MATCHED_META_COLS))
  for (om in names(feature_cols)) {
    validate_analytic_frame(
      full[[om]], id_col = "ID", outcome = "status",
      cov_cols = MATCHED_COV_COLS, feature_cols = feature_cols[[om]],
      label = paste0("matched / ", om)
    )
    if (length(feature_cols[[om]]) != MATCHED_EXPECTED_FEATURES[[om]]) {
      stop(sprintf("%s: expected %d feature columns, found %d",
                   om, MATCHED_EXPECTED_FEATURES[[om]], length(feature_cols[[om]])))
    }
  }

  aligned <- lapply(full, function(d) d[match(ids, as.character(d$ID)), , drop = FALSE])

  # The clinical covariates and the outcome must agree across layers for the
  # same men; if they do not, the layers are not describing the same cohort.
  ref_cov <- aligned[[1]][, MATCHED_COV_COLS, drop = FALSE]
  ref_status <- as.character(aligned[[1]]$status)
  same_column <- function(a, b) {
    if (is.factor(a) || is.factor(b)) {
      identical(as.character(a), as.character(b))
    } else {
      isTRUE(all.equal(a, b, check.attributes = FALSE))
    }
  }
  for (om in names(aligned)[-1]) {
    cov_same <- vapply(MATCHED_COV_COLS, function(col) {
      same_column(ref_cov[[col]], aligned[[om]][[col]])
    }, logical(1))
    if (!all(cov_same)) {
      stop("clinical covariates differ across layers for the matched participants")
    }
    if (!identical(ref_status, as.character(aligned[[om]]$status))) {
      stop("vital status differs across layers for the matched participants")
    }
  }

  status <- factor(ref_status, levels = c("Deceased", "Active"))
  if (anyNA(status) || sum(status == "Deceased") == 0L || sum(status == "Active") == 0L) {
    stop("Matched cohort status is missing or has only one class")
  }
  clinical <- tibble::tibble(ID = ids, status = status) |>
    dplyr::bind_cols(ref_cov)
  if (anyNA(aligned[[1]]$age_v4) || any(!is.finite(aligned[[1]]$age_v4))) {
    stop("Matched cohort age_v4 contains missing or infinite values")
  }

  features <- lapply(names(aligned), function(om) {
    tibble::tibble(ID = ids, status = status) |>
      dplyr::bind_cols(aligned[[om]][, feature_cols[[om]], drop = FALSE])
  })
  names(features) <- names(aligned)

  list(
    ids = ids,
    n = length(ids),
    n_deaths = sum(status == "Deceased"),
    age = aligned[[1]]$age_v4,
    status = status,
    clinical = clinical,
    features = features,
    n_features = vapply(feature_cols, length, integer(1))
  )
}

# --- per-fold feature screening -------------------------------------------
#
# limma on the same median-imputed modelling scale as the primary analysis,
# Deceased versus Active, BH-adjusted P below `fdr`, and run on the analysis set
# only so the assessment fold plays no part in choosing features. The fallback
# to the `fallback_n` smallest adjusted P values is also shared with the primary
# implementation.
#
# The number retained varies substantially across folds in this cohort (roughly
# 4 to 94 proteins), because the matched participants carry a weaker univariate
# proteomic signal than the full proteomics sample. That variation is a real
# property of the threshold rule and is recorded per fold in the `tuning` element
# of the saved output.
screen_limma <- function(df, feature_cols,
                         fdr = MATCHED$limma_fdr,
                         fallback_n = MATCHED$limma_fallback_n) {
  limma_select_features_train(
    df_train = df, outcome = "status", feature_cols = feature_cols,
    limma_fdr = fdr, fallback_n = fallback_n
  )
}

# --- model specifications --------------------------------------------------
#
# Recipes, model families and tuning grids reproduce the primary analysis. The
# one deliberate change is the event level: the outcome factor is ordered with
# Deceased first, so tidymodels treats death as the event by default and the
# predicted probability of death is read directly from .pred_Deceased. The
# primary and matched analyses both store death risk directly.
matched_recipe <- function(df, kind, screen = FALSE, feature_cols = NULL) {
  rec <- recipes::recipe(status ~ ., data = df) |>
    recipes::update_role(ID, new_role = "id")
  if (isTRUE(screen)) {
    rec <- step_limma_select(
      rec, feature_cols = feature_cols, outcome = "status",
      limma_fdr = MATCHED$limma_fdr,
      fallback_n = MATCHED$limma_fallback_n
    )
  }
  rec <- rec |> recipes::step_impute_median(recipes::all_numeric_predictors())
  if (kind == "clinical") {
    rec <- rec |> recipes::step_impute_mode(recipes::all_nominal_predictors())
  }
  rec <- rec |> recipes::step_zv(recipes::all_predictors())
  if (kind == "en") {
    rec <- rec |> recipes::step_normalize(recipes::all_numeric_predictors())
  }
  rec
}

matched_spec <- function(kind) {
  switch(
    kind,
    clinical = parsnip::rand_forest(trees = tune(), mtry = tune(), min_n = tune()) |>
      parsnip::set_engine("ranger", probability = TRUE, num.threads = 1L) |>
      parsnip::set_mode("classification"),
    en = parsnip::logistic_reg(penalty = tune(), mixture = tune()) |>
      parsnip::set_engine("glmnet") |>
      parsnip::set_mode("classification"),
    xgb = parsnip::boost_tree(
      trees = tune(), tree_depth = tune(), learn_rate = tune(), mtry = tune(),
      min_n = tune(), loss_reduction = tune(), sample_size = tune(), stop_iter = 30
    ) |>
      parsnip::set_engine("xgboost", eval_metric = "auc", nthread = 1L) |>
      parsnip::set_mode("classification"),
    stop("unknown model kind: ", kind)
  )
}

matched_grid <- function(kind, rec, df, mtry_cap = NULL, n_analysis = NULL) {
  if (kind == "en") {
    return(dials::grid_regular(
      dials::penalty(range = c(-10, 0)),
      dials::mixture(range = c(0, 1)),
      levels = c(50, 11)
    ))
  }
  X <- recipes::bake(recipes::prep(rec, training = df), new_data = df) |>
    dplyr::select(-dplyr::any_of(c("ID", "status")))
  # mtry is bounded by the smallest post-recipe predictor count observed across
  # the inner analysis folds, so zero-variance removal cannot invalidate a grid.
  n_pred <- ncol(X)
  if (!is.null(mtry_cap)) n_pred <- min(n_pred, as.integer(mtry_cap))
  mtry_hi <- max(1L, n_pred)
  params <- if (kind == "clinical") {
    dials::parameters(
      dials::trees(range = c(300L, 2000L)),
      dials::mtry(range = c(1L, mtry_hi)),
      dials::min_n(range = c(2L, 40L))
    )
  } else {
    dials::parameters(
      dials::trees(range = c(500L, 4000L)),
      dials::tree_depth(range = c(2L, 4L)),
      dials::learn_rate(range = c(-2.5, -0.5)),
      dials::mtry(range = c(min(10L, mtry_hi), min(140L, mtry_hi))),
      # Bounded by the analysis set the model is tuned on. See
      # analysis_xgb_min_n_range() in 00_setup.R: min_n is xgboost's
      # min_child_weight, a sum of Hessian weights, so a value above roughly
      # n/16 makes every split impossible and the model predicts a constant.
      dials::min_n(range = analysis_xgb_min_n_range(
        if (is.null(n_analysis)) nrow(df) else n_analysis
      )),
      dials::loss_reduction(range = c(-4, 1)),
      sample_size = dials::sample_prop(range = c(0.5, 0.9))
    )
  }
  dials::grid_space_filling(params, size = nrow(params) * 10)
}

# --- one outer fold, one model --------------------------------------------
#
# Tune on an inner 3-fold split of the analysis set, fit the winning
# configuration on the whole analysis set, and predict the assessment fold. The
# assessment fold is touched exactly once, at the end.
tune_fit_predict <- function(analysis_df, assessment_df, kind, screen = FALSE,
                             feature_cols = NULL, resample_seed = MATCHED$seed) {
  set.seed(resample_seed)
  inner <- rsample::vfold_cv(analysis_df, v = MATCHED$inner_v, strata = status)
  rec <- matched_recipe(analysis_df, kind, screen = screen,
                        feature_cols = feature_cols)
  mtry_cap <- NULL
  if (kind %in% c("clinical", "xgb")) {
    mtry_cap <- min(vapply(inner$splits, function(sp) {
      inner_data <- rsample::analysis(sp)
      baked <- recipes::bake(
        recipes::prep(rec, training = inner_data), new_data = inner_data
      )
      ncol(dplyr::select(baked, -dplyr::any_of(c("ID", "status"))))
    }, integer(1)))
  }
  spec <- matched_spec(kind)
  set.seed(resample_seed + 1L)
  # The xgboost grid is bounded by the smallest inner analysis fold, because
  # that is the data each candidate is actually fitted on during tuning.
  inner_min_n <- min(vapply(inner$splits, function(sp) {
    nrow(rsample::analysis(sp))
  }, integer(1)))
  grid <- matched_grid(kind, rec, analysis_df, mtry_cap = mtry_cap,
                       n_analysis = inner_min_n)
  wf <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)

  tuned <- tune::tune_grid(
    wf, resamples = inner, grid = grid,
    metrics = yardstick::metric_set(yardstick::roc_auc),
    control = analysis_tune_control(save_pred = TRUE, verbose = FALSE)
  )
  best <- tune::select_best(tuned, metric = "roc_auc")
  inner_pred <- tune::collect_predictions(tuned, parameters = best)
  if (!all(c(".row", ".pred_Deceased") %in% names(inner_pred))) {
    stop(kind, ": inner tuning did not retain the predictions needed for stacking")
  }
  if (nrow(inner_pred) != nrow(analysis_df) || anyDuplicated(inner_pred$.row) ||
      !setequal(inner_pred$.row, seq_len(nrow(analysis_df)))) {
    stop(kind, ": inner out-of-fold predictions are incomplete")
  }
  inner_pred <- inner_pred[order(inner_pred$.row), , drop = FALSE]
  analysis_oof <- tibble::tibble(
    ID = as.character(analysis_df$ID[inner_pred$.row]),
    truth = as.character(analysis_df$status[inner_pred$.row]),
    risk = as.numeric(inner_pred$.pred_Deceased)
  )
  if (anyNA(analysis_oof$risk) || any(!is.finite(analysis_oof$risk)) ||
      any(analysis_oof$risk < 0 | analysis_oof$risk > 1)) {
    stop(kind, ": invalid inner out-of-fold probabilities")
  }
  final <- tune::finalize_workflow(wf, best) |> parsnip::fit(data = analysis_df)

  preds <- stats::predict(final, new_data = assessment_df, type = "prob")
  risk <- as.numeric(preds$.pred_Deceased)
  if (anyNA(risk) || any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
    stop(kind, ": invalid assessment-set predicted probabilities")
  }
  list(
    predictions = tibble::tibble(
      ID = as.character(assessment_df$ID),
      truth = as.character(assessment_df$status),
      risk = risk
    ),
    analysis_oof = analysis_oof,
    best = best,
    n_features = ncol(recipes::bake(recipes::prep(rec, training = analysis_df),
                                    new_data = analysis_df)) - 2L
  )
}

# --- one model across the whole nested design ------------------------------
run_matched_model <- function(dat, kind, label, omics, screen = FALSE,
                              feature_cols = NULL) {
  set.seed(MATCHED$seed)
  outer <- rsample::vfold_cv(dat, v = MATCHED$outer_v,
                             repeats = MATCHED$n_repeats, strata = status)
  preds <- vector("list", nrow(outer))
  stack_training <- vector("list", nrow(outer))
  meta <- vector("list", nrow(outer))

  for (i in seq_len(nrow(outer))) {
    sp <- outer$splits[[i]]
    an <- rsample::analysis(sp)
    as_ <- rsample::assessment(sp)

    if (screen) sel <- screen_limma(an, feature_cols)

    # rsample only adds `id2` when repeats > 1; with a single repeat the fold
    # label lives in `id`.
    repeated <- "id2" %in% names(outer)
    rep_lab <- if (repeated) outer$id[[i]] else "Repeat1"
    fold_lab <- if (repeated) outer$id2[[i]] else outer$id[[i]]

    out <- tune_fit_predict(
      an, as_, kind, screen = screen, feature_cols = feature_cols,
      resample_seed = MATCHED$seed + i
    )
    preds[[i]] <- out$predictions |>
      dplyr::mutate(repeat_id = rep_lab, fold_id = fold_lab)
    stack_training[[i]] <- out$analysis_oof |>
      dplyr::mutate(repeat_id = rep_lab, fold_id = fold_lab)
    meta[[i]] <- tibble::tibble(
      repeat_id = rep_lab, fold_id = fold_lab,
      n_features = if (screen) length(sel) else out$n_features,
      selected_features = list(if (screen) sel else character())
    ) |> dplyr::bind_cols(out$best |> dplyr::select(-dplyr::any_of(".config")))
    message(sprintf("    %-22s %s / %s  (%d features)", label,
                    rep_lab, fold_lab, meta[[i]]$n_features))
  }

  by_fold <- dplyr::bind_rows(preds)
  if (nrow(by_fold) != nrow(dat) * MATCHED$n_repeats) {
    stop(label, ": repeated outer-fold predictions are incomplete")
  }
  if (anyDuplicated(paste(by_fold$repeat_id, by_fold$ID, sep = "|"))) {
    stop(label, ": duplicate participant within a repeat")
  }
  expected_ids <- sort(as.character(dat$ID))
  for (rp in unique(by_fold$repeat_id)) {
    got <- sort(as.character(by_fold$ID[by_fold$repeat_id == rp]))
    if (!identical(got, expected_ids)) stop(label, ": a repeat does not cover every participant")
  }
  # Each repeat covers all participants exactly once, so a repeat-level
  # prediction vector is complete and averaging across repeats is well defined.
  averaged <- by_fold |>
    dplyr::group_by(ID, truth) |>
    dplyr::summarise(risk = mean(risk), .groups = "drop") |>
    dplyr::mutate(model = label, omics = omics)
  if (nrow(averaged) != nrow(dat) || anyDuplicated(as.character(averaged$ID))) {
    stop(label, ": averaged predictions do not contain one row per participant")
  }

  list(
    predictions = averaged,
    predictions_by_repeat = by_fold |> dplyr::mutate(model = label, omics = omics),
    stack_training = dplyr::bind_rows(stack_training) |>
      dplyr::mutate(model = label, omics = omics),
    tuning = dplyr::bind_rows(meta) |> dplyr::mutate(model = label, omics = omics)
  )
}

# Age-only logistic regression evaluated on the same repeated outer folds as
# the other matched-cohort models. The fitted probability of death is stored,
# so the LR label and the [0, 1] risk contract are both literal.
run_matched_age_logistic <- function(dat) {
  set.seed(MATCHED$seed)
  outer <- rsample::vfold_cv(
    dat, v = MATCHED$outer_v, repeats = MATCHED$n_repeats, strata = status
  )
  preds <- vector("list", nrow(outer))
  tuning <- vector("list", nrow(outer))
  for (i in seq_len(nrow(outer))) {
    an <- rsample::analysis(outer$splits[[i]]) |>
      dplyr::mutate(death = as.integer(status == "Deceased"))
    as_ <- rsample::assessment(outer$splits[[i]])
    fit <- stats::glm(death ~ age_v4, data = an, family = stats::binomial())
    risk <- as.numeric(stats::predict(fit, newdata = as_, type = "response"))
    if (anyNA(risk) || any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
      stop("LR age-only: invalid outer-assessment probabilities")
    }
    repeated <- "id2" %in% names(outer)
    rep_lab <- if (repeated) outer$id[[i]] else "Repeat1"
    fold_lab <- if (repeated) outer$id2[[i]] else outer$id[[i]]
    preds[[i]] <- tibble::tibble(
      ID = as.character(as_$ID), truth = as.character(as_$status), risk = risk,
      repeat_id = rep_lab, fold_id = fold_lab
    )
    cf <- stats::coef(fit)
    tuning[[i]] <- tibble::tibble(
      repeat_id = rep_lab, fold_id = fold_lab,
      intercept = unname(cf[[1L]]), age_coefficient = unname(cf[[2L]])
    )
  }
  by_repeat <- dplyr::bind_rows(preds) |>
    dplyr::mutate(model = "LR age-only", omics = "clinical")
  if (nrow(by_repeat) != nrow(dat) * MATCHED$n_repeats ||
      anyDuplicated(paste(by_repeat$repeat_id, by_repeat$ID, sep = "|"))) {
    stop("LR age-only: repeated outer-fold predictions are incomplete")
  }
  averaged <- by_repeat |>
    dplyr::group_by(ID, truth) |>
    dplyr::summarise(risk = mean(risk), .groups = "drop") |>
    dplyr::mutate(model = "LR age-only", omics = "clinical")
  list(
    predictions = averaged,
    predictions_by_repeat = by_repeat,
    tuning = dplyr::bind_rows(tuning) |>
      dplyr::mutate(model = "LR age-only", omics = "clinical")
  )
}

# Fit a ridge-logistic stack inside every matched outer analysis set. The
# meta-training predictors are inner out-of-fold probabilities produced using
# only that outer analysis set. The corresponding outer assessment fold is
# used once, after the base models and ridge penalty have been selected.
run_matched_stack <- function(clinical_result, omics_result, label, omics) {
  outer_clin <- clinical_result$predictions_by_repeat
  outer_omics <- omics_result$predictions_by_repeat
  inner_clin <- clinical_result$stack_training
  inner_omics <- omics_result$stack_training
  fold_keys <- unique(outer_clin[c("repeat_id", "fold_id")])
  preds <- vector("list", nrow(fold_keys))
  tuning <- vector("list", nrow(fold_keys))

  for (i in seq_len(nrow(fold_keys))) {
    rp <- fold_keys$repeat_id[[i]]
    fd <- fold_keys$fold_id[[i]]
    in_fold <- function(d) d$repeat_id == rp & d$fold_id == fd

    tr_clin <- inner_clin[in_fold(inner_clin), c("ID", "truth", "risk")]
    tr_omics <- inner_omics[in_fold(inner_omics), c("ID", "truth", "risk")]
    meta_train <- dplyr::inner_join(
      tr_clin, tr_omics, by = "ID", suffix = c("_clinical", "_omics")
    )
    if (nrow(meta_train) != nrow(tr_clin) || nrow(meta_train) != nrow(tr_omics) ||
        anyDuplicated(meta_train$ID) ||
        !identical(meta_train$truth_clinical, meta_train$truth_omics)) {
      stop(label, " / ", omics, ": invalid outer-analysis stacking data")
    }
    meta_train <- meta_train |>
      dplyr::transmute(
        truth = factor(truth_clinical, levels = c("Deceased", "Active")),
        z_clinical = matched_logit(risk_clinical),
        z_omics = matched_logit(risk_omics)
      )

    as_clin <- outer_clin[in_fold(outer_clin), c("ID", "truth", "risk")]
    as_omics <- outer_omics[in_fold(outer_omics), c("ID", "truth", "risk")]
    meta_assess <- dplyr::inner_join(
      as_clin, as_omics, by = "ID", suffix = c("_clinical", "_omics")
    )
    if (nrow(meta_assess) != nrow(as_clin) || nrow(meta_assess) != nrow(as_omics) ||
        anyDuplicated(meta_assess$ID) ||
        !identical(meta_assess$truth_clinical, meta_assess$truth_omics)) {
      stop(label, " / ", omics, ": invalid outer-assessment stacking data")
    }
    assess_x <- meta_assess |>
      dplyr::transmute(
        z_clinical = matched_logit(risk_clinical),
        z_omics = matched_logit(risk_omics)
      )

    spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 0) |>
      parsnip::set_engine("glmnet") |>
      parsnip::set_mode("classification")
    wf <- workflows::workflow() |>
      workflows::add_model(spec) |>
      workflows::add_formula(truth ~ z_clinical + z_omics)
    set.seed(MATCHED$seed + i)
    inner <- rsample::vfold_cv(meta_train, v = MATCHED$inner_v, strata = truth)
    tuned <- tune::tune_grid(
      wf, resamples = inner,
      grid = dials::grid_regular(dials::penalty(range = c(-10, 2)), levels = 60),
      metrics = yardstick::metric_set(yardstick::roc_auc),
      control = analysis_tune_control(save_pred = FALSE, verbose = FALSE)
    )
    best <- tune::select_best(tuned, metric = "roc_auc")
    fit <- tune::finalize_workflow(wf, best) |> parsnip::fit(data = meta_train)
    risk <- as.numeric(stats::predict(fit, new_data = assess_x, type = "prob")$.pred_Deceased)
    if (anyNA(risk) || any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
      stop(label, " / ", omics, ": invalid outer-assessment stack probabilities")
    }
    preds[[i]] <- tibble::tibble(
      ID = as.character(meta_assess$ID),
      truth = as.character(meta_assess$truth_clinical), risk = risk,
      repeat_id = rp, fold_id = fd, model = label, omics = omics
    )
    tuning[[i]] <- best |>
      dplyr::select(-dplyr::any_of(".config")) |>
      dplyr::mutate(repeat_id = rp, fold_id = fd, model = label, omics = omics)
  }

  by_repeat <- dplyr::bind_rows(preds)
  expected_n <- nrow(clinical_result$predictions) * MATCHED$n_repeats
  if (nrow(by_repeat) != expected_n ||
      anyDuplicated(paste(by_repeat$repeat_id, by_repeat$ID, sep = "|"))) {
    stop(label, " / ", omics, ": repeated outer-fold stack predictions are incomplete")
  }
  averaged <- by_repeat |>
    dplyr::group_by(ID, truth) |>
    dplyr::summarise(risk = mean(risk), .groups = "drop") |>
    dplyr::mutate(model = label, omics = omics)
  list(
    predictions = averaged,
    predictions_by_repeat = by_repeat,
    stack_training = NULL,
    tuning = dplyr::bind_rows(tuning)
  )
}

# --- driver ----------------------------------------------------------------
assert_packages(matched_packages, "matched-cohort refit")
message("Building matched cohort ...")
# The re-fitted layers, not the primary ones: they carry the provenance
# manifests the matched stage records against, and they store risk as the
# predicted death probability throughout. See 00_setup.R on the two paths.
model_outputs <- load_prediction_model_outputs(
  paths,
  dir = paths$model_outputs_refit,
  require_provenance = TRUE,
  expect_score = "deceased"
)
primary_source_files <- attr(model_outputs, "source_files")
cohort <- build_matched_cohort(model_outputs)
rm(model_outputs); invisible(gc())
message(sprintf("  n = %d, deaths = %d, features = %s",
                cohort$n, cohort$n_deaths,
                paste(sprintf("%s %d", names(cohort$n_features), cohort$n_features),
                      collapse = ", ")))

# Forked workers, not multisession: the proteomics models screen features
# through step_limma_select(), whose S3 methods are defined in the global
# environment and cannot reach a separate R process. See
# analysis_parallel_plan() in 00_setup.R.
if (requireNamespace("future", quietly = TRUE)) {
  old_future_plan <- future::plan()
  matched_plan_label <- analysis_parallel_plan(MATCHED$n_workers)
  # Attach before any future is created, so forked workers inherit these rather
  # than resolving them from the shared filesystem. See 00_setup.R.
  matched_preloaded <- analysis_preload_worker_packages()
  message("[PARALLEL] plan = ", matched_plan_label,
          " | workers = ", MATCHED$n_workers,
          " | preloaded = ", paste(matched_preloaded, collapse = ", "))
}

matched_code_files <- file.path(project_root, c(
  "code/single_omics/00_setup.R", "code/single_omics/analysis_provenance.R",
  "code/single_omics/feature_selection.R", "code/single_omics/helpers.R",
  "code/single_omics/06_fit_matched_cohort_models.R",
  "code/single_omics/07_make_matched_cohort_tables.R"
))
matched_provenance <- make_provenance(
  stage = "matched_cohort",
  inputs = primary_source_files,
  code_files = matched_code_files,
  config = MATCHED,
  seed = MATCHED$seed,
  packages = matched_packages
)
write_manifest_csv(matched_provenance,
                   file.path(matched_dir, "matched_cohort_provenance_manifest.csv"))

results <- list()

# --- checkpointing ---------------------------------------------------------
#
# Each job is saved as soon as it finishes, so a run that dies partway through
# -- or one moved to a cluster -- can pick up where it stopped. The resampling
# design is stored alongside the result and checked on reload: a checkpoint
# written under a different number of repeats or folds is not interchangeable
# with one written under the current settings, and silently reusing it would
# mix two designs in one table.
checkpoint_path <- function(key) file.path(matched_dir, sprintf("matched_%s.rds", key))

save_checkpoint <- function(key, obj) {
  obj$design <- MATCHED[c("seed", "outer_v", "inner_v", "n_repeats", "limma_fdr",
                         "limma_fallback_n", "screened_layers")]
  obj$provenance <- matched_provenance
  atomic_save_rds(obj, checkpoint_path(key))
  obj
}

load_checkpoint <- function(key) {
  f <- checkpoint_path(key)
  if (!isTRUE(MATCHED$resume) || !file.exists(f)) return(NULL)
  obj <- tryCatch(readRDS(f), error = function(e) NULL)
  if (is.null(obj)) return(NULL)
  want <- MATCHED[c("seed", "outer_v", "inner_v", "n_repeats", "limma_fdr",
                    "limma_fallback_n", "screened_layers")]
  if (!identical(obj$design, want) || !provenance_matches(obj, matched_provenance)) {
    message(sprintf("  checkpoint for %s was written under a different design; refitting", key))
    return(NULL)
  }
  message(sprintf("  reusing checkpoint for %s", key))
  obj
}

# glmnet and xgboost emit convergence and tuning warnings routinely. They are
# usually benign, but a public repository should record them rather than let
# thousands scroll past unread, so they are collected and saved with the output.
matched_warnings <- character()
collect_warnings <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    matched_warnings <<- c(matched_warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
}

# The clinical model is fitted once, not once per layer. This is the point of
# the redesign: a single clinical benchmark with a single AUC.
message("Fitting clinical model (RF, 14 covariates) ...")
results$clinical <- load_checkpoint("clinical")
if (is.null(results$clinical)) {
  t0 <- Sys.time()
  results$clinical <- collect_warnings(
    run_matched_model(cohort$clinical, "clinical", "RF clinical", "clinical")
  )
  message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  results$clinical <- save_checkpoint("clinical", results$clinical)
}

# Fit age-only logistic regression in every outer analysis set and generate
# repeated out-of-fold death probabilities on the matched cohort.
results$age_only <- run_matched_age_logistic(
  cohort$clinical |> dplyr::select(ID, status, age_v4)
)

# One job per layer and algorithm. Only proteomics is screened.
for (om in names(cohort$features)) {
  dat <- cohort$features[[om]]
  feat <- setdiff(names(dat), c("ID", "status"))
  screen <- om %in% MATCHED$screened_layers
  for (kind in c("en", "xgb")) {
    label <- if (kind == "en") "EN omics" else "XGB omics"
    key <- paste(om, kind, sep = "_")
    message(sprintf("Fitting %s / %s ...", om, label))
    results[[key]] <- load_checkpoint(key)
    if (is.null(results[[key]])) {
      t0 <- Sys.time()
      results[[key]] <- collect_warnings(run_matched_model(
        dat, kind, label, om, screen = screen, feature_cols = feat
      ))
      message(sprintf("  done in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
      results[[key]] <- save_checkpoint(key, results[[key]])
    }
  }
}

# Construct the two prespecified stacks for each omics layer. These are fitted
# here—not during table rendering—because stacking is part of the statistical
# training procedure and must respect the matched outer folds.
for (om in names(cohort$features)) {
  results[[paste(om, "stack_en", sep = "_")]] <- run_matched_stack(
    results$clinical, results[[paste(om, "en", sep = "_")]],
    label = "RF clinical + EN", omics = om
  )
  results[[paste(om, "stack_xgb", sep = "_")]] <- run_matched_stack(
    results$clinical, results[[paste(om, "xgb", sep = "_")]],
    label = "RF clinical + XGB", omics = om
  )
}

# --- repeat-level AUC ------------------------------------------------------
#
# `predictions` averages each participant's risk across the five repeats, and
# the headline AUC is computed from that averaged vector. That is the AUC of a
# repeat-averaged predictor: averaging cancels fold-assignment noise, so it sits
# systematically above the mean of the per-repeat AUCs and is not what a reader
# usually understands by "cross-validated AUC".
#
# Both quantities are therefore recorded. The averaged AUC (computed downstream
# in 07) is unchanged; these additions report the AUC of each repeat separately
# plus their mean and spread, so the gap between the two is visible rather than
# implicit. Purely additive -- no existing output changes.
matched_auc_by_repeat <- function(by_repeat) {
  if (is.null(by_repeat) || !nrow(by_repeat)) {
    return(tibble::tibble(
      model = character(), omics = character(), repeat_id = character(),
      n = integer(), n_deaths = integer(), roc_auc = numeric()
    ))
  }
  by_repeat |>
    dplyr::mutate(truth = factor(truth, levels = c("Deceased", "Active"))) |>
    dplyr::group_by(model, omics, repeat_id) |>
    dplyr::summarise(
      n = dplyr::n(),
      n_deaths = sum(truth == "Deceased"),
      roc_auc = as.numeric(
        yardstick::roc_auc_vec(truth = truth, estimate = risk,
                               event_level = "first")
      ),
      .groups = "drop"
    )
}

matched_auc_repeat_table <- matched_auc_by_repeat(
  dplyr::bind_rows(lapply(results, `[[`, "predictions_by_repeat"))
)
matched_auc_repeat_summary <- matched_auc_repeat_table |>
  dplyr::group_by(model, omics) |>
  dplyr::summarise(
    n_repeats = dplyr::n(),
    mean_roc_auc = mean(roc_auc),
    sd_roc_auc = stats::sd(roc_auc),
    min_roc_auc = min(roc_auc),
    max_roc_auc = max(roc_auc),
    .groups = "drop"
  )

matched_cohort_outputs <- list(
  cohort = cohort[c("ids", "n", "n_deaths", "status", "n_features")],
  design = MATCHED,
  predictions = dplyr::bind_rows(lapply(results, `[[`, "predictions")),
  predictions_by_repeat = dplyr::bind_rows(lapply(results, `[[`, "predictions_by_repeat")),
  auc_by_repeat = matched_auc_repeat_table,
  auc_by_repeat_summary = matched_auc_repeat_summary,
  tuning = dplyr::bind_rows(lapply(results, `[[`, "tuning")),
  warnings = if (length(matched_warnings)) {
    as.data.frame(table(matched_warnings), stringsAsFactors = FALSE) |>
      stats::setNames(c("warning", "n")) |>
      dplyr::arrange(dplyr::desc(n))
  } else {
    data.frame(warning = character(), n = integer())
  },
  session_info = utils::sessionInfo()
)

matched_cohort_outputs$provenance <- matched_provenance
atomic_save_rds(
  matched_cohort_outputs,
  file.path(matched_dir, "matched_cohort_model_outputs.rds")
)
if (exists("old_future_plan", inherits = FALSE)) future::plan(old_future_plan)
message("Saved -> ", file.path(matched_dir, "matched_cohort_model_outputs.rds"))
