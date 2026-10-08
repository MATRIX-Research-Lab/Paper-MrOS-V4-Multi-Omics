# Single-omics prediction pipeline (RF / EN / XGB / stacks + FS)
# Formerly omics_pipeline_withFS_XGB100_v2.R

logit <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

library(tidyverse)
library(tidymodels)
library(kernelshap)
library(future)
library(limma)
library(mRMRe)
library(xgboost)

run_omics_ensemble_pipeline = function(
    df,
    outcome,
    id_col = "ID",
    event_level = "second",
    event_class = NULL,
    event_prob_col = NULL,
    n_workers = 10,
    seed = 123,
    cov_cols = NULL,
    feature_cols = NULL,
    omics_label = "", # "Microbiome", "Metabolomics", "Proteomics"
    fs_method = c("none", "limma", "mrmr", "sis"),   
    limma_fdr = 0.05,
    mrmr_top = 50,
    sis_top = NULL,              # default uses n/log(n) inside each training split
    sis_cap = 200,               # optional upper cap, if sis_top is too large
    sis_method = c("spearman", "pearson")  # "spearman" or "pearson"
) {
  
  # ---- Basic checks ----
  if (!id_col %in% names(df)) stop("id_col not found in df: ", id_col)
  if (!outcome %in% names(df)) stop("outcome not found in df: ", outcome)
  
  if (!is.factor(df[[outcome]])) df[[outcome]] = as.factor(df[[outcome]])
  if (nlevels(df[[outcome]]) < 2) stop("Outcome must have at least 2 levels.")
  
  if (is.null(event_class)) {
    event_level = "second"
    event_class = levels(df[[outcome]])[2]
  }
  if (!event_class %in% levels(df[[outcome]])) stop("event_class not in levels(df[[outcome]]).")
  
  if (is.null(event_prob_col)) event_prob_col = paste0(".pred_", event_class)
  
  if (is.null(cov_cols)) {
    if (exists("cov_cols", envir = parent.frame(), inherits = TRUE)) {
      cov_cols <- get("cov_cols", envir = parent.frame(), inherits = TRUE)
    } else {
      stop("`cov_cols` not provided and not found in the calling environment.")
    }
  }
  if (is.null(feature_cols)) {
    if (exists("feature_cols", envir = parent.frame(), inherits = TRUE)) {
      feature_cols <- get("feature_cols", envir = parent.frame(), inherits = TRUE)
    } else {
      stop("`feature_cols` not provided and not found in the calling environment.")
    }
  }
  
  # ---- Seed and parallel plan ----
  set.seed(seed)
  old_plan = future::plan()
  on.exit({
    future::plan(old_plan)
  }, add = TRUE)
  
  future::plan(future::multisession, workers = n_workers)
  
  # Def emit
  emit <- function(...) {
    cat(..., "\n", sep = "")
    flush.console()
  }
  
  # ----------------------------
  # 1) Train/test split + CV folds
  # ----------------------------
  split = rsample::initial_split(df, prop = 0.8, strata = all_of(outcome))
  train = rsample::training(split)
  test  = rsample::testing(split)
  
  # NEW: create stable row index inside train for OOF join
  # NEW: tidymodels tune_grid() returns ".row" which is row position in training set
  # NEW: we will join OOF by ".train_row" (same meaning)
  train <- train %>% dplyr::mutate(.train_row = dplyr::row_number())
  
  folds = rsample::vfold_cv(train, v = 5, strata = all_of(outcome))
  auc_set = yardstick::metric_set(yardstick::roc_auc)
  
  # --------------------------
  # 2) MODEL 1: Covariates RF
  # --------------------------
  train_cov = train %>% dplyr::select(all_of(c(id_col, outcome, cov_cols)))
  test_cov  = test  %>% dplyr::select(all_of(c(id_col, outcome, cov_cols)))
  
  rec_rf_cov = recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = train_cov) %>%
    recipes::update_role(all_of(id_col), new_role = "id") %>%
    recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
    recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    step_dummy(all_nominal_predictors())
  
  rf_cov_spec = parsnip::rand_forest(
    trees = tune(),
    mtry  = tune(),
    min_n = tune()
  ) %>%
    parsnip::set_engine("ranger", probability = TRUE) %>%
    parsnip::set_mode("classification")
  
  wf_rf_cov = workflows::workflow() %>% workflows::add_recipe(rec_rf_cov) %>% workflows::add_model(rf_cov_spec)
  
  prep_cov  = recipes::prep(rec_rf_cov, training = train_cov, verbose = FALSE)
  baked_cov = recipes::bake(prep_cov, new_data = train_cov)
  X_cov     = baked_cov %>% dplyr::select(-all_of(c(id_col, outcome)))
  
  params_cov = dials::parameters(
    dials::trees(range = c(300L, 2000L)),
    dials::finalize(dials::mtry(), X_cov),
    dials::min_n(range = c(2L, 40L))
  )
  
  grid_cov = dials::grid_space_filling(params_cov, size = nrow(params_cov) * 10)
  
  set.seed(seed)
  tuned_cov = tune::tune_grid(
    wf_rf_cov,
    resamples = folds,
    grid = grid_cov,
    metrics = auc_set,
    control = tune::control_grid(save_pred = TRUE)
  )
  
  best_cov = tune::select_best(tuned_cov, metric = "roc_auc")
  wf_cov_final = tune::finalize_workflow(wf_rf_cov, best_cov)
  
  # OOF predictions for stacking (Model 1)
  pred_oof_cov = tune::collect_predictions(tuned_cov, parameters = best_cov) %>%
    dplyr::transmute(.train_row = .row, truth = .data[[outcome]], p1 = .data[[event_prob_col]])  # NEW: keep stable key
  
  # Fit final on full training and predict test (Model 1)
  fit_cov = parsnip::fit(wf_cov_final, data = train_cov)
  
  prob_test_cov = predict(fit_cov, new_data = test_cov, type = "prob") %>%
    dplyr::bind_cols(test_cov %>% dplyr::select(all_of(c(id_col, outcome))) %>% dplyr::rename(truth = all_of(outcome))) %>%
    dplyr::mutate(p1 = .data[[event_prob_col]])
  
  auc_test_1 = yardstick::roc_auc(prob_test_cov, truth = truth, p1, event_level = event_level)
  
  # --------------------------
  # 3) MODEL 2: Microbiome Elastic Net (glmnet)
  # --------------------------
  train_micro = train %>% dplyr::select(all_of(c(id_col, outcome, feature_cols, ".train_row")))  # NEW: keep .train_row
  test_micro  = test  %>% dplyr::select(all_of(c(id_col, outcome, feature_cols)))
  
  # --------------------------
  # Optional feature selection on OMICS only (train-only)
  # --------------------------
  fs_method <- match.arg(fs_method)
  sis_method <- match.arg(sis_method)
  
  # NEW: helper functions for FS (limma / mrmr/ SIS) to be called inside each CV fold
  select_features_limma <- function(df_train, outcome, feature_cols, limma_fdr) {
    y <- droplevels(df_train[[outcome]])
    if (!is.factor(y) || nlevels(y) != 2) {
      stop("limma feature selection currently supports binary factor outcome only.")
    }
    
    X <- df_train %>% dplyr::select(all_of(feature_cols))
    
    # median impute
    X_imp <- X
    for (j in seq_len(ncol(X_imp))) {
      v <- X_imp[[j]]
      if (anyNA(v)) X_imp[[j]][is.na(v)] <- median(v, na.rm = TRUE)
    }
    
    mat <- t(as.matrix(X_imp))  # features x samples
    design <- stats::model.matrix(~ y)
    fit <- limma::lmFit(mat, design)
    fit <- limma::eBayes(fit)
    tab <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "P")
    
    limma_keep <- rownames(tab)[tab$adj.P.Val < limma_fdr]
    if (length(limma_keep) == 0) {
      # fallback: keep top 50 by adj.P.Val
      limma_keep <- rownames(tab)[seq_len(min(50, nrow(tab)))]
    }
    limma_keep
  }
  
  select_features_mrmr <- function(df_train, outcome, feature_cols, mrmr_top) {
    y <- droplevels(df_train[[outcome]])
    if (!is.factor(y) || nlevels(y) != 2) {
      stop("mrmr feature selection currently supports binary factor outcome only.")
    }
    
    X <- df_train %>% dplyr::select(all_of(feature_cols))
    
    # median impute
    X_imp <- X
    for (j in seq_len(ncol(X_imp))) {
      v <- X_imp[[j]]
      if (anyNA(v)) X_imp[[j]][is.na(v)] <- median(v, na.rm = TRUE)
    }
    
    # outcome must be numeric for mRMRe
    y_num <- as.numeric(y)
    
    k_target <- min(mrmr_top, length(feature_cols))
    if (k_target <= 0) stop("mrmr_top must be >= 1")
    if (length(feature_cols) <= k_target) return(feature_cols)
    
    # IMPORTANT: map sanitized names back to original names
    name_map <- setNames(feature_cols, make.names(feature_cols))
    
    df_mrmr <- data.frame(y = y_num, X_imp)  # may sanitize names internally
    dd <- mRMRe::mRMR.data(data = df_mrmr)
    
    res <- mRMRe::mRMR.classic(data = dd, target_indices = 1, feature_count = k_target)
    idx <- mRMRe::solutions(res)[[1]]
    
    selected_sanitized <- setdiff(mRMRe::featureNames(dd)[idx], "y")
    selected_feature_cols <- unname(name_map[selected_sanitized])
    selected_feature_cols <- selected_feature_cols[!is.na(selected_feature_cols)]
    
    if (length(selected_feature_cols) == 0) feature_cols else selected_feature_cols
  }
  
  select_features_sis <- function(df_train, outcome, feature_cols,
                                  sis_top = NULL, sis_cap = 200,
                                  sis_method = "spearman") {
    y <- droplevels(df_train[[outcome]])
    if (!is.factor(y) || nlevels(y) != 2) {
      stop("SIS feature selection currently supports binary factor outcome only.")
    }
    
    n <- nrow(df_train)
    if (is.null(sis_top)) sis_top <- floor(n / log(n))
    if (!is.null(sis_cap) && is.finite(sis_cap)) sis_top <- min(sis_top, sis_cap)
    sis_top <- max(1, sis_top)
    k_target <- min(sis_top, length(feature_cols))
    
    y01 <- as.numeric(y == levels(y)[2])
    
    X <- df_train %>% dplyr::select(all_of(feature_cols))
    
    # median impute
    X_imp <- X
    for (j in seq_len(ncol(X_imp))) {
      v <- X_imp[[j]]
      if (anyNA(v)) X_imp[[j]][is.na(v)] <- median(v, na.rm = TRUE)
    }
    
    scores <- vapply(X_imp, function(x) {
      suppressWarnings(abs(stats::cor(x, y01, method = sis_method, use = "complete.obs")))
    }, numeric(1))
    
    scores[is.na(scores)] <- 0
    
    ord <- order(scores, decreasing = TRUE)
    feature_cols[ord][seq_len(k_target)]
  }
  
  # NEW: Route B implementation (outer CV loop with FS inside each fold)
  # NEW: inside each outer fold, do inner CV tuning on analysis only, then predict outer assessment
  inner_v <- 3
  
  pred_oof_micro_en_list <- list()  # NEW: store OOF per fold
  pred_oof_micro_xgb_list <- list() # NEW: store OOF per fold
  selected_features_by_fold <- list() # NEW: record selected features per outer fold (optional debug)
  
  for (k in seq_along(folds$splits)) {
    sp <- folds$splits[[k]]
    
    analysis_df <- rsample::analysis(sp)
    assess_df   <- rsample::assessment(sp)
    
    # NEW: FS on OUTER analysis only
    selected_feature_cols <- feature_cols
    if (fs_method == "limma") {
      selected_feature_cols <- select_features_limma(
        df_train = analysis_df,
        outcome = outcome,
        feature_cols = feature_cols,
        limma_fdr = limma_fdr
      )
    }
    if (fs_method == "mrmr") {
      selected_feature_cols <- select_features_mrmr(
        df_train = analysis_df,
        outcome = outcome,
        feature_cols = feature_cols,
        mrmr_top = mrmr_top
      )
    }
    if (fs_method == "sis") {
      selected_feature_cols <- select_features_sis(
        df_train = analysis_df,
        outcome = outcome,
        feature_cols = feature_cols,
        sis_top = sis_top,
        sis_cap = sis_cap,
        sis_method = sis_method
      )
    }
    # FS list for each fold
    selected_features_by_fold[[k]] <- list(
      fold = k,
      analysis_n = nrow(analysis_df),
      assess_n = nrow(assess_df),
      n_features = length(selected_feature_cols),
      features = selected_feature_cols
    )
    
    emit(
      "[FS] outer fold ", k,
      " | analysis n=", nrow(analysis_df),
      " | assess n=", nrow(assess_df),
      " | selected features=", length(selected_feature_cols)
    )
    
    # NEW: build fold-specific train/val with selected features
    analysis_micro <- analysis_df %>%
      dplyr::select(all_of(c(id_col, outcome, ".train_row", selected_feature_cols)))
    assess_micro <- assess_df %>%
      dplyr::select(all_of(c(id_col, outcome, ".train_row", selected_feature_cols)))
    
    # NEW: inner folds for tuning (analysis only)
    inner_folds <- rsample::vfold_cv(analysis_micro, v = inner_v, strata = all_of(outcome))
    
    # NEW: Model 2 (EN) tuned on inner folds, fit on analysis_micro, predict assess_micro
    rec_en_micro = recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = analysis_micro) %>%
      recipes::update_role(all_of(id_col), new_role = "id") %>%
      recipes::update_role(all_of(".train_row"), new_role = "id") %>%  # NEW: prevent leakage into predictors
      recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
      recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
      recipes::step_zv(recipes::all_nominal_predictors()) %>%
      recipes::step_dummy(recipes::all_nominal_predictors()) %>%
      recipes::step_zv(recipes::all_predictors()) %>%
      recipes::step_normalize(recipes::all_numeric_predictors())
    
    en_spec = parsnip::logistic_reg(
      penalty = tune(),
      mixture = tune()
    ) %>%
      parsnip::set_engine("glmnet") %>%
      parsnip::set_mode("classification")
    
    wf_micro_en = workflows::workflow() %>% workflows::add_recipe(rec_en_micro) %>% workflows::add_model(en_spec)
    
    en_grid = dials::grid_regular(
      dials::penalty(range = c(-10, 0)),
      dials::mixture(range = c(0, 1)),
      levels = c(50, 11)
    )
    
    set.seed(seed + k)  # NEW: stable but fold-different seed
    tuned_micro_en_k = tune::tune_grid(
      wf_micro_en,
      resamples = inner_folds,
      grid = en_grid,
      metrics = auc_set,
      control = tune::control_grid(save_pred = FALSE)
    )
    
    notes_en_k <- tune::collect_notes(tuned_micro_en_k)
    if (nrow(notes_en_k) > 0) {
      emit(
        "[EN tune notes] outer fold ", k,
        " | showing up to 25 rows (see .err for full tibble if needed)"
      )
      print(utils::head(notes_en_k, 25))
    }
    
    metrics_en_k <- tune::collect_metrics(tuned_micro_en_k)
    if (!any(metrics_en_k$.metric == "roc_auc" & is.finite(metrics_en_k$mean))) {
      emit("[EN tune failure] outer fold ", k, " | no finite roc_auc; fallback to fixed EN hyperparameters")
      if (nrow(notes_en_k) > 0) {
        emit("[EN tune failure] outer fold ", k, " | full notes:")
        print(notes_en_k)
      }
      best_micro_en_k <- tibble::tibble(penalty = 0.01, mixture = 0.5)
    } else {
      best_micro_en_k = tune::select_best(tuned_micro_en_k, metric = "roc_auc")
    }
    wf_micro_en_final_k = tune::finalize_workflow(wf_micro_en, best_micro_en_k)
    
    fit_micro_en_k = parsnip::fit(wf_micro_en_final_k, data = analysis_micro)
    
    prob_assess_en_k = predict(fit_micro_en_k, new_data = assess_micro, type = "prob") %>%
      dplyr::bind_cols(assess_micro %>% dplyr::select(all_of(".train_row"))) %>%  # NEW: only key
      dplyr::mutate(p2 = .data[[event_prob_col]]) %>%
      dplyr::select(.train_row, p2) 
    
    pred_oof_micro_en_list[[k]] <- prob_assess_en_k
    
    # NEW: Model 3 (XGB) tuned on inner folds, fit on analysis_micro, predict assess_micro
    rec_xgb_micro = recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = analysis_micro) %>%
      recipes::update_role(all_of(id_col), new_role = "id") %>%
      recipes::update_role(all_of(".train_row"), new_role = "id") %>%  # NEW
      recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
      recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
      recipes::step_zv(recipes::all_nominal_predictors()) %>%
      recipes::step_dummy(recipes::all_nominal_predictors()) %>%
      recipes::step_zv(recipes::all_predictors())
    
    xgb_micro_spec = parsnip::boost_tree(
      trees = tune(),
      tree_depth = tune(),
      learn_rate = tune(),
      mtry = tune(),
      min_n = tune(),
      loss_reduction = tune(),
      sample_size = tune(),
      stop_iter = 30
    ) %>%
      parsnip::set_engine("xgboost", eval_metric = "auc", event_level = event_level) %>%
      parsnip::set_mode("classification")
    
    wf_micro_xgb = workflows::workflow() %>% workflows::add_recipe(rec_xgb_micro) %>% workflows::add_model(xgb_micro_spec)
    
    # NEW: finalize mtry based on analysis_micro after recipe prep (fold-specific)
    prep_micro_k  = recipes::prep(rec_xgb_micro, training = analysis_micro, verbose = FALSE)
    baked_micro_k = recipes::bake(prep_micro_k, new_data = analysis_micro)
    X_micro_k     = baked_micro_k %>% dplyr::select(-all_of(c(id_col, ".train_row", outcome)))
    
    params_micro_k = dials::parameters(
      dials::trees(range = c(500L, 4000L)),
      dials::tree_depth(range = c(2L, 4L)),
      dials::learn_rate(range = c(-2.5, -0.5)),
      dials::finalize(dials::mtry(range = c(10L, 140L)), X_micro_k),
      dials::min_n(range = c(10L, 120L)),
      dials::loss_reduction(range = c(-4, 1)),
      dials::sample_prop(range = c(0.5, 0.9))
    )
    
    grid_micro_k = dials::grid_space_filling(params_micro_k, size = nrow(params_micro_k) * 10)
    
    set.seed(seed + 1000 + k)  # NEW
    tuned_micro_xgb_k = tune::tune_grid(
      wf_micro_xgb,
      resamples = inner_folds,
      grid = grid_micro_k,
      metrics = auc_set,
      control = tune::control_grid(save_pred = FALSE)
    )
    
    best_micro_xgb_k = tune::select_best(tuned_micro_xgb_k, metric = "roc_auc")
    wf_micro_xgb_final_k = tune::finalize_workflow(wf_micro_xgb, best_micro_xgb_k)
    
    fit_micro_xgb_k = parsnip::fit(wf_micro_xgb_final_k, data = analysis_micro)
    
    prob_assess_xgb_k = predict(fit_micro_xgb_k, new_data = assess_micro, type = "prob") %>%
      dplyr::bind_cols(assess_micro %>% dplyr::select(all_of(".train_row"))) %>%  # NEW
      dplyr::mutate(p3 = .data[[event_prob_col]]) %>%
      dplyr::select(.train_row, p3)  
    
    pred_oof_micro_xgb_list[[k]] <- prob_assess_xgb_k
  }
  
  # NEW: combine OOF predictions across outer folds
  pred_oof_micro_en <- dplyr::bind_rows(pred_oof_micro_en_list) %>%
    dplyr::arrange(.train_row)
  
  pred_oof_micro_xgb <- dplyr::bind_rows(pred_oof_micro_xgb_list) %>%
    dplyr::arrange(.train_row)
  
  # NEW: OOF completeness checks (critical for stacking correctness)
  n_train <- nrow(train)
  
  check_oof <- function(oof_df, name) {
    if (!(".train_row" %in% names(oof_df))) stop(name, ": missing .train_row")
    if (nrow(oof_df) != n_train) stop(name, ": nrow(oof) != n_train (", nrow(oof_df), " vs ", n_train, ")")
    if (anyDuplicated(oof_df$.train_row) > 0) stop(name, ": duplicated .train_row detected")
    if (!all(sort(oof_df$.train_row) == seq_len(n_train))) stop(name, ": .train_row does not cover 1:n_train")
    invisible(TRUE)
  }
  
  check_oof(pred_oof_cov, "OOF cov (p1)")
  check_oof(pred_oof_micro_en, "OOF EN (p2)")
  check_oof(pred_oof_micro_xgb, "OOF XGB (p3)")
  
  ## FS summary table
  fs_summary <- purrr::map_dfr(selected_features_by_fold, function(x) {
    tibble::tibble(
      fold = x$fold,
      feature = x$features
    )
  })
  
  fs_counts_by_fold <- purrr::map_dfr(selected_features_by_fold, function(x) {
    tibble::tibble(
      fold = x$fold,
      n_features = x$n_features,
      analysis_n = x$analysis_n,
      assess_n = x$assess_n
    )
  })
  
  fs_freq <- fs_summary %>%
    dplyr::count(feature, name = "n_folds") %>%
    dplyr::arrange(dplyr::desc(n_folds), feature)
  
  # --------------------------
  # Prediction correlation diagnostics (OOF)
  # --------------------------
  corr_rf_en <- cor(pred_oof_cov$p1, pred_oof_micro_en$p2, use = "complete.obs")
  corr_rf_xgb <- cor(pred_oof_cov$p1, pred_oof_micro_xgb$p3, use = "complete.obs")
  
  corr_rf_en_spearman <- cor(pred_oof_cov$p1, pred_oof_micro_en$p2, method = "spearman", use = "complete.obs")
  corr_rf_xgb_spearman <- cor(pred_oof_cov$p1, pred_oof_micro_xgb$p3, method = "spearman", use = "complete.obs")
  
  emit(sprintf("[CORR OOF] RF vs EN  (pearson)=%.3f | (spearman)=%.3f", corr_rf_en, corr_rf_en_spearman))
  emit(sprintf("[CORR OOF] RF vs XGB (pearson)=%.3f | (spearman)=%.3f", corr_rf_xgb, corr_rf_xgb_spearman))
  
  # NEW: Fit final OMICS models on FULL TRAIN for TEST prediction
  # NEW: Here we do FS on full training only once to define final feature set (for final model + SHAP)
  selected_feature_cols <- feature_cols
  if (fs_method == "limma") {
    selected_feature_cols <- select_features_limma(
      df_train = train_micro,
      outcome = outcome,
      feature_cols = feature_cols,
      limma_fdr = limma_fdr
    )
  }
  if (fs_method == "mrmr") {
    selected_feature_cols <- select_features_mrmr(
      df_train = train_micro,
      outcome = outcome,
      feature_cols = feature_cols,
      mrmr_top = mrmr_top
    )
  }
  if (fs_method == "sis") {
    selected_feature_cols <- select_features_sis(
      df_train = train_micro,   
      outcome = outcome,
      feature_cols = feature_cols,
      sis_top = sis_top,
      sis_cap = sis_cap,
      sis_method = sis_method
    )
  }
  final_selected_features <- selected_feature_cols
  
  # Apply selected features to omics datasets
  train_micro <- train_micro %>% dplyr::select(all_of(c(id_col, outcome, selected_feature_cols)))
  test_micro  <- test_micro  %>% dplyr::select(all_of(c(id_col, outcome, selected_feature_cols)))
  
  # Update feature_cols for downstream (XGB recipe, SHAP, stacking)
  feature_cols <- selected_feature_cols
  
  # NEW: tune final EN on full training (no FS inside; FS already fixed above)
  rec_en_micro = recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = train_micro) %>%
    recipes::update_role(all_of(id_col), new_role = "id") %>%
    recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
    recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_predictors()) %>%
    recipes::step_normalize(recipes::all_numeric_predictors())
  
  en_spec = parsnip::logistic_reg(
    penalty = tune(),
    mixture = tune()
  ) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")
  
  wf_micro_en = workflows::workflow() %>% workflows::add_recipe(rec_en_micro) %>% workflows::add_model(en_spec)
  
  en_grid = dials::grid_regular(
    dials::penalty(range = c(-10, 0)),
    dials::mixture(range = c(0, 1)),
    levels = c(50, 11)
  )
  
  # NEW: use the original outer folds for tuning on full training (performance estimate is not taken from this step)
  set.seed(seed)
  tuned_micro_en = tune::tune_grid(
    wf_micro_en,
    resamples = folds,
    grid = en_grid,
    metrics = auc_set,
    control = tune::control_grid(save_pred = TRUE)
  )
  
  notes_en <- tune::collect_notes(tuned_micro_en)
  if (nrow(notes_en) > 0) {
    emit("[EN tune notes] final training | showing up to 25 rows")
    print(utils::head(notes_en, 25))
  }
  
  metrics_en <- tune::collect_metrics(tuned_micro_en)
  if (!any(metrics_en$.metric == "roc_auc" & is.finite(metrics_en$mean))) {
    emit("[EN tune failure] final training | no finite roc_auc; fallback to fixed EN hyperparameters")
    if (nrow(notes_en) > 0) {
      emit("[EN tune failure] final training | full notes:")
      print(notes_en)
    }
    best_micro_en <- tibble::tibble(penalty = 0.01, mixture = 0.5)
  } else {
    best_micro_en = tune::select_best(tuned_micro_en, metric = "roc_auc")
  }
  wf_micro_en_final = tune::finalize_workflow(wf_micro_en, best_micro_en)
  
  # Fit final on full training and predict test (Model 2)
  fit_micro_en = parsnip::fit(wf_micro_en_final, data = train_micro)
  
  prob_test_micro_en = predict(fit_micro_en, new_data = test_micro, type = "prob") %>%
    dplyr::bind_cols(test_micro %>% dplyr::select(all_of(c(id_col, outcome))) %>% dplyr::rename(truth = all_of(outcome))) %>%
    dplyr::mutate(p2 = .data[[event_prob_col]])
  
  auc_test_2 = yardstick::roc_auc(prob_test_micro_en, truth = truth, p2, event_level = event_level)
  
  # --------------------------
  # 4) MODEL 3: Microbiome XGB (strongly regularized)
  # --------------------------
  rec_xgb_micro = recipes::recipe(stats::as.formula(paste(outcome, "~ .")), data = train_micro) %>%
    recipes::update_role(all_of(id_col), new_role = "id") %>%
    recipes::step_impute_median(recipes::all_numeric_predictors()) %>%
    recipes::step_impute_mode(recipes::all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_nominal_predictors()) %>%
    recipes::step_dummy(recipes::all_nominal_predictors()) %>%
    recipes::step_zv(recipes::all_predictors())
  
  xgb_micro_spec = parsnip::boost_tree(
    trees = tune(),
    tree_depth = tune(),
    learn_rate = tune(),
    mtry = tune(),
    min_n = tune(),
    loss_reduction = tune(),
    sample_size = tune(),
    stop_iter = 30
  ) %>%
    parsnip::set_engine("xgboost", eval_metric = "auc", event_level = event_level) %>%
    parsnip::set_mode("classification")
  
  wf_micro_xgb = workflows::workflow() %>% workflows::add_recipe(rec_xgb_micro) %>% workflows::add_model(xgb_micro_spec)
  
  prep_micro  = recipes::prep(rec_xgb_micro, training = train_micro, verbose = FALSE)
  baked_micro = recipes::bake(prep_micro, new_data = train_micro)
  X_micro     = baked_micro %>% dplyr::select(-all_of(c(id_col, outcome)))
  
  params_micro = dials::parameters(
    dials::trees(range = c(500L, 4000L)),
    dials::tree_depth(range = c(2L, 4L)),
    dials::learn_rate(range = c(-2.5, -0.5)),
    dials::finalize(dials::mtry(range = c(10L, 140L)), X_micro),
    dials::min_n(range = c(10L, 120L)),
    dials::loss_reduction(range = c(-4, 1)),
    dials::sample_prop(range = c(0.5, 0.9))
  )
  
  grid_micro = dials::grid_space_filling(params_micro, size = nrow(params_micro) * 10)
  
  set.seed(seed)
  tuned_micro_xgb = tune::tune_grid(
    wf_micro_xgb,
    resamples = folds,
    grid = grid_micro,
    metrics = auc_set,
    control = tune::control_grid(save_pred = TRUE)
  )
  
  best_micro_xgb = tune::select_best(tuned_micro_xgb, metric = "roc_auc")
  wf_micro_xgb_final = tune::finalize_workflow(wf_micro_xgb, best_micro_xgb)
  
  # Fit final on full training and predict test (Model 3)
  fit_micro_xgb = parsnip::fit(wf_micro_xgb_final, data = train_micro)
  
  prob_test_micro = predict(fit_micro_xgb, new_data = test_micro, type = "prob") %>%
    dplyr::bind_cols(test_micro %>% dplyr::select(all_of(c(id_col, outcome))) %>% dplyr::rename(truth = all_of(outcome))) %>%
    dplyr::mutate(p3 = .data[[event_prob_col]])
  
  auc_test_3 = yardstick::roc_auc(prob_test_micro, truth = truth, p3, event_level = event_level)
  
  # --------------------------
  # Model 3 SHAP (XGBoost)
  # --------------------------
  
  rec_prep_m3 <- workflows::extract_recipe(fit_micro_xgb)
  
  X_test_baked_m3 <- recipes::bake(rec_prep_m3, new_data = test_micro)
  
  X_explain_m3 <- dplyr::select(
    X_test_baked_m3,
    -dplyr::all_of(c(id_col, outcome))
  )
  
  obj_m3 <- list(
    booster = workflows::extract_fit_parsnip(fit_micro_xgb)$fit
  )
  
  pred_fun_m3 <- function(object, X) {
    as.numeric(
      predict(
        object$booster,
        newdata = as.matrix(X)
      )
    )
  }
  
  ks_m3 <- kernelshap::kernelshap(
    object = obj_m3,
    pred_fun = pred_fun_m3,
    X = X_explain_m3,
    parallel = FALSE
  )
  
  sv_m3 <- shapviz::shapviz(
    ks_m3,
    X = X_explain_m3
  )
  
  n_feats_m3 <- ncol(X_explain_m3)
  
  p_shap3_bar <- shapviz::sv_importance(
    sv_m3,
    kind = "bar",
    max_display = min(100L, n_feats_m3)
  ) +
    ggplot2::labs(
      title = paste0("Model 3: XGBoost (", omics_label, ")")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  p_shap3_bee <- shapviz::sv_importance(
    sv_m3,
    kind = "beeswarm",
    max_display = min(100L, n_feats_m3)
  ) +
    ggplot2::labs(
      title = paste0("Model 3: XGBoost (", omics_label, ")")
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    )
  
  
  # --------------------------
  # 5) Ensemble Model 4: Model 1 + Model 2 (ridge logistic on OOF preds)
  # --------------------------
  meta_12 = pred_oof_cov %>%
    dplyr::inner_join(pred_oof_micro_en, by = ".train_row") %>%  # NEW: join key changed
    dplyr::mutate(
      truth = as.factor(truth),
      z1 = logit(p1),
      z2 = logit(p2)
    ) %>%
    dplyr::select(truth, z1, z2)
  
  meta_folds = rsample::vfold_cv(meta_12, v = 5, strata = truth)
  
  meta12_spec = parsnip::logistic_reg(penalty = tune(), mixture = 0) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")
  
  meta12_wf = workflows::workflow() %>%
    workflows::add_model(meta12_spec) %>%
    workflows::add_formula(truth ~ z1 + z2)
  
  meta12_grid = dials::grid_regular(dials::penalty(range = c(-10, 2)), levels = 60)
  
  set.seed(seed)
  meta12_tuned = tune::tune_grid(
    meta12_wf,
    resamples = meta_folds,
    grid = meta12_grid,
    metrics = auc_set,
    control = tune::control_grid(save_pred = FALSE)
  )
  
  meta12_best = tune::select_best(meta12_tuned, metric = "roc_auc")
  meta12_fit  = tune::finalize_workflow(meta12_wf, meta12_best) %>% parsnip::fit(meta_12)
  
  # NEW: extract meta12 coefficients (b0, b1, b2) 
  coef_meta12_tbl <- broom::tidy(workflows::extract_fit_parsnip(meta12_fit)) %>%
    dplyr::group_by(term) %>% dplyr::slice(1) %>% dplyr::ungroup()
  
  b0 <- coef_meta12_tbl$estimate[coef_meta12_tbl$term == "(Intercept)"]
  b1 <- coef_meta12_tbl$estimate[coef_meta12_tbl$term == "z1"]
  b2 <- coef_meta12_tbl$estimate[coef_meta12_tbl$term == "z2"]
  
  if (length(b0) == 0) b0 <- 0
  if (length(b1) == 0) stop("meta12: coefficient for z1 not found")
  if (length(b2) == 0) stop("meta12: coefficient for z2 not found")
  
  emit(sprintf("[META12 COEF] b0=%.4f | b1(RF)=%.4f | b2(EN)=%.4f", b0, b1, b2))
  emit(sprintf("[META12 WEIGHT RATIO] |b2|/|b1| = %.4f", abs(b2)/abs(b1)))
  
  meta_test_12 = prob_test_cov %>%
    dplyr::select(all_of(c(id_col)), truth, p1) %>%
    dplyr::inner_join(prob_test_micro_en %>% dplyr::select(all_of(c(id_col)), p2), by = id_col) %>%
    dplyr::mutate(z1 = logit(p1), z2 = logit(p2))
  
  prob_test_4 = predict(meta12_fit, new_data = meta_test_12, type = "prob") %>%
    dplyr::bind_cols(meta_test_12 %>% dplyr::select(truth)) %>%
    dplyr::mutate(p4 = .data[[event_prob_col]])
  
  auc_test_4 = yardstick::roc_auc(prob_test_4, truth = truth, p4, event_level = event_level)
  
  # --------------------------
  # Model 4 SHAP (stacked RF + EN) -- match Model 5 style
  # --------------------------
  obj_m4 <- list(
    logit = logit,
    fit_cov = fit_cov,
    fit_micro = fit_micro_en,
    meta12_fit = meta12_fit,
    cov_cols = cov_cols,
    feature_cols = feature_cols,
    id_col = id_col,
    id_values = test[[id_col]],
    event_class = event_class
  )
  
  pred_fun_stack_prob_m4 <- function(object, X) {
    requireNamespace("workflows"); requireNamespace("parsnip")
    requireNamespace("recipes"); requireNamespace("hardhat")
    
    id_col <- object$id_col
    id_values <- object$id_values
    pcol <- paste0(".pred_", object$event_class)
    logit <- object$logit
    
    X2 <- X
    X2[[id_col]] <- id_values
    
    nd_cov <- X2 %>% dplyr::select(all_of(c(id_col, object$cov_cols)))
    nd_micro <- X2 %>% dplyr::select(all_of(c(id_col, object$feature_cols)))
    
    p1_tbl <- predict(object$fit_cov, new_data = nd_cov, type = "prob")
    p2_tbl <- predict(object$fit_micro, new_data = nd_micro, type = "prob")
    
    p1 <- p1_tbl[[pcol]]
    p2 <- p2_tbl[[pcol]]
    
    meta_df <- tibble::tibble(z1 = logit(p1), z2 = logit(p2))
    p4_tbl <- predict(object$meta12_fit, new_data = meta_df, type = "prob")
    p4_tbl[[pcol]]
  }
  
  X_explain_m4 <- test %>% dplyr::select(all_of(c(cov_cols, feature_cols)))
  
  ks_m4 <- kernelshap::kernelshap(
    object = obj_m4,
    X = X_explain_m4,
    pred_fun = pred_fun_stack_prob_m4,
    parallel = FALSE
  )
  
  sv_m4 <- shapviz::shapviz(ks_m4, X = X_explain_m4)
  
  n_feats_m4 <- ncol(X_explain_m4)
  
  p_shap4_bar <- shapviz::sv_importance(sv_m4, kind = "bar", max_display = min(100L, n_feats_m4), show_numbers = FALSE) +
    ggplot2::labs(title = "Model 4: Ensemble (1+2)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  p_shap4_bee <- shapviz::sv_importance(sv_m4, kind = "beeswarm", max_display = min(100L, n_feats_m4), show_numbers = FALSE) +
    ggplot2::labs(title = "Model 4: Ensemble (1+2)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  # --------------------------
  # 6) Ensemble Model 5:  Model 1 + Model 3 (ridge logistic on OOF preds)
  # --------------------------
  meta_13 = pred_oof_cov %>%
    dplyr::inner_join(pred_oof_micro_xgb, by = ".train_row") %>%  # NEW: join key changed
    dplyr::mutate(
      truth = as.factor(truth),
      z1 = logit(p1),
      z3 = logit(p3)
    ) %>%
    dplyr::select(truth, z1, z3)
  
  meta13_folds = rsample::vfold_cv(meta_13, v = 5, strata = truth)
  
  meta13_spec = parsnip::logistic_reg(penalty = tune(), mixture = 0) %>%
    parsnip::set_engine("glmnet") %>%
    parsnip::set_mode("classification")
  
  meta13_wf = workflows::workflow() %>%
    workflows::add_model(meta13_spec) %>%
    workflows::add_formula(truth ~ z1 + z3)
  
  meta13_grid = dials::grid_regular(dials::penalty(range = c(-10, 2)), levels = 60)
  
  set.seed(seed)
  meta13_tuned = tune::tune_grid(
    meta13_wf,
    resamples = meta13_folds,
    grid = meta13_grid,
    metrics = auc_set,
    control = tune::control_grid(save_pred = FALSE)
  )
  
  meta13_best = tune::select_best(meta13_tuned, metric = "roc_auc")
  meta13_fit  = tune::finalize_workflow(meta13_wf, meta13_best) %>% parsnip::fit(meta_13)
  
  # NEW: extract meta13 coefficients (c0, c1, c3) for RF+XGB
  coef_meta13_tbl <- broom::tidy(workflows::extract_fit_parsnip(meta13_fit)) %>%
    dplyr::group_by(term) %>% dplyr::slice(1) %>% dplyr::ungroup()
  
  c0 <- coef_meta13_tbl$estimate[coef_meta13_tbl$term == "(Intercept)"]
  c1 <- coef_meta13_tbl$estimate[coef_meta13_tbl$term == "z1"]
  c3 <- coef_meta13_tbl$estimate[coef_meta13_tbl$term == "z3"]
  
  if (length(c0) == 0) c0 <- 0
  if (length(c1) == 0) stop("meta13: coefficient for z1 not found")
  if (length(c3) == 0) stop("meta13: coefficient for z3 not found")
  
  emit(sprintf("[META13 COEF] c0=%.4f | c1(RF)=%.4f | c3(XGB)=%.4f", c0, c1, c3))
  emit(sprintf("[META13 WEIGHT RATIO] |c3|/|c1| = %.4f", abs(c3)/abs(c1)))
  

  meta_test_13 = prob_test_cov %>%
    dplyr::select(all_of(c(id_col)), truth, p1) %>%
    dplyr::inner_join(prob_test_micro %>% dplyr::select(all_of(c(id_col)), p3), by = id_col) %>%
    dplyr::mutate(z1 = logit(p1), z3 = logit(p3))
  
  prob_test_5 = predict(meta13_fit, new_data = meta_test_13, type = "prob") %>%
    dplyr::bind_cols(meta_test_13 %>% dplyr::select(truth)) %>%
    dplyr::mutate(p5 = .data[[event_prob_col]])
  
  auc_test_5 = yardstick::roc_auc(prob_test_5, truth = truth, p5, event_level = event_level)
  
  # --------------------------
  # 7) ROC plots
  # --------------------------
  roc_df <- function(df, truth_col, prob_col, model_label, event_level = "second") {
    truth <- rlang::sym(truth_col)
    prob  <- rlang::sym(prob_col)
    
    yardstick::roc_curve(df, truth = !!truth, !!prob, event_level = event_level) %>%
      dplyr::mutate(model = model_label)
  }
  
  roc_124 = dplyr::bind_rows(
    roc_df(prob_test_cov %>% dplyr::transmute(truth, p1), "truth", "p1",
           paste0("Model 1: RF cov (AUC = ", round(auc_test_1$.estimate, 3), ")"),
           event_level = event_level),
    roc_df(prob_test_micro_en %>% dplyr::transmute(truth, p2), "truth", "p2",
           paste0("Model 2: EN ", omics_label, " (AUC = ", round(auc_test_2$.estimate, 3), ")"),
           event_level = event_level),
    roc_df(prob_test_4 %>% dplyr::transmute(truth, p4), "truth", "p4",
           paste0("Model 4: Ensemble (1+2) (AUC = ", round(auc_test_4$.estimate, 3), ")"),
           event_level = event_level)
  )
  
  p_roc_124 = ggplot2::ggplot(roc_124, ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", color = NULL, title = "ROC: Ensemble 1+2") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.7, 0.2)
    )
  
  roc_135 = dplyr::bind_rows(
    roc_df(prob_test_cov %>% dplyr::transmute(truth, p1), "truth", "p1",
           paste0("Model 1: RF cov (AUC = ", round(auc_test_1$.estimate, 3), ")"),
           event_level = event_level),
    roc_df(prob_test_micro %>% dplyr::transmute(truth, p3), "truth", "p3",
           paste0("Model 3: XGB ", omics_label, " (AUC = ", round(auc_test_3$.estimate, 3), ")"),
           event_level = event_level),
    roc_df(prob_test_5 %>% dplyr::transmute(truth, p5), "truth", "p5",
           paste0("Model 5: Ensemble (1+3) (AUC = ", round(auc_test_5$.estimate, 3), ")"),
           event_level = event_level)
  )
  
  p_roc_135 = ggplot2::ggplot(roc_135, ggplot2::aes(x = 1 - specificity, y = sensitivity, color = model)) +
    ggplot2::geom_path(linewidth = 1) +
    ggplot2::geom_abline(linetype = 2) +
    ggplot2::coord_equal() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "inside",
      legend.position.inside = c(0.7, 0.2)
    )
  
  # --------------------------
  # 8) Model 1 SHAP + Model 2 coef plot
  # --------------------------
  rec_prep = workflows::extract_recipe(fit_cov)
  X_test_baked = recipes::bake(rec_prep, new_data = test_cov)
  X_explain_m1 = dplyr::select(X_test_baked, -dplyr::all_of(c(id_col, outcome)))
  
  obj_m1 = list(
    fit_rf = workflows::extract_fit_parsnip(fit_cov)$fit,
    event_class = event_class
  )
  
  pred_fun_m1 = function(object, X) {
    requireNamespace("ranger")
    
    fit_rf = object$fit_rf
    event_class = object$event_class
    p1_matrix = predict(fit_rf, data = X, type = "response")$predictions
    p1_matrix[, event_class]
  }
  
  ks_m1 = kernelshap::kernelshap(
    object = obj_m1,
    pred_fun = pred_fun_m1,
    X = X_explain_m1,
    parallel = FALSE
  )
  
  sv_m1 = shapviz::shapviz(ks_m1, X = X_explain_m1)
  
  p_shap1_bar = shapviz::sv_importance(sv_m1, kind = "bar", max_display = 10) +
    ggplot2::labs(title = "Model 1: RF Covariates") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  p_shap1_bee = shapviz::sv_importance(sv_m1, kind = "beeswarm", max_display = 10) +
    ggplot2::labs(title = "Model 1: RF Covariates") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  # Extract nonzero coefficients from Model 2
  coef_en = broom::tidy(fit_micro_en) %>%
    dplyr::filter(term != "(Intercept)", estimate != 0) %>%
    dplyr::arrange(dplyr::desc(abs(estimate)))
  
  coef_en_top = coef_en %>% dplyr::slice_head(n = 20)
  
  p_en_coef = coef_en_top %>%
    dplyr::mutate(sign = dplyr::if_else(estimate >= 0, "Positive", "Negative")) %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(term, abs(estimate)), y = estimate, fill = sign)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste0("Model 2: Elastic Net (", omics_label, ")"),
                  x = NULL, y = "Coefficient", fill = NULL) +
    ggplot2::theme(legend.position = "none", panel.grid.minor = ggplot2::element_blank())
  
  # --------------------------
  # 9) Model 5 SHAP (stacked)
  # --------------------------
  obj_m5 = list(
    logit = logit,
    fit_cov = fit_cov,
    fit_micro = fit_micro_xgb,
    meta13_fit = meta13_fit,
    cov_cols = cov_cols,
    feature_cols = feature_cols,
    id_col = id_col,
    id_values = test[[id_col]],
    event_class = event_class
  )
  
  pred_fun_stack_prob = function(object, X) {
    requireNamespace("workflows"); requireNamespace("parsnip")
    requireNamespace("recipes"); requireNamespace("hardhat")
    
    id_col = object$id_col
    id_values = object$id_values
    pcol = paste0(".pred_", object$event_class)
    logit = object$logit
    
    X2 = X
    X2[[id_col]] = id_values
    
    nd_cov   = X2 %>% dplyr::select(all_of(c(id_col, object$cov_cols)))
    nd_micro = X2 %>% dplyr::select(all_of(c(id_col, object$feature_cols)))
    
    p1_tbl = predict(object$fit_cov, new_data = nd_cov, type = "prob")
    p3_tbl = predict(object$fit_micro, new_data = nd_micro, type = "prob")
    
    p1 = p1_tbl[[pcol]]
    p3 = p3_tbl[[pcol]]
    
    meta_df = tibble::tibble(z1 = logit(p1), z3 = logit(p3))
    p5_tbl = predict(object$meta13_fit, new_data = meta_df, type = "prob")
    p5_tbl[[pcol]]
  }
  
  X_explain_m5 = test %>% dplyr::select(all_of(c(cov_cols, feature_cols)))
  
  ks_m5 = kernelshap::kernelshap(
    object = obj_m5,
    X = X_explain_m5,
    pred_fun = pred_fun_stack_prob,
    parallel = FALSE
  )
  
  sv_m5 = shapviz::shapviz(ks_m5, X = X_explain_m5)
  
  n_feats_m5 <- ncol(X_explain_m5)
  
  p_shap5_bar = shapviz::sv_importance(sv_m5, kind = "bar", max_display = min(100L, n_feats_m5), show_numbers = FALSE) +
    ggplot2::labs(title = "Model 5: Ensemble (1+3)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  p_shap5_bee = shapviz::sv_importance(sv_m5, kind = "beeswarm", max_display = min(100L, n_feats_m5), show_numbers = FALSE) +
    ggplot2::labs(title = "Model 5: Ensemble (1+3)") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  
  # --------------------------
  # Export 
  # --------------------------
  list(
    split = split,
    feature_selection = list(
      per_fold = selected_features_by_fold,
      counts_by_fold = fs_counts_by_fold,
      freq_table = fs_freq,
      final_selected_features = final_selected_features
    ),
    # NEW: helpful debug output
    selected_features_by_fold = selected_features_by_fold,
    auc = list(
      model1_rf_cov = auc_test_1,
      model2_en_micro = auc_test_2,
      model3_xgb_micro = auc_test_3,
      model4_ens_12 = auc_test_4,
      model5_ens_13 = auc_test_5
    ),
    fits = list(
      fit_cov = fit_cov,
      fit_micro_en = fit_micro_en,
      fit_micro_xgb = fit_micro_xgb,
      meta12_fit = meta12_fit,
      meta13_fit = meta13_fit
    ),
    predictions = list(
      prob_test_cov = prob_test_cov,
      prob_test_micro_en = prob_test_micro_en,
      prob_test_micro = prob_test_micro,
      prob_test_4 = prob_test_4,
      prob_test_5 = prob_test_5,
      # NEW: export OOF preds (strict FS-in-CV)
      pred_oof_cov = pred_oof_cov,
      pred_oof_micro_en = pred_oof_micro_en,
      pred_oof_micro_xgb = pred_oof_micro_xgb
    ),
    plots = list(
      p_roc_124 = p_roc_124,
      p_roc_135 = p_roc_135,
      p_shap1_bar = p_shap1_bar,
      p_shap1_bee = p_shap1_bee,
      p_shap3_bar = p_shap3_bar, # NEW: XGB
      p_shap3_bee = p_shap3_bee, # NEW: XGB
      p_en_coef = p_en_coef,
      p_shap4_bar = p_shap4_bar, # NEW: RF+EN
      p_shap4_bee = p_shap4_bee,  # NEW: RF+EN
      p_shap5_bar = p_shap5_bar,
      p_shap5_bee = p_shap5_bee
    ),
    coef_en = coef_en,
    # NEW：cor
    diagnostics = list(
      corr_rf_en = corr_rf_en,
      corr_rf_xgb = corr_rf_xgb,
      corr_rf_en_spearman = corr_rf_en_spearman,
      corr_rf_xgb_spearman = corr_rf_xgb_spearman
    ),
    meta_weights = list(
      # RF + EN
      meta12 = tibble::tibble(
        intercept = b0,
        rf_weight = b1,
        en_weight = b2,
        ratio_abs_en_over_rf = abs(b2) / abs(b1)
      ),
       # RF + XGB
      meta13 = tibble::tibble(
        intercept = c0,
        rf_weight = c1,
        xgb_weight = c3,
        ratio_abs_xgb_over_rf = abs(c3) / abs(c1)
      )
      
    )
  )
}