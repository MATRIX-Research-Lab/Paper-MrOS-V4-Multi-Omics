# Feature-selection routines shared by the primary and matched analyses.
# Keeping the limma implementation in one file prevents the two analyses from
# silently using different transformations or missing-value rules.

limma_select_features_train <- function(df_train, outcome, feature_cols,
                                        limma_fdr = 0.05, fallback_n = 50L) {
  if (!outcome %in% names(df_train)) stop("Missing outcome column: ", outcome)
  if (!length(feature_cols)) stop("No feature columns supplied to limma screening")
  if (!is.finite(limma_fdr) || limma_fdr <= 0 || limma_fdr > 1) {
    stop("limma_fdr must be in (0, 1]")
  }
  if (!is.finite(fallback_n) || fallback_n < 1) stop("fallback_n must be positive")
  y <- droplevels(as.factor(df_train[[outcome]]))
  if (nlevels(y) != 2L) stop("limma screening requires a binary outcome")
  X <- df_train[, feature_cols, drop = FALSE]
  if (any(!vapply(X, is.numeric, logical(1)))) {
    stop("limma screening requires numeric feature columns")
  }
  if (any(vapply(X, function(x) any(is.infinite(x)), logical(1)))) {
    stop("limma screening does not accept infinite feature values")
  }
  X_imp <- X
  for (j in seq_along(X_imp)) {
    v <- X_imp[[j]]
    if (anyNA(v)) {
      med <- stats::median(v, na.rm = TRUE)
      if (!is.finite(med)) med <- 0
      X_imp[[j]][is.na(v)] <- med
    }
  }
  mat <- t(as.matrix(X_imp))
  design <- stats::model.matrix(~ y)
  fit <- limma::eBayes(limma::lmFit(mat, design))
  tab <- limma::topTable(fit, coef = 2, number = Inf,
                         adjust.method = "BH", sort.by = "P")
  selected <- rownames(tab)[tab$adj.P.Val < limma_fdr]
  # The fallback triggers below two features, not below one. A screen that
  # returns a single predictor is not merely uninformative: glmnet rejects a
  # one-column design ("x should be a matrix with 2 or more columns"), and the
  # error escapes tune's per-configuration handling when the winning workflow is
  # refitted, killing the whole layer. Falling back to the top `fallback_n` by
  # adjusted P keeps the screen's contract -- always a usable feature set.
  if (length(selected) < 2L) {
    selected <- rownames(tab)[seq_len(min(as.integer(fallback_n), nrow(tab)))]
  }
  selected
}

# A supervised recipes step that learns the limma screen from the data supplied
# to prep(). During tune_grid(), recipes prep each resample independently, so an
# inner assessment fold cannot influence the features used to predict that fold.
step_limma_select <- function(recipe, feature_cols, outcome,
                              limma_fdr = 0.05, fallback_n = 50L,
                              role = NA, trained = FALSE,
                              selected_features = NULL, skip = FALSE,
                              id = recipes::rand_id("limma_select")) {
  recipes::add_step(
    recipe,
    step_limma_select_new(
      feature_cols = feature_cols, outcome = outcome,
      limma_fdr = limma_fdr, fallback_n = fallback_n,
      role = role, trained = trained,
      selected_features = selected_features,
      skip = skip, id = id
    )
  )
}

step_limma_select_new <- function(feature_cols, outcome, limma_fdr, fallback_n,
                                  role, trained, selected_features, skip, id) {
  structure(
    list(
      feature_cols = as.character(feature_cols),
      outcome = as.character(outcome),
      limma_fdr = limma_fdr,
      fallback_n = fallback_n,
      role = role,
      trained = trained,
      selected_features = selected_features,
      skip = skip,
      id = id
    ),
    class = c("step_limma_select", "step")
  )
}

prep.step_limma_select <- function(x, training, info = NULL, ...) {
  missing <- setdiff(c(x$outcome, x$feature_cols), names(training))
  if (length(missing)) {
    stop("limma recipe step is missing columns: ", paste(missing, collapse = ", "))
  }
  selected <- limma_select_features_train(
    df_train = training, outcome = x$outcome, feature_cols = x$feature_cols,
    limma_fdr = x$limma_fdr, fallback_n = x$fallback_n
  )
  step_limma_select_new(
    feature_cols = x$feature_cols, outcome = x$outcome,
    limma_fdr = x$limma_fdr, fallback_n = x$fallback_n,
    role = x$role, trained = TRUE, selected_features = selected,
    skip = x$skip, id = x$id
  )
}

bake.step_limma_select <- function(object, new_data, ...) {
  if (!isTRUE(object$trained)) stop("limma recipe step has not been trained")
  drop <- setdiff(object$feature_cols, object$selected_features)
  new_data[, setdiff(names(new_data), drop), drop = FALSE]
}

print.step_limma_select <- function(x, width = max(20, options()$width - 30), ...) {
  cat("Outcome-guided limma feature screening for ", length(x$feature_cols),
      " candidate predictors", if (isTRUE(x$trained)) " [trained]" else "", "\n",
      sep = "")
  invisible(x)
}

tidy.step_limma_select <- function(x, ...) {
  if (isTRUE(x$trained)) {
    tibble::tibble(terms = x$selected_features, id = x$id)
  } else {
    tibble::tibble(terms = x$feature_cols, id = x$id)
  }
}

required_pkgs.step_limma_select <- function(x, ...) c("limma")
