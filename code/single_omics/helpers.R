processed_files <- function(paths) {
  c(
    microbiome = file.path(paths$processed, "mb_asv_data_v4.rds"),
    metabolomics_matrix = file.path(paths$processed, "mbx_expr_data_v4.rds"),
    metabolomics_metadata = file.path(paths$processed, "mbx_sample_metadata_v4.rds"),
    proteomics_matrix = file.path(paths$processed, "prot_feature_v4.rds"),
    proteomics_metadata = file.path(paths$processed, "prot_sample_metadata_merge_v4.rds")
  )
}

check_processed_inputs <- function(paths) {
  files <- processed_files(paths)
  data.frame(
    input = names(files),
    path = file.path("data", "processed_data", basename(files)),
    exists = file.exists(files),
    stringsAsFactors = FALSE
  )
}

check_model_dependencies <- function() {
  pkgs <- c("tidymodels", "xgboost", "limma", "kernelshap", "shapviz")
  data.frame(
    package = pkgs,
    installed = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    stringsAsFactors = FALSE
  )
}

model_refit_inputs_available <- function(paths) {
  inputs <- check_processed_inputs(paths)
  deps <- check_model_dependencies()
  all(inputs$exists) && all(deps$installed)
}

# Both readers below take an explicit directory. There are two sets of
# layer-specific model outputs and they are not interchangeable:
#
#   paths$model_outputs        the primary analysis (01_-05_ read this)
#   paths$model_outputs_refit  the re-fitted layers that 06_ builds the matched
#                              cohort from
#
# See the comment on both paths in 00_setup.R for why they are kept apart.
model_output_files <- function(paths, dir = paths$model_outputs) {
  # Only the three canonical objects are consumed. Anything else in the
  # directory -- the generated "*_age_baseline.rds" files, a stray export -- is
  # never picked up by accident.
  pick <- function(omics) {
    canonical <- file.path(dir, paste0(omics, "_model_outputs.rds"))
    if (!file.exists(canonical)) {
      stop("Missing canonical ", omics, " model output: ", canonical,
           ". Run the prespecified primary runner.")
    }
    canonical
  }
  c(
    microbiome = pick("microbiome"),
    metabolomics = pick("metabolomics"),
    proteomics = pick("proteomics")
  )
}

# `expect_score` and `require_provenance` are both defaulted for the primary
# analysis, whose saved objects predate the provenance manifests and store the
# survival-oriented score. 06_ passes the strict values, because the objects it
# refits from are the hardened ones. See validate_prediction_output().
load_prediction_model_outputs <- function(paths, dir = paths$model_outputs,
                                          require_provenance = FALSE,
                                          expect_score = "auto") {
  files <- model_output_files(paths, dir)
  out <- lapply(files, function(f) {
    x <- readRDS(f)
    if (exists("validate_prediction_output", mode = "function")) {
      validate_prediction_output(x, basename(f), expect_score = expect_score)
    }
    if (isTRUE(require_provenance) &&
        (is.null(x$provenance) || is.null(x$provenance$signature))) {
      stop("Model output lacks provenance: ", f,
           ". Rebuild the canonical primary outputs with the hardened runner.")
    }
    x
  })
  names(out) <- names(files)
  attr(out, "source_files") <- unname(files)
  out
}

manuscript_table1 <- function() {
  data.frame(
    variable = c(
      "Race", "  African American", "  Asian", "  Others", "  White",
      "Education", "  <= High School", "  College", "  Graduate School",
      "Self-Rated Overall Health", "  Good/Excellent", "  Very Poor/Poor/Fair",
      "Marital Status", "  Married", "  Others", "  Widowed",
      "Smoking Status", "  Missing", "  Non-Smoker", "  Past or Current Smoker",
      "Diabetes", "High Blood Pressure", "Cancer", "BMI", "  Missing",
      "Teng 3MS Score", "  Missing", "Geriatric Depression Scale", "  Missing",
      "PASE Score", "Total Medications", "Age"
    ),
    overall = c(
      "", "25 (2.8%)", "32 (3.6%)", "26 (3.0%)", "796 (90.6%)",
      "", "158 (18.0%)", "332 (37.8%)", "389 (44.3%)",
      "", "872 (99.2%)", "7 (0.8%)",
      "", "633 (72.0%)", "65 (7.4%)", "181 (20.6%)",
      "", "91 (10.4%)", "367 (41.8%)", "421 (47.9%)",
      "132 (15.0%)", "451 (51.3%)", "416 (47.3%)", "26.90 (3.70)", "1",
      "92.94 (6.67)", "6", "1.68 (1.82)", "1",
      "123.38 (66.05)", "9.00 (4.68)", "84.01 (3.86)"
    ),
    deceased = c(
      "", "12 (2.5%)", "13 (2.7%)", "11 (2.3%)", "449 (92.6%)",
      "", "94 (19.4%)", "190 (39.2%)", "201 (41.4%)",
      "", "479 (98.8%)", "6 (1.2%)",
      "", "338 (69.7%)", "30 (6.2%)", "117 (24.1%)",
      "", "58 (12.0%)", "183 (37.7%)", "244 (50.3%)",
      "87 (17.9%)", "261 (53.8%)", "241 (49.7%)", "26.80 (3.82)", "1",
      "91.93 (7.46)", "4", "2.02 (1.99)", "1",
      "114.52 (68.62)", "9.48 (4.62)", "85.08 (4.16)"
    ),
    active = c(
      "", "13 (3.3%)", "19 (4.8%)", "15 (3.8%)", "347 (88.1%)",
      "", "64 (16.2%)", "142 (36.0%)", "188 (47.7%)",
      "", "393 (99.7%)", "1 (0.3%)",
      "", "295 (74.9%)", "35 (8.9%)", "64 (16.2%)",
      "", "33 (8.4%)", "184 (46.7%)", "177 (44.9%)",
      "45 (11.4%)", "190 (48.2%)", "175 (44.4%)", "27.02 (3.54)", "0",
      "94.18 (5.28)", "2", "1.27 (1.49)", "0",
      "134.30 (61.09)", "8.41 (4.70)", "82.69 (2.96)"
    ),
    p_value = c(
      "0.14", "", "", "", "",
      "0.2", "", "", "",
      "0.14", "", "",
      "0.009", "", "", "",
      "0.017", "", "", "",
      "0.007", "0.10", "0.12", "0.3", "",
      "<0.001", "", "<0.001", "",
      "<0.001", "<0.001", "<0.001"
    ),
    stringsAsFactors = FALSE
  )
}

write_table1_pdf <- function(table1, output_pdf) {
  grDevices::pdf(output_pdf, width = 11, height = 8.5)
  on.exit(grDevices::dev.off(), add = TRUE)

  # Body only: no title and no footnote. The number, title, legend and the note
  # on which tests produced the P values all live in main_tables.tex, so a table
  # cannot carry one caption while the manuscript carries another. This is the
  # same convention the supplementary tables already follow.
  grid::grid.newpage()

  x <- c(0.04, 0.38, 0.56, 0.73, 0.90)
  headers <- c("Variable", "Overall N = 879", "Deceased N = 485", "Active N = 394", "p-value")
  y0 <- 0.90
  row_h <- 0.024
  for (j in seq_along(headers)) {
    grid::grid.text(headers[j], x = x[j], y = y0, just = "left",
                    gp = grid::gpar(fontsize = 8, fontface = "bold"))
  }
  grid::grid.lines(x = c(0.04, 0.97), y = c(y0 - 0.012, y0 - 0.012))

  for (i in seq_len(nrow(table1))) {
    y <- y0 - i * row_h
    is_group <- table1$overall[i] == "" && !grepl("^  ", table1$variable[i])
    face <- if (is_group) "bold" else "plain"
    grid::grid.text(table1$variable[i], x = x[1], y = y, just = "left",
                    gp = grid::gpar(fontsize = 7.2, fontface = face))
    grid::grid.text(table1$overall[i], x = x[2], y = y, just = "left", gp = grid::gpar(fontsize = 7.2))
    grid::grid.text(table1$deceased[i], x = x[3], y = y, just = "left", gp = grid::gpar(fontsize = 7.2))
    grid::grid.text(table1$active[i], x = x[4], y = y, just = "left", gp = grid::gpar(fontsize = 7.2))
    grid::grid.text(table1$p_value[i], x = x[5], y = y, just = "left", gp = grid::gpar(fontsize = 7.2))
  }
  TRUE
}

# Figure 1: participant flow (panel a) and the analysis pipeline (panel b).
#
# Counts are taken from the manuscript's analytic samples so that the figure and
# Table 1 cannot drift apart. `flowchart_counts()` is the single place they are
# defined; `participant_flow_counts()` recomputes them from the saved model
# outputs so the hard-coded values can be checked.
# Counts drawn in Figure 1. Declared here and asserted against the data by
# participant_flow_counts(), which 04_ calls before drawing anything -- the same
# declare-then-verify contract 03_ applies to Table 1.
#
# The feature counts come from the analysis contract in 00_setup.R so all three
# boxes report the same quantity: the block that actually enters the models.
# Figure 1 previously drew "463 genera", a number from the superseded Greengenes
# annotation. The SILVA release this analysis uses contains 429 distinct
# genus-level lineages in total (351 distinct genus labels), so 463 cannot arise
# at any stage -- filtering only removes taxa -- and 147 is what survives the
# ANCOM-BC2 prevalence filter and is modelled. The manuscript Methods already
# say 147.
flowchart_counts <- function() {
  list(
    total = 879L, deaths = 485L, alive = 394L,
    microbiome = 765L, metabolomics = 447L, proteomics = 446L,
    all_three = 332L,
    genera = unname(analysis_contract$microbiome[["features"]]),
    metabolites = unname(analysis_contract$metabolomics[["features"]]),
    proteins = unname(analysis_contract$proteomics[["features"]])
  )
}

# Recompute every count in Figure 1 from the saved model outputs. Returns a data
# frame comparing the recomputed value with the value drawn in the figure, so a
# mismatch is visible rather than silent.
#
# Feature counts are included, not just participants. The stale genus count sat
# in the figure precisely because this check covered participants only -- and
# because nothing called it.
participant_flow_counts <- function(model_outputs) {
  ids <- lapply(model_outputs, function(o) as.character(o$split$data$ID))
  cnt <- flowchart_counts()
  meta_cols <- unique(c(analysis_non_feature_columns, analysis_covariates))
  feats <- vapply(model_outputs,
                  function(o) length(setdiff(names(o$split$data), meta_cols)),
                  integer(1))
  status <- do.call(rbind, lapply(model_outputs, function(o) {
    data.frame(ID = as.character(o$split$data$ID),
               status = as.character(o$split$data$status),
               stringsAsFactors = FALSE)
  }))
  status <- status[!duplicated(status$ID), , drop = FALSE]
  data.frame(
    quantity = c(names(ids), "all_three", "total", "deaths", "alive",
                 paste0(names(feats), "_features")),
    recomputed = c(vapply(ids, function(x) length(unique(x)), integer(1)),
                   length(Reduce(intersect, ids)),
                   length(unique(unlist(ids))),
                   sum(status$status == "Deceased"),
                   sum(status$status == "Active"),
                   unname(feats)),
    in_figure = c(cnt$microbiome, cnt$metabolomics, cnt$proteomics,
                  cnt$all_three, cnt$total, cnt$deaths, cnt$alive,
                  cnt$genera, cnt$metabolites, cnt$proteins),
    stringsAsFactors = FALSE
  )
}

write_flowchart_pdf <- function(output_file) {
  n <- flowchart_counts()
  # cairo_pdf keeps the UTF-8 glyphs (>=, en dash, middle dot); the base pdf
  # device transliterates them.
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(output_file, width = 10, height = 5.6)
  } else {
    grDevices::pdf(output_file, width = 10, height = 5.6)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport())

  fill_data <- "#EAF1FA"; fill_step <- "#EFEAF7"; fill_out <- "#F2F2F2"
  edge <- "#33475B"

  box <- function(label, x, y, w, h, fill = fill_data, fontsize = 7.4, face = "plain") {
    grid::grid.roundrect(x, y, w, h, r = grid::unit(0.035, "snpc"),
                         gp = grid::gpar(fill = fill, col = edge, lwd = 0.9))
    grid::grid.text(label, x, y, gp = grid::gpar(fontsize = fontsize, fontface = face,
                                                 lineheight = 1.15))
  }
  down <- function(x, y0, y1) {
    grid::grid.segments(x, y0, x, y1,
                        arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed"),
                        gp = grid::gpar(lwd = 0.9, fill = edge, col = edge))
  }
  elbow <- function(x0, y0, x1, y1) {
    ym <- (y0 + y1) / 2
    grid::grid.lines(c(x0, x0, x1, x1), c(y0, ym, ym, y1),
                     arrow = grid::arrow(length = grid::unit(0.06, "inches"), type = "closed"),
                     gp = grid::gpar(lwd = 0.9, fill = edge, col = edge))
  }
  panel <- function(letter, x) {
    grid::grid.text(letter, x, 0.955, just = "left",
                    gp = grid::gpar(fontsize = 13, fontface = "bold"))
  }

  # ---- panel a: participant flow -------------------------------------------
  panel("a", 0.025)
  ax <- 0.245
  box(sprintf("MrOS Visit 4 (2014–2016)\nclinical phenotyping, stool and serum biospecimens"),
      ax, 0.855, 0.42, 0.115)
  down(ax, 0.795, 0.755)
  box(sprintf("Participants with ≥ 1 omics layer and\nascertained vital status\nn = %d", n$total),
      ax, 0.685, 0.42, 0.125, face = "bold")
  bx <- c(0.095, 0.245, 0.395); by <- 0.485
  for (i in seq_along(bx)) elbow(ax, 0.622, bx[i], by + 0.062)
  box(sprintf("Microbiome\nn = %d\n%d genera", n$microbiome, n$genera), bx[1], by, 0.135, 0.125)
  box(sprintf("Metabolome\nn = %d\n%s features", n$metabolomics, format(n$metabolites, big.mark = ",")),
      bx[2], by, 0.135, 0.125)
  box(sprintf("Proteome\nn = %d\n%s features", n$proteomics, format(n$proteins, big.mark = ",")),
      bx[3], by, 0.135, 0.125)
  for (i in seq_along(bx)) elbow(bx[i], 0.422, ax, 0.365)
  box(sprintf("All three omics layers available\nn = %d", n$all_three), ax, 0.300, 0.42, 0.105)
  down(ax, 0.248, 0.205)
  box(sprintf("Mortality follow-up through September 2024\n%d deaths (%.1f%%), %d alive",
              n$deaths, 100 * n$deaths / n$total, n$alive),
      ax, 0.140, 0.42, 0.115, fill = fill_out)

  # ---- panel b: analysis pipeline ------------------------------------------
  panel("b", 0.525)
  px <- 0.755
  box("Clinical covariates (14)\nand omics features", px, 0.855, 0.42, 0.105)
  down(px, 0.803, 0.762)
  box(paste("Preprocessing and feature screening",
            "microbiome: CLR / DESeq2 / ANCOM-BC2",
            "metabolome and proteome: limma / mRMR / SIS",
            sep = "\n"),
      px, 0.690, 0.42, 0.135, fill = fill_step)
  for (i in c(0.625, 0.755, 0.885)) elbow(px, 0.622, i, 0.545)
  box("Random forest\nclinical", 0.625, 0.485, 0.115, 0.105, fill = fill_step)
  box("Elastic net\nomics", 0.755, 0.485, 0.115, 0.105, fill = fill_step)
  box("XGBoost\nomics", 0.885, 0.485, 0.115, 0.105, fill = fill_step)
  for (i in c(0.625, 0.755, 0.885)) elbow(i, 0.432, px, 0.372)
  box("Ensemble (stacking)\nridge meta-learner: RF + EN, RF + XGB", px, 0.310, 0.42, 0.105,
      fill = fill_step, face = "bold")
  down(px, 0.257, 0.213)
  box(paste("Nested cross-validation (outer 5-fold, inner 3-fold) and 20% held-out test set",
            "AUC-ROC · SHAP feature attribution · cross-omics correlation networks",
            sep = "\n"),
      px, 0.150, 0.42, 0.115, fill = fill_out)

  grid::grid.lines(c(0.495, 0.495), c(0.06, 0.94),
                   gp = grid::gpar(col = "grey80", lwd = 0.7))
  TRUE
}

write_ggplot_panel_figure <- function(output_file, plots, labels, ncol, width, height) {
  stopifnot(length(plots) == length(labels))
  nrow <- ceiling(length(plots) / ncol)
  grDevices::pdf(output_file, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  layout <- grid::grid.layout(nrow, ncol)
  grid::pushViewport(grid::viewport(layout = layout))

  for (i in seq_along(plots)) {
    row <- ceiling(i / ncol)
    col <- ((i - 1) %% ncol) + 1
    grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
    print(plots[[i]] + ggplot2::labs(title = NULL), newpage = FALSE)
    grid::grid.text(
      labels[i],
      x = 0.03,
      y = 0.98,
      just = c("left", "top"),
      gp = grid::gpar(fontsize = 10, fontface = "bold")
    )
    grid::popViewport()
  }
  grid::popViewport()
}

hide_legend <- function(plot) {
  plot + ggplot2::theme(legend.position = "none")
}

# "age" is listed first so the LR age-only baseline appears at the top of the
# ROC legend; the remaining five keep their original order/colours.
roc_palette <- c(
  age = "#000000",
  clinical = "#4D4D4D",
  en = "#0072B2",
  xgb = "#D55E00",
  rf_en = "#009E73",
  rf_xgb = "#CC79A7"
)

# A fitted glm keeps its formula/terms environment, which for a formula built
# inside this function captures the whole model-output object and bloats the
# saved RDS to hundreds of MB. Detach that environment and drop the bulky
# per-row components so only the small, useful parts remain (coefficients, etc.).
strip_glm_for_storage <- function(fit) {
  fit$data <- NULL
  fit$model <- NULL
  fit$residuals <- NULL
  fit$fitted.values <- NULL
  fit$effects <- NULL
  fit$linear.predictors <- NULL
  fit$weights <- NULL
  fit$prior.weights <- NULL
  fit$y <- NULL
  fit$qr <- NULL
  attr(fit$terms, ".Environment") <- globalenv()
  if (!is.null(fit$formula)) environment(fit$formula) <- globalenv()
  fit
}

# Age-only baseline: logistic regression (status ~ age_v4) fit on the SAME
# train/test split stored in the model-output RDS, so it is directly comparable
# to the other five models within each omics panel and needs no model refitting.
# Returns predictions in the same tibble shape as the pipeline's prob_test_*,
# an AUC tibble matching the stored yardstick format, the DeLong 95% CI, and the
# ROC curve coordinates for plotting.
fit_age_baseline <- function(model_output) {
  sp <- model_output$split
  if (is.null(sp) || is.null(sp$data) || !("age_v4" %in% names(sp$data))) {
    return(NULL)
  }
  dat <- sp$data
  train <- dat[sp$in_id, , drop = FALSE]
  test <- dat[-sp$in_id, , drop = FALSE]
  train$status <- factor(train$status)
  lv <- levels(train$status)
  test$status <- factor(test$status, levels = lv)
  fit <- stats::glm(status ~ age_v4, data = train, family = stats::binomial())
  p_second <- stats::predict(fit, newdata = test, type = "response")
  pred <- tibble::tibble(
    .pred_Deceased = 1 - p_second,
    .pred_Active = p_second,
    ID = as.character(test$ID),
    truth = test$status,
    p_age = p_second
  )
  # roc_death() pins cases = Deceased and direction = "<". Do not pass `lv`
  # here: the modelling frames carry status as levels c("Deceased", "Active"),
  # so handing those to pROC makes survival the event and yields the ROC of the
  # complementary task -- identical AUC, mirrored curve. `direction = "auto"`
  # hides that, because it always returns the orientation with AUC >= 0.5.
  r <- roc_death(pred$truth, pred$.pred_Deceased)
  ci <- as.numeric(pROC::ci.auc(r, conf.level = 0.95))
  co <- pROC::coords(r, "all", ret = c("specificity", "sensitivity"), transpose = FALSE)
  fit <- strip_glm_for_storage(fit)
  list(
    fit = fit,
    predictions = list(prob_test_age = pred),
    auc = list(age_only = tibble::tibble(
      .metric = "roc_auc", .estimator = "binary", .estimate = ci[2]
    )),
    ci = c(lower = ci[1], estimate = ci[2], upper = ci[3]),
    roc_curve = data.frame(specificity = co$specificity, sensitivity = co$sensitivity)
  )
}

# 95% DeLong CI for each of the five stored models, computed from the saved
# test-set predictions (prob_test_*). The point AUC reproduces the stored value.
auc_ci_by_model <- function(model_output, conf_level = 0.95) {
  res <- list()
  preds <- model_output$predictions
  if (is.null(preds)) return(res)
  for (nm in names(preds)) {
    d <- preds[[nm]]
    if (!is.data.frame(d)) next
    pcol <- grep("^p[1-5]$", names(d), value = TRUE)
    if (length(pcol) != 1) next
    if (!all(c("truth", ".pred_Deceased") %in% names(d))) next
    tr <- d[["truth"]]
    if (!is.factor(tr)) tr <- factor(tr)
    if (length(levels(droplevels(tr))) != 2) next
    mnum <- as.integer(sub("^p", "", pcol))
    # roc_death() pins cases = Deceased and direction = "<". Do not use
    # direction = "auto" here: it returns whichever orientation gives
    # AUC >= 0.5, so a genuinely below-chance model is reported as its
    # complement and an inverted score never surfaces. That masking is what
    # kept the survival-oriented ROC curves in Figure 2 invisible for so long.
    ci <- tryCatch({
      r <- roc_death(tr, d[[".pred_Deceased"]])
      as.numeric(pROC::ci.auc(r, conf.level = conf_level))
    }, error = function(e) rep(NA_real_, 3))
    res[[as.character(mnum)]] <- c(auc = ci[2], lower = ci[1], upper = ci[3])
  }
  res
}

# Curve coordinates are recomputed from the stored test-set predictions; only
# the legend text is taken from the saved plot objects.
#
# The saved ROC curves cannot be used directly. They were built with survival as
# the event -- the modelling frames carry status as levels
# c("Deceased", "Active"), and the first level was taken as the event -- so each
# saved curve is the ROC of predicting *Active*. That has the same AUC as the
# intended curve but a mirrored shape, and at a fixed specificity the sensitivity
# it reports is off by up to 0.27. Recomputing through roc_death(), which pins
# cases = Deceased and direction = "<", puts every line back on the death task.
#
# The printed AUCs are still parsed from the saved labels so the legend continues
# to match Table 2 exactly. AUC is invariant to the swap, so those numbers were
# never affected.
roc_bind_five <- function(model_output, omics_label) {
  d124 <- as.data.frame(model_output$plots$p_roc_124$data)
  d135 <- as.data.frame(model_output$plots$p_roc_135$data)
  d135 <- d135[!as.character(d135$model) %in% intersect(as.character(d124$model), as.character(d135$model)), , drop = FALSE]
  saved <- dplyr::bind_rows(d124, d135)
  raw <- as.character(saved$model)
  model_number <- as.integer(sub("^Model ([0-9]+):.*", "\\1", raw))
  auc_raw <- suppressWarnings(as.numeric(sub(".*AUC = ([0-9.]+).*", "\\1", raw)))
  auc_of <- tapply(auc_raw, model_number, function(z) z[[1]])

  # Model number -> the stored test-set prediction frame it was drawn from.
  pred_of <- c("1" = "prob_test_cov", "2" = "prob_test_micro_en",
               "3" = "prob_test_micro", "4" = "prob_test_4", "5" = "prob_test_5")
  group_of <- c("1" = "clinical", "2" = "en", "3" = "xgb",
                "4" = "rf_en", "5" = "rf_xgb")
  # The figure legend carries the 95% DeLong interval alongside each AUC, so the
  # legend matches Table 2 and the reader can see how wide the test-set
  # intervals are without leaving the figure.
  ci_map <- tryCatch(auc_ci_by_model(model_output), error = function(e) list())

  parts <- list()
  for (m in sort(unique(model_number))) {
    key <- as.character(m)
    pred <- model_output$predictions[[pred_of[[key]]]]
    if (is.null(pred) || !all(c("truth", ".pred_Deceased") %in% names(pred))) next
    co <- pROC::coords(roc_death(pred$truth, pred$.pred_Deceased), "all",
                       ret = c("specificity", "sensitivity"), transpose = FALSE)
    v <- ci_map[[key]]
    ci_txt <- if (is.null(v) || any(is.na(v))) {
      ""
    } else {
      sprintf(" (%.2f–%.2f)", v[["lower"]], v[["upper"]])
    }
    a <- auc_of[[key]]
    label <- switch(
      key,
      "1" = sprintf("RF clinical  AUC = %.3f%s", a, ci_txt),
      "2" = sprintf("EN %s  AUC = %.3f%s", omics_label, a, ci_txt),
      "3" = sprintf("XGB %s  AUC = %.3f%s", omics_label, a, ci_txt),
      "4" = sprintf("RF clinical + EN %s  AUC = %.3f%s", omics_label, a, ci_txt),
      "5" = sprintf("RF clinical + XGB %s  AUC = %.3f%s", omics_label, a, ci_txt)
    )
    parts[[key]] <- data.frame(
      specificity = co$specificity,
      sensitivity = co$sensitivity,
      group = group_of[[key]],
      label = label,
      stringsAsFactors = FALSE
    )
  }
  roc <- do.call(rbind, parts)
  rownames(roc) <- NULL
  ab <- tryCatch(fit_age_baseline(model_output), error = function(e) NULL)
  if (!is.null(ab) && nrow(ab$roc_curve) > 0) {
    age_df <- data.frame(
      specificity = ab$roc_curve$specificity,
      sensitivity = ab$roc_curve$sensitivity,
      group = "age",
      label = sprintf("LR age-only  AUC = %.3f (%.2f–%.2f)",
                      ab$ci[["estimate"]], ab$ci[["lower"]], ab$ci[["upper"]]),
      stringsAsFactors = FALSE
    )
    roc <- rbind(roc, age_df)
  }
  roc$group <- factor(roc$group, levels = names(roc_palette))
  roc
}

theme_roc <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.text = ggplot2::element_text(colour = "black", size = base_size),
      axis.title = ggplot2::element_text(colour = "black", size = base_size + 1),
      plot.title = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size, colour = "black"),
      legend.key.height = grid::unit(9, "pt"),
      legend.key.width = grid::unit(14, "pt"),
      legend.background = ggplot2::element_rect(
        fill = scales::alpha("white", 0.9), colour = "grey70", linewidth = 0.2
      ),
      legend.margin = ggplot2::margin(2, 3, 2, 3),
      legend.spacing.y = grid::unit(0, "pt"),
      plot.margin = ggplot2::margin(4, 6, 4, 4)
    )
}

make_roc_panel <- function(model_output, title, omics_label) {
  roc <- roc_bind_five(model_output, omics_label)
  groups <- levels(roc$group)[levels(roc$group) %in% unique(as.character(roc$group))]
  label_for_group <- tapply(as.character(roc$label), as.character(roc$group), function(x) x[[1]])[groups]

  roc_linetype <- c(
    clinical = "solid", en = "solid", xgb = "solid",
    rf_en = "solid", rf_xgb = "solid", age = "22"
  )
  ggplot2::ggplot(
    roc,
    ggplot2::aes(x = 1 - specificity, y = sensitivity, group = group, colour = group, linetype = group)
  ) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey55") +
    ggplot2::geom_path(linewidth = 0.55, lineend = "round") +
    ggplot2::scale_colour_manual(values = roc_palette[groups], breaks = groups, labels = unname(label_for_group)) +
    ggplot2::scale_linetype_manual(values = roc_linetype[groups], breaks = groups, labels = unname(label_for_group)) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    ggplot2::scale_y_continuous(limits = c(0, 1), expand = c(0, 0), breaks = seq(0, 1, 0.25)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "1 \u2212 Specificity", y = "Sensitivity", title = title, colour = NULL, linetype = NULL) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(ncol = 1, byrow = TRUE, override.aes = list(linewidth = 0.8)),
      linetype = ggplot2::guide_legend(ncol = 1, byrow = TRUE)
    ) +
    theme_roc(base_size = 7) +
    ggplot2::theme(
      legend.position = c(0.98, 0.02),
      legend.justification = c(1, 0),
      legend.direction = "vertical"
    )
}

write_figure2_pdf <- function(model_outputs, output_file) {
  fig2 <- ggpubr::ggarrange(
    make_roc_panel(model_outputs$microbiome, "Microbiome", "microbiome"),
    make_roc_panel(model_outputs$metabolomics, "Metabolomics", "metabolomics"),
    make_roc_panel(model_outputs$proteomics, "Proteomics", "proteomics"),
    labels = c("a", "b", "c"),
    ncol = 3,
    nrow = 1,
    common.legend = FALSE,
    font.label = list(face = "bold", size = 10, family = "Helvetica"),
    hjust = -0.2,
    vjust = 1.4
  )
  ggplot2::ggsave(output_file, plot = fig2, width = 88 * 3, height = 88, units = "mm")
}

# Build a tidy AUC + 95% CI table (DeLong) for the five stored models plus the
# age-only logistic baseline, across all omics, and write it to an .xlsx file.
build_auc_ci_table <- function(model_outputs) {
  model_lab <- c(
    "1" = "RF clinical", "2" = "EN omics", "3" = "XGB omics",
    "4" = "RF clinical + EN", "5" = "RF clinical + XGB"
  )
  rows <- list()
  for (om in names(model_outputs)) {
    o <- model_outputs[[om]]
    ci_map <- auc_ci_by_model(o)
    for (mn in names(ci_map)) {
      v <- ci_map[[mn]]
      rows[[length(rows) + 1]] <- data.frame(
        omics = om, model = unname(model_lab[[mn]]),
        auc = v[["auc"]], ci_lower = v[["lower"]], ci_upper = v[["upper"]],
        stringsAsFactors = FALSE
      )
    }
    ab <- tryCatch(fit_age_baseline(o), error = function(e) NULL)
    if (!is.null(ab)) {
      rows[[length(rows) + 1]] <- data.frame(
        omics = om, model = "Age-only (logistic)",
        auc = ab$ci[["estimate"]], ci_lower = ab$ci[["lower"]], ci_upper = ab$ci[["upper"]],
        stringsAsFactors = FALSE
      )
    }
  }
  tab <- do.call(rbind, rows)
  tab$`AUC (95% CI)` <- sprintf("%.3f (%.3f\u2013%.3f)", tab$auc, tab$ci_lower, tab$ci_upper)
  tab
}

write_auc_ci_table <- function(model_outputs, output_file) {
  tab <- build_auc_ci_table(model_outputs)
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  openxlsx::write.xlsx(tab, output_file, overwrite = TRUE)
  invisible(tab)
}

# Fit the LR age-only baseline for every omics from the saved split and write one
# RDS per omics (same directory / naming convention as the model-output files).
# Returns the fitted baselines invisibly. Does not refit the five main models.
fit_and_save_age_baselines <- function(model_outputs, paths) {
  dir.create(paths$model_outputs, recursive = TRUE, showWarnings = FALSE)
  res <- list()
  for (om in names(model_outputs)) {
    ab <- fit_age_baseline(model_outputs[[om]])
    if (is.null(ab)) next
    saveRDS(ab, file.path(paths$model_outputs, paste0(om, "_age_baseline.rds")))
    res[[om]] <- ab
  }
  invisible(res)
}

# Wide AUC (95% CI) matrix: one row per model (LR age-only first), one column per
# omics. Cell values are "AUC (lower-upper)".
build_auc_ci_wide <- function(model_outputs) {
  model_lab <- c(
    "1" = "RF clinical", "2" = "EN omics", "3" = "XGB omics",
    "4" = "RF clinical + EN", "5" = "RF clinical + XGB"
  )
  model_order <- c("LR age-only", unname(model_lab))
  omics_order <- names(model_outputs)
  cell <- function(v) sprintf("%.3f (%.3f\u2013%.3f)", v[["auc"]], v[["lower"]], v[["upper"]])
  mat <- matrix("", nrow = length(model_order), ncol = length(omics_order),
                dimnames = list(model_order, omics_order))
  for (om in omics_order) {
    o <- model_outputs[[om]]
    cm <- auc_ci_by_model(o)
    for (mn in names(cm)) {
      mat[unname(model_lab[[mn]]), om] <- cell(cm[[mn]])
    }
    ab <- tryCatch(fit_age_baseline(o), error = function(e) NULL)
    if (!is.null(ab)) {
      mat["LR age-only", om] <- cell(c(
        auc = ab$ci[["estimate"]], lower = ab$ci[["lower"]], upper = ab$ci[["upper"]]
      ))
    }
  }
  df <- data.frame(Model = rownames(mat), mat, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- NULL
  names(df)[-1] <- tools::toTitleCase(names(df)[-1])
  df
}

# Draw a booktabs-style three-line table (top rule / header rule / bottom rule)
# of the AUC (95% CI) wide table onto the currently active graphics device.
# Body only: the number, title and legend live in main_tables.tex. `title`
# stays as an argument so a caller can still label an ad hoc rendering, but
# the released table is drawn without one.
draw_auc_ci_three_line <- function(df, title = NULL) {
  n_row <- nrow(df)
  n_data <- ncol(df) - 1L
  grid::grid.newpage()

  left <- 0.06
  right <- 0.96
  data_x <- seq(0.40, 0.84, length.out = n_data)
  y_top <- 0.80
  y_head <- 0.75
  y_mid <- 0.71
  row_h <- 0.095

  if (!is.null(title)) {
    grid::grid.text(title, x = left, y = 0.90, just = "left",
                    gp = grid::gpar(fontsize = 12, fontface = "bold"))
  }
  grid::grid.lines(x = c(left, right), y = c(y_top, y_top), gp = grid::gpar(lwd = 1.4))
  grid::grid.text("Model", x = left, y = y_head, just = "left",
                  gp = grid::gpar(fontsize = 9, fontface = "bold"))
  for (j in seq_len(n_data)) {
    grid::grid.text(names(df)[j + 1L], x = data_x[j], y = y_head, just = "centre",
                    gp = grid::gpar(fontsize = 9, fontface = "bold"))
  }
  grid::grid.lines(x = c(left, right), y = c(y_mid, y_mid), gp = grid::gpar(lwd = 0.8))

  for (i in seq_len(n_row)) {
    yy <- y_mid - i * row_h
    grid::grid.text(df$Model[i], x = left, y = yy, just = "left",
                    gp = grid::gpar(fontsize = 8.5))
    for (j in seq_len(n_data)) {
      grid::grid.text(df[i, j + 1L], x = data_x[j], y = yy, just = "centre",
                      gp = grid::gpar(fontsize = 8.5))
    }
  }
  y_bottom <- y_mid - (n_row + 0.5) * row_h
  grid::grid.lines(x = c(left, right), y = c(y_bottom, y_bottom), gp = grid::gpar(lwd = 1.4))
  grid::grid.text(
    "AUC with 95% DeLong confidence interval on the held-out test set. LR age-only = logistic regression on age only.",
    x = left, y = y_bottom - 0.06, just = "left", gp = grid::gpar(fontsize = 7)
  )
  invisible(TRUE)
}

# Booktabs-style three-line table of the test-set AUC with 95% CI, one column per
# omics, written to a PDF file.
write_auc_ci_three_line_pdf <- function(model_outputs, output_pdf,
                                        title = NULL) {
  df <- build_auc_ci_wide(model_outputs)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf, width = 9, height = 4.6)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_auc_ci_three_line(df, title = title)
  TRUE
}

format_feature_label <- function(x, context = c("bubble", "supplement")) {
  context <- match.arg(context)
  x <- as.character(x)
  clinical_bubble <- c(
    age_v4 = "age", race = "Race", edu = "Education", ol_health = "Self-rated health",
    bmi = "BMI", mstat = "Marital status", smoke = "Smoking status", diab = "Diabetes",
    hbp = "Hypertension", cancer = "Cancer", tmm_score = "Teng 3MS score",
    gds = "GDS-15", pase = "PASE", total_meds = "Total medications",
    site = "Clinical site", batch = "Batch"
  )
  clinical_supplement <- c(
    age_v4 = "Age", race = "Race", edu = "Education", ol_health = "Self-rated health",
    bmi = "BMI", mstat = "Marital status", smoke = "Smoking status", diab = "Diabetes",
    hbp = "Hypertension", cancer = "Cancer", tmm_score = "Teng mMMSE score",
    gds = "Depressive symptoms (GDS)", pase = "Physical activity (PASE)",
    total_meds = "Total medications", site = "Clinical site", batch = "Batch"
  )
  clinical <- if (context == "bubble") clinical_bubble else clinical_supplement
  out <- clinical[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

feature_class <- function(x) {
  clinical <- c(
    "age_v4", "race", "edu", "ol_health", "bmi", "mstat", "smoke", "diab",
    "hbp", "cancer", "tmm_score", "gds", "pase", "total_meds", "site", "batch"
  )
  ifelse(x %in% clinical, "Clinical", "Omics")
}

select_top_shap_features <- function(shap_data, n = 20) {
  shap_data |>
    dplyr::group_by(feature) |>
    dplyr::summarise(mean_abs = mean(abs(value), na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mean_abs)) |>
    dplyr::slice_head(n = n) |>
    dplyr::pull(feature)
}

shap_summary <- function(shap_data, model_name, top_n = 20) {
  shap_data |>
    dplyr::group_by(feature) |>
    dplyr::summarise(
      importance = mean(abs(value), na.rm = TRUE),
      direction_value = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(importance)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::mutate(
      model = model_name,
      direction = ifelse(direction_value >= 0, "Higher risk", "Lower risk"),
      class = feature_class(feature)
    )
}

en_summary <- function(model_output, model_name = "EN", top_n = 20) {
  model_output$coef_en |>
    dplyr::filter(term != "(Intercept)", estimate != 0) |>
    dplyr::arrange(dplyr::desc(abs(estimate))) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::transmute(
      feature = term,
      importance = abs(estimate),
      direction_value = estimate,
      model = model_name,
      direction = ifelse(estimate >= 0, "Higher risk", "Lower risk"),
      class = "Omics"
    )
}

en_coefficient_data <- function(model_output) {
  if (!is.null(model_output$plots$p_en_coef$data)) {
    return(model_output$plots$p_en_coef$data)
  }
  model_output$coef_en
}

saved_shap_data <- function(model_output, plot_name) {
  plot_obj <- model_output$plots[[plot_name]]
  if (is.null(plot_obj) || is.null(plot_obj$data)) {
    stop("Missing saved plot object: plots$", plot_name)
  }
  required <- c("value", "feature", "color")
  missing <- setdiff(required, names(plot_obj$data))
  if (length(missing)) {
    stop("Saved plot object ", plot_name, " lacks columns: ", paste(missing, collapse = ", "))
  }
  plot_obj$data
}

top_bar_features <- function(model_output, plot_name, top_n = 20) {
  bar <- as.data.frame(model_output$plots[[plot_name]]$data)
  bar$feature <- as.character(bar$feature)
  bar <- bar[order(bar$value, decreasing = TRUE), , drop = FALSE]
  bar <- bar[seq_len(min(top_n, nrow(bar))), , drop = FALSE]
  bar$rank <- seq_len(nrow(bar))
  bar
}

top_en_features <- function(model_output, top_n = 20) {
  en <- model_output$coef_en
  en$term <- as.character(en$term)
  en <- en[en$term != "(Intercept)" & en$estimate != 0, , drop = FALSE]
  en <- en[order(abs(en$estimate), decreasing = TRUE), , drop = FALSE]
  en <- en[seq_len(min(top_n, nrow(en))), , drop = FALSE]
  data.frame(feature = en$term, value = abs(en$estimate), signed = en$estimate, rank = seq_len(nrow(en)))
}

shap_direction <- function(model_output, plot_name, features, flip = FALSE) {
  bee <- saved_shap_data(model_output, plot_name)
  bee$feature <- as.character(bee$feature)
  features <- as.character(features)
  out <- bee |>
    dplyr::filter(feature %in% features) |>
    dplyr::group_by(feature) |>
    dplyr::summarise(direction_value = mean(value, na.rm = TRUE), .groups = "drop")
  if (flip) {
    out$direction_value <- -out$direction_value
  }
  out
}

overlap_data <- function(model_output, left, right, top_n = 20) {
  if (left == "EN") {
    left_tbl <- top_en_features(model_output, top_n)
  } else {
    left_tbl <- top_bar_features(model_output, left, top_n)
    names(left_tbl)[names(left_tbl) == "value"] <- "importance"
    left_tbl <- data.frame(feature = left_tbl$feature, value = left_tbl$importance, signed = left_tbl$importance, rank = left_tbl$rank)
  }
  right_tbl <- top_bar_features(model_output, right, top_n)
  shared <- intersect(left_tbl$feature, right_tbl$feature)
  left_tbl <- left_tbl[left_tbl$feature %in% shared, , drop = FALSE]
  right_tbl <- right_tbl[right_tbl$feature %in% shared, , drop = FALSE]
  right_dir <- shap_direction(model_output, sub("_bar$", "_bee", right), shared, flip = identical(right, "p_shap3_bar"))
  left_dir <- if (left == "EN") {
    data.frame(feature = left_tbl$feature, direction_value = left_tbl$signed)
  } else {
    shap_direction(model_output, sub("_bar$", "_bee", left), shared)
  }
  left_plot <- data.frame(
    feature = left_tbl$feature,
    model = if (left == "EN") "EN" else "RF+EN",
    importance = left_tbl$value,
    direction_value = left_dir$direction_value[match(left_tbl$feature, left_dir$feature)]
  )
  right_plot <- data.frame(
    feature = right_tbl$feature,
    model = if (right == "p_shap3_bar") "XGB" else "RF+XGB",
    importance = right_tbl$value,
    direction_value = right_dir$direction_value[match(right_tbl$feature, right_dir$feature)]
  )
  direction_by_feature <- left_plot$direction_value[match(right_plot$feature, left_plot$feature)]
  right_plot$direction_value <- direction_by_feature
  rank_tbl <- merge(
    left_tbl[, c("feature", "rank"), drop = FALSE],
    right_tbl[, c("feature", "rank"), drop = FALSE],
    by = "feature",
    suffixes = c("_left", "_right")
  )
  rank_tbl$mean_rank <- rowMeans(rank_tbl[, c("rank_left", "rank_right")])
  rank_tbl <- rank_tbl[order(rank_tbl$mean_rank, rank_tbl$rank_left), , drop = FALSE]
  feature_order <- as.character(rank_tbl$feature)
  list(data = dplyr::bind_rows(left_plot, right_plot), features = feature_order)
}

make_bubble_panel <- function(overlap, title, show_class = FALSE, context = "bubble") {
  plot_data <- overlap$data
  plot_data$feature <- as.character(plot_data$feature)
  overlap$features <- as.character(overlap$features)
  plot_data$direction <- ifelse(plot_data$direction_value < 0, "Higher mortality risk (adverse)", "Lower mortality risk (protective)")
  plot_data$class <- feature_class(plot_data$feature)
  feature_labels <- format_feature_label(overlap$features, context = context)
  if (show_class && requireNamespace("ggtext", quietly = TRUE)) {
    feature_classes <- feature_class(overlap$features)
    feature_labels <- ifelse(
      feature_classes == "Omics",
      paste0("<span style='color:#008b45'>", feature_labels, "</span>"),
      feature_labels
    )
    label_lookup <- stats::setNames(feature_labels, overlap$features)
  } else {
    label_lookup <- stats::setNames(feature_labels, overlap$features)
  }
  plot_data$feature_label <- factor(label_lookup[plot_data$feature], levels = rev(feature_labels))

  base <- ggplot2::ggplot(plot_data, ggplot2::aes(x = model, y = feature_label, size = importance, fill = direction))
  if (show_class) {
    base <- base +
      ggplot2::geom_point(ggplot2::aes(shape = class), color = "grey25", alpha = 0.95, stroke = 0.25) +
      ggplot2::scale_shape_manual(
        values = c(Clinical = 21, Omics = 22),
        labels = c(Clinical = "Clinical", Omics = "Omics (microbiome/\nmetabolome/proteome)")
      ) +
      ggplot2::labs(shape = "Feature class")
  } else {
    base <- base + ggplot2::geom_point(shape = 21, color = "grey25", alpha = 0.95, stroke = 0.25)
  }

  axis_text_y <- if (show_class && requireNamespace("ggtext", quietly = TRUE)) {
    ggtext::element_markdown(size = 8, colour = "black")
  } else {
    ggplot2::element_text(size = 8, colour = "black")
  }

  base +
    ggplot2::scale_fill_manual(values = c("Higher mortality risk (adverse)" = "#B64B4B", "Lower mortality risk (protective)" = "#3B6EA8")) +
    ggplot2::scale_size_continuous(
      range = c(2.3, 8.5),
      breaks = function(x) seq(x[1], x[2], length.out = 4),
      labels = function(x) c("Low", "Medium", "High", "Very high")[seq_along(x)]
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL,
      fill = if (show_class) "Mortality risk" else "Direction",
      size = "Relative feature importance"
    ) +
    ggplot2::theme_classic(base_size = 9, base_family = "Helvetica") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 10),
      axis.text.x = ggplot2::element_text(face = "bold", size = 9),
      axis.text.y = axis_text_y,
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(size = 7),
      legend.title = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(4, 8, 4, 8)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(order = 1, nrow = 2, override.aes = list(size = 4)),
      shape = ggplot2::guide_legend(order = 2, nrow = 2, override.aes = list(size = 4)),
      size = ggplot2::guide_legend(order = 3, nrow = 1)
    )
}

# =============================================================================
# Bubble figures (Fig3 / Fig4) ported verbatim from the figures_new logic:
#   figures_new/all_omics_bubble_shared.R  +  figures_new/fig3_xgb_beeswarm_visual.R
# Panels are rebuilt from the raw data frames in the RDS (coef_en / p_shapX_bar /
# p_shapX_bee) so output matches Figure3_all_omics_bubble.Rmd /
# Figure4_all_omics_bubble.Rmd exactly. Clinical labels/classification reuse the
# helpers defined further below (clinical_id_labels, is_clinical_id, ...).
# =============================================================================

col_adv <- "#E64B35"
col_pro <- "#3C5488"
omics_label_green <- "#2E7D32"
col_omics_label <- omics_label_green

canonicalize_preprocess_feature_id <- function(ch) {
  ch <- as.character(ch)
  ifelse(ch == "age", "age_v4", ch)
}
canonicalize_preprocess_feature_col <- function(v) {
  if (is.factor(v)) {
    levels(v) <- canonicalize_preprocess_feature_id(levels(v))
    v
  } else {
    canonicalize_preprocess_feature_id(v)
  }
}
canonicalize_clinical_ids_df <- function(d) {
  if (is.null(d) || nrow(d) == 0) return(d)
  d <- as.data.frame(d)
  for (nm in intersect(names(d), c("feature", "term", "variable"))) {
    d[[nm]] <- canonicalize_preprocess_feature_col(d[[nm]])
  }
  d
}

sanitize_pipeline_obj <- function(obj) {
  if (!is.null(obj$coef_en)) obj$coef_en <- canonicalize_clinical_ids_df(obj$coef_en)
  if (!is.null(obj$plots)) {
    for (nm in names(obj$plots)) {
      p <- obj$plots[[nm]]
      if (inherits(p, "ggplot") && !is.null(p$data)) {
        obj$plots[[nm]]$data <- canonicalize_clinical_ids_df(p$data)
      }
    }
  }
  obj
}

pretty_feature_label <- function(x, max_len) {
  x <- as.character(x)
  max_len <- as.integer(max_len)
  out <- vapply(seq_along(x), function(i) {
    xi <- x[[i]]
    if (is_preprocess_clinical_id(xi)) pretty_preprocess_feature(xi, max_len) else xi
  }, character(1))
  long <- nchar(out) > max_len
  out[long] <- paste0(substr(out[long], 1L, max_len - 1L), "\u2026")
  out
}

shap_risk_dir <- function(v) ifelse(v >= 0, "Protective", "Adverse")
assoc_sign_dir <- function(v) ifelse(v >= 0, "Positive", "Negative")

quartile_shap_direction <- function(fv, sh, label_style) {
  label_style <- match.arg(label_style, c("mortality", "association"))
  ok <- stats::complete.cases(fv, sh)
  fv <- fv[ok]; sh <- sh[ok]
  m_fallback <- mean(sh, na.rm = TRUE)
  if (length(fv) < 4L) {
    return(if (label_style == "mortality") shap_risk_dir(m_fallback) else assoc_sign_dir(m_fallback))
  }
  q_low <- stats::quantile(fv, 0.25, na.rm = TRUE, names = FALSE)
  q_high <- stats::quantile(fv, 0.75, na.rm = TRUE, names = FALSE)
  if (!is.finite(q_low) || !is.finite(q_high) || q_low >= q_high) {
    return(if (label_style == "mortality") shap_risk_dir(m_fallback) else assoc_sign_dir(m_fallback))
  }
  low_shap <- mean(sh[fv <= q_low], na.rm = TRUE)
  high_shap <- mean(sh[fv >= q_high], na.rm = TRUE)
  if (!is.finite(low_shap) || !is.finite(high_shap)) {
    return(if (label_style == "mortality") shap_risk_dir(m_fallback) else assoc_sign_dir(m_fallback))
  }
  delta <- high_shap - low_shap
  if (label_style == "association") {
    if (delta > 0) return("Positive")
    return("Negative")
  }
  if (delta > 0) return("Protective")
  "Adverse"
}

beeswarm_feature_value_col <- function(bee_df, allow_shapviz_color = TRUE) {
  d <- as.data.frame(bee_df)
  nm <- names(d)
  color_names <- c("color", "colour", "Color", "Colour")
  cand_priority <- if (allow_shapviz_color) c(color_names, "feature_value") else c("feature_value")
  for (cand in cand_priority) {
    if (cand %in% nm) return(cand)
  }
  reserved <- c("feature", "value", "Var1", "Var2")
  rest <- setdiff(nm, reserved)
  rest <- rest[!grepl("^\\.\\.\\.", rest) & !grepl("^\\.", rest)]
  if (!allow_shapviz_color) rest <- rest[!tolower(rest) %in% tolower(color_names)]
  hits <- character(0)
  for (col in rest) {
    v <- d[[col]]
    if (!is.numeric(v) && !is.integer(v)) next
    if (all(is.na(v))) next
    rng <- range(as.numeric(v), na.rm = TRUE)
    if (is.finite(rng[1]) && is.finite(rng[2]) && rng[1] >= -1e-6 && rng[2] <= 1 + 1e-6) {
      hits <- c(hits, col)
    }
  }
  if (length(hits) == 1L) return(hits[[1]])
  NULL
}

direction_from_shapviz_bee <- function(bee, label_style = c("mortality", "association"),
                                       allow_shapviz_color = TRUE) {
  label_style <- match.arg(label_style)
  bee <- as.data.frame(bee)
  bee$feature <- as.character(bee$feature)
  fv_col <- beeswarm_feature_value_col(bee, allow_shapviz_color = allow_shapviz_color)
  if (is.null(fv_col)) {
    ms <- tapply(bee$value, bee$feature, mean, na.rm = TRUE)
    dlab <- if (label_style == "mortality") {
      shap_risk_dir(as.numeric(unname(ms)))
    } else {
      assoc_sign_dir(as.numeric(unname(ms)))
    }
    return(tibble::tibble(feature = names(ms), direction = dlab))
  }
  spl <- split(bee, bee$feature)
  feat_names <- names(spl)
  direction <- vapply(spl, function(sub) {
    quartile_shap_direction(sub[[fv_col]], sub$value, label_style)
  }, FUN.VALUE = character(1))
  tibble::tibble(feature = feat_names, direction = unname(direction))
}

# Fig.3 XGB direction from the *rendered* beeswarm (ggplot2 point colours).
xgb_direction_from_beeswarm_visual <- function(p_bee, xgb_threshold = 0.002) {
  if (!inherits(p_bee, "ggplot")) {
    stop("p_bee must be a ggplot object (e.g. obj$plots$p_shap3_bee).", call. = FALSE)
  }
  gb <- ggplot2::ggplot_build(p_bee)
  ld <- NULL
  for (i in seq_along(gb$data)) {
    d <- gb$data[[i]]
    if (is.null(d) || nrow(d) < 2L) next
    if (all(c("x", "colour") %in% names(d))) { ld <- d; break }
  }
  if (is.null(ld)) {
    stop("xgb_direction_from_beeswarm_visual: no point layer with x and colour in ggplot_build().",
         call. = FALSE)
  }
  pd <- as.data.frame(p_bee$data)
  if (nrow(pd) != nrow(ld)) {
    stop("xgb_direction_from_beeswarm_visual: nrow(plot$data) (", nrow(pd),
         ") != built layer (", nrow(ld), ").", call. = FALSE)
  }
  pd$shap_value <- ld$x
  pd$colour_hex <- as.character(ld$colour)
  pd$feature <- as.character(pd$feature)
  red_score <- rep(NA_real_, nrow(pd))
  okc <- !is.na(pd$colour_hex) & nzchar(pd$colour_hex)
  if (any(okc)) {
    rgb_mat <- grDevices::col2rgb(pd$colour_hex[okc], alpha = FALSE)
    red_score[okc] <- rgb_mat[1, ] - rgb_mat[3, ]
  }
  dd <- dplyr::mutate(pd, red_score = red_score)
  dd |>
    dplyr::group_by(feature) |>
    dplyr::summarise(
      mean_shap_red = mean(shap_value[red_score >= stats::quantile(red_score, 0.75, na.rm = TRUE)], na.rm = TRUE),
      mean_shap_blue = mean(shap_value[red_score <= stats::quantile(red_score, 0.25, na.rm = TRUE)], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      delta = mean_shap_red - mean_shap_blue,
      direction_xgb = dplyr::case_when(
        is.na(delta) ~ "Unknown",
        delta > xgb_threshold ~ "Protective",
        delta < -xgb_threshold ~ "Adverse",
        TRUE ~ "Neutral"
      )
    )
}

extract_fig3 <- function(obj, top_n = 20, flip_model3_xgb_shap_sign = FALSE) {
  obj <- sanitize_pipeline_obj(obj)
  en <- obj$coef_en |>
    dplyr::filter(estimate != 0) |>
    dplyr::mutate(
      abs_est   = abs(estimate),
      direction = ifelse(estimate > 0, "Protective", "Adverse")
    ) |>
    dplyr::arrange(dplyr::desc(abs_est)) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::select(feature = term, estimate, abs_est, direction_en = direction, rank_en = rank)

  bar <- as.data.frame(obj$plots$p_shap3_bar$data) |>
    dplyr::mutate(feature = as.character(feature))
  p_bee <- obj$plots$p_shap3_bee
  if (isTRUE(flip_model3_xgb_shap_sign)) {
    bd <- as.data.frame(p_bee$data)
    if (!"value" %in% names(bd)) {
      stop("extract_fig3: p_shap3_bee$data has no `value` column (SHAP).", call. = FALSE)
    }
    bd$value <- -as.numeric(bd$value)
    # ggplot2 >= 4.0: replace plot data with `+ data.frame` ( `%+%` deprecated ).
    p_bee <- p_bee + bd
  }
  sign_xgb <- xgb_direction_from_beeswarm_visual(p_bee, xgb_threshold = 0.002) |>
    dplyr::select(feature, direction_xgb)
  xgb <- bar |>
    dplyr::rename(mean_abs_shap = value) |>
    dplyr::arrange(dplyr::desc(mean_abs_shap)) |>
    dplyr::mutate(rank_xgb = dplyr::row_number()) |>
    dplyr::left_join(sign_xgb, by = "feature")

  list(en = en, xgb = xgb, top_n = top_n)
}

extract_fig4 <- function(obj, top_n = 20) {
  obj <- sanitize_pipeline_obj(obj)
  get_shap <- function(bar_key, bee_key) {
    bar <- as.data.frame(obj$plots[[bar_key]]$data) |>
      dplyr::mutate(feature = as.character(feature))
    bee <- as.data.frame(obj$plots[[bee_key]]$data) |>
      dplyr::mutate(feature = as.character(feature))
    sign <- direction_from_shapviz_bee(bee, label_style = "mortality")
    bar |>
      dplyr::rename(mean_abs_shap = value) |>
      dplyr::arrange(dplyr::desc(mean_abs_shap)) |>
      dplyr::mutate(rank = dplyr::row_number()) |>
      dplyr::left_join(sign |> dplyr::select(feature, direction), by = "feature")
  }
  m4 <- get_shap("p_shap4_bar", "p_shap4_bee")
  m5 <- get_shap("p_shap5_bar", "p_shap5_bee")
  list(m4 = m4, m5 = m5, top_n = top_n)
}

bubble_fig3 <- function(d, omics_label, label_max_len) {
  inter <- intersect(d$en$feature[d$en$rank_en <= d$top_n],
                     d$xgb$feature[d$xgb$rank_xgb <= d$top_n])
  if (length(inter) == 0) return(NULL)

  en_f  <- d$en  |> dplyr::filter(feature %in% inter)
  xgb_f <- d$xgb |> dplyr::filter(feature %in% inter)

  avg_r <- dplyr::full_join(
    en_f  |> dplyr::select(feature, r1 = rank_en),
    xgb_f |> dplyr::select(feature, r2 = rank_xgb), by = "feature"
  ) |> dplyr::mutate(avg_r = (r1 + r2) / 2)

  ord <- avg_r |> dplyr::arrange(dplyr::desc(avg_r))
  labs <- make.unique(pretty_feature_label(ord$feature, label_max_len), sep = " ")
  lab_map <- stats::setNames(labs, ord$feature)

  en_long <- en_f |>
    dplyr::transmute(feature, model = "EN", importance = abs_est,
                     direction = as.character(direction_en), rank = rank_en)
  xgb_long <- xgb_f |>
    dplyr::transmute(
      feature, model = "XGB", importance = mean_abs_shap,
      direction = dplyr::coalesce(as.character(direction_xgb), "Unknown"),
      rank = rank_xgb
    )

  dplyr::bind_rows(en_long, xgb_long) |>
    dplyr::mutate(
      label     = factor(lab_map[feature], levels = labs),
      model     = factor(model, levels = c("EN", "XGB")),
      direction = factor(direction, levels = c("Adverse", "Protective", "Neutral", "Unknown")),
      omics     = omics_label
    ) |>
    dplyr::group_by(model) |>
    dplyr::mutate(imp_norm = importance / max(importance, na.rm = TRUE)) |>
    dplyr::ungroup()
}

bubble_fig4 <- function(d, omics_label, label_max_len) {
  inter <- intersect(d$m4$feature[d$m4$rank <= d$top_n],
                     d$m5$feature[d$m5$rank <= d$top_n])
  if (length(inter) == 0) return(NULL)

  feat_type <- tibble::tibble(
    feature   = inter,
    feat_type = ifelse(is_clinical_id(inter), "Clinical", "Omics")
  )

  m4f <- d$m4 |> dplyr::filter(feature %in% inter)
  m5f <- d$m5 |> dplyr::filter(feature %in% inter)

  avg_r <- dplyr::full_join(
    m4f |> dplyr::select(feature, r1 = rank),
    m5f |> dplyr::select(feature, r2 = rank), by = "feature"
  ) |> dplyr::mutate(avg_r = (r1 + r2) / 2) |>
    dplyr::left_join(feat_type, by = "feature") |>
    dplyr::mutate(feat_type = factor(feat_type, levels = c("Omics", "Clinical"))) |>
    dplyr::arrange(feat_type, dplyr::desc(avg_r))

  labs <- make.unique(pretty_feature_label(avg_r$feature, label_max_len), sep = " ")
  lab_map <- stats::setNames(labs, avg_r$feature)
  lab_lev <- labs

  m4_long <- m4f |>
    dplyr::transmute(feature, model = "RF+EN",  importance = mean_abs_shap, direction, rank)
  m5_long <- m5f |>
    dplyr::transmute(feature, model = "RF+XGB", importance = mean_abs_shap, direction, rank)

  dplyr::bind_rows(m4_long, m5_long) |>
    dplyr::mutate(
      label     = factor(lab_map[feature], levels = lab_lev),
      model     = factor(model, levels = c("RF+EN", "RF+XGB")),
      direction = factor(direction, levels = c("Adverse", "Protective")),
      omics     = omics_label,
      feat_type = ifelse(is_clinical_id(feature), "Clinical", "Omics")
    ) |>
    dplyr::group_by(model) |>
    dplyr::mutate(imp_norm = importance / max(importance, na.rm = TRUE)) |>
    dplyr::ungroup()
}

draw_panel <- function(df, col_a, col_b, title_str, base_sz = 11,
                       show_omics_label = TRUE, omics_label_txt = NULL,
                       margin_left = 14, margin_right = NULL,
                       omics_y_axis_green = TRUE,
                       direction_scale_style = c("fig4", "fig3")) {
  direction_scale_style <- match.arg(direction_scale_style)
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot2::ggplot() +
             ggplot2::labs(title = title_str, subtitle = "Data not available") +
             ggplot2::theme_void(base_size = base_sz) +
             ggplot2::theme(plot.title = ggplot2::element_text(
               face = "bold", size = base_sz + 2, hjust = 0.5)))
  }

  lab_ft <- df |>
    dplyr::mutate(label_chr = as.character(label)) |>
    dplyr::distinct(label_chr, feat_type)
  lvl <- levels(df$label)
  ft_lvl <- lab_ft$feat_type[match(lvl, lab_ft$label_chr)]
  ft_lvl[is.na(ft_lvl)] <- "Clinical"
  axis_y_colour <- ifelse(ft_lvl == "Omics" & isTRUE(omics_y_axis_green),
                          omics_label_green, "black")
  axis_y_face   <- ifelse(ft_lvl == "Omics", "italic", "plain")

  all_labs  <- levels(df$label)
  omics_pos <- which(all_labs %in% as.character(unique(df$label[df$feat_type == "Omics"])))
  clin_pos  <- which(all_labs %in% as.character(unique(df$label[df$feat_type == "Clinical"])))
  div_y     <- if (length(omics_pos) > 0 && length(clin_pos) > 0) {
    (max(omics_pos) + min(clin_pos)) / 2
  } else {
    NULL
  }

  right_margin <- if (show_omics_label) 26 else 10
  right_margin <- if (!is.null(margin_right)) margin_right else right_margin

  fill_scale <- if (direction_scale_style == "fig3") {
    ggplot2::scale_fill_manual(
      values = c(Adverse = col_adv, Protective = col_pro, Neutral = "#BDBDBD", Unknown = "#9E9E9E"),
      name = "Direction",
      labels = c("Higher mortality risk (adverse)", "Lower mortality risk (protective)",
                 "No clear direction", "Unknown"),
      drop = FALSE
    )
  } else {
    ggplot2::scale_fill_manual(
      values = c(Adverse = col_adv, Protective = col_pro, Unknown = "#9E9E9E"),
      name   = "Direction (mortality risk)",
      labels = c("Higher mortality risk (adverse)", "Lower mortality risk (protective)", "Unavailable"),
      drop = TRUE
    )
  }
  col_scale <- if (direction_scale_style == "fig3") {
    ggplot2::scale_color_manual(
      values = c(Adverse = col_adv, Protective = col_pro, Neutral = "#BDBDBD", Unknown = "#9E9E9E"),
      guide = "none", drop = FALSE
    )
  } else {
    ggplot2::scale_color_manual(
      values = c(Adverse = col_adv, Protective = col_pro, Unknown = "#9E9E9E"),
      guide = "none", drop = TRUE
    )
  }

  ggplot2::ggplot(df,
    ggplot2::aes(x = model, y = label, size = imp_norm, fill = direction, color = direction)) +
    ggplot2::geom_point(shape = 21, stroke = 0.45, alpha = 0.90) +
    {if (!is.null(div_y))
      ggplot2::geom_hline(yintercept = div_y, linetype = "dashed",
                          colour = "grey65", linewidth = 0.35)} +
    {if (show_omics_label && length(clin_pos) > 0)
      ggplot2::annotate("text", x = 2.62, y = mean(range(clin_pos)),
        label = "Clinical", angle = -90, size = 3.55,
        colour = "grey35", fontface = "italic")} +
    {if (show_omics_label && length(omics_pos) > 0 && !is.null(omics_label_txt))
      ggplot2::annotate("text", x = 2.62, y = mean(range(omics_pos)),
        label = omics_label_txt, angle = -90, size = 3.55,
        colour = col_omics_label, fontface = "italic")} +
    fill_scale +
    col_scale +
    ggplot2::scale_size_continuous(
      range  = c(2.1, 8.5),
      name   = "Relative feature importance",
      breaks = c(0.25, 0.5, 0.75, 1.0),
      labels = c("Low", "Medium", "High", "Very high")
    ) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(title = title_str, x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_classic(base_size = base_sz) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", size = base_sz + 2.5, hjust = 0.5),
      axis.line    = ggplot2::element_blank(),
      axis.ticks   = ggplot2::element_blank(),
      axis.text.x.top = ggplot2::element_text(face = "bold", size = base_sz + 2),
      axis.text.y  = ggplot2::element_text(size = base_sz + 0.15, colour = axis_y_colour, face = axis_y_face),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90", linewidth = 0.25),
      panel.grid.major.x = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box  = "horizontal",
      legend.title = ggplot2::element_text(size = base_sz + 0.35, face = "bold"),
      legend.text  = ggplot2::element_text(size = base_sz - 0.35),
      legend.key   = ggplot2::element_rect(fill = NA, colour = NA),
      plot.margin  = ggplot2::margin(6, right_margin, 6, margin_left)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(override.aes = list(size = 5), order = 1, nrow = 1),
      size = ggplot2::guide_legend(order = 2, nrow = 1)
    )
}

# The direction key lists only the categories the panels actually use.
#
# It used to hard-code all four levels, so "No clear direction" and "Unknown"
# were advertised even when -- as in the published figure -- every shared top-20
# feature has a definite direction and no grey bubble is ever drawn. A key
# describing swatches that do not appear reads as a rendering fault. Callers pass
# the levels present in the plotted data; the swatches keep their original pitch,
# so a figure that genuinely used all four is unchanged.
make_legend_fig3_standalone <- function(directions_present = c("Adverse", "Protective",
                                                              "Neutral", "Unknown")) {
  spec <- tibble::tibble(
    level = c("Adverse", "Protective", "Neutral", "Unknown"),
    lab = c(
      "Higher mortality risk\n(adverse)",
      "Lower mortality risk\n(protective)",
      "No clear direction",
      "Unknown"
    ),
    fill_col = c(col_adv, col_pro, "#BDBDBD", "#9E9E9E")
  )
  spec <- spec[spec$level %in% directions_present, , drop = FALSE]
  if (!nrow(spec)) spec <- spec[1, , drop = FALSE]
  y_dir <- 1.55
  dir3_df <- tibble::tibble(
    x = 0.52 + 0.86 * (seq_len(nrow(spec)) - 1),
    y = rep(y_dir, nrow(spec)),
    lab = spec$lab,
    fill_col = spec$fill_col
  )
  imp_df <- tibble::tibble(
    x = c(4.85, 5.63, 6.41, 7.19), y = rep(1.42, 4),
    lab = c("Low", "Medium", "High", "Very high"), sz = c(1.2, 2.5, 4.2, 6.5)
  )
  ggplot2::ggplot() +
    ggplot2::geom_point(data = dir3_df, ggplot2::aes(x = x, y = y), fill = dir3_df$fill_col,
      shape = 21, size = 4.5, colour = "white", stroke = 0.2) +
    ggplot2::geom_text(data = dir3_df, ggplot2::aes(x = x + 0.20, y = y, label = lab),
      hjust = 0, size = 3.15, family = "Helvetica") +
    ggplot2::geom_point(data = imp_df, ggplot2::aes(x = x, y = y, size = sz),
      shape = 21, fill = "grey40", colour = "grey40", stroke = 0.2, alpha = 0.85) +
    ggplot2::geom_text(data = imp_df, ggplot2::aes(x = x, y = y - 0.42, label = lab),
      hjust = 0.5, size = 3.55, family = "Helvetica") +
    ggplot2::scale_size_identity() +
    ggplot2::annotate("text", x = mean(dir3_df$x) + 0.04, y = 3.0, label = "Direction",
                      hjust = 0.5, size = 3.65, fontface = "bold") +
    ggplot2::annotate("text", x = 6.0, y = 3.0, label = "Relative feature importance",
                      hjust = 0.5, size = 3.65, fontface = "bold") +
    ggplot2::xlim(0.35, 8.0) + ggplot2::ylim(0.2, 3.5) +
    ggplot2::theme_void()
}

make_legend_fig4_standalone <- function() {
  feat_df <- tibble::tibble(
    x = c(0.88, 0.88), y = c(2, 1),
    lab = c("Clinical", "Omics (microbiome/\nmetabolome/proteome)"),
    fill_col = c("grey80", col_omics_label)
  )
  dir4_df <- tibble::tibble(
    x = c(2.5, 2.5), y = c(2, 1),
    lab = c("Lower mortality risk\n(protective)", "Higher mortality risk\n(adverse)"),
    fill_col = c(col_pro, col_adv)
  )
  imp_df <- tibble::tibble(
    x = c(4.6, 5.38, 6.16, 6.94), y = rep(1.42, 4),
    lab = c("Low", "Medium", "High", "Very high"), sz = c(1.2, 2.5, 4.2, 6.5)
  )
  ggplot2::ggplot() +
    ggplot2::geom_tile(data = feat_df, ggplot2::aes(x = x, y = y), fill = feat_df$fill_col,
      width = 0.28, height = 0.28, colour = "grey50", linewidth = 0.3) +
    ggplot2::geom_text(data = feat_df, ggplot2::aes(x = x + 0.20, y = y, label = lab),
      hjust = 0, size = 3.55, family = "Helvetica") +
    ggplot2::geom_point(data = dir4_df, ggplot2::aes(x = x, y = y), fill = dir4_df$fill_col,
      shape = 21, size = 4.5, colour = "white", stroke = 0.2) +
    ggplot2::geom_text(data = dir4_df, ggplot2::aes(x = x + 0.20, y = y, label = lab),
      hjust = 0, size = 3.55, family = "Helvetica") +
    ggplot2::geom_point(data = imp_df, ggplot2::aes(x = x, y = y, size = sz),
      shape = 21, fill = "grey40", colour = "grey40", stroke = 0.2, alpha = 0.85) +
    ggplot2::geom_text(data = imp_df, ggplot2::aes(x = x, y = y - 0.42, label = lab),
      hjust = 0.5, size = 3.55, family = "Helvetica") +
    ggplot2::scale_size_identity() +
    ggplot2::annotate("text", x = 0.95, y = 3.0, label = "Feature class",
                      hjust = 0.5, size = 3.65, fontface = "bold") +
    ggplot2::annotate("text", x = 2.5, y = 3.0, label = "Mortality risk",
                      hjust = 0.5, size = 3.65, fontface = "bold") +
    ggplot2::annotate("text", x = 5.75, y = 3.0, label = "Relative feature importance",
                      hjust = 0.5, size = 3.65, fontface = "bold") +
    ggplot2::xlim(0.45, 7.45) + ggplot2::ylim(0.2, 3.5) +
    ggplot2::theme_void()
}

write_bubble_figure <- function(model_outputs, output_file, figure = c("fig3", "fig4"),
                                flip_model3_xgb_shap_sign = TRUE) {
  # Directions match figures_new bubble logic; manuscript PDFs are panels + legend
  # only (no figure title / subtitle block).
  figure <- match.arg(figure)
  top_n <- 20L
  omics_order <- c("microbiome", "metabolomics", "proteomics")
  labels_omics <- c(microbiome = "Microbiome", metabolomics = "Metabolomics",
                    proteomics = "Proteomics")
  label_max <- c(microbiome = 72L, metabolomics = 72L, proteomics = 96L)
  panel_rel_widths <- c(1, 1, 1.48)
  # Panel letters are drawn by cowplot in the top-left corner of each panel
  # (Nature style: lower case, bold) instead of being prefixed to the title.
  panel_letters <- c("a", "b", "c")
  panel_letter_size <- 15

  if (figure == "fig3") {
    titles <- c(microbiome = "Microbiome", metabolomics = "Metabolomics",
                proteomics = "Proteomics")
    data_list <- lapply(omics_order, function(nm)
      extract_fig3(model_outputs[[nm]], top_n,
                   flip_model3_xgb_shap_sign = isTRUE(flip_model3_xgb_shap_sign)))
    names(data_list) <- omics_order

    bubble_data <- lapply(omics_order, function(nm) {
      df <- bubble_fig3(data_list[[nm]], labels_omics[[nm]], label_max[[nm]])
      if (!is.null(df)) df$feat_type <- "Omics"
      df
    })
    names(bubble_data) <- omics_order
    # Legend key is built from what the panels draw, so an unused category (a
    # grey "No clear direction" or "Unknown" swatch) is never advertised.
    directions_present <- unique(as.character(unlist(
      lapply(bubble_data, function(df) if (is.null(df)) NULL else df$direction)
    )))
    directions_present <- directions_present[!is.na(directions_present)]

    panels <- lapply(omics_order, function(nm) {
      df <- bubble_data[[nm]]
      n <- if (!is.null(df)) length(unique(df$label)) else 0
      ttl <- if (n > 0) sprintf("%s  (n = %d)", titles[[nm]], n) else titles[[nm]]
      ml <- if (nm == "proteomics") 20L else 14L
      draw_panel(df, "EN", "XGB", ttl, show_omics_label = FALSE,
                 margin_left = ml,
                 margin_right = if (nm == "proteomics") 22 else NULL,
                 omics_y_axis_green = FALSE, direction_scale_style = "fig3") +
        ggplot2::theme(legend.position = "none")
    })
    n_feat <- sapply(omics_order, function(nm) {
      d <- data_list[[nm]]
      length(intersect(d$en$feature[d$en$rank_en <= top_n],
                       d$xgb$feature[d$xgb$rank_xgb <= top_n]))
    })
    row_h <- max(n_feat) * 0.40 + 1.25
    row_panels <- cowplot::plot_grid(plotlist = panels, ncol = 3, align = "h",
                                     rel_widths = panel_rel_widths,
                                     labels = panel_letters,
                                     label_fontface = "bold",
                                     label_size = panel_letter_size,
                                     label_x = 0, label_y = 1,
                                     hjust = -0.4, vjust = 1.2)
    fig <- cowplot::plot_grid(
      row_panels, make_legend_fig3_standalone(directions_present), ncol = 1,
      rel_heights = c(row_h, 1.45)
    )
    total_h <- row_h + 1.85
  } else {
    titles <- c(microbiome = "Microbiome-integrated",
                metabolomics = "Metabolome-integrated",
                proteomics = "Proteome-integrated")
    data_list <- lapply(omics_order, function(nm) extract_fig4(model_outputs[[nm]], top_n))
    names(data_list) <- omics_order

    panels <- lapply(omics_order, function(nm) {
      df <- bubble_fig4(data_list[[nm]], labels_omics[[nm]], label_max[[nm]])
      n <- if (!is.null(df)) length(unique(df$label)) else 0
      ttl <- if (n > 0) sprintf("%s  (n = %d)", titles[[nm]], n) else titles[[nm]]
      ml <- if (nm == "proteomics") 20L else 14L
      draw_panel(df, "RF+EN", "RF+XGB", ttl,
                 show_omics_label = TRUE, omics_label_txt = labels_omics[[nm]],
                 margin_left = ml,
                 margin_right = if (nm == "proteomics") 36 else NULL) +
        ggplot2::theme(legend.position = "none")
    })
    n_feat <- sapply(omics_order, function(nm) {
      d <- data_list[[nm]]
      length(intersect(d$m4$feature[d$m4$rank <= top_n],
                       d$m5$feature[d$m5$rank <= top_n]))
    })
    row_h <- max(n_feat) * 0.40 + 1.45
    row_panels <- cowplot::plot_grid(plotlist = panels, ncol = 3, align = "h",
                                     rel_widths = panel_rel_widths,
                                     labels = panel_letters,
                                     label_fontface = "bold",
                                     label_size = panel_letter_size,
                                     label_x = 0, label_y = 1,
                                     hjust = -0.4, vjust = 1.2)
    fig <- cowplot::plot_grid(
      row_panels, make_legend_fig4_standalone(), ncol = 1,
      rel_heights = c(row_h, 1.55)
    )
    total_h <- row_h + 2.1
  }

  cowplot::save_plot(output_file, fig, base_width = 17.8, base_height = total_h)
}

# =============================================================================
# Supplementary figures ported from the figures_new logic
# (Figure3_Importance.Rmd for supp1, Figure4_Ensemble.Rmd for supp2). Panels are
# rebuilt from the raw data frames in the RDS (coef_en / plots$p_shapX_bee$data)
# so output is independent of the ggplot2 version and matches the figures_new
# reference exactly. Clinical-covariate labels/classification are ported from
# figures_new/clinical_feature_labels.R.
# =============================================================================

clinical_id_labels <- c(
  age_v4 = "Age", age_v5 = "Age at visit 5 (y)", race = "Race", edu = "Education",
  site = "Site", ol_health = "Self-rated health", bmi = "BMI", mstat = "Marital status",
  smoke = "Smoking status", diab = "Diabetes", hbp = "Hypertension", cancer = "Cancer",
  tmm_score = "Teng 3MS score", gds = "GDS-15", pase = "PASE", total_meds = "Total medications",
  abx = "Oral antibiotic (2 wk)", cr_cmm = "Muscle mass (creatine %)", grip = "Grip strength (kg)",
  chair = "Chair stands (per 10 s)", gait400m = "Gait speed 400 m (m/s)", gait6m = "Gait speed 6 m (m/s)",
  b4thd = "Total hip BMD", b4fnd = "Femoral neck BMD", b4lsd = "Lumbar spine BMD",
  faprev4 = "Prior fracture before V4", hqdrfefl = "HR-pQCT radius failure load",
  hqdtfefl = "HR-pQCT tibia failure load", hqptfefl = "HR-pQCT diaph. tibia failure load",
  status = "Vital status", fr_chs4 = "Frailty (CHS)"
)

# Supplementary Figure 2 uses the manuscript's expanded clinical labels.
clinical_id_labels_supp2 <- clinical_id_labels
clinical_id_labels_supp2[c("tmm_score", "gds", "pase")] <- c(
  "Teng mMMSE score", "Depressive symptoms (GDS)", "Physical activity (PASE)"
)

clinical_covariate_ids <- function() names(clinical_id_labels)

is_preprocess_clinical_id <- function(x) {
  x <- as.character(x)
  stem <- sub("_.*$", "", x)
  ids <- clinical_covariate_ids()
  x %in% ids | stem %in% ids
}

pretty_preprocess_feature <- function(x, max_len = NA_integer_) {
  x <- as.character(x)
  stem <- sub("_.*$", "", x)
  m1 <- match(x, names(clinical_id_labels), nomatch = 0L)
  m2 <- match(stem, names(clinical_id_labels), nomatch = 0L)
  idx <- ifelse(m1 > 0L, m1, ifelse(m2 > 0L, m2, NA_integer_))
  out <- ifelse(!is.na(idx), unname(clinical_id_labels)[idx], x)
  if (!is.na(max_len) && isTRUE(max_len > 0L)) {
    long <- nchar(out) > max_len
    out[long] <- paste0(substr(out[long], 1L, max_len - 1L), "\u2026")
  }
  out
}

pretty_plot_feature_label <- function(x, max_len = 60L) {
  x <- as.character(x)
  out <- vapply(seq_along(x), function(i) {
    xi <- x[[i]]
    if (is_preprocess_clinical_id(xi)) pretty_preprocess_feature(xi, NA_integer_) else xi
  }, character(1))
  long <- nchar(out) > max_len
  out[long] <- paste0(substr(out[long], 1L, max_len - 1L), "\u2026")
  out
}

CLIN_IDS <- clinical_covariate_ids()

is_clinical_id <- function(x) {
  x <- as.character(x)
  stem <- sub("_.*$", "", x)
  x %in% CLIN_IDS | stem %in% CLIN_IDS
}

pretty_name <- function(x) pretty_plot_feature_label(x, 60L)

pretty_supp2_name <- function(x, max_len = 60L) {
  x <- as.character(x)
  stem <- sub("_.*$", "", x)
  key <- ifelse(x %in% names(clinical_id_labels_supp2), x, stem)
  out <- x
  matched <- !is.na(key) & key %in% names(clinical_id_labels_supp2)
  out[matched] <- unname(clinical_id_labels_supp2[key[matched]])
  long <- nchar(out) > max_len
  out[long] <- paste0(substr(out[long], 1L, max_len - 1L), "\u2026")
  out
}

# --- Elastic-net coefficient panels (figures_new Figure3 rows a-c) ------------
supp_col_pos <- "#E64B35"
supp_col_neg <- "#3C5488"

theme_supp_en <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.line         = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.text         = ggplot2::element_text(colour = "black", size = base_size),
      axis.text.y       = ggplot2::element_text(size = base_size - 0.5),
      axis.title        = ggplot2::element_text(colour = "black", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5,
                                                margin = ggplot2::margin(b = 3)),
      plot.title        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = base_size, colour = "black"),
      legend.key        = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.height = grid::unit(9, "pt"),
      legend.key.width  = grid::unit(12, "pt"),
      legend.background = ggplot2::element_blank(),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      plot.margin       = ggplot2::margin(4, 10, 4, 4)
    )
}

top_en <- function(obj, n = 20) {
  d <- obj$coef_en |>
    dplyr::filter(estimate != 0) |>
    dplyr::mutate(abs_est = abs(estimate)) |>
    dplyr::slice_max(abs_est, n = n, with_ties = FALSE) |>
    dplyr::arrange(estimate) |>
    dplyr::mutate(
      direction = ifelse(estimate > 0, "Positive", "Negative"),
      direction = factor(direction, levels = c("Positive", "Negative")),
      label = pretty_plot_feature_label(as.character(term), 34L)
    )
  d$label <- factor(make.unique(d$label, sep = " "), levels = make.unique(d$label, sep = " "))
  d
}

plot_en_panel <- function(df, subtitle) {
  ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = label, fill = direction)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, colour = "black") +
    ggplot2::scale_fill_manual(values = c(Positive = supp_col_pos, Negative = supp_col_neg), drop = FALSE) +
    ggplot2::labs(x = "EN coefficient", y = NULL, subtitle = subtitle, fill = NULL) +
    theme_supp_en(base_size = 7) +
    ggplot2::theme(
      legend.position        = "inside",
      legend.position.inside = c(0.98, 0.02),
      legend.justification   = c(1, 0),
      legend.direction       = "vertical",
      legend.background      = ggplot2::element_rect(fill = scales::alpha("white", 0.9),
                                                     colour = "grey70", linewidth = 0.2),
      legend.margin          = ggplot2::margin(2, 3, 2, 3)
    )
}

# --- SHAP beeswarm panels (figures_new Figure3 rows d-f / Figure4) ------------
supp_shap_gradient <- c("#3C5488", "#CFCFCF", "#E64B35")

theme_supp_bee <- function(base_size = 7) {
  ggplot2::theme_classic(base_size = base_size, base_family = "Helvetica") +
    ggplot2::theme(
      axis.line         = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks        = ggplot2::element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(2, "pt"),
      axis.text         = ggplot2::element_text(colour = "black", size = base_size),
      axis.text.y       = ggplot2::element_text(size = base_size - 0.5),
      axis.title        = ggplot2::element_text(colour = "black", size = base_size + 1),
      plot.subtitle     = ggplot2::element_text(face = "bold", size = base_size + 2, hjust = 0.5,
                                                margin = ggplot2::margin(b = 3)),
      plot.title        = ggplot2::element_blank(),
      legend.title      = ggplot2::element_blank(),
      legend.text       = ggplot2::element_text(size = base_size, colour = "black"),
      legend.key        = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key.height = grid::unit(9, "pt"),
      legend.key.width  = grid::unit(12, "pt"),
      legend.background = ggplot2::element_blank(),
      legend.margin     = ggplot2::margin(0, 0, 0, 0),
      plot.margin       = ggplot2::margin(4, 16, 4, 4)
    )
}

# Top-N SHAP beeswarm rows by mean |SHAP| for a given saved plot (p_shapX_bee).
top_shap_bee <- function(obj, plot_name, n = 20, label_style = c("default", "supp2")) {
  label_style <- match.arg(label_style)
  if (is.null(obj$plots[[plot_name]])) stop("No ", plot_name, " in RDS")
  d <- as.data.frame(obj$plots[[plot_name]]$data)
  rank_tbl <- d |>
    dplyr::group_by(feature) |>
    dplyr::summarise(mean_abs = mean(abs(value)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(mean_abs)) |>
    dplyr::slice_head(n = n)
  keep <- as.character(rank_tbl$feature)
  d <- d |> dplyr::filter(as.character(feature) %in% keep)
  label_fn <- if (label_style == "supp2") pretty_supp2_name else pretty_name
  lab_map <- stats::setNames(label_fn(keep), keep)
  lab_map <- make.unique(lab_map, sep = " ")
  names(lab_map) <- keep
  d$label <- factor(lab_map[as.character(d$feature)], levels = rev(unname(lab_map)))
  d$is_clinical <- is_clinical_id(d$feature)
  d
}

# Negate SHAP for legacy Model-3 XGB (matches figures_new flip_model3_xgb_shap_sign).
maybe_flip_shap <- function(bee_df, flip) {
  if (!isTRUE(flip)) return(bee_df)
  bee_df$value <- -as.numeric(bee_df$value)
  bee_df
}

# omics label colour = green (figures_new Figure4 convention).
supp_omics_label_colour <- "#2E7D32"

plot_shap_bee_panel <- function(df, subtitle, colour_omics = FALSE) {
  lvl <- levels(df$label)
  tag <- df |> dplyr::distinct(label, is_clinical)
  m <- match(lvl, tag$label)
  is_clin_lvl <- ifelse(is.na(m), TRUE, tag$is_clinical[m])

  axis_text_colour <- if (colour_omics) {
    ifelse(is_clin_lvl, "black", supp_omics_label_colour)
  } else {
    "black"
  }

  ggplot2::ggplot(df, ggplot2::aes(x = value, y = label, colour = color)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.3, colour = "grey55") +
    ggbeeswarm::geom_quasirandom(groupOnX = FALSE, size = 0.45, alpha = 0.75, shape = 16, width = 0.38) +
    ggplot2::scale_colour_gradientn(colours = supp_shap_gradient, limits = c(0, 1),
                                    breaks = c(0, 1), labels = c("Low", "High"), name = "Feature value") +
    ggplot2::labs(x = "SHAP value", y = NULL, subtitle = subtitle) +
    ggplot2::guides(colour = ggplot2::guide_colourbar(
      title.position = "right", title.hjust = 0.5, title.vjust = 0.5,
      barwidth = grid::unit(3.5, "pt"), barheight = grid::unit(55, "pt"), ticks = FALSE
    )) +
    theme_supp_bee(base_size = 7) +
    ggplot2::theme(
      legend.position    = "right",
      legend.title       = ggplot2::element_text(angle = 90, size = 6.5, colour = "black"),
      legend.text        = ggplot2::element_text(size = 6, colour = "black"),
      legend.margin      = ggplot2::margin(0, 0, 0, 2),
      legend.box.spacing = grid::unit(2, "pt"),
      axis.text.y = ggplot2::element_text(
        face = ifelse(is_clin_lvl, "plain", "italic"),
        size = 6.5, colour = axis_text_colour
      )
    )
}

write_supplementary_importance_pdf <- function(model_outputs, output_file, figure = c("supp1", "supp2")) {
  figure <- match.arg(figure)
  omics <- c("microbiome", "metabolomics", "proteomics")
  omics_titles <- c(microbiome = "Microbiome", metabolomics = "Metabolomics", proteomics = "Proteomics")

  if (figure == "supp1") {
    en_panels <- lapply(omics, function(nm) {
      plot_en_panel(top_en(model_outputs[[nm]], 20), paste0(omics_titles[[nm]], " - EN"))
    })
    xgb_panels <- lapply(omics, function(nm) {
      bee <- maybe_flip_shap(top_shap_bee(model_outputs[[nm]], "p_shap3_bee", 20), flip = TRUE)
      # Figure3 convention: XGB omics labels stay black (italic), not green.
      plot_shap_bee_panel(bee, paste0(omics_titles[[nm]], " \u2013 XGBoost"), colour_omics = FALSE)
    })
    plots <- c(en_panels, xgb_panels)
  } else {
    rfen_panels <- lapply(omics, function(nm) {
      bee <- top_shap_bee(model_outputs[[nm]], "p_shap4_bee", 20, label_style = "supp2")
      plot_shap_bee_panel(bee, paste0(omics_titles[[nm]], " \u2013 RF + EN"), colour_omics = TRUE)
    })
    rfxgb_panels <- lapply(omics, function(nm) {
      bee <- top_shap_bee(model_outputs[[nm]], "p_shap5_bee", 20, label_style = "supp2")
      plot_shap_bee_panel(bee, paste0(omics_titles[[nm]], " \u2013 RF + XGBoost"), colour_omics = TRUE)
    })
    plots <- c(rfen_panels, rfxgb_panels)
  }

  fig <- ggpubr::ggarrange(
    plotlist = plots,
    labels = c("a", "b", "c", "d", "e", "f"),
    ncol = 3,
    nrow = 2,
    font.label = list(face = "bold", size = 10, family = "Helvetica"),
    hjust = -0.2,
    vjust = 1.4
  )
  ggplot2::ggsave(output_file, fig, width = 120 * 3, height = 110 * 2, units = "mm")
}

# =============================================================================
# Model comparison tables (Supplementary Tables 9 and 10)
#
# Two questions are answered here, both from the *saved* model outputs so that
# no model is refitted:
#
#   (1) Within an omics layer, does adding omics information to the clinical
#       model improve discrimination? Because every model in a layer is scored
#       on the same held-out test subjects, this is a paired comparison and is
#       tested with the paired DeLong test.
#
#   (2) Across omics layers, which layer carries the most mortality
#       information? Test-set AUCs are NOT comparable across layers because the
#       layers were sampled from different participants (the identical clinical
#       model scores 0.775 / 0.709 / 0.647 in the three subsets). The matched
#       comparison below scores all three layers on the same complete-case
#       participants instead.
#
# Stored probabilities `p1`..`p5` are P(Active), i.e. survival; risk of death is
# therefore 1 - p. `.pred_Deceased` already carries that orientation.
# =============================================================================

# pROC ROC object with a fixed, explicit orientation (higher risk -> Deceased),
# so a mis-signed predictor shows up as AUC < 0.5 rather than being silently
# flipped by pROC's direction = "auto".
roc_death <- function(truth, risk) {
  keep <- !is.na(truth) & !is.na(risk)
  pROC::roc(
    response  = factor(as.character(truth)[keep], levels = c("Active", "Deceased")),
    predictor = as.numeric(risk)[keep],
    levels    = c("Active", "Deceased"),
    direction = "<",
    quiet     = TRUE
  )
}

# Every analytic participant of one omics layer with an out-of-sample predicted
# risk: out-of-fold predictions for the training portion (mapped back to IDs via
# `.train_row`) and held-out predictions for the test portion. Returns one tibble
# per model type. Ensemble models are excluded because no out-of-fold ensemble
# predictions are stored.
out_of_sample_predictions <- function(model_output) {
  tr <- rsample::training(model_output$split)
  tr_id <- as.character(tr$ID)
  tr_status <- as.character(tr$status)
  preds <- model_output$predictions
  spec <- list(
    clinical = c("pred_oof_cov", "prob_test_cov", "p1"),
    en       = c("pred_oof_micro_en", "prob_test_micro_en", "p2"),
    xgb      = c("pred_oof_micro_xgb", "prob_test_micro", "p3")
  )
  lapply(spec, function(k) {
    oof <- preds[[k[[1]]]]
    tst <- preds[[k[[2]]]]
    if (is.null(oof) || is.null(tst)) return(NULL)
    # Resolve the stored score convention rather than assuming it.
    #
    # The p1..p3 column is survival-oriented in the primary outputs
    # (p == .pred_Active) and death-oriented in the re-fitted ones
    # (p == .pred_Deceased). The out-of-fold frames carry only the score, so the
    # convention is read off the paired test frame, which carries both, and the
    # same transform is applied. Hard-coding `1 - p` silently inverts every
    # out-of-fold risk on the other generation -- it turns a clinical AUC of
    # 0.669 into 0.331 -- and nothing downstream would flag it.
    score_col <- k[[3]]
    if (!all(c(score_col, ".pred_Deceased") %in% names(tst))) {
      stop("out_of_sample_predictions: cannot resolve the score convention for ", score_col)
    }
    p_test <- as.numeric(tst[[score_col]])
    p_death <- as.numeric(tst$.pred_Deceased)
    is_death <- isTRUE(all.equal(p_test, p_death, tolerance = 1e-12, check.attributes = FALSE))
    is_alive <- isTRUE(all.equal(p_test, 1 - p_death, tolerance = 1e-12, check.attributes = FALSE))
    if (!is_death && !is_alive) {
      stop("out_of_sample_predictions: stored ", score_col,
           " is neither the predicted death risk nor its complement")
    }
    oof_risk <- as.numeric(oof[[score_col]])
    if (!is_death) oof_risk <- 1 - oof_risk
    dplyr::bind_rows(
      tibble::tibble(
        ID    = tr_id[oof$.train_row],
        truth = tr_status[oof$.train_row],
        risk  = oof_risk,
        source = "out_of_fold"
      ),
      tibble::tibble(
        ID    = as.character(tst$ID),
        truth = as.character(tst$truth),
        risk  = tst$.pred_Deceased,
        source = "held_out_test"
      )
    )
  })
}

# Percentile bootstrap CI for the difference in paired AUCs. pROC's roc.test
# gives the P value but not an interval for the difference itself.
#
# Note that the two therefore come from different procedures: every delta-AUC
# row in the supplementary tables pairs a *bootstrap percentile* interval with a
# *DeLong* P value. They agree in all but one published row -- metabolomics,
# EN omics vs LR age-only, where the interval is [+0.0005, +0.2959] and
# P = 0.054 -- because a percentile interval and a DeLong test are not
# constrained to draw the same boundary on a marginal comparison. Both are
# valid; the caption should say which is which. If they must agree by
# construction, derive the interval from the DeLong statistic instead
# (SE = delta / Z from roc.test, interval = delta +/- 1.96 SE), but note that
# changes every published interval.
delta_auc_ci <- function(truth, risk_a, risk_b, n_boot = 2000L, conf = 0.95,
                         seed = NULL) {
  if (length(n_boot) != 1L || !is.finite(n_boot) || n_boot < 1 ||
      n_boot != as.integer(n_boot)) {
    stop("n_boot must be a positive integer")
  }
  if (length(conf) != 1L || !is.finite(conf) || conf <= 0 || conf >= 1) {
    stop("conf must be between 0 and 1")
  }
  if (!is.null(seed)) {
    if (length(seed) != 1L || !is.finite(seed)) stop("seed must be finite")
    set.seed(as.integer(seed))
  }
  if (length(truth) != length(risk_a) || length(truth) != length(risk_b) ||
      length(truth) < 2L) {
    stop("truth and both risk vectors must have the same length of at least 2")
  }
  # Both callers pass death risk: build_delta_auc_table() reads .pred_Deceased
  # directly, and the matched builders store risk that way. The saved primary
  # objects also carry a survival-oriented p1..p5 column, which is never fed in
  # here without being flipped first -- see out_of_sample_predictions().
  if (anyNA(truth) || anyNA(risk_a) || anyNA(risk_b) ||
      any(!is.finite(risk_a)) || any(!is.finite(risk_b)) ||
      any(risk_a < 0 | risk_a > 1) || any(risk_b < 0 | risk_b > 1) ||
      length(unique(as.character(truth))) != 2L) {
    stop("truth and risks must be complete, binary, and finite probabilities")
  }
  n <- length(truth)
  d <- vapply(seq_len(n_boot), function(i) {
    idx <- sample.int(n, n, replace = TRUE)
    if (length(unique(truth[idx])) < 2L) return(NA_real_)
    as.numeric(pROC::auc(roc_death(truth[idx], risk_a[idx]))) -
      as.numeric(pROC::auc(roc_death(truth[idx], risk_b[idx])))
  }, numeric(1))
  stats::quantile(d, c((1 - conf) / 2, 1 - (1 - conf) / 2), na.rm = TRUE)
}

# Supplementary Table 9: within-layer paired comparisons of each model against
# the clinical benchmark and against the age-only floor, on the held-out test
# set. BH adjustment is applied within each comparison family (each reference
# model), across the 12 model x omics combinations.
# Each reference model defines a comparison family that answers one question,
# and P values are adjusted within the family:
#
#   vs the clinical model  - do omics markers alone match the clinical panel
#                            (EN/XGB omics), and does adding omics improve on it
#                            (RF+EN, RF+XGB)?
#   vs the age-only model  - does the clinical panel beat age alone, and do
#                            omics markers alone beat age alone?
#
# Stacked models are deliberately NOT compared with the age-only model: such a
# comparison bundles "does the rest of the clinical panel beat age?" with "does
# omics add to clinical?" and cannot separate them. Both questions are answered
# individually by the rows above.
build_delta_auc_table <- function(model_outputs, n_boot = 2000L, seed = 123L) {
  set.seed(seed)
  model_lab <- c(
    prob_test_micro_en = "EN omics",
    prob_test_micro    = "XGB omics",
    prob_test_4        = "RF clinical + EN",
    prob_test_5        = "RF clinical + XGB"
  )
  # which models are compared with which reference
  family_models <- list(
    `RF clinical` = c("EN omics", "XGB omics", "RF clinical + EN", "RF clinical + XGB"),
    `LR age-only` = c("RF clinical", "EN omics", "XGB omics")
  )
  rows <- list()
  for (om in names(model_outputs)) {
    preds <- model_outputs[[om]]$predictions
    base <- preds$prob_test_cov
    age <- fit_age_baseline(model_outputs[[om]])
    truth <- as.character(base$truth)
    r_clin <- roc_death(truth, base$.pred_Deceased)
    r_age <- roc_death(truth, age$predictions$prob_test_age$.pred_Deceased)
    # the clinical model itself is a comparand in the age-only family
    candidates <- c(as.list(model_lab), list(`RF clinical` = "RF clinical"))
    for (k in names(candidates)) {
      d <- if (k == "RF clinical") base else preds[[k]]
      if (is.null(d)) next
      this_model <- if (k == "RF clinical") "RF clinical" else unname(model_lab[[k]])
      r <- roc_death(truth, d$.pred_Deceased)
      for (ref in c("RF clinical", "LR age-only")) {
        if (!(this_model %in% family_models[[ref]])) next
        r_ref <- if (ref == "RF clinical") r_clin else r_age
        risk_ref <- if (ref == "RF clinical") {
          base$.pred_Deceased
        } else {
          age$predictions$prob_test_age$.pred_Deceased
        }
        ci <- delta_auc_ci(truth, d$.pred_Deceased, risk_ref, n_boot = n_boot)
        rows[[length(rows) + 1L]] <- tibble::tibble(
          omics = om, model = this_model, reference = ref,
          n_test = nrow(d), n_deaths = sum(truth == "Deceased"),
          auc = as.numeric(pROC::auc(r)),
          auc_ref = as.numeric(pROC::auc(r_ref)),
          delta = as.numeric(pROC::auc(r)) - as.numeric(pROC::auc(r_ref)),
          delta_lo = ci[[1]], delta_hi = ci[[2]],
          p_value = pROC::roc.test(r, r_ref, method = "delong", paired = TRUE)$p.value
        )
      }
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::group_by(reference) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()
}

# Generic booktabs-style three-line table renderer, used for the supplementary
# comparison tables. Column 1 is left-aligned, the rest centred.
draw_three_line_table <- function(df, title, footnote = NULL, fontsize = 8) {
  n_row <- nrow(df)
  n_data <- ncol(df) - 1L
  grid::grid.newpage()
  left <- 0.04
  right <- 0.98
  data_x <- seq(0.34, 0.94, length.out = n_data)
  row_h <- min(0.045, 0.62 / max(n_row + 2L, 1L))
  y_top <- 0.88
  y_head <- y_top - row_h
  y_mid <- y_head - 0.6 * row_h

  grid::grid.text(title, x = left, y = 0.95, just = "left",
                  gp = grid::gpar(fontsize = fontsize + 3, fontface = "bold"))
  grid::grid.lines(x = c(left, right), y = c(y_top, y_top), gp = grid::gpar(lwd = 1.4))
  grid::grid.text(names(df)[1], x = left, y = y_head, just = "left",
                  gp = grid::gpar(fontsize = fontsize, fontface = "bold"))
  for (j in seq_len(n_data)) {
    grid::grid.text(names(df)[j + 1L], x = data_x[j], y = y_head, just = "centre",
                    gp = grid::gpar(fontsize = fontsize, fontface = "bold"))
  }
  grid::grid.lines(x = c(left, right), y = c(y_mid, y_mid), gp = grid::gpar(lwd = 0.8))
  for (i in seq_len(n_row)) {
    yy <- y_mid - i * row_h
    grid::grid.text(as.character(df[[1]][i]), x = left, y = yy, just = "left",
                    gp = grid::gpar(fontsize = fontsize))
    for (j in seq_len(n_data)) {
      grid::grid.text(as.character(df[[j + 1L]][i]), x = data_x[j], y = yy, just = "centre",
                      gp = grid::gpar(fontsize = fontsize))
    }
  }
  y_bottom <- y_mid - (n_row + 0.5) * row_h
  grid::grid.lines(x = c(left, right), y = c(y_bottom, y_bottom), gp = grid::gpar(lwd = 1.4))
  if (!is.null(footnote)) {
    for (i in seq_along(footnote)) {
      grid::grid.text(footnote[[i]], x = left, y = y_bottom - 0.035 * i, just = "left",
                      gp = grid::gpar(fontsize = fontsize - 1.5))
    }
  }
  invisible(TRUE)
}

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))

write_delta_auc_pdf <- function(tab, output_pdf) {
  df <- tab |>
    dplyr::arrange(reference, omics, model) |>
    dplyr::transmute(
      Reference = reference,
      Omics = tools::toTitleCase(omics),
      Model = model,
      `Test n (deaths)` = sprintf("%d (%d)", n_test, n_deaths),
      AUC = sprintf("%.3f", auc),
      `AUC difference (95% CI)` = sprintf("%+.3f (%+.3f to %+.3f)", delta, delta_lo, delta_hi),
      P = fmt_p(p_value),
      `P (BH)` = fmt_p(p_adj)
    )
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf, width = 11, height = 8)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_three_line_table(
    df,
    title = "Supplementary Table 9. Paired within-layer model comparisons (held-out test set)",
    footnote = c(
      "Each model is compared with the reference model on the same held-out test participants using the paired DeLong test.",
      "AUC difference 95% CI is a 2,000-replicate percentile bootstrap. P (BH) is Benjamini-Hochberg adjusted within each reference family (12 comparisons).",
      "Computed from the stored test-set predictions; no model was refitted."
    )
  )
  TRUE
}

# Calibration of the held-out test-set predictions.
#
# Discrimination (AUC) says whether the model ranks participants correctly;
# calibration says whether the predicted probabilities are numerically right.
# Both are needed before a model can be read as a risk calculator, and TRIPOD
# asks for both. Two standard summaries are reported:
#
#   calibration-in-the-large (intercept): fit y ~ offset(logit(p)); 0 means the
#     average predicted risk matches the observed event rate.
#   calibration slope: fit y ~ logit(p); 1 is ideal. Below 1 means predictions
#     are too extreme (the usual sign of overfitting); above 1 means they are
#     too tightly clustered around the base rate.
#
# The Brier score is the mean squared difference between predicted probability
# and outcome; it is compared with the Brier score of a model that predicts the
# observed event rate for everyone, which is the relevant "no information" floor.
build_calibration_table <- function(model_outputs) {
  model_lab <- c(
    prob_test_cov      = "RF clinical",
    prob_test_micro_en = "EN omics",
    prob_test_micro    = "XGB omics",
    prob_test_4        = "RF clinical + EN",
    prob_test_5        = "RF clinical + XGB"
  )
  rows <- list()
  for (om in names(model_outputs)) {
    preds <- model_outputs[[om]]$predictions
    for (k in names(model_lab)) {
      d <- preds[[k]]
      if (is.null(d)) next
      y <- as.integer(as.character(d$truth) == "Deceased")
      p <- pmin(pmax(d$.pred_Deceased, 1e-6), 1 - 1e-6)
      lp <- stats::qlogis(p)
      slope <- tryCatch(
        unname(stats::coef(stats::glm(y ~ lp, family = stats::binomial()))[2]),
        error = function(e) NA_real_
      )
      inter <- tryCatch(
        unname(stats::coef(stats::glm(y ~ stats::offset(lp), family = stats::binomial()))[1]),
        error = function(e) NA_real_
      )
      rows[[length(rows) + 1L]] <- tibble::tibble(
        omics = om, model = unname(model_lab[[k]]), n = length(y), events = sum(y),
        mean_predicted = mean(p), observed = mean(y),
        brier = mean((p - y)^2), brier_null = mean(y) * (1 - mean(y)),
        calibration_intercept = inter, calibration_slope = slope
      )
    }
  }
  dplyr::bind_rows(rows)
}

write_calibration_pdf <- function(tab, output_pdf) {
  df <- tab |>
    dplyr::transmute(
      Omics = tools::toTitleCase(omics), Model = model,
      `n (events)` = sprintf("%d (%d)", n, events),
      `Predicted / observed` = sprintf("%.2f / %.2f", mean_predicted, observed),
      `Brier (null)` = sprintf("%.3f (%.3f)", brier, brier_null),
      Intercept = sprintf("%+.2f", calibration_intercept),
      Slope = sprintf("%.2f", calibration_slope)
    )
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf, width = 10, height = 7)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_three_line_table(
    df,
    title = "Supplementary Table 13. Calibration of the held-out test-set predictions",
    footnote = c(
      "Calibration-in-the-large (intercept) from y ~ offset(logit(p)); 0 is ideal. Calibration slope from y ~ logit(p); 1 is ideal.",
      "Brier score is the mean squared error of the predicted probability; the null value is that of predicting the observed event rate for everyone.",
      "Test sets contain 90-154 participants, so the slopes in particular are estimated imprecisely and should be read as indicative."
    )
  )
  TRUE
}

# Table 1 as a Word file.
#
# The submission copy of Table 1 was previously maintained by hand, which meant
# the manuscript folder held a table no script could regenerate. This writes the
# same data frame that feeds the PDF version into a .docx, so the submitted table
# and the repository stay in step.
write_table1_docx <- function(table1, output_docx) {
  if (!requireNamespace("officer", quietly = TRUE) ||
      !requireNamespace("flextable", quietly = TRUE)) {
    return(FALSE)
  }
  df <- as.data.frame(table1, stringsAsFactors = FALSE)
  names(df) <- c("Variable", "Overall", "Deceased", "Alive", "P-value")
  n <- attr(table1, "group_sizes")
  if (is.null(n)) n <- c(overall = 879L, deceased = 485L, alive = 394L)
  ft <- flextable::flextable(df)
  ft <- flextable::set_header_labels(
    ft,
    Overall  = sprintf("Overall\nN = %d", n[["overall"]]),
    Deceased = sprintf("Deceased\nN = %d", n[["deceased"]]),
    Alive    = sprintf("Alive\nN = %d", n[["alive"]])
  )
  ft <- flextable::add_footer_lines(ft, c(
    "Values are n (%) for categorical variables and mean (SD) for continuous variables.",
    "P-values from Pearson's chi-squared test, Fisher's exact test, or the Wilcoxon rank sum test."
  ))
  ft <- flextable::font(ft, fontname = "Times New Roman", part = "all")
  ft <- flextable::fontsize(ft, size = 9, part = "all")
  ft <- flextable::fontsize(ft, size = 8, part = "footer")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, j = 2:5, align = "center", part = "all")
  ft <- flextable::autofit(ft)
  # Three-line (booktabs) rules, matching the other manuscript tables.
  rule <- officer::fp_border(color = "black", width = 1.2)
  thin <- officer::fp_border(color = "black", width = 0.6)
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline_top(ft, border = rule, part = "header")
  ft <- flextable::hline_bottom(ft, border = thin, part = "header")
  ft <- flextable::hline_bottom(ft, border = rule, part = "body")

  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, "Table 1. Participant Characteristics by Mortality Status",
                               style = "Normal")
  doc <- flextable::body_add_flextable(doc, ft)
  dir.create(dirname(output_docx), recursive = TRUE, showWarnings = FALSE)
  print(doc, target = output_docx)
  TRUE
}

# Table 1 computed from the analytic sample rather than transcribed.
#
# `manuscript_table1()` holds the published values as a literal data frame, which
# meant Table 1 could not be regenerated from data. This computes the same table
# from the union of the three analytic samples, and is used in preference when
# the model-output objects are available. `compare_table1()` checks the computed
# table against the stored one so any drift is visible rather than silent.
# The phenotype columns Table 1 and its per-layer counterparts describe. Held in
# one place so the union table and the layer tables cannot drift apart.
table1_columns <- c(
  "ID", "status", "age_v4", "bmi", "tmm_score", "gds", "pase", "total_meds",
  "race", "edu", "ol_health", "mstat", "smoke", "diab", "hbp", "cancer"
)

table1_analytic_sample <- function(model_outputs) {
  d <- lapply(model_outputs, function(o) o$split$data)
  d <- dplyr::bind_rows(lapply(d, function(x) {
    x[, intersect(names(x), table1_columns), drop = FALSE]
  }))
  dplyr::distinct(d, ID, .keep_all = TRUE)
}

# Participant characteristics inside each layer's own analytic sample.
#
# The same builder as Table 1, applied per layer instead of to the union of the
# three. Table 1 describes the 879 men with at least one layer; the supplement
# reports one of these per layer, because the layers do not share a sample.
#
# Returns the rendered table plus the counts, so the caller can state n / deaths
# / alive in the recorded note. The table body carries no N row: by the
# convention used throughout results/tables, the body is the table and the
# caption in supp.tex supplies the header counts.
build_layer_characteristics_tables <- function(model_outputs) {
  lapply(model_outputs, function(o) {
    d <- o$split$data
    d <- d[, intersect(names(d), table1_columns), drop = FALSE]
    d <- dplyr::distinct(d, ID, .keep_all = TRUE)
    status <- as.character(d$status)
    list(
      table = table1_core(d),
      n = nrow(d),
      deaths = sum(status == "Deceased"),
      alive = sum(status == "Active")
    )
  })
}

compute_table1 <- function(model_outputs) {
  table1_core(table1_analytic_sample(model_outputs))
}

# Shared builder. `group` is a two-level character vector; the first level fills
# the "deceased" column and the second the "active" column, so the same code
# produces Table 1 (by vital status) and the complete-case membership table.
#
# `include_status` prepends a vital-status block. It is opt-in because vital
# status is itself the grouping variable for Table 1 and the per-layer tables,
# where the row would read 100% / 0% and say nothing. It is a real comparison
# only when the grouping is something else -- complete-case membership -- and
# there it carries the mortality contrast the manuscript cites.
table1_core <- function(d, group = as.character(d$status),
                        levels_ab = c("Deceased", "Active"),
                        include_status = FALSE) {
  grp <- as.character(group)
  dec <- grp == levels_ab[[1]]
  alv <- grp == levels_ab[[2]]
  # Chi-squared without continuity correction and the Wilcoxon rank sum test,
  # which is what the published Table 1 used; displayed to two significant
  # figures, or three decimals below 0.1.
  fmt_p <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("<0.001")
    if (p < 0.1) return(sprintf("%.3f", p))
    sprintf("%.2f", p)
  }
  rows <- list()
  add <- function(variable, overall, deceased, active, p_value = "") {
    rows[[length(rows) + 1L]] <<- data.frame(
      variable = variable, overall = overall, deceased = deceased,
      active = active, p_value = p_value, stringsAsFactors = FALSE
    )
  }
  pct <- function(k, n) sprintf("%d (%.1f%%)", k, 100 * k / n)

  categorical <- function(var, label) {
    x <- as.character(d[[var]])
    tab <- table(x, grp)
    p <- tryCatch(stats::chisq.test(tab, correct = FALSE)$p.value, warning = function(w)
      tryCatch(stats::fisher.test(tab, simulate.p.value = TRUE)$p.value,
               error = function(e) NA_real_), error = function(e) NA_real_)
    add(label, "", "", "", fmt_p(p))
    for (lv in sort(unique(x[!is.na(x)]))) {
      add(paste0("  ", lv), pct(sum(x == lv, na.rm = TRUE), length(x)),
          pct(sum(x == lv & dec, na.rm = TRUE), sum(dec)),
          pct(sum(x == lv & alv, na.rm = TRUE), sum(alv)))
    }
  }
  binary_yes <- function(var, label) {
    x <- as.character(d[[var]])
    p <- tryCatch(stats::chisq.test(table(x, grp), correct = FALSE)$p.value,
                  warning = function(w) NA_real_, error = function(e) NA_real_)
    add(label, pct(sum(x == "Yes", na.rm = TRUE), length(x)),
        pct(sum(x == "Yes" & dec, na.rm = TRUE), sum(dec)),
        pct(sum(x == "Yes" & alv, na.rm = TRUE), sum(alv)), fmt_p(p))
  }
  continuous <- function(var, label) {
    x <- suppressWarnings(as.numeric(d[[var]]))
    ms <- function(i) sprintf("%.2f (%.2f)", mean(x[i], na.rm = TRUE), stats::sd(x[i], na.rm = TRUE))
    p <- tryCatch(stats::wilcox.test(x[dec], x[alv])$p.value,
                  warning = function(w) NA_real_, error = function(e) NA_real_)
    add(label, ms(rep(TRUE, length(x))), ms(dec), ms(alv), fmt_p(p))
    if (any(is.na(x))) {
      add("  Missing", as.character(sum(is.na(x))),
          as.character(sum(is.na(x) & dec)), as.character(sum(is.na(x) & alv)))
    }
  }

  # Drawn first: it is the outcome, not a covariate.
  if (isTRUE(include_status)) {
    st <- as.character(d$status)
    p_st <- tryCatch(stats::chisq.test(table(st, grp), correct = FALSE)$p.value,
                     warning = function(w) NA_real_, error = function(e) NA_real_)
    add("Vital Status", "", "", "", fmt_p(p_st))
    for (lv in c("Deceased", "Active")) {
      add(if (identical(lv, "Active")) "  Alive" else "  Deceased",
          pct(sum(st == lv, na.rm = TRUE), length(st)),
          pct(sum(st == lv & dec, na.rm = TRUE), sum(dec)),
          pct(sum(st == lv & alv, na.rm = TRUE), sum(alv)))
    }
  }

  categorical("race", "Race")
  categorical("edu", "Education")
  categorical("ol_health", "Self-Rated Overall Health")
  categorical("mstat", "Marital Status")
  categorical("smoke", "Smoking Status")
  binary_yes("diab", "Diabetes")
  binary_yes("hbp", "High Blood Pressure")
  binary_yes("cancer", "Cancer")
  continuous("bmi", "BMI")
  continuous("tmm_score", "Teng 3MS Score")
  continuous("gds", "Geriatric Depression Scale")
  continuous("pase", "PASE Score")
  continuous("total_meds", "Total Medications")
  continuous("age_v4", "Age")

  out <- do.call(rbind, rows)
  attr(out, "group_sizes") <- c(overall = nrow(d), deceased = sum(dec), alive = sum(alv))
  out
}

# Non-fatal check that the computed table still reproduces the published values.
compare_table1 <- function(computed, stored = manuscript_table1()) {
  key <- function(x) paste(trimws(x$variable), x$overall, x$deceased, x$active)
  setdiff_rows <- setdiff(key(stored), key(computed))
  if (length(setdiff_rows)) {
    record_issue("03_make_manuscript_tables", "table1 drift",
                 paste("Rows in the stored Table 1 not reproduced by computation:",
                       paste(utils::head(setdiff_rows, 5), collapse = " | ")))
  }
  length(setdiff_rows) == 0L
}

# Is the omics signal simply an age proxy?
#
# The plasma proteome is a well-known correlate of chronological age, so a
# reasonable objection to a proteomic mortality signal is that it re-expresses
# age. This assembles the evidence needed to answer that in one place: how
# strongly each omics risk score tracks age, how age alone performs, and how
# much the omics score adds over age.
build_age_proxy_table <- function(model_outputs) {
  model_lab <- c(prob_test_micro_en = "EN omics", prob_test_micro = "XGB omics")
  rows <- list()
  for (om in names(model_outputs)) {
    o <- model_outputs[[om]]
    preds <- o$predictions
    dat <- o$split$data
    age_ref <- fit_age_baseline(o)
    ids <- as.character(preds$prob_test_cov$ID)
    age <- suppressWarnings(as.numeric(dat$age_v4[match(ids, as.character(dat$ID))]))
    truth <- as.character(preds$prob_test_cov$truth)
    r_age <- roc_death(truth, age_ref$predictions$prob_test_age$.pred_Deceased)
    for (k in names(model_lab)) {
      d <- preds[[k]]
      if (is.null(d)) next
      risk <- d$.pred_Deceased
      r_om <- roc_death(truth, risk)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        omics = om, model = unname(model_lab[[k]]),
        rho_age = stats::cor(risk, age, method = "spearman", use = "complete.obs"),
        auc_omics = as.numeric(pROC::auc(r_om)),
        auc_age_only = as.numeric(pROC::auc(r_age)),
        delta = as.numeric(pROC::auc(r_om)) - as.numeric(pROC::auc(r_age)),
        p_value = pROC::roc.test(r_om, r_age, method = "delong", paired = TRUE)$p.value
      )
    }
  }
  # The clinical model is shown for reference: it contains age explicitly, so
  # its correlation with age is the natural yardstick for the omics scores.
  ref <- lapply(names(model_outputs), function(om) {
    o <- model_outputs[[om]]
    d <- o$predictions$prob_test_cov
    dat <- o$split$data
    age <- suppressWarnings(as.numeric(dat$age_v4[match(as.character(d$ID), as.character(dat$ID))]))
    tibble::tibble(omics = om, model = "RF clinical (reference)",
                   rho_age = stats::cor(d$.pred_Deceased, age, method = "spearman", use = "complete.obs"),
                   auc_omics = NA_real_, auc_age_only = NA_real_,
                   delta = NA_real_, p_value = NA_real_)
  })
  dplyr::bind_rows(dplyr::bind_rows(rows), dplyr::bind_rows(ref)) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"))
}

write_age_proxy_pdf <- function(tab, output_pdf) {
  df <- tab |>
    dplyr::arrange(omics, model) |>
    dplyr::transmute(
      Omics = tools::toTitleCase(omics), Model = model,
      `Spearman rho with age` = sprintf("%+.2f", rho_age),
      AUC = ifelse(is.na(auc_omics), "-", sprintf("%.3f", auc_omics)),
      `AUC age only` = ifelse(is.na(auc_age_only), "-", sprintf("%.3f", auc_age_only)),
      `AUC difference` = ifelse(is.na(delta), "-", sprintf("%+.3f", delta)),
      `P (BH)` = ifelse(is.na(p_adj), "-", fmt_p(p_adj))
    )
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  grDevices::pdf(output_pdf, width = 10.5, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_three_line_table(
    df,
    title = "Supplementary Table 15. Relationship between omics risk scores and chronological age",
    footnote = c(
      "Spearman correlation between the predicted risk and chronological age on the held-out test set.",
      "The clinical model contains age as a covariate and is shown as the reference correlation.",
      "AUC difference compares each omics model with the logistic age-only model on the same test participants (paired DeLong)."
    )
  )
  TRUE
}

# Supplementary Table 16: the complete-case subset.
#
# The cross-layer comparison is made in the participants who have all three
# omics layers, and those participants are not a random sample of the cohort:
# assay availability was driven by sub-study membership. This describes that
# subset by vital status (panel a) and quantifies how it differs from the rest
# of the analytic sample (panel b), so the reader can judge how far the matched
# result generalises.
# Rank overlap and directional concordance between two feature rankings.
#
# Overlap: features appearing in both models' top N. Concordance: of those, how
# many carry the same adverse-versus-protective classification. N is 20, 50 and
# 100, as reported in the supplement.
#
# The two published comparisons differ only in where the rankings come from, so
# one builder serves both, taking a function that yields the pair. Each
# extractor returns complete ranked lists rather than truncated ones, so a
# single call per layer covers every N.
#
# The elastic net contributes only its non-zero coefficients, so its list can be
# shorter than N -- 50, 68 and 20 features for microbiome, metabolomics and
# proteomics. The overlap is then bounded by that length rather than by N, while
# the percentage column stays a percentage of N. That is why the proteomic rows
# at N = 50 and N = 100 describe 20 features.
build_rank_overlap_table <- function(model_outputs, ranks, ns = c(20L, 50L, 100L)) {
  rows <- list()
  for (om in names(model_outputs)) {
    pair <- ranks(model_outputs[[om]])
    for (n in ns) {
      a <- pair$a[pair$a$rank <= n, , drop = FALSE]
      b <- pair$b[pair$b$rank <= n, , drop = FALSE]
      inter <- intersect(a$feature, b$feature)
      da <- a$direction[match(inter, a$feature)]
      db <- b$direction[match(inter, b$feature)]
      concordant <- sum(!is.na(da) & !is.na(db) & da == db)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        omics = om,
        top_n = as.integer(n),
        overlap = length(inter),
        overlap_pct = 100 * length(inter) / n,
        concordant = concordant,
        concordant_pct = if (length(inter)) 100 * concordant / length(inter) else NA_real_
      )
    }
  }
  dplyr::bind_rows(rows)
}

# Omics-only elastic net versus XGBoost: elastic net ranked by absolute
# coefficient, XGBoost by mean absolute SHAP, directions as in Figure 3.
#
# `flip_model3_xgb_shap_sign = TRUE` is not optional here: it is the convention
# write_bubble_figure() uses to draw Figure 3, and this table is the numeric
# statement of that same figure. Taking extract_fig3()'s FALSE default instead
# silently inverted every XGB direction relative to the figure, which drove the
# concordance column to ~0% while the overlap counts stayed right -- the tell
# that only the sign, not the ranking, was wrong. With the flip the nine cells
# reproduce the published values exactly (100/76.5/56.5, 100/100/100,
# 100/100/86.4). If the figure's convention is ever revisited, both call sites
# must move together.
ranks_en_vs_xgb <- function(obj) {
  d <- extract_fig3(obj, flip_model3_xgb_shap_sign = TRUE)
  list(
    a = dplyr::transmute(d$en, feature, rank = rank_en, direction = direction_en),
    b = dplyr::transmute(d$xgb, feature, rank = rank_xgb, direction = direction_xgb)
  )
}

# The two stacked ensembles. Both columns are SHAP-based, so the ensembles are
# scored on the same footing and disagreement reflects genuine instability of a
# feature's direction rather than a difference in how direction was defined.
ranks_stacked_ensembles <- function(obj) {
  d <- extract_fig4(obj)
  list(
    a = dplyr::transmute(d$m4, feature, rank, direction),
    b = dplyr::transmute(d$m5, feature, rank, direction)
  )
}

# The 332 men with all three omics layers, compared with the remaining 547.
#
# Only the membership contrast is built. A second panel describing the 332 by
# vital status was dropped: the manuscript never cited it, and the characteristics
# of that subset by vital status are already covered by the three per-layer
# tables. The mortality contrast the manuscript does cite -- 50.9% versus 57.8%
# -- is now a row of this table rather than a footnote, via `include_status`.
build_complete_case_tables <- function(model_outputs) {
  d <- table1_analytic_sample(model_outputs)
  ids <- Reduce(intersect, lapply(model_outputs, function(o) as.character(o$split$data$ID)))
  cc <- as.character(d$ID) %in% ids
  membership <- table1_core(
    d,
    group = ifelse(cc, "CompleteCase", "Other"),
    levels_ab = c("CompleteCase", "Other"),
    include_status = TRUE
  )
  list(
    membership = membership,
    n_cc = sum(cc), n_other = sum(!cc),
    deaths_cc = sum(cc & d$status == "Deceased"),
    deaths_other = sum(!cc & d$status == "Deceased")
  )
}

# =============================================================================
# Supplementary tables for supp.tex
#
# These are embedded with \includegraphics and captioned in LaTeX, so the PDFs
# carry the table body only: no title, no footnotes, and a page sized to the
# content rather than a fixed sheet with whitespace around it. One panel per
# file, so each has its own \caption.
# =============================================================================

# Two column layouts, one renderer. Everything below the layout block -- rules,
# alignment, row pitch -- is shared; only the column centres and the font size
# are computed differently, and the two rules do not produce the same PDF:
#
#   "fixed"    data columns evenly spaced over a reserved band, label column
#              given the rest. What the primary supplementary tables were
#              rendered with.
#   "autofit"  column widths measured from the rendered strings, with the font
#              shrunk if the table cannot fit. What the matched-cohort tables
#              were rendered with; it handles their long model names.
#
# Both are kept because the published PDFs must stay byte-identical to what
# produced them. Pick per table; do not "upgrade" existing callers.
# A column name containing "\n" is drawn as a stacked header, and its width is
# measured from the widest line rather than the whole string. That is what keeps
# a header like "Concordant directions\n(% of overlap)" from forcing the autofit
# layout to shrink the font until the table is unreadable. With no "\n" anywhere
# the geometry below reduces exactly to the single-line case, so tables that were
# already rendered are unaffected.
draw_bare_table <- function(df, fontsize = 8, first_col_x = 0.02,
                            layout = c("fixed", "autofit")) {
  layout <- match.arg(layout)
  n_row <- nrow(df)
  n_col <- ncol(df)
  hdr_lines <- strsplit(names(df), "\n", fixed = TRUE)
  n_hdr <- max(lengths(hdr_lines))
  grid::grid.newpage()
  left <- first_col_x
  right <- 0.99

  if (identical(layout, "fixed")) {
    fs <- fontsize
    # centres[1] is unused: column 1 is left-aligned at `left`.
    centres <- c(NA_real_, seq(0.36, 0.955, length.out = n_col - 1L))
  } else {
    gap <- 0.012
    avail <- right - left

    # Width of the widest entry per column, headers measured bold as they are drawn.
    widths_at <- function(fs) {
      vapply(seq_len(n_col), function(j) {
        vp_h <- grid::viewport(gp = grid::gpar(fontsize = fs, fontface = "bold"))
        grid::pushViewport(vp_h)
        wh <- max(grid::convertWidth(grid::stringWidth(hdr_lines[[j]]), "npc", valueOnly = TRUE))
        grid::popViewport()
        cells <- as.character(df[[j]])
        wc <- 0
        if (length(cells)) {
          vp_c <- grid::viewport(gp = grid::gpar(fontsize = fs))
          grid::pushViewport(vp_c)
          wc <- max(grid::convertWidth(grid::stringWidth(cells), "npc", valueOnly = TRUE))
          grid::popViewport()
        }
        max(wh, wc)
      }, numeric(1))
    }

    fs <- fontsize
    cw <- widths_at(fs)
    need <- sum(cw) + gap * (n_col - 1)
    if (need > avail) {
      fs <- max(4.5, fs * (avail / need) * 0.98)
      cw <- widths_at(fs)
      need <- sum(cw) + gap * (n_col - 1)
    }
    cw <- cw + max(0, avail - need) / n_col   # spread any slack

    # Band edges; column 1 is left-aligned, the rest centred inside their band.
    edges <- left + c(0, cumsum(cw + gap))
    centres <- edges[-1] - gap - cw / 2
  }

  top <- 0.955
  # Extra header lines claim vertical space so the body is not squeezed.
  row_h <- (top - 0.03) / (n_row + 2.2 + 0.8 * (n_hdr - 1))
  line_h <- 0.8 * row_h
  y_head <- top - 1.1 * row_h                       # first (topmost) header line
  y_mid <- y_head - (n_hdr - 1) * line_h - 0.62 * row_h

  # Header lines are drawn top-down and the block is bottom-aligned on the rule,
  # so a one-line header in a table that has two-line headers elsewhere sits on
  # the same baseline as their last line.
  draw_header <- function(lines, x, just) {
    pad <- n_hdr - length(lines)
    for (k in seq_along(lines)) {
      grid::grid.text(lines[[k]], x = x, y = y_head - (pad + k - 1) * line_h,
                      just = just, gp = grid::gpar(fontsize = fs, fontface = "bold"))
    }
  }

  grid::grid.lines(x = c(left, right), y = c(top, top), gp = grid::gpar(lwd = 1.3))
  draw_header(hdr_lines[[1]], left, "left")
  for (j in seq_len(n_col - 1L)) {
    draw_header(hdr_lines[[j + 1L]], centres[j + 1L], "centre")
  }
  grid::grid.lines(x = c(left, right), y = c(y_mid, y_mid), gp = grid::gpar(lwd = 0.7))
  for (i in seq_len(n_row)) {
    yy <- y_mid - i * row_h
    grid::grid.text(as.character(df[[1]][i]), x = left, y = yy, just = "left",
                    gp = grid::gpar(fontsize = fs))
    for (j in seq_len(n_col - 1L)) {
      grid::grid.text(as.character(df[[j + 1L]][i]), x = centres[j + 1L], y = yy, just = "centre",
                      gp = grid::gpar(fontsize = fs))
    }
  }
  y_bottom <- y_mid - (n_row + 0.55) * row_h
  grid::grid.lines(x = c(left, right), y = c(y_bottom, y_bottom), gp = grid::gpar(lwd = 1.3))
  invisible(TRUE)
}

# Page height is set from the row count so the PDF crops to the table.
write_bare_table <- function(df, output_pdf, width = 9, fontsize = 8,
                             row_in = 0.20, pad_in = 0.45,
                             layout = c("fixed", "autofit")) {
  layout <- match.arg(layout)
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  n_hdr <- max(lengths(strsplit(names(df), "\n", fixed = TRUE)))
  height <- pad_in + (nrow(df) + 2 + 0.8 * (n_hdr - 1)) * row_in
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(output_pdf, width = width, height = height)
  } else {
    grDevices::pdf(output_pdf, width = width, height = height)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  draw_bare_table(df, fontsize = fontsize, layout = layout)
  TRUE
}

# --- panel builders: data frame per supplementary table ----------------------

supp_delta_auc_panel <- function(tab) {
  tab |>
    dplyr::arrange(reference, omics, model) |>
    dplyr::transmute(
      Reference = reference, Omics = tools::toTitleCase(omics), Model = model,
      `Test n (deaths)` = sprintf("%d (%d)", n_test, n_deaths),
      AUC = sprintf("%.3f", auc),
      `AUC difference (95% CI)` = sprintf("%+.3f (%+.3f to %+.3f)", delta, delta_lo, delta_hi),
      P = fmt_p(p_value), `P (BH)` = fmt_p(p_adj)
    )
}

supp_calibration_panel <- function(tab) {
  tab |>
    dplyr::transmute(
      Omics = tools::toTitleCase(omics), Model = model,
      `n (events)` = sprintf("%d (%d)", n, events),
      `Predicted / observed` = sprintf("%.2f / %.2f", mean_predicted, observed),
      `Brier (null)` = sprintf("%.3f (%.3f)", brier, brier_null),
      Intercept = sprintf("%+.2f", calibration_intercept),
      Slope = sprintf("%.2f", calibration_slope)
    )
}

supp_age_proxy_panel <- function(tab) {
  tab |>
    dplyr::arrange(omics, model) |>
    dplyr::transmute(
      Omics = tools::toTitleCase(omics), Model = model,
      `Spearman rho with age` = sprintf("%+.2f", rho_age),
      AUC = ifelse(is.na(auc_omics), "-", sprintf("%.3f", auc_omics)),
      `AUC age only` = ifelse(is.na(auc_age_only), "-", sprintf("%.3f", auc_age_only)),
      `AUC difference` = ifelse(is.na(delta), "-", sprintf("%+.3f", delta)),
      `P (BH)` = ifelse(is.na(p_adj), "-", fmt_p(p_adj))
    )
}

supp_table1_panel <- function(tab, labels) {
  out <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(out) <- labels
  out
}

# Overlap percentages print as whole numbers and concordance to one decimal,
# matching how these two tables have always been reported.
supp_rank_overlap_panel <- function(tab) {
  data.frame(
    `Omics` = tools::toTitleCase(tab$omics),
    `Top N` = as.character(tab$top_n),
    `Overlap count` = as.character(tab$overlap),
    `Overlap\n(% of N)` = sprintf("%.0f", tab$overlap_pct),
    `Concordant directions\n(n)` = as.character(tab$concordant),
    `Concordant directions\n(% of overlap)` =
      ifelse(is.na(tab$concordant_pct), "", sprintf("%.1f", tab$concordant_pct)),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Matched-cohort refit (06_fit_matched_cohort_models.R)
#
# Unlike the layer-specific analysis, every model here is fitted on the same 332
# men using the same folds, so there is a single clinical benchmark with a single
# AUC and the cross-layer comparison needs no provenance caveat. These builders
# read the saved base-model and stacked outer-fold predictions produced by 06_.
# =============================================================================

load_matched_cohort_outputs <- function(paths, require_provenance = TRUE) {
  f <- file.path(paths$data, "model_outputs", "matched_cohort",
                 "matched_cohort_model_outputs.rds")
  if (!file.exists(f)) {
    stop("Missing matched-cohort output; run 06_fit_matched_cohort_models.R first: ", f)
  }
  x <- readRDS(f)
  if (isTRUE(require_provenance) &&
      (is.null(x$provenance) || is.null(x$provenance$signature))) {
    stop("Matched-cohort output lacks provenance; rerun 06_fit_matched_cohort_models.R")
  }
  required <- c("cohort", "design", "predictions", "predictions_by_repeat", "tuning")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop("Matched-cohort output is incomplete: missing ", paste(missing, collapse = ", "))
  if (!identical(as.integer(x$cohort$n), 332L) ||
      !identical(as.integer(x$cohort$n_deaths), 169L)) {
    stop("Matched-cohort output does not describe the prespecified 332-men, 169-death cohort")
  }
  p <- x$predictions
  if (!is.data.frame(p) || !all(c("ID", "truth", "risk", "model", "omics") %in% names(p))) {
    stop("Matched-cohort predictions lack the required ID/truth/risk/model/omics columns")
  }
  if (anyNA(p$risk) || any(!is.finite(p$risk)) || any(p$risk < 0 | p$risk > 1)) {
    stop("Matched-cohort predictions contain invalid risks")
  }
  if (any(!as.character(p$truth) %in% c("Active", "Deceased"))) {
    stop("Matched-cohort predictions contain an unknown vital-status label")
  }
  cohort_ids <- as.character(x$cohort$ids)
  if (length(cohort_ids) != x$cohort$n || anyNA(cohort_ids) ||
      anyDuplicated(cohort_ids)) {
    stop("Matched-cohort participant IDs are missing or duplicated")
  }
  cohort_truth <- as.character(x$cohort$status)
  if (length(cohort_truth) != length(cohort_ids)) {
    stop("Matched-cohort status vector does not match participant IDs")
  }
  expected_keys <- c(
    "RF clinical|clinical", "LR age-only|clinical",
    paste(rep(c("EN omics", "XGB omics", "RF clinical + EN", "RF clinical + XGB"),
              each = 3),
          rep(c("microbiome", "metabolomics", "proteomics"), 4), sep = "|")
  )
  keys <- paste(p$model, p$omics, sep = "|")
  if (!all(expected_keys %in% unique(keys))) {
    stop("Matched-cohort predictions are missing one or more prespecified model/layer combinations")
  }
  for (key in unique(keys)) {
    d <- p[keys == key, , drop = FALSE]
    if (nrow(d) != x$cohort$n || anyDuplicated(as.character(d$ID)) ||
        !setequal(as.character(d$ID), cohort_ids)) {
      stop("Matched-cohort predictions must contain one row per participant for ", key)
    }
    d <- d[match(cohort_ids, as.character(d$ID)), , drop = FALSE]
    if (!identical(as.character(d$truth), cohort_truth)) {
      stop("Matched-cohort truth labels are not aligned for ", key)
    }
  }
  if (is.data.frame(x$predictions_by_repeat) && nrow(x$predictions_by_repeat)) {
    br <- x$predictions_by_repeat
    if (!all(c("ID", "truth", "risk", "model", "omics", "repeat_id", "fold_id") %in% names(br))) {
      stop("Matched repeated predictions lack required resampling columns")
    }
    br_keys <- paste(br$model, br$omics, sep = "|")
    expected_repeat_keys <- expected_keys
    if (!all(expected_repeat_keys %in% unique(br_keys))) {
      stop("Matched repeated predictions are missing one or more prespecified model/layer combinations")
    }
    for (key in unique(br_keys)) {
      d <- br[br_keys == key, , drop = FALSE]
      if (anyDuplicated(paste(d$repeat_id, d$ID, sep = "|")) ||
          nrow(d) != x$cohort$n * x$design$n_repeats) {
        stop("Matched repeated predictions are incomplete for ", key)
      }
      repeats <- unique(d$repeat_id)
      if (length(repeats) != x$design$n_repeats ||
          any(vapply(repeats, function(rp) {
            ids <- as.character(d$ID[d$repeat_id == rp])
            length(ids) == x$cohort$n && setequal(ids, cohort_ids)
          }, logical(1)) == FALSE)) {
        stop("Matched repeated predictions do not cover the cohort in every repeat for ", key)
      }
      if (anyNA(d$risk) || any(!is.finite(d$risk)) || any(d$risk < 0 | d$risk > 1)) {
        stop("Matched repeated predictions contain invalid risks for ", key)
      }
    }
    ref <- br[br_keys == "RF clinical|clinical",
              c("repeat_id", "ID", "fold_id"), drop = FALSE]
    ref <- ref[order(ref$repeat_id, ref$ID), , drop = FALSE]
    for (key in expected_repeat_keys) {
      d <- br[br_keys == key, c("repeat_id", "ID", "fold_id"), drop = FALSE]
      d <- d[order(d$repeat_id, d$ID), , drop = FALSE]
      if (!identical(as.character(d$repeat_id), as.character(ref$repeat_id)) ||
          !identical(as.character(d$ID), as.character(ref$ID)) ||
          !identical(as.character(d$fold_id), as.character(ref$fold_id))) {
        stop("Matched models do not share identical outer-fold assignments for ", key)
      }
    }
  } else {
    stop("Matched-cohort output has no repeated out-of-fold predictions")
  }
  tuning <- x$tuning
  if (!is.data.frame(tuning) ||
      !all(c("model", "omics", "repeat_id", "fold_id", "n_features",
             "selected_features") %in% names(tuning))) {
    stop("Matched-cohort tuning output lacks fold-level feature-selection records")
  }
  prot <- tuning[tuning$omics == "proteomics" &
                   tuning$model %in% c("EN omics", "XGB omics"), , drop = FALSE]
  expected_prot_rows <- 2L * x$design$outer_v * x$design$n_repeats
  if (nrow(prot) != expected_prot_rows ||
      any(lengths(prot$selected_features) != prot$n_features) ||
      any(lengths(prot$selected_features) < 1L)) {
    stop("Matched proteomics output lacks complete outer-fold limma selections")
  }
  x
}

matched_layer_label <- function(x) tools::toTitleCase(x)

# The three omics layers, in reporting order.
matched_layers <- function() c("microbiome", "metabolomics", "proteomics")

# Layer-level models present in the matched nested-CV output, in reporting order.
matched_omics_models <- function(mc) {
  order_m <- c("EN omics", "XGB omics", "RF clinical + EN", "RF clinical + XGB")
  order_m[order_m %in% unique(mc$predictions$model)]
}

matched_risk <- function(mc, model, omics) {
  d <- mc$predictions[mc$predictions$model == model & mc$predictions$omics == omics, ]
  if (!nrow(d)) stop("no predictions for ", model, " / ", omics)
  d[order(d$ID), ]
}

# --- Table: discrimination in the matched cohort ---------------------------
build_matched_auc_table <- function(mc, ...) {
  rows <- list()
  add <- function(model, omics) {
    d <- matched_risk(mc, model, omics)
    ci <- as.numeric(pROC::ci.auc(roc_death(d$truth, d$risk)))
    rows[[length(rows) + 1L]] <<- tibble::tibble(
      model = model, omics = omics, layer = matched_layer_label(omics),
      n = nrow(d), n_deaths = sum(d$truth == "Deceased"),
      auc = ci[[2]], lower = ci[[1]], upper = ci[[3]]
    )
  }
  add("RF clinical", "clinical")
  add("LR age-only", "clinical")
  for (om in matched_layers()) for (m in matched_omics_models(mc)) add(m, om)
  dplyr::bind_rows(rows)
}

# --- Table: pairwise cross-layer comparisons -------------------------------
#
# Only the omics models appear here. The clinical model is fitted once on the
# matched cohort, so it has a single AUC and nothing to compare across layers --
# which is the entire reason for the refit.
build_matched_crosslayer_tests <- function(mc, n_boot = 2000L, seed = 123L,
                                           ...) {
  set.seed(seed)
  layers <- matched_layers()
  rows <- list()
  for (m in matched_omics_models(mc)) {
    d <- lapply(layers, function(om) matched_risk(mc, m, om))
    names(d) <- layers
    truth <- d[[1]]$truth
    pairs <- utils::combn(layers, 2L)
    for (j in seq_len(ncol(pairs))) {
      a <- pairs[1L, j]; b <- pairs[2L, j]
      ra <- roc_death(truth, d[[a]]$risk); rb <- roc_death(truth, d[[b]]$risk)
      ci <- delta_auc_ci(truth, d[[a]]$risk, d[[b]]$risk, n_boot = n_boot)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        model = m,
        comparison = paste(matched_layer_label(a), "vs", matched_layer_label(b)),
        auc_a = as.numeric(pROC::auc(ra)), auc_b = as.numeric(pROC::auc(rb)),
        delta = as.numeric(pROC::auc(ra)) - as.numeric(pROC::auc(rb)),
        delta_lo = ci[[1]], delta_hi = ci[[2]],
        p_value = pROC::roc.test(ra, rb, method = "delong", paired = TRUE)$p.value
      )
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::group_by(model) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) |>
    dplyr::ungroup()
}

# --- Table: omics panels against the single clinical model -----------------
build_matched_vs_clinical <- function(mc, n_boot = 2000L, seed = 123L,
                                      ...) {
  set.seed(seed)
  clin <- matched_risk(mc, "RF clinical", "clinical")
  r_clin <- roc_death(clin$truth, clin$risk)
  rows <- list()
  for (om in matched_layers()) {
    for (m in matched_omics_models(mc)) {
      d <- matched_risk(mc, m, om)
      r_om <- roc_death(d$truth, d$risk)
      ci <- delta_auc_ci(clin$truth, d$risk, clin$risk, n_boot = n_boot)
      rows[[length(rows) + 1L]] <- tibble::tibble(
        omics = om, layer = matched_layer_label(om), model = m, n = nrow(d),
        auc_omics = as.numeric(pROC::auc(r_om)),
        auc_clinical = as.numeric(pROC::auc(r_clin)),
        delta = as.numeric(pROC::auc(r_om)) - as.numeric(pROC::auc(r_clin)),
        delta_lo = ci[[1]], delta_hi = ci[[2]],
        p_value = pROC::roc.test(r_om, r_clin, method = "delong", paired = TRUE)$p.value
      )
    }
  }
  dplyr::bind_rows(rows) |>
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"))
}

supp_matched_auc_panel <- function(tab) {
  tab |>
    dplyr::transmute(
      Model = model, Data = layer,
      `n (deaths)` = sprintf("%d (%d)", n, n_deaths),
      `AUC (95% CI)` = sprintf("%.3f (%.3f-%.3f)", auc, lower, upper)
    )
}

supp_matched_test_panel <- function(tab) {
  tab |>
    dplyr::arrange(model, comparison) |>
    dplyr::transmute(
      Model = model, Comparison = comparison,
      `AUC difference (95% CI)` = sprintf("%+.3f (%+.3f to %+.3f)", delta, delta_lo, delta_hi),
      P = fmt_p(p_value), `P (BH)` = fmt_p(p_adj)
    )
}

supp_matched_vs_clinical_panel <- function(tab) {
  tab |>
    dplyr::transmute(
      Data = layer, Model = model,
      `AUC omics` = sprintf("%.3f", auc_omics),
      `AUC clinical` = sprintf("%.3f", auc_clinical),
      `AUC difference (95% CI)` = sprintf("%+.3f (%+.3f to %+.3f)", delta, delta_lo, delta_hi),
      P = fmt_p(p_value), `P (BH)` = fmt_p(p_adj)
    )
}
matched_logit <- function(p, eps = 1e-6) {
  p <- pmin(pmax(p, eps), 1 - eps)
  log(p / (1 - p))
}

# Cross-layer contrasts and omics-versus-clinical contrasts share a column
# structure, so they are reported as one table with a Comparison column rather
# than as lettered sub-panels. Benjamini-Hochberg adjustment stays within its
# own family -- the two ask different questions -- and the family is named in
# the table so the reader can see which P values were adjusted together.
build_matched_comparisons <- function(mc, n_boot = 2000L, seed = 123L,
                                      ...) {
  cross <- build_matched_crosslayer_tests(mc, n_boot = n_boot, seed = seed) |>
    dplyr::transmute(family = "Between omics layers", model, comparison,
                     auc_a, auc_b, delta, delta_lo, delta_hi, p_value, p_adj)
  vsc <- build_matched_vs_clinical(mc, n_boot = n_boot, seed = seed) |>
    dplyr::transmute(family = "Versus the clinical model", model,
                     comparison = paste(layer, "vs Clinical"),
                     auc_a = auc_omics, auc_b = auc_clinical,
                     delta, delta_lo, delta_hi, p_value, p_adj)
  dplyr::bind_rows(cross, vsc)
}

# The two AUCs being differenced are printed alongside the difference. Without
# them the table states that two layers are indistinguishable but not what either
# scored, so the point estimates quoted in the Results had no table to land in.
# They are given in the order named by the Comparison cell.
supp_matched_comparisons_panel <- function(tab) {
  tab |>
    dplyr::transmute(
      Family = family, Model = model, Comparison = comparison,
      `AUC (first vs second)` = sprintf("%.3f vs %.3f", auc_a, auc_b),
      `AUC difference (95% CI)` = sprintf("%+.3f (%+.3f to %+.3f)", delta, delta_lo, delta_hi),
      P = fmt_p(p_value), `P (BH)` = fmt_p(p_adj)
    )
}
