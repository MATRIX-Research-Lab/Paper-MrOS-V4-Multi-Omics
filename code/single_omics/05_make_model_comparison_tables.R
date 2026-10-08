# Supplementary tables: paired model comparisons, calibration, and the
# composition of the complete-case subset.
#
# All are computed from the saved model-output RDS files; no model is refitted,
# so the AUCs reproduce those in Table 2 exactly.
#
# Output convention: one PDF per table panel, containing the table body only.
# Captions live in manuscript/supp.tex, which includes these files with
# \includegraphics. Filenames describe content rather than a supplementary
# number, so tables can be reordered without renaming anything.
# CSV copies are written alongside for inspection.

model_outputs <- load_prediction_model_outputs(paths)

emit <- function(artifact, expr) {
  unlink(file.path(paths$tables, artifact))
  tryCatch(
    force(expr),
    error = function(err) {
      record_artifact(artifact, FALSE, "generation_failed", conditionMessage(err))
      record_issue("05_make_model_comparison_tables", "generation failed",
                   paste(artifact, conditionMessage(err)))
    }
  )
}

tbl <- function(name) file.path(paths$tables, name)
note <- function(artifact, detail) {
  record_artifact(artifact, file.exists(tbl(artifact)),
                  "generated_from_model_outputs", detail)
}

# --- within-layer paired model comparisons ---------------------------------
emit("SuppTable_model_comparisons.pdf", {
  delta_tab <- build_delta_auc_table(model_outputs)
  write_bare_table(supp_delta_auc_panel(delta_tab), tbl("SuppTable_model_comparisons.pdf"),
                   width = 10, fontsize = 7.5)
  utils::write.csv(delta_tab, tbl("SuppTable_model_comparisons.csv"), row.names = FALSE)
  note("SuppTable_model_comparisons.pdf", "Paired DeLong comparisons against the clinical and age-only reference models.")
})

# The cross-layer comparison built here from the primary layer-specific fits
# (SuppTable_crosslayer_*) and the pooled out-of-fold AUC table were removed on
# 2026-08-12. The cross-layer question is now answered by the matched-cohort
# refit in 06_/07_, which scores every layer on one shared set of folds instead
# of splicing together predictions of mixed provenance; the pooled out-of-fold
# AUCs were never reported in the manuscript. Same reason for
# SuppTable_omics_vs_clinical, whose comparison is now a family within
# SuppTable_matched_comparisons.

# --- calibration ------------------------------------
emit("SuppTable_calibration.pdf", {
  calib <- build_calibration_table(model_outputs)
  write_bare_table(supp_calibration_panel(calib), tbl("SuppTable_calibration.pdf"), width = 10, fontsize = 7.5)
  utils::write.csv(calib, tbl("SuppTable_calibration.csv"), row.names = FALSE)
  note("SuppTable_calibration.pdf", "Calibration intercept, slope and Brier score on the held-out test set.")
})

# --- age proxy --------------------------------------
emit("SuppTable_age_proxy.pdf", {
  ap <- build_age_proxy_table(model_outputs)
  write_bare_table(supp_age_proxy_panel(ap), tbl("SuppTable_age_proxy.pdf"), width = 10, fontsize = 7.5)
  utils::write.csv(ap, tbl("SuppTable_age_proxy.csv"), row.names = FALSE)
  note("SuppTable_age_proxy.pdf", "Correlation of each omics risk score with chronological age.")
})

# --- participant characteristics, one table per layer ----------------------
#
# Table 1 describes the union of the three analytic samples; the layers do not
# share a sample, so the supplement reports each layer separately. Same builder,
# applied per layer.
emit("SuppTable_characteristics_microbiome.pdf", {
  ch <- build_layer_characteristics_tables(model_outputs)
  for (om in names(ch)) {
    art <- sprintf("SuppTable_characteristics_%s.pdf", om)
    write_bare_table(
      supp_table1_panel(ch[[om]]$table,
                        c("Variable", "Overall", "Deceased", "Active", "P-value")),
      tbl(art), width = 8, fontsize = 7.5)
    utils::write.csv(ch[[om]]$table,
                     tbl(sprintf("SuppTable_characteristics_%s.csv", om)),
                     row.names = FALSE)
    note(art, sprintf("Characteristics of the %d %s participants (%d deceased, %d alive).",
                      ch[[om]]$n, om, ch[[om]]$deaths, ch[[om]]$alive))
  }
})

# --- feature-ranking overlap and directional concordance --------------------
#
# Both tables were hand-entered in the supplement and went stale; see
# manuscript/bug_fix.md section 3. They are computed here from the same
# extractors that build Figures 3 and 4, so they cannot drift from the figures.
emit("SuppTable_en_xgb_overlap.pdf", {
  tab <- build_rank_overlap_table(model_outputs, ranks_en_vs_xgb)
  # layout = "autofit": the fixed layout spaces columns evenly and centres the
  # last one at 0.955, so a wide header ("Concordant directions") ran off the
  # right edge and was clipped. Autofit measures the rendered strings instead.
  write_bare_table(supp_rank_overlap_panel(tab), tbl("SuppTable_en_xgb_overlap.pdf"),
                   width = 9, fontsize = 7.5, layout = "autofit")
  utils::write.csv(tab, tbl("SuppTable_en_xgb_overlap.csv"), row.names = FALSE)
  note("SuppTable_en_xgb_overlap.pdf",
       "Elastic net versus XGBoost, omics-only: rank overlap and directional concordance at the top 20/50/100.")
})

emit("SuppTable_ensemble_overlap.pdf", {
  tab <- build_rank_overlap_table(model_outputs, ranks_stacked_ensembles)
  write_bare_table(supp_rank_overlap_panel(tab), tbl("SuppTable_ensemble_overlap.pdf"),
                   width = 9, fontsize = 7.5, layout = "autofit")
  utils::write.csv(tab, tbl("SuppTable_ensemble_overlap.csv"), row.names = FALSE)
  note("SuppTable_ensemble_overlap.pdf",
       "RF + elastic net versus RF + XGBoost stacked ensembles: rank overlap and directional concordance.")
})

# --- the complete-case subset versus the remainder ---
#
# The companion panel describing the 332 by vital status (SuppTable_complete_case)
# was removed on 2026-08-12: it was never cited. The mortality contrast the
# manuscript does cite is now a Vital Status row of the table below.
emit("SuppTable_complete_case_vs_rest.pdf", {
  cc <- build_complete_case_tables(model_outputs)
  # Group sizes come from the table itself rather than being restated here, so
  # the headers cannot drift from the counts the rows were computed on. They go
  # on a second header line to avoid widening the columns.
  n <- attr(cc$membership, "group_sizes")
  write_bare_table(
    supp_table1_panel(cc$membership,
                      c("Variable",
                        sprintf("Overall\n(N = %d)", n[["overall"]]),
                        sprintf("All three layers\n(N = %d)", n[["deceased"]]),
                        sprintf("Remainder\n(N = %d)", n[["alive"]]),
                        "P-value")),
    tbl("SuppTable_complete_case_vs_rest.pdf"), width = 8.5, fontsize = 7.5)
  utils::write.csv(cc$membership, tbl("SuppTable_complete_case_vs_rest.csv"), row.names = FALSE)
  note("SuppTable_complete_case_vs_rest.pdf",
       sprintf("The %d participants with all three layers versus the remaining %d, including the mortality contrast.",
               cc$n_cc, cc$n_other))
})
