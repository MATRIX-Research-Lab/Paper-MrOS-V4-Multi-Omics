# Supplementary tables for the matched-cohort refit.
#
# All base-model and stacked predictions are computed during the matched nested
# CV stage and read from its provenance-validated output. Table generation does
# not fit or tune models. Output follows 05_: one PDF per table panel
# containing the table body only and no caption, filenames describing content
# rather than a supplementary number so tables can be reordered without
# renaming. CSV copies are written alongside.

assert_packages(c("dplyr", "tidymodels", "glmnet", "pROC"),
                "matched-cohort table generation")
mc <- load_matched_cohort_outputs(paths, require_provenance = TRUE)

emit <- function(artifact, expr) {
  tryCatch(
    force(expr),
    error = function(err) {
      record_artifact(artifact, FALSE, "generation_failed", conditionMessage(err))
      record_issue("07_make_matched_cohort_tables", "generation failed",
                   paste(artifact, conditionMessage(err)), severity = "fatal")
      stop(err)
    }
  )
}

tbl <- function(name) file.path(paths$tables, name)
note <- function(artifact, detail) {
  record_artifact(artifact, file.exists(tbl(artifact)), "generated_from_matched_refit", detail)
}

message(sprintf("Matched cohort: n = %d, deaths = %d, %d repeats of %d-fold outer CV",
                mc$cohort$n, mc$cohort$n_deaths, mc$design$n_repeats, mc$design$outer_v))

# --- discrimination by layer, single clinical benchmark --------------------
emit("SuppTable_matched_auc.pdf", {
  tab <- build_matched_auc_table(mc)
  # layout = "autofit": these tables carry long model names, so column widths
  # are measured from the rendered strings rather than evenly spaced. See
  # draw_bare_table() -- the primary supplementary tables use the fixed layout.
  write_bare_table(supp_matched_auc_panel(tab), tbl("SuppTable_matched_auc.pdf"),
                   width = 8, layout = "autofit")
  utils::write.csv(tab, tbl("SuppTable_matched_auc.csv"), row.names = FALSE)
  note("SuppTable_matched_auc.pdf",
       sprintf("Refit on the %d matched participants; one clinical model, identical folds.",
               mc$cohort$n))
})

# --- all paired comparisons: between layers, and versus the clinical model ---
emit("SuppTable_matched_comparisons.pdf", {
  tab <- build_matched_comparisons(mc)
  write_bare_table(supp_matched_comparisons_panel(tab),
                   tbl("SuppTable_matched_comparisons.pdf"), width = 10, fontsize = 7.5,
                   layout = "autofit")
  utils::write.csv(tab, tbl("SuppTable_matched_comparisons.csv"), row.names = FALSE)
  note("SuppTable_matched_comparisons.pdf",
       "Paired DeLong contrasts between omics layers and against the single clinical model.")
})
