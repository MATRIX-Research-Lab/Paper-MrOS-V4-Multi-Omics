table_pdf <- artifact_output_file("Table1.pdf")

# Table 1 is computed from the union of the three analytic samples. There is no
# fallback to the stored literal and compare_table1()'s verdict is acted on, not
# merely logged: silently substituting the hard-coded values would make a failed
# reconstruction look like a successful run, which is how Table 1 came to be
# maintained by hand in the first place.
table1 <- compute_table1(load_prediction_model_outputs(paths))
if (!compare_table1(table1)) {
  stop("Computed Table 1 does not match the stored reference; inspect the data harmonisation before publishing.")
}

made_table <- write_table1_pdf(table1, table_pdf)
made_docx <- tryCatch(
  write_table1_docx(table1, file.path(paths$tables, "Table1.docx")),
  error = function(e) FALSE
)
record_artifact(
  "Table1.docx", isTRUE(made_docx),
  if (isTRUE(made_docx)) "generated" else "generation_failed",
  "Word version of Table 1 for submission, generated from the same computed table as the PDF."
)

if (made_table && file.exists(table_pdf)) {
  record_artifact("Table1.pdf", TRUE, "generated", "Regenerated from processed metadata.")
} else {
  record_artifact("Table1.pdf", FALSE, "missing_source", "Processed metadata was not sufficient to regenerate Table 1.")
  record_issue(
    "03_make_manuscript_tables",
    "missing source code",
    "Processed metadata was not sufficient to regenerate Table 1."
  )
}

# Table 2: three-line (booktabs) table of test-set AUC with 95% CI for all
# models (including the LR age-only baseline), one column per omics. Also written
# to Excel for convenience. Computed from the saved model-output predictions.
auc_ci_pdf <- artifact_output_file("Table2_AUC_CI.pdf")
made_auc_ci <- tryCatch(
  {
    model_outputs <- load_prediction_model_outputs(paths)
    write_auc_ci_three_line_pdf(model_outputs, auc_ci_pdf)
    write_auc_ci_table(model_outputs, file.path(paths$tables, "AUC_CI_table.xlsx"))
    TRUE
  },
  error = function(err) {
    record_issue("03_make_manuscript_tables", "auc ci table failed", conditionMessage(err))
    FALSE
  }
)

if (made_auc_ci && file.exists(auc_ci_pdf)) {
  record_artifact(
    "Table2_AUC_CI.pdf", TRUE, "generated_from_model_outputs",
    "AUC with 95% DeLong CI for all models incl LR age-only baseline; also exported to AUC_CI_table.xlsx."
  )
} else {
  record_artifact(
    "Table2_AUC_CI.pdf", FALSE, "missing_source",
    "Model-output RDS files were not available to build the AUC-CI table."
  )
}
