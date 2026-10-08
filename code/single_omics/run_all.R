# Render every manuscript figure and table from the saved model outputs.
#
#   Rscript code/single_omics/run_all.R
#
# Two analyses are rendered here, from two different saved objects:
#
#   01_-05_  the primary analysis -- each omics layer fitted in its own analytic
#            sample, read from data/model_outputs/single_omics/
#   07_      the matched-cohort refit -- all layers on the same 332 men with
#            shared folds, read from data/model_outputs/matched_cohort/
#
# Nothing here fits a model. The two expensive stages are deliberately outside
# this workflow and are run on their own:
#
#   code/single_omics/preprocess_data.R              rebuilds the frozen analytic frames from
#                                       the restricted raw data
#   code/single_omics/run_omics_pred.R               refits the primary layers
#   code/single_omics/06_fit_matched_cohort_models.R refits the matched cohort (over an hour)
#
# Keeping them out means a clean checkout renders from the frozen inputs and
# cannot silently regenerate them from data it may not even have.

find_workflow_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  starts <- c(getwd(), if (length(fa)) dirname(fa[[1]]) else NULL)
  for (s in unique(starts)) {
    d <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(20L)) {
      if (file.exists(file.path(d, "code", "single_omics", "00_setup.R"))) return(d)
      p <- dirname(d)
      if (identical(p, d)) break
      d <- p
    }
  }
  stop("Could not locate the repository root")
}

workflow_root <- find_workflow_root()
setwd(workflow_root)
source(file.path(workflow_root, "code", "single_omics", "00_setup.R"))

run_step <- function(label, file) {
  message(label)
  tryCatch(
    {
      source(file, local = new.env(parent = globalenv()))
      TRUE
    },
    error = function(err) {
      record_issue(label, "execution_error", conditionMessage(err), severity = "fatal")
      FALSE
    }
  )
}

steps <- list(
  c("01 prepare data", "code/single_omics/01_prepare_data.R"),
  c("02 fit prediction models", "code/single_omics/02_fit_prediction_models.R"),
  c("03 make manuscript tables", "code/single_omics/03_make_manuscript_tables.R"),
  c("04 make manuscript figures", "code/single_omics/04_make_manuscript_figures.R"),
  c("05 make model comparison tables", "code/single_omics/05_make_model_comparison_tables.R"),
  c("07 make matched cohort tables", "code/single_omics/07_make_matched_cohort_tables.R")
)

for (x in steps) {
  ok <- run_step(x[[1]], x[[2]])
  if (!isTRUE(ok)) {
    message("Stopping after fatal step: ", x[[1]])
    break
  }
}

# Figures produced by the other entry points -- code/cross_omics for the
# heat maps, LaTeX for the workflow schematics. Confirm they survived the render
# rather than claiming credit for them: results/ used to lose them on every run.
for (a in names(external_artifact_source)) {
  f <- artifact_output_file(a)
  present <- file.exists(f)
  producer <- external_artifact_source[[a]]
  record_artifact(
    a, present,
    if (present) "present_not_regenerated" else "missing_external_output",
    paste0("Produced by ", producer, "; run_all.R only verifies it is still present.")
  )
}

# Explicit generation summary (so a successful run is visible in the console).
message("")
message("=== Generated artifacts ===")
for (dir_label in c("figures", "tables")) {
  d <- paths[[dir_label]]
  files <- list.files(d, full.names = TRUE)
  message(dir_label, " -> ", d)
  if (!length(files)) {
    message("  (empty)")
    next
  }
  for (f in sort(files)) {
    message(sprintf("  %-55s %8.1f KB", basename(f), file.size(f) / 1024))
  }
}
message("")
if (exists("issues", envir = globalenv()) && nrow(get("issues", envir = globalenv())) > 0) {
  message("=== Issues / notes ===")
  print(get("issues", envir = globalenv()))
}
if (exists("artifact_status", envir = globalenv())) {
  message("=== Artifact status ===")
  print(get("artifact_status", envir = globalenv())[, c("artifact", "generated", "status")])
}

utils::write.csv(issues, file.path(paths$results, "workflow_issues.csv"), row.names = FALSE)
utils::write.csv(artifact_status, file.path(paths$results, "workflow_artifact_status.csv"), row.names = FALSE)
fatal_issues <- issues[issues$severity == "fatal", , drop = FALSE]
missing_artifacts <- artifact_status[!artifact_status$generated, , drop = FALSE]
if (nrow(fatal_issues) || nrow(missing_artifacts)) {
  stop("Workflow did not complete: ", nrow(fatal_issues), " fatal issue(s) and ",
       nrow(missing_artifacts), " missing artifact(s). See results/workflow_issues.csv and workflow_artifact_status.csv.")
}
message("Done. Open: ", paths$results)
