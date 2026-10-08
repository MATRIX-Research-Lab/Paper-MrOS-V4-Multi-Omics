ready <- model_refit_inputs_available(paths)

if (!ready) {
  record_issue(
    "02_fit_prediction_models",
    "data-path issue",
    "Model refitting was not started because processed inputs or model packages are missing."
  )
} else {
  record_issue(
    "02_fit_prediction_models",
    "model refit not included",
    paste(
      "Saved model-output RDS files are used as the default reproducibility inputs.",
      "Full model refitting is not part of this cleaned public workflow."
    )
  )
}

# Age-only baseline: logistic regression (status ~ age_v4). This is a lightweight
# fit that reuses the train/test split stored inside each model-output RDS, so it
# does NOT refit any of the five main models and needs no processed inputs. The
# fitted baselines are saved next to the model-output files for reuse by Fig2 and
# the AUC-CI tables.
tryCatch(
  {
    age_model_outputs <- load_prediction_model_outputs(paths)
    saved <- fit_and_save_age_baselines(age_model_outputs, paths)
    record_issue(
      "02_fit_prediction_models",
      "age-only baseline fit",
      paste0(
        "Fitted and saved LR age-only baseline (status ~ age_v4) for: ",
        paste(names(saved), collapse = ", "), "."
      )
    )
  },
  error = function(err) {
    record_issue("02_fit_prediction_models", "age baseline failed", conditionMessage(err))
  }
)
