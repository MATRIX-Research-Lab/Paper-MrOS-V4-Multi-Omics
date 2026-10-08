input_status <- check_processed_inputs(paths)

raw_required_rel <- file.path("data", "raw_data", c(
  file.path("V1FEB24", "v1feb24.sas7bdat"),
  file.path("V4FEB24", "v4feb24.sas7bdat"),
  file.path("OMICS", "Microbiome", "DADA_Silva_Mros_data_release.rds"),
  file.path("OMICS", "Metabolomics", "CPMC-0101-23PHML+ DATA TABLES.xlsx"),
  file.path("OMICS", "Proteomics", "SS-2338045_v4.1_Serum.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
))
raw_required <- file.path(project_root, raw_required_rel)

raw_status <- data.frame(
  path = raw_required_rel,
  exists = file.exists(raw_required),
  stringsAsFactors = FALSE
)

if (!all(raw_status$exists)) {
  record_issue(
    "01_prepare_data",
    "data-path issue",
    "Raw-data reconstruction cannot run because one or more raw inputs are missing."
  )
}

if (!requireNamespace("SomaDataIO", quietly = TRUE)) {
  record_issue(
    "01_prepare_data",
    "package/version change",
    "SomaDataIO is not installed; proteomics ADAT reconstruction cannot run in this R environment."
  )
}

rebuild_proteomics_inputs <- function(paths) {
  needed <- c(
    file.path(paths$processed, "prot_feature_v4.rds"),
    file.path(paths$processed, "prot_sample_metadata_merge_v4.rds")
  )
  if (all(file.exists(needed))) {
    return(FALSE)
  }
  if (!all(raw_status$exists) || !requireNamespace("SomaDataIO", quietly = TRUE)) {
    return(FALSE)
  }

  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
  })

  q_sas <- function(x, p) stats::quantile(x, probs = p, na.rm = TRUE, type = 2, names = FALSE)
  is_tag <- function(x, tag) haven::is_tagged_na(x) & haven::na_tag(x) == tag

  v1 <- haven::read_sas(file.path(paths$raw, "V1FEB24", "v1feb24.sas7bdat"))
  df_base <- v1 %>%
    transmute(
      ID,
      age_enroll = if_else(is.na(GIAGE1), median(GIAGE1, na.rm = TRUE), GIAGE1),
      race = case_when(
        is.na(GIRACE) ~ "Missing",
        GIRACE == 1 ~ "White",
        GIRACE == 2 ~ "African American",
        GIRACE == 3 ~ "Asian",
        TRUE ~ "Others"
      ),
      edu = case_when(
        is.na(GIEDUC) ~ "Missing",
        GIEDUC %in% 1:4 ~ "<= High School",
        GIEDUC %in% 5:6 ~ "College",
        GIEDUC %in% 7:8 ~ "Graduate School",
        TRUE ~ "Others"
      )
    )

  v4 <- haven::read_sas(file.path(paths$raw, "V4FEB24", "v4feb24.sas7bdat"))
  df_v4 <- v4 %>%
    transmute(
      ID,
      type = if_else(V4TYPE == "1", "Yes", "No, SAQ only"),
      site = case_when(
        is.na(SITE) ~ "Missing",
        SITE == "BI" ~ "Birmingham",
        SITE == "MN" ~ "Minneapolis",
        SITE == "PA" ~ "Palo Alto",
        SITE == "PI" ~ "Pittsburgh",
        SITE == "PO" ~ "Portland",
        TRUE ~ "San Diego"
      ),
      ol_health = case_when(
        is.na(QLHEALTH) ~ "Missing",
        QLHEALTH %in% 1:3 ~ "Good/Excellent",
        TRUE ~ "Very Poor/Poor/Fair"
      ),
      hgt = HWHGT,
      bmi = HWBMI,
      mstat = case_when(is.na(GIMSTAT) ~ "Missing", GIMSTAT == 1 ~ "Married", GIMSTAT == 2 ~ "Widowed", TRUE ~ "Others"),
      smoke = case_when(is.na(TURSMOKE) ~ "Missing", TURSMOKE == 0 ~ "Non-Smoker", TRUE ~ "Past or Current Smoker"),
      diab = case_when(is.na(MHDIAB) ~ "Missing", MHDIAB == 0 ~ "No", TRUE ~ "Yes"),
      hbp = case_when(is.na(MHBP) ~ "Missing", MHBP == 0 ~ "No", TRUE ~ "Yes"),
      cancer = case_when(is.na(MHCANCER) ~ "Missing", MHCANCER == 0 ~ "No", TRUE ~ "Yes"),
      tmm_score = TMMSCORE,
      gds = DPGDS15,
      pase = PASCORE,
      grip = GSGRPMAX,
      chair = NFCHAIR10,
      gait400m = NF4WLKSPD,
      gait6m = NFWLKSPD,
      unable_grip = is_tag(GSGRPMAX, "u") | is_tag(GSGRPMAX, "r"),
      unable_walk = (NFWLKNA1 == 2 & NFWLKNA2 == 2),
      wgt34change = HW34WPC,
      intend_wgt_loss = case_when(is.na(MHWTLOSS) ~ "Missing", MHWTLOSS == 0 ~ "No", TRUE ~ "Yes"),
      exhaust = case_when(is.na(DPENER) ~ "Missing", DPENER == 0 ~ "No", TRUE ~ "Yes")
    )

  cr <- haven::read_sas(file.path(paths$raw, "CR4AUG17", "CR4AUG17.SAS7BDAT")) %>%
    transmute(ID, cr_cmm = if_else(is.na(CRPERCMM), median(CRPERCMM, na.rm = TRUE), CRPERCMM))

  med <- haven::read_sas(file.path(paths$raw, "M4AUG16", "M4AUG16.SAS7BDAT")) %>%
    transmute(
      ID,
      abx = case_when(is.na(M1ANTIB2) ~ "Missing", M1ANTIB2 == 0 ~ "No", TRUE ~ "Yes"),
      total_meds = M1MEDSIN
    )

  ef <- haven::read_sas(file.path(paths$raw, "EFAUG24", "efaug24.sas7bdat")) %>%
    transmute(
      ID,
      status = case_when(
        is.na(EFSTATUS) ~ "Missing",
        EFSTATUS == 0 ~ "Active",
        EFSTATUS == 1 ~ "Deceased",
        EFSTATUS == 2 ~ "Terminated",
        TRUE ~ "Postcard Only"
      ),
      fu_yt_base = FUCYTIME,
      fu_yt_v4 = FUV4YT,
      fu_yt_v5 = FUV5YT
    )

  sample_metadata <- list(df_base, df_v4, cr, med, ef) %>%
    purrr::reduce(left_join, by = "ID") %>%
    mutate(
      age_v4 = age_enroll + fu_yt_base - fu_yt_v4,
      age_v5 = age_enroll + fu_yt_base - fu_yt_v5
    ) %>%
    select(
      ID, status, age_v4, age_v5, race, edu, site, ol_health, bmi, mstat,
      smoke, diab, hbp, cancer, tmm_score, gds, pase, total_meds,
      grip, chair, gait400m, gait6m, cr_cmm, abx
    )

  prot_flags <- readr::read_csv(
    file.path(paths$raw, "OMICS", "Proteomics", "sample_flags.csv"),
    show_col_types = FALSE
  )
  prot_metadata <- prot_flags %>%
    left_join(sample_metadata, by = "ID") %>%
    mutate(
      ATTEND_V5 = as.factor(ATTEND_V5),
      D3CR_MB_COHORT = as.factor(D3CR_MB_COHORT)
    ) %>%
    arrange(ID)

  adat <- SomaDataIO::read_adat(
    file.path(paths$raw, "OMICS", "Proteomics", "SS-2338045_v4.1_Serum.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat")
  )
  prot_feature <- as.data.frame(adat) %>%
    filter(SampleType == "Sample") %>%
    select(SampleId, matches("^seq\\.")) %>%
    rename(ID = SampleId) %>%
    arrange(ID)

  saveRDS(prot_feature, file.path(paths$processed, "prot_feature_v4.rds"))
  saveRDS(prot_metadata, file.path(paths$processed, "prot_sample_metadata_merge_v4.rds"))
  TRUE
}

rebuilt_proteomics <- rebuild_proteomics_inputs(paths)
if (rebuilt_proteomics) {
  record_issue(
    "01_prepare_data",
    "processed data regenerated",
    "Missing proteomics processed inputs were rebuilt from the raw ADAT and clinical metadata."
  )
}

input_status <- check_processed_inputs(paths)
missing_inputs <- input_status$input[!input_status$exists]
if (length(missing_inputs)) {
  record_issue(
    "01_prepare_data",
    "data-path issue",
    paste("Missing processed inputs:", paste(missing_inputs, collapse = ", "))
  )
} else {
  message("All expected processed inputs are present under data/processed_data.")
}
