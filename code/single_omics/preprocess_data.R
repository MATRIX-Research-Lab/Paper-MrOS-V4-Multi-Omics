#!/usr/bin/env Rscript

# MrOS preprocessing and analytic-frame construction
#
# This is the only local preprocessing entry point.  It reads the restricted
# MrOS source files under data/raw_data/, constructs the shared clinical
# metadata, creates the lower-level processed omics objects, and writes the
# three complete modeling frames consumed by the cluster jobs.
#
# Run from the repository root:
#   Rscript code/single_omics/preprocess_data.R --overwrite
#
# The script is intentionally a plain R script rather than an interactive Rmd:
# every input path is resolved from the repository, exploratory print-only
# chunks are omitted, and the expected analytic-sample contracts are checked
# before any modeling frame is written.

`%||%` <- function(x, y) if (!is.null(x)) x else y

find_project_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  starts <- c(
    getwd(),
    if (length(file_arg)) dirname(file_arg[[1]]) else character(),
    tryCatch(dirname(sys.frame(1)$ofile), error = function(e) character())
  )
  starts <- unique(starts[nzchar(starts) & !is.na(starts)])
  for (start in starts) {
    d <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in seq_len(20L)) {
      if (file.exists(file.path(d, "code", "single_omics", "00_setup.R"))) {
        return(d)
      }
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Could not locate the repository root", call. = FALSE)
}

parse_args <- function(args) {
  out <- list(
    overwrite = FALSE,
    seed = as.integer(Sys.getenv("MROS_PREPROCESS_SEED", "123")),
    # Single-process is the portable default; parallel sockets are optional.
    ancom_workers = as.integer(Sys.getenv("MROS_ANCOM_WORKERS", "1"))
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("--overwrite", "-f")) {
      out$overwrite <- TRUE
      i <- i + 1L
    } else if (arg == "--seed" && i < length(args)) {
      out$seed <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (arg == "--ancom-workers" && i < length(args)) {
      out$ancom_workers <- as.integer(args[[i + 1L]])
      i <- i + 2L
    } else if (arg %in% c("--help", "-h")) {
      cat("Usage: Rscript code/single_omics/preprocess_data.R [--overwrite] [--seed N] [--ancom-workers N]\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", arg, call. = FALSE)
    }
  }
  if (is.na(out$seed) || out$seed < 1L) stop("--seed must be positive")
  if (is.na(out$ancom_workers) || out$ancom_workers < 1L) {
    stop("--ancom-workers must be positive")
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
project_root <- find_project_root()
setwd(project_root)
source(file.path(project_root, "code", "single_omics", "00_setup.R"))

preprocess_packages <- c(
  "dplyr", "tidyr", "purrr", "haven", "readr", "openxlsx", "readxl",
  "SomaDataIO", "phyloseq", "microbiome", "ANCOMBC"
)
assert_packages(preprocess_packages, "MrOS preprocessing")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

as_id <- function(x) trimws(as.character(x))

q_sas <- function(x, p) {
  stats::quantile(as.numeric(x), probs = p, na.rm = TRUE, type = 2, names = FALSE)
}

is_tag <- function(x, tag) {
  haven::is_tagged_na(x) & haven::na_tag(x) == tag
}

impute_median <- function(x, label) {
  x_num <- as.numeric(x)
  med <- stats::median(x_num, na.rm = TRUE)
  if (!is.finite(med)) stop(label, " has no finite value for median imputation")
  ifelse(is.na(x_num), med, x_num)
}

raw_rel <- c(
  file.path("V1FEB24", "v1feb24.sas7bdat"),
  file.path("V4FEB24", "v4feb24.sas7bdat"),
  file.path("CR4AUG17", "CR4AUG17.SAS7BDAT"),
  file.path("M4AUG16", "M4AUG16.SAS7BDAT"),
  file.path("B4JUN19", "b4jun19.sas7bdat"),
  file.path("FAAUG24", "faaug24.sas7bdat"),
  file.path("HQ4FEB18", "HQ4FEB18.SAS7BDAT"),
  file.path("EFAUG24", "efaug24.sas7bdat"),
  file.path("OMICS", "Microbiome", "DADA_Silva_Mros_data_release.rds"),
  file.path("OMICS", "Metabolomics", "sample_flags.csv"),
  file.path("OMICS", "Metabolomics", "CPMC-0101-23PHML+ DATA TABLES.xlsx"),
  file.path("OMICS", "Proteomics", "sample_flags.csv"),
  file.path("OMICS", "Proteomics", "SomaScan_7K_Annotated_Content.xlsx"),
  file.path(
    "OMICS", "Proteomics",
    "SS-2338045_v4.1_Serum.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
  )
)
raw_files <- file.path(paths$raw, raw_rel)
missing_raw <- raw_files[!file.exists(raw_files)]
if (length(missing_raw)) {
  stop(
    "Missing raw MrOS inputs:\n  ",
    paste(file.path("data", "raw_data", raw_rel[!file.exists(raw_files)]), collapse = "\n  ")
  )
}

output_files <- c(
  "mb_asv_data_v4.rds",
  "mbx_expr_data_v4.rds",
  "mbx_feature_metadata_v4.rds",
  "mbx_sample_metadata_v4.rds",
  "prot_feature_v4.rds",
  "prot_sample_metadata_merge_v4.rds",
  "df_status_micro_ancombc.rds",
  "df_status_metab.rds",
  "df_status_prot.rds"
)
output_paths <- file.path(paths$processed, output_files)
if (any(file.exists(output_paths)) && !isTRUE(args$overwrite)) {
  stop(
    "Processed outputs already exist. Use --overwrite only after confirming that a fresh rebuild is intended:\n  ",
    paste(output_paths[file.exists(output_paths)], collapse = "\n  ")
  )
}

covariates <- analysis_covariates
non_feature_columns <- analysis_non_feature_columns

# ---------------------------------------------------------------------------
# Shared clinical metadata
# ---------------------------------------------------------------------------
message("Reading clinical source files ...")

v1feb24 <- haven::read_sas(file.path(paths$raw, "V1FEB24", "v1feb24.sas7bdat"))
df_base <- v1feb24 %>%
  transmute(
    ID = as_id(ID),
    age_enroll = impute_median(GIAGE1, "GIAGE1"),
    race = case_when(
      is.na(GIRACE) ~ "Missing",
      as.numeric(GIRACE) == 1 ~ "White",
      as.numeric(GIRACE) == 2 ~ "African American",
      as.numeric(GIRACE) == 3 ~ "Asian",
      TRUE ~ "Others"
    ),
    edu = case_when(
      is.na(GIEDUC) ~ "Missing",
      as.numeric(GIEDUC) %in% 1:4 ~ "<= High School",
      as.numeric(GIEDUC) %in% 5:6 ~ "College",
      as.numeric(GIEDUC) %in% 7:8 ~ "Graduate School",
      TRUE ~ "Others"
    )
  )
assert_unique_ids(df_base, "ID", "baseline clinical metadata")

v4feb24 <- haven::read_sas(file.path(paths$raw, "V4FEB24", "v4feb24.sas7bdat"))
df_v4 <- v4feb24 %>%
  transmute(
    ID = as_id(ID),
    type = case_when(
      is.na(V4TYPE) ~ "Missing",
      as.character(V4TYPE) == "1" ~ "Yes",
      TRUE ~ "No, SAQ only"
    ),
    site = case_when(
      is.na(SITE) ~ "Missing",
      as.character(SITE) == "BI" ~ "Birmingham",
      as.character(SITE) == "MN" ~ "Minneapolis",
      as.character(SITE) == "PA" ~ "Palo Alto",
      as.character(SITE) == "PI" ~ "Pittsburgh",
      as.character(SITE) == "PO" ~ "Portland",
      TRUE ~ "San Diego"
    ),
    ol_health = case_when(
      is.na(QLHEALTH) ~ "Missing",
      as.numeric(QLHEALTH) %in% 1:3 ~ "Good/Excellent",
      TRUE ~ "Very Poor/Poor/Fair"
    ),
    hgt = as.numeric(HWHGT),
    bmi = as.numeric(HWBMI),
    mstat = case_when(
      is.na(GIMSTAT) ~ "Missing",
      as.numeric(GIMSTAT) == 1 ~ "Married",
      as.numeric(GIMSTAT) == 2 ~ "Widowed",
      TRUE ~ "Others"
    ),
    smoke = case_when(
      is.na(TURSMOKE) ~ "Missing",
      as.numeric(TURSMOKE) == 0 ~ "Non-Smoker",
      TRUE ~ "Past or Current Smoker"
    ),
    diab = case_when(
      is.na(MHDIAB) ~ "Missing",
      as.numeric(MHDIAB) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    hbp = case_when(
      is.na(MHBP) ~ "Missing",
      as.numeric(MHBP) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    cancer = case_when(
      is.na(MHCANCER) ~ "Missing",
      as.numeric(MHCANCER) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    intend_wgt_loss = case_when(
      is.na(MHWTLOSS) ~ "Missing",
      as.numeric(MHWTLOSS) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    exhaust = case_when(
      is.na(DPENER) ~ "Missing",
      as.numeric(DPENER) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    tmm_score = as.numeric(TMMSCORE),
    gds = as.numeric(DPGDS15),
    pase = as.numeric(PASCORE),
    unable_grip = is_tag(GSGRPMAX, "u") | is_tag(GSGRPMAX, "r"),
    grip = as.numeric(GSGRPMAX),
    chair = as.numeric(NFCHAIR10),
    gait400m = as.numeric(NF4WLKSPD),
    unable_walk = as.numeric(NFWLKNA1) == 2 & as.numeric(NFWLKNA2) == 2,
    gait6m = as.numeric(NFWLKSPD),
    wgt34change = as.numeric(HW34WPC)
  )
assert_unique_ids(df_v4, "ID", "visit-4 clinical metadata")

df_v4_cr <- haven::read_sas(file.path(paths$raw, "CR4AUG17", "CR4AUG17.SAS7BDAT")) %>%
  transmute(ID = as_id(ID), cr_cmm = impute_median(CRPERCMM, "CRPERCMM"))
df_v4_med <- haven::read_sas(file.path(paths$raw, "M4AUG16", "M4AUG16.SAS7BDAT")) %>%
  transmute(
    ID = as_id(ID),
    abx = case_when(
      is.na(M1ANTIB2) ~ "Missing",
      as.numeric(M1ANTIB2) == 0 ~ "No",
      TRUE ~ "Yes"
    ),
    total_meds = as.numeric(M1MEDSIN)
  )
df_v4_dax <- haven::read_sas(file.path(paths$raw, "B4JUN19", "b4jun19.sas7bdat")) %>%
  transmute(
    ID = as_id(ID),
    b4thd = as.numeric(B4THD),
    b4fnd = as.numeric(B4FND),
    b4lsd = as.numeric(B4LSD)
  )
df_v4_fx <- haven::read_sas(file.path(paths$raw, "FAAUG24", "faaug24.sas7bdat")) %>%
  transmute(
    ID = as_id(ID),
    faprev4 = case_when(
      is.na(FAPREV4) ~ "Missing",
      as.numeric(FAPREV4) == 0 ~ "No",
      TRUE ~ "Yes"
    )
  )
df_v4_hq <- haven::read_sas(file.path(paths$raw, "HQ4FEB18", "HQ4FEB18.SAS7BDAT")) %>%
  transmute(
    ID = as_id(ID),
    hqdrfefl = as.numeric(HQDRFEFL),
    hqdtfefl = as.numeric(HQDTFEFL),
    hqptfefl = as.numeric(HQPTFEFL)
  )
df_v4_status <- haven::read_sas(file.path(paths$raw, "EFAUG24", "efaug24.sas7bdat")) %>%
  transmute(
    ID = as_id(ID),
    status = case_when(
      is.na(EFSTATUS) ~ "Missing",
      as.numeric(EFSTATUS) == 0 ~ "Active",
      as.numeric(EFSTATUS) == 1 ~ "Deceased",
      as.numeric(EFSTATUS) == 2 ~ "Terminated",
      TRUE ~ "Postcard Only"
    ),
    fu_yt_base = as.numeric(FUCYTIME),
    fu_yt_v4 = as.numeric(FUV4YT),
    fu_yt_v5 = as.numeric(FUV5YT)
  )

clinical_for_join <- list(df_base, df_v4, df_v4_cr, df_v4_med,
                          df_v4_dax, df_v4_fx, df_v4_hq, df_v4_status)
for (i in seq_along(clinical_for_join)) {
  assert_unique_ids(clinical_for_join[[i]], "ID", paste0("clinical input ", i))
}

# Frailty variables are retained for the metadata and matched-cohort contract.
# The original Rmd used the 55th percentile as the second BMI cutpoint and
# compared the categorical weight-loss variable to numeric 0/1.  Both were
# transcription errors; fixing them does not alter any of the 14 prespecified
# mortality-model covariates.
cut_dat <- df_v4 %>% filter(type == "Yes")
eligible_grip <- cut_dat %>% filter(!is.na(grip) | unable_grip)
bmi25 <- q_sas(df_v4$bmi, 0.25)
bmi50 <- q_sas(df_v4$bmi, 0.50)
bmi75 <- q_sas(df_v4$bmi, 0.75)
grip_cuts <- eligible_grip %>%
  mutate(
    bmiqrt = case_when(
      !is.na(bmi) & bmi <= bmi25 ~ 1,
      !is.na(bmi) & bmi <= bmi50 ~ 2,
      !is.na(bmi) & bmi <= bmi75 ~ 3,
      !is.na(bmi) & bmi > bmi75 ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  group_by(bmiqrt) %>%
  summarize(grip20 = q_sas(grip, 0.20), .groups = "drop")
eligible_walk <- cut_dat %>% filter(!is.na(gait6m) | unable_walk)
hgt50 <- q_sas(df_v4$hgt, 0.50)
wlk_cuts <- eligible_walk %>%
  mutate(hgtmed = case_when(
    !is.na(hgt) & hgt < hgt50 ~ 1,
    !is.na(hgt) & hgt >= hgt50 ~ 2,
    TRUE ~ NA_real_
  )) %>%
  group_by(hgtmed) %>%
  summarize(gait6m20 = q_sas(gait6m, 0.20), .groups = "drop")
gait6m_lm <- wlk_cuts$gait6m20[wlk_cuts$hgtmed == 1]
gait6m_gem <- wlk_cuts$gait6m20[wlk_cuts$hgtmed == 2]
pase20 <- q_sas(df_v4$pase, 0.20)

df_v4_frail <- df_v4 %>%
  mutate(
    fr_sh4 = case_when(
      !is.na(wgt34change) & wgt34change <= -5 & intend_wgt_loss == "No" ~ 1,
      !is.na(wgt34change) & wgt34change <= -5 & intend_wgt_loss == "Yes" ~ 0,
      !is.na(wgt34change) & wgt34change > -5 ~ 0,
      TRUE ~ NA_real_
    ),
    fr_ex4 = case_when(exhaust == "No" ~ 1, exhaust == "Yes" ~ 0, TRUE ~ NA_real_),
    fr_ac4 = case_when(
      !is.na(pase) & pase <= pase20 ~ 1,
      !is.na(pase) & pase > pase20 ~ 0,
      TRUE ~ NA_real_
    ),
    fr_sl4 = case_when(
      !is.na(hgt) & hgt < hgt50 &
        (unable_walk | (!is.na(gait6m) & gait6m <= gait6m_lm)) ~ 1,
      !is.na(hgt) & hgt < hgt50 & !is.na(gait6m) & gait6m > gait6m_lm ~ 0,
      !is.na(hgt) & hgt >= hgt50 &
        (unable_walk | (!is.na(gait6m) & gait6m <= gait6m_gem)) ~ 1,
      !is.na(hgt) & hgt >= hgt50 & !is.na(gait6m) & gait6m > gait6m_gem ~ 0,
      TRUE ~ NA_real_
    ),
    bmiqrt = case_when(
      !is.na(bmi) & bmi <= bmi25 ~ 1,
      !is.na(bmi) & bmi <= bmi50 ~ 2,
      !is.na(bmi) & bmi <= bmi75 ~ 3,
      !is.na(bmi) & bmi > bmi75 ~ 4,
      TRUE ~ NA_real_
    )
  ) %>%
  left_join(grip_cuts, by = "bmiqrt") %>%
  mutate(
    fr_wk4 = case_when(
      !is.na(grip20) & (unable_grip | (!is.na(grip) & grip <= grip20)) ~ 1,
      !is.na(grip20) & !is.na(grip) & grip > grip20 ~ 0,
      TRUE ~ NA_real_
    )
  ) %>%
  rowwise() %>%
  mutate(
    fr_cnt4 = sum(!is.na(c(fr_sh4, fr_wk4, fr_ex4, fr_ac4, fr_sl4))),
    fr_chsb4 = {
      x <- c(fr_sh4, fr_wk4, fr_ex4, fr_ac4, fr_sl4)
      if (type != "Yes") NA_real_ else if (fr_cnt4 <= 2) NA_real_ else
        if (fr_cnt4 %in% c(3, 4)) 5 * mean(x, na.rm = TRUE) else
          if (fr_cnt4 == 5) sum(x) else NA_real_
    },
    fr_chsn4 = if_else(!is.na(fr_chsb4), round(fr_chsb4, 0), NA_real_),
    fr_chs4 = case_when(
      is.na(fr_chsn4) ~ "Missing",
      fr_chsn4 == 0 ~ "Robust",
      fr_chsn4 %in% c(1, 2) ~ "Intermediate",
      fr_chsn4 >= 3 ~ "Frail"
    )
  ) %>%
  ungroup() %>%
  select(ID, type, fr_sh4, fr_wk4, fr_ex4, fr_ac4, fr_sl4,
         fr_cnt4, fr_chsb4, fr_chsn4, fr_chs4)

df_sample_metadata_merge <- list(
  df_base,
  df_v4 %>% select(-type, -hgt, -unable_grip, -unable_walk,
                   -wgt34change, -intend_wgt_loss, -exhaust),
  df_v4_cr, df_v4_med, df_v4_dax, df_v4_fx, df_v4_hq, df_v4_status,
  df_v4_frail %>% select(ID, fr_chs4)
) %>%
  reduce(left_join, by = "ID") %>%
  mutate(
    age_v4 = age_enroll + fu_yt_base - fu_yt_v4,
    age_v5 = age_enroll + fu_yt_base - fu_yt_v5
  ) %>%
  select(-age_enroll, -fu_yt_base, -fu_yt_v4, -fu_yt_v5) %>%
  arrange(ID)
assert_unique_ids(df_sample_metadata_merge, "ID", "combined clinical metadata")

# ---------------------------------------------------------------------------
# Microbiome: lower-level object and ANCOM-BC2 bias-corrected frame
# ---------------------------------------------------------------------------
message("Reading and processing microbiome data ...")
mb_v4 <- readRDS(file.path(paths$raw, "OMICS", "Microbiome", "DADA_Silva_Mros_data_release.rds"))
df_v4_asv <- data.frame(microbiome::abundances(mb_v4), check.names = FALSE)
mb_sample_metadata_v4 <- microbiome::meta(mb_v4) %>%
  transmute(ID = as_id(SampleID), batchno = batchno, batch = batch)
mb_feature_metadata_v4 <- as.matrix(phyloseq::tax_table(mb_v4))
assert_unique_ids(mb_sample_metadata_v4, "ID", "microbiome sample metadata")

mb_sample_metadata_merge_v4 <- mb_sample_metadata_v4 %>%
  left_join(df_sample_metadata_merge, by = "ID") %>%
  select(-batchno)
assert_unique_ids(mb_sample_metadata_merge_v4, "ID", "merged microbiome metadata")

ASV <- phyloseq::otu_table(df_v4_asv, taxa_are_rows = TRUE)
META <- phyloseq::sample_data(mb_sample_metadata_merge_v4)
phyloseq::sample_names(META) <- mb_sample_metadata_merge_v4$ID
TAX <- phyloseq::tax_table(mb_feature_metadata_v4)
asv_data <- phyloseq::phyloseq(ASV, TAX, META)
asv_data <- phyloseq::subset_samples(asv_data, !is.na(abx) & abx == "No")
asv_status <- phyloseq::subset_samples(asv_data, status %in% c("Active", "Deceased"))

clean_taxonomic_name <- function(name) {
  parts <- strsplit(trimws(as.character(name)), "_ ", fixed = TRUE)[[1]]
  genus_idx <- grep("^g__", parts)
  if (length(genus_idx)) {
    genus <- sub("^g__", "", parts[[genus_idx[[1]]]])
    if (nzchar(genus) && !is.na(genus)) return(genus)
  }
  for (level in c("f__", "o__", "c__", "p__", "k__")) {
    idx <- grep(paste0("^", level), parts)
    if (length(idx)) {
      value <- sub(paste0("^", level), "", parts[[idx[[1]]]])
      if (nzchar(value) && !is.na(value)) return(paste0(level, value))
    }
  }
  as.character(name)
}

set.seed(args$seed)
ancombc_result <- ANCOMBC::ancombc2(
  data = asv_status,
  tax_level = "Genus",
  fix_formula = "batch + site",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo_sens = TRUE,
  prv_cut = 0.10,
  lib_cut = 1000,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = args$ancom_workers,
  verbose = TRUE
)
if (is.null(ancombc_result$bias_correct_log_table)) {
  stop("ANCOM-BC2 did not return bias_correct_log_table")
}

df_genus_bc <- exp(as.matrix(ancombc_result$bias_correct_log_table))
df_genus_bc[is.na(df_genus_bc)] <- 0
df_genus_bc <- as.data.frame(t(df_genus_bc), check.names = FALSE)
df_genus_bc <- tibble::rownames_to_column(df_genus_bc, "ID")
raw_taxa_names <- setdiff(names(df_genus_bc), "ID")
clean_taxa_names <- vapply(raw_taxa_names, clean_taxonomic_name, character(1))
clean_taxa_names[is.na(clean_taxa_names) | !nzchar(clean_taxa_names)] <-
  raw_taxa_names[is.na(clean_taxa_names) | !nzchar(clean_taxa_names)]
names(df_genus_bc) <- c("ID", make.unique(clean_taxa_names, sep = "."))
df_genus_bc$ID <- as_id(df_genus_bc$ID)

micro_meta_status <- as(phyloseq::sample_data(asv_status), "data.frame")
micro_meta_status$ID <- as_id(micro_meta_status$ID)
df_status_micro_ancombc <- micro_meta_status %>%
  left_join(df_genus_bc, by = "ID") %>%
  filter(status %in% c("Active", "Deceased")) %>%
  mutate(status = factor(status, levels = c("Deceased", "Active")))

rm(ancombc_result, asv_status, df_genus_bc, micro_meta_status)

# ---------------------------------------------------------------------------
# Metabolomics: matrix, metadata, and final modeling frame
# ---------------------------------------------------------------------------
message("Reading and processing metabolomics data ...")
metab_dir <- file.path(paths$raw, "OMICS", "Metabolomics")
mbx_sample_flag_v4 <- readr::read_csv(file.path(metab_dir, "sample_flags.csv"), show_col_types = FALSE) %>%
  mutate(ID = as_id(ID))
mbx_sample_metadata_v4 <- readxl::read_excel(
  file.path(metab_dir, "CPMC-0101-23PHML+ DATA TABLES.xlsx"), sheet = 3
) %>%
  transmute(ID = as_id(CLIENT_IDENTIFIER), sample_name = as.character(PARENT_SAMPLE_NAME)) %>%
  left_join(mbx_sample_flag_v4, by = "ID")
assert_unique_ids(mbx_sample_metadata_v4, "ID", "metabolomics sample metadata")

df_v4_mbx <- openxlsx::read.xlsx(
  file.path(metab_dir, "CPMC-0101-23PHML+ DATA TABLES.xlsx"), sheet = 6
)
names(df_v4_mbx)[1] <- "sample_name"
df_v4_mbx <- df_v4_mbx %>%
  mutate(sample_name = as.character(sample_name)) %>%
  left_join(mbx_sample_metadata_v4 %>% select(ID, sample_name), by = "sample_name") %>%
  select(ID, everything(), -sample_name) %>%
  mutate(ID = as_id(ID)) %>%
  arrange(ID)
assert_unique_ids(df_v4_mbx, "ID", "metabolomics feature matrix")

mbx_feature_metadata_v4 <- openxlsx::read.xlsx(
  file.path(metab_dir, "CPMC-0101-23PHML+ DATA TABLES.xlsx"), sheet = 2
)
if (!all(c("CHEM_ID", "PLOT_NAME") %in% names(mbx_feature_metadata_v4))) {
  stop("Metabolomics feature metadata must contain CHEM_ID and PLOT_NAME")
}
metab_map <- setNames(
  make.unique(ifelse(
    is.na(mbx_feature_metadata_v4$PLOT_NAME) | !nzchar(as.character(mbx_feature_metadata_v4$PLOT_NAME)),
    paste0("CHEM_", mbx_feature_metadata_v4$CHEM_ID),
    as.character(mbx_feature_metadata_v4$PLOT_NAME)
  )),
  as.character(mbx_feature_metadata_v4$CHEM_ID)
)
old_names <- names(df_v4_mbx)
new_names <- vapply(old_names, function(x) {
  if (x %in% names(metab_map)) unname(metab_map[[x]]) else x
}, character(1))
names(df_v4_mbx) <- new_names

mbx_sample_metadata_merge_v4 <- mbx_sample_metadata_v4 %>%
  left_join(df_sample_metadata_merge, by = "ID") %>%
  select(-sample_name) %>%
  mutate(
    ATTEND_V5 = as.factor(ATTEND_V5),
    D3CR_MB_COHORT = as.factor(D3CR_MB_COHORT)
  ) %>%
  arrange(ID) %>%
  filter(abx == "No")
assert_unique_ids(mbx_sample_metadata_merge_v4, "ID", "merged metabolomics metadata")

df_status_metab <- mbx_sample_metadata_merge_v4 %>%
  filter(status %in% c("Active", "Deceased")) %>%
  left_join(df_v4_mbx, by = "ID") %>%
  drop_na(status) %>%
  mutate(status = factor(status, levels = c("Deceased", "Active")))

# ---------------------------------------------------------------------------
# Proteomics: ADAT, annotation mapping, metadata, and final modeling frame
# ---------------------------------------------------------------------------
message("Reading and processing proteomics data ...")
prot_dir <- file.path(paths$raw, "OMICS", "Proteomics")
prot_sample_flag_v4 <- readr::read_csv(file.path(prot_dir, "sample_flags.csv"), show_col_types = FALSE) %>%
  mutate(ID = as_id(ID))
assert_unique_ids(prot_sample_flag_v4, "ID", "proteomics sample flags")

prot_sample_metadata_merge_v4 <- prot_sample_flag_v4 %>%
  left_join(df_sample_metadata_merge, by = "ID") %>%
  mutate(
    ATTEND_V5 = as.factor(ATTEND_V5),
    D3CR_MB_COHORT = as.factor(D3CR_MB_COHORT)
  ) %>%
  arrange(ID)
assert_unique_ids(prot_sample_metadata_merge_v4, "ID", "merged proteomics metadata")

adat <- SomaDataIO::read_adat(file.path(
  prot_dir,
  "SS-2338045_v4.1_Serum.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
))
prot_feature_v4 <- as.data.frame(adat) %>%
  filter(SampleType == "Sample") %>%
  select(SampleId, matches("^seq[.]")) %>%
  rename(ID = SampleId) %>%
  mutate(ID = as_id(ID)) %>%
  arrange(ID)
assert_unique_ids(prot_feature_v4, "ID", "proteomics feature matrix")

# The archived preprocessing specification removed the first ADAT sample row
# (`df_v4_prot[-1, ]`).  In the current release that row is BI0840.  Encode the
# same data-freeze rule by ID so a change in ADAT row order cannot remove a
# different participant without stopping the rebuild.
proteomics_excluded_ids <- "BI0840"
if (!all(proteomics_excluded_ids %in% prot_feature_v4$ID)) {
  stop("Expected historical proteomics exclusion is absent: ",
       paste(setdiff(proteomics_excluded_ids, prot_feature_v4$ID), collapse = ", "))
}
prot_feature_v4 <- prot_feature_v4 %>%
  filter(!ID %in% proteomics_excluded_ids)

# The current annotation workbook contains a short legal/preamble section, so
# locate the SeqId header instead of relying on a version-specific `skip` value.
annotation_raw <- readxl::read_excel(
  file.path(prot_dir, "SomaScan_7K_Annotated_Content.xlsx"),
  sheet = "Annotations", col_names = FALSE, .name_repair = "minimal"
)
header_row <- which(trimws(as.character(annotation_raw[[1]])) == "SeqId")[1]
if (is.na(header_row)) stop("Could not find SeqId header in SomaScan annotation workbook")
annotation_names <- as.character(unlist(annotation_raw[header_row, ], use.names = FALSE))
annotation_names[is.na(annotation_names) | !nzchar(annotation_names)] <-
  paste0("annotation_", which(is.na(annotation_names) | !nzchar(annotation_names)))
annotation_names <- make.unique(annotation_names)
prot_annotation <- annotation_raw[(header_row + 1L):nrow(annotation_raw), , drop = FALSE]
names(prot_annotation) <- annotation_names
seq_col <- match("SeqId", names(prot_annotation))
target_col <- match("Target Full Name", names(prot_annotation))
if (is.na(target_col)) target_col <- match("Target Name", names(prot_annotation))
if (is.na(seq_col) || is.na(target_col)) {
  stop("SomaScan annotation workbook must contain SeqId and a target-name column")
}
annotation_seq <- as.character(prot_annotation[[seq_col]])
annotation_target <- as.character(prot_annotation[[target_col]])
annotation_target[is.na(annotation_target) | !nzchar(trimws(annotation_target))] <-
  paste0("SeqId_", annotation_seq[is.na(annotation_target) | !nzchar(trimws(annotation_target))])
annotation_target <- make.unique(trimws(annotation_target))
prot_map <- setNames(annotation_target, paste0("seq.", gsub("-", ".", annotation_seq)))
prot_map <- prot_map[!is.na(names(prot_map)) & nzchar(names(prot_map))]
rename_spec <- setNames(intersect(names(prot_feature_v4), names(prot_map)),
                        unname(prot_map[intersect(names(prot_feature_v4), names(prot_map))]))
if (length(rename_spec)) {
  prot_feature_v4 <- dplyr::rename(prot_feature_v4, !!!rename_spec)
}

df_status_prot <- prot_feature_v4 %>%
  left_join(
    prot_sample_metadata_merge_v4 %>%
      filter(status %in% c("Active", "Deceased"), abx == "No"),
    by = "ID"
  ) %>%
  filter(!is.na(status)) %>%
  mutate(status = factor(status, levels = c("Deceased", "Active")))

# ---------------------------------------------------------------------------
# Contracts, writes, and provenance
# ---------------------------------------------------------------------------
expected <- analysis_contract

check_model_frame <- function(df, layer) {
  if (!is.data.frame(df)) stop(layer, ": expected a data.frame")
  df <- as.data.frame(df)
  df$ID <- as_id(df$ID)
  assert_unique_ids(df, "ID", paste0(layer, " modeling frame"))
  feature_cols <- setdiff(names(df), unique(c(covariates, non_feature_columns)))
  validate_analytic_frame(
    df, id_col = "ID", outcome = "status", cov_cols = covariates,
    feature_cols = feature_cols, label = paste0("preprocessing / ", layer)
  )
  n <- nrow(df)
  deaths <- sum(as.character(df$status) == "Deceased")
  exp <- expected[[layer]]
  if (n != exp[["n"]] || deaths != exp[["deaths"]] || length(feature_cols) != exp[["features"]]) {
    stop(
      sprintf(
        "%s contract failed: observed n=%d, deaths=%d, features=%d; expected n=%d, deaths=%d, features=%d",
        layer, n, deaths, length(feature_cols), exp[["n"]], exp[["deaths"]], exp[["features"]]
      )
    )
  }
  message(sprintf("%s: n=%d, deaths=%d, features=%d", layer, n, deaths, length(feature_cols)))
  invisible(df)
}

df_status_micro_ancombc <- check_model_frame(df_status_micro_ancombc, "microbiome")
df_status_metab <- check_model_frame(df_status_metab, "metabolomics")
df_status_prot <- check_model_frame(df_status_prot, "proteomics")

objects_to_write <- list(
  mb_asv_data_v4 = asv_data,
  mbx_expr_data_v4 = df_v4_mbx,
  mbx_feature_metadata_v4 = mbx_feature_metadata_v4,
  mbx_sample_metadata_v4 = mbx_sample_metadata_merge_v4,
  prot_feature_v4 = prot_feature_v4,
  prot_sample_metadata_merge_v4 = prot_sample_metadata_merge_v4,
  df_status_micro_ancombc = df_status_micro_ancombc,
  df_status_metab = df_status_metab,
  df_status_prot = df_status_prot
)
for (nm in names(objects_to_write)) {
  target <- file.path(paths$processed, paste0(nm, ".rds"))
  atomic_save_rds(objects_to_write[[nm]], target)
  message("Wrote ", file.path("data", "processed_data", basename(target)))
}

dir.create(paths$manifests, recursive = TRUE, showWarnings = FALSE)
raw_manifest <- data.frame(
  path = file.path("data", "raw_data", raw_rel),
  bytes = vapply(raw_files, function(f) as.numeric(file.info(f)$size), numeric(1)),
  md5 = unname(as.character(tools::md5sum(raw_files))),
  stringsAsFactors = FALSE
)
utils::write.csv(raw_manifest, file.path(paths$manifests, "preprocessing_raw_manifest.csv"), row.names = FALSE)
preprocess_code_files <- file.path(
  project_root, "code", "single_omics",
  c("00_setup.R", "preprocess_data.R", "analysis_provenance.R",
    "helpers.R", "feature_selection.R")
)
preprocess_code_manifest <- data.frame(
  path = file.path("code", "single_omics", basename(preprocess_code_files)),
  md5 = unname(as.character(tools::md5sum(preprocess_code_files))),
  stringsAsFactors = FALSE
)
utils::write.csv(preprocess_code_manifest,
                 file.path(paths$manifests, "preprocessing_code_manifest.csv"),
                 row.names = FALSE)
writeLines(capture.output(utils::sessionInfo()), file.path(paths$manifests, "preprocessing_session_info.txt"))

message("Preprocessing complete. Raw-input manifest and session information were written under results/manifests/.")
