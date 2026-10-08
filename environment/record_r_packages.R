#!/usr/bin/env Rscript
#
# Verify the analysis environment against the contracts in 00_setup.R, then
# record every installed R package and the full session information.
#
# Invoked by environment/capture_environment.sh, which writes the conda lock
# first. This script reads that lock to separate conda-provided packages from
# the ones installed from CRAN, so the CRAN set never has to be maintained by
# hand.

find_repository_root <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  starts <- c(getwd(), if (length(file_arg)) dirname(file_arg[[1]]) else NULL)
  for (s in unique(starts[nzchar(starts)])) {
    d <- normalizePath(s, winslash = "/", mustWork = FALSE)
    for (i in seq_len(20L)) {
      if (file.exists(file.path(d, "code", "single_omics", "00_setup.R"))) return(d)
      parent <- dirname(d)
      if (identical(parent, d)) break
      d <- parent
    }
  }
  stop("Could not locate the repository root")
}

project_root <- find_repository_root()
setwd(project_root)
out_dir <- file.path(project_root, "environment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Verify -------------------------------------------------------------
# Sourcing 00_setup.R is what makes this authoritative: the package lists and
# the API check live with the analysis code, so this cannot drift from what the
# pipeline actually requires.
source(file.path(project_root, "code", "single_omics", "00_setup.R"))

assert_packages(render_packages, "render stage")
assert_packages(primary_packages, "primary stage")
assert_packages(matched_packages, "matched-cohort stage")
assert_modeling_api()
message("contracts satisfied: render, primary, matched, modeling API")

# ---- 1b. Library integrity -------------------------------------------------
# The contracts only load the packages they name. A damaged *transitive*
# dependency -- digest, Rcpp, codetools -- passes every contract and then kills
# a job hours later with "there is no package called 'digest'".
#
# conda is no help: its metadata records the package as installed, so
# `conda install` reports "All requested packages already installed" and does
# nothing, whatever the state of the files on disk. Check the library directly.
# Driven by installed.packages() rather than by listing directories: an R
# library also holds non-package directories (translations, and build
# leftovers), which a raw directory scan reports as damaged.
inventory_all <- utils::installed.packages()
broken <- character()
for (i in seq_len(nrow(inventory_all))) {
  pkg <- rownames(inventory_all)[[i]]
  pkg_dir <- file.path(inventory_all[i, "LibPath"], pkg)
  if (!file.exists(file.path(pkg_dir, "Meta", "package.rds")) ||
      !file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    broken <- c(broken, pkg)
  }
}
if (length(broken)) {
  stop(
    "Incomplete package installations (", length(broken), "): ",
    paste(utils::head(broken, 20), collapse = ", "),
    if (length(broken) > 20) ", ..." else "", "\n",
    "conda reports these as installed and will not repair them. Rebuild at a ",
    "fresh prefix:\n",
    "  bash environment/restore_environment.sh <new-prefix>"
  )
}
message("library integrity: ", nrow(inventory_all), " packages complete")

# ---- 2. Inventory ----------------------------------------------------------
installed <- utils::installed.packages(
  fields = c("Package", "Version", "Priority", "Built")
)
installed <- as.data.frame(installed, stringsAsFactors = FALSE)
installed <- installed[order(tolower(installed$Package)), , drop = FALSE]

# Packages shipped with R itself are recorded by the R version, not separately.
is_base <- !is.na(installed$Priority) &
  installed$Priority %in% c("base", "recommended")
inventory <- installed[!is_base, c("Package", "Version", "Built"), drop = FALSE]

# Which of these did conda provide? Test membership against the explicit lock
# rather than a hand-maintained list. conda names R packages r-<lowercase> and
# Bioconductor packages bioconductor-<lowercase>.
lock_path <- file.path(out_dir, "conda-linux-64-r-4.4.3.lock")
if (!file.exists(lock_path)) {
  stop("Missing ", lock_path, ". Run environment/capture_environment.sh, ",
       "which writes the conda lock before calling this script.")
}
lock_lines <- readLines(lock_path, warn = FALSE)

from_conda <- vapply(inventory$Package, function(pkg) {
  patterns <- paste0("/", c("r-", "bioconductor-"), tolower(pkg), "-")
  any(vapply(patterns, function(p) any(grepl(p, lock_lines, fixed = TRUE)),
             logical(1)))
}, logical(1))

inventory$Source <- ifelse(from_conda, "conda", "CRAN")

utils::write.csv(
  inventory,
  file.path(out_dir, "r-packages.csv"),
  row.names = FALSE
)

# The CRAN-installed subset is what the conda lock cannot rebuild, so it gets
# its own file and its own restore step.
#
# Restricted to packages the cluster stages actually require. An environment
# that also carries the local preprocessing stack would otherwise record
# ANCOMBC, CVXR, gmp and their dependencies here, and a later restore would try
# to rebuild a toolchain the cluster never loads. Everything installed is still
# recorded in full in r-packages.csv; this file is specifically the restore
# input.
cran_only <- inventory[inventory$Source == "CRAN" &
                         inventory$Package %in% render_packages,
                       c("Package", "Version"), drop = FALSE]
utils::write.csv(
  cran_only,
  file.path(out_dir, "cran-packages.csv"),
  row.names = FALSE
)

# ---- 3. Session information ------------------------------------------------
session_path <- file.path(out_dir, "r-session-info.txt")
con <- file(session_path, open = "wt")
sink(con, type = "output")
# sessionInfo() already reports the BLAS and LAPACK libraries on Linux.
print(utils::sessionInfo())
cat("\n--- Relevant environment variables ---\n")
vars <- c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
          "R_LIBS_USER", "R_LIBS_SITE", "CONDA_PREFIX")
for (v in vars) cat(sprintf("%-22s %s\n", v, Sys.getenv(v, unset = "<unset>")))
sink(type = "output")
close(con)

message("recorded ", nrow(inventory), " packages (",
        sum(from_conda), " conda, ", sum(!from_conda), " CRAN)")
