#!/usr/bin/env Rscript
#
# Install the CRAN-only packages at the versions recorded in
# environment/cran-packages.csv.
#
# conda rebuilds everything in conda-linux-64-r-4.4.3.lock exactly. These few
# packages are not available from conda-forge or bioconda, so they are pinned
# by version here instead. Current sources are tried first, then the CRAN
# Archive, which is where a version goes once it is superseded.

# Cap BLAS/OpenMP threads for the R CMD INSTALL subprocesses this script
# spawns. Conda's OpenBLAS starts one thread per core -- 128 on a Zaratan node
# -- and login nodes cap RLIMIT_NPROC at 256, so installs die with
# "blas_thread_init: pthread_create failed". The parent process has already
# loaded BLAS and is unaffected, but children inherit these.
for (v in c("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")) {
  if (!nzchar(Sys.getenv(v))) do.call(Sys.setenv, stats::setNames(list("1"), v))
}

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grepl("^--file=", args)])
script_dir <- if (length(script_path)) {
  dirname(normalizePath(script_path[[1]]))
} else {
  file.path(getwd(), "environment")
}

csv_path <- file.path(script_dir, "cran-packages.csv")
if (!file.exists(csv_path)) {
  stop("Missing ", csv_path)
}

pinned <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
if (!nrow(pinned)) {
  message("No CRAN-only packages recorded; nothing to do.")
  quit(save = "no", status = 0)
}

# Install only what the cluster stages actually require.
#
# cran-packages.csv is derived from what happened to be installed when the
# environment was captured, so a record made on a machine that also carried the
# local preprocessing stack lists ANCOMBC, CVXR, gmp and friends. Those are
# preprocessing-only; installing them here would rebuild a toolchain the
# cluster never loads, and a failure in any of them would abort a restore that
# was otherwise complete. render_packages in 00_setup.R is the authoritative
# list, so intersect against it.
project_root <- dirname(script_dir)
setup_file <- file.path(project_root, "code", "single_omics", "00_setup.R")
required <- tryCatch({
  env <- new.env(parent = globalenv())
  sys.source(setup_file, envir = env)
  get("render_packages", envir = env)
}, error = function(e) {
  message("NOTE: could not read render_packages from 00_setup.R (",
          conditionMessage(e), "); installing every recorded package.")
  NULL
})

if (!is.null(required)) {
  skipped <- setdiff(pinned$Package, required)
  pinned <- pinned[pinned$Package %in% required, , drop = FALSE]
  if (length(skipped)) {
    message("skipping ", length(skipped),
            " recorded package(s) the cluster stages do not require:")
    message("  ", paste(skipped, collapse = ", "))
  }
  if (!nrow(pinned)) {
    message("Nothing required from CRAN; conda provides everything.")
    quit(save = "no", status = 0)
  }
}

# Not every package outside the conda lock comes from CRAN. Bioconductor
# packages are not on CRAN at all, so a CRAN-only restore fails on them. Both
# repository families are configured, and the tarball fallbacks below try each.
BIOC_VERSION <- "3.20"

options(
  repos = c(
    CRAN     = "https://cloud.r-project.org",
    BioCsoft = sprintf("https://bioconductor.org/packages/%s/bioc", BIOC_VERSION),
    BioCann  = sprintf("https://bioconductor.org/packages/%s/data/annotation", BIOC_VERSION),
    BioCexp  = sprintf("https://bioconductor.org/packages/%s/data/experiment", BIOC_VERSION)
  ),
  timeout = 600,
  Ncpus = 1L
)

# Installing a pinned tarball by URL uses repos = NULL, which switches off
# dependency resolution entirely: R installs that one file and nothing else. So
# resolve and install the package's missing dependencies from the configured
# repositories first. Without this, a package is attempted before the packages
# it needs -- cran-packages.csv is ordered alphabetically, not topologically --
# and fails with "dependencies ... are not available".
ensure_dependencies <- function(package) {
  db <- tryCatch(utils::available.packages(), error = function(e) NULL)
  if (is.null(db) || !package %in% rownames(db)) return(invisible(FALSE))
  deps <- tools::package_dependencies(
    package, db = db, which = c("Depends", "Imports", "LinkingTo"),
    recursive = TRUE
  )[[package]]
  if (!length(deps)) return(invisible(TRUE))
  installed <- rownames(utils::installed.packages())
  missing <- setdiff(deps, installed)
  if (!length(missing)) return(invisible(TRUE))
  message("  installing ", length(missing), " missing dependencies of ",
          package, " ...")
  utils::install.packages(missing)
  # install.packages() only warns on a failed dependency, and the package that
  # needs it then fails to load for a reason the log has scrolled past. Name
  # them here instead.
  still_missing <- setdiff(missing, rownames(utils::installed.packages()))
  if (length(still_missing)) {
    stop(package, ": these dependencies failed to install: ",
         paste(still_missing, collapse = ", "),
         ".\nScroll up for the first compiler or download error. If the log ",
         "shows 'blas_thread_init: pthread_create failed', set ",
         "OPENBLAS_NUM_THREADS=1 (and OMP_NUM_THREADS=1, MKL_NUM_THREADS=1) ",
         "in the shell and re-run.")
  }
  invisible(TRUE)
}

source_urls <- function(package, version) {
  tarball <- sprintf("%s_%s.tar.gz", package, version)
  cran <- "https://cloud.r-project.org/src/contrib"
  bioc <- sprintf("https://bioconductor.org/packages/%s/bioc/src/contrib",
                  BIOC_VERSION)
  bioc_annot <- sprintf("https://bioconductor.org/packages/%s/data/annotation/src/contrib",
                        BIOC_VERSION)
  bioc_exp <- sprintf("https://bioconductor.org/packages/%s/data/experiment/src/contrib",
                      BIOC_VERSION)
  c(
    file.path(cran, tarball),
    file.path(cran, "Archive", package, tarball),
    file.path(bioc, tarball),
    file.path(bioc, "Archive", package, tarball),
    file.path(bioc_annot, tarball),
    file.path(bioc_exp, tarball)
  )
}

# Read the installed version from DESCRIPTION without loading the package.
#
# requireNamespace() would be the obvious check, but it loads the namespace and
# therefore every dependency, so a package that installed perfectly well is
# reported as failed when one of its imports is broken. kernelshap imports
# doFuture: with doFuture missing, requireNamespace("kernelshap") is FALSE even
# though kernelshap itself is on disk at the right version.
installed_version <- function(package) {
  tryCatch(as.character(utils::packageVersion(package)),
           error = function(e) NA_character_)
}

install_pinned <- function(package, version) {
  if (identical(installed_version(package), version)) {
    message("already installed at recorded version: ", package, " ", version)
    return(invisible(TRUE))
  }
  ensure_dependencies(package)
  for (url in source_urls(package, version)) {
    # Warnings are not failures here: install.packages() warns about all sorts
    # of benign things. Catching them as failures made the loop fall through to
    # the next URL after a successful install.
    tryCatch(
      utils::install.packages(url, repos = NULL, type = "source"),
      error = function(e) message("  ", conditionMessage(e))
    )
    if (identical(installed_version(package), version)) {
      message("installed: ", package, " ", version)
      return(invisible(TRUE))
    }
  }
  stop("Could not install ", package, " ", version,
       " from CRAN or Bioconductor ", BIOC_VERSION, ".\n",
       "If this package is not actually needed by the cluster stages, the ",
       "environment record is stale: remove it from the environment and re-run ",
       "environment/capture_environment.sh so cran-packages.csv lists only what ",
       "the analysis uses. Otherwise download the tarball and install it with ",
       "R CMD INSTALL.")
}

for (i in seq_len(nrow(pinned))) {
  install_pinned(pinned$Package[[i]], pinned$Version[[i]])
}

message("all recorded CRAN packages installed")
