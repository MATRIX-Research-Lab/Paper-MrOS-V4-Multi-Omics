# Copy every generated figure and table into the submission folder and verify
# the copies are byte-identical, so "figures and tables" can never drift from
# what the code produces.
#
#   Rscript code/single_omics/sync_submission_assets.R "../manuscript/figures and tables"
#   Rscript code/single_omics/sync_submission_assets.R "../manuscript/figures and tables" --dry-run
#
# Two ordering rules matter here, because this script deletes files in a folder
# that is not the repository's own:
#
#   * an empty results/ is a failure, never an instruction to delete everything;
#   * stale artefacts are removed only after every copy has been verified, so a
#     failed sync cannot leave the destination stripped of what it came in with.

args <- commandArgs(trailingOnly = TRUE)
dry_run <- "--dry-run" %in% args
args <- args[args != "--dry-run"]
dest_arg <- if (length(args)) args[[1]] else file.path("..", "manuscript", "figures and tables")

# Resolve the destination against the caller's working directory *before*
# 00_setup.R runs: it calls setwd(project_root), which would otherwise silently
# retarget a relative path to a different folder than the one the user typed.
dest <- normalizePath(dest_arg, winslash = "/", mustWork = FALSE)

find_repository_root <- function() {
  d <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
  for (i in seq_len(20L)) {
    if (file.exists(file.path(d, "code", "single_omics", "00_setup.R"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate the repository root; run this from inside the checkout.")
}

source(file.path(find_repository_root(), "code", "single_omics", "00_setup.R"))

src <- c(list.files(paths$figures, full.names = TRUE),
         list.files(paths$tables, full.names = TRUE))
src <- src[!grepl("^\\.", basename(src))]

# Nothing to copy means the render has not run, or stopped early. Continuing
# would classify every file in the destination as stale and delete the whole
# submission folder.
if (!length(src)) {
  stop("No generated artefacts under results/figures or results/tables. ",
       "Run code/single_omics/run_all.R first; refusing to sync an empty results/ into ", dest)
}

dir.create(dest, recursive = TRUE, showWarnings = FALSE)

# Anything in the destination the pipeline does not produce. Directories are
# left alone: unlink() would not remove them anyway, and this script only owns
# the artefact files.
stale <- setdiff(list.files(dest), basename(src))
stale <- stale[!grepl("^\\.", stale)]
stale <- stale[!dir.exists(file.path(dest, stale))]

if (dry_run) {
  message(sprintf("[dry run] would copy %d artefacts into %s", length(src), dest))
  if (length(stale)) {
    message("[dry run] would remove: ", paste(stale, collapse = ", "))
  } else {
    message("[dry run] nothing to remove")
  }
  quit(save = "no")
}

invisible(file.copy(src, dest, overwrite = TRUE))

ok <- vapply(src, function(f) {
  d <- file.path(dest, basename(f))
  file.exists(d) && identical(tools::md5sum(f)[[1]], tools::md5sum(d)[[1]])
}, logical(1))

message(sprintf("%d/%d artefacts synced and verified identical", sum(ok), length(ok)))
if (!all(ok)) {
  stop("Mismatch after sync; nothing was removed from ", dest, ": ",
       paste(basename(src)[!ok], collapse = ", "))
}

# Only now: remove what the pipeline no longer produces, so stale hand-made
# artefacts cannot survive there unnoticed.
if (length(stale)) {
  message("Removing artefacts not produced by the pipeline: ",
          paste(stale, collapse = ", "))
  unlink(file.path(dest, stale))
}
