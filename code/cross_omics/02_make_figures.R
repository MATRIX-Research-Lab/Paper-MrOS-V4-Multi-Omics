## 02_make_figures.R -- cross-omics figures for the blocks from 01_find_blocks.R.
##
##   Rscript code/cross_omics/02_make_figures.R      # or source() it from R / RStudio
##
## Draws the cross-correlation square only: no within-omics triangles, a tick mark
## on every row and column, block boxes numbered (1) and (2), and either feature
## names or "Number of ..." counts on the axes (the `names` field in JOBS).
##
## Paper PDFs, written under results/figures/:
##   Fig5_micro-metab_correlation_square.pdf    (taxon and metabolite names on both axes)
##   Supplementary_Fig3_micro_prot_square.pdf   (numeric positions)
##   Supplementary_Fig4_meta_prot_square.pdf    (numeric positions)
## PNG copies go to results/figures/cross_omics/. PNG needs the ragg package;
## PDF uses cairo_pdf when X11 is available and falls back to pdf() otherwise.

## Work inside code/cross_omics however the script was started: Rscript, source()
## in R, the Source button in RStudio, or Run / pasting the code into the console.
locate_script_dir <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)                       # Rscript
  f <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
  if (length(f)) return(dirname(normalizePath(f[[1]])))
  for (i in rev(seq_len(sys.nframe()))) {                        # source()
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(dirname(normalizePath(of)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&          # RStudio Run / console
      isTRUE(try(rstudioapi::isAvailable(), silent = TRUE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) NULL)
    if (!is.null(p) && nzchar(p)) return(dirname(normalizePath(p)))
  }
  NULL
}
script_dir <- locate_script_dir()
if (!is.null(script_dir)) setwd(script_dir)
if (!file.exists("common.R")) {
  stop("common.R was not found in the working directory:\n  ", getwd(),
       "\nSet the working directory to the folder holding this script, then run again:\n",
       '  setwd("<path to the repository>/code/cross_omics")')
}

source("common.R")
source("plot3in1_v1.R")        # still required: matlab_jet(), map_col(), unique_keep_last(), MATLAB_GREY
## source("plot3in1_v3.R")     # only needed by the commented-out 3-in-1 block at the bottom
source("plot_square_v4.R")

blocks_file <- file.path(tables_dir(), "cross_omics_blocks.rds")
if (!file.exists(blocks_file)) {
  stop("Missing ", blocks_file, ". Run 01_find_blocks.R first.")
}
blocks <- readRDS(blocks_file)
extra <- file.path(figures_dir(), "cross_omics")
dir.create(extra, recursive = TRUE, showWarnings = FALSE)

LEGEND   <- FALSE   # TRUE adds a block-size legend under the panels
NAME_PS  <- 14      # point size of feature names used as tick labels
ROW_FONT <- 1       # font for the row (vertical-axis) names: 1 plain, 3 italic (genus names)
COL_SRT  <- 45      # angle of the column (horizontal-axis) names: 90 vertical, 45 or 30 tilted
PAGE_W   <- 20      # inches, as in the 3-in-1 version

## One entry per figure.
##   rows / cols : omics on the vertical / horizontal axis (rows = first omics of the pair)
##   names       : "both", "rows", "cols" or "none" -- which axes show feature names
##                 instead of numeric positions
##   xlab / ylab : axis titles
##   step        : numeric tick-label step c(x, y), NULL = automatic (ignored on a named axis)
##   x_unit / y_unit : words for the optional legend
JOBS <- list(
  list(pair = "micro_metab", rows = "micro", cols = "metab", names = "both",
       xlab = "Metabolome", ylab = "Microbiome", step = NULL,
       x_unit = "metabolites", y_unit = "microbes",
       pdf = "Fig5_micro-metab_correlation_square.pdf"),
  list(pair = "micro_prot",  rows = "micro", cols = "prot",  names = "none",
       xlab = "Number of proteins", ylab = "Number of microbes", step = NULL,
       x_unit = "proteins", y_unit = "microbes",
       pdf = "Supplementary_Fig3_micro_prot_square.pdf"),
  list(pair = "metab_prot",  rows = "metab", cols = "prot",  names = "none",
       xlab = "Number of proteins", ylab = "Number of metabolites", step = NULL,
       x_unit = "proteins", y_unit = "metabolites",
       pdf = "Supplementary_Fig4_meta_prot_square.pdf"))

## Page height: 8.6 in for numeric axes; with feature names the margins grow to fit
## the longest name, and the page grows with them so the heatmap keeps its size.
page_size <- function(row_labels, col_labels) {
  leg <- if (LEGEND) 0.7 else 0
  if (is.null(row_labels) && is.null(col_labels)) return(c(PAGE_W, 8.6 + leg))
  mar_b <- if (is.null(col_labels)) 1.3 else col_label_depth_in(col_labels, NAME_PS, COL_SRT) + 0.25 + 0.55
  mar_l <- if (is.null(row_labels)) 1.5 else label_extent_in(row_labels, NAME_PS, ROW_FONT) + 0.25 + 0.55
  bw <- 0.47 * PAGE_W - mar_l - 0.3           # heatmap width inside one panel
  c(PAGE_W, 0.6 + bw / 1.2 + mar_b + 0.2 + leg)
}

for (j in JOBS) {
  S  <- load_pair(j$rows, j$cols, within = FALSE)   # only the cross matrices are needed
  bl <- blocks[[j$pair]]
  rows <- unique_keep_last(unlist(lapply(bl, `[[`, "X")))
  cols <- unique_keep_last(unlist(lapply(bl, `[[`, "Y")))
  rl <- if (j$names %in% c("both", "rows")) S$names_x[rows] else NULL
  cl <- if (j$names %in% c("both", "cols")) S$names_y[cols] else NULL
  fs <- page_size(rl, cl)
  for (f in c(file.path(figures_dir(), j$pdf), file.path(extra, paste0(j$pair, "_square.png")))) {
    open_square_device(f, figsize = fs)
    plot_square_v4(S$R_xyA, S$R_xyD, bl, step = j$step, bar = TRUE, legend = LEGEND,
                   xlabel_str = j$xlab, ylabel_str = j$ylab, x_unit = j$x_unit, y_unit = j$y_unit,
                   row_labels = rl, col_labels = cl, name_ps = NAME_PS, row_font = ROW_FONT,
                   col_srt = COL_SRT, figsize = fs)
    dev.off()
    message("wrote ", f)
  }
  rm(S); gc()
}


## ===========================================================================
## ORIGINAL 3-IN-1 FIGURE CODE -- COMMENTED OUT, NOTHING BELOW THIS LINE RUNS
## ===========================================================================
## This is the previous content of 02_make_figures.R, kept verbatim for
## reference. It drew each panel as a square plus the two within-omics
## correlation triangles, and wrote these files:
##   Fig5_micro-metab_correlation.pdf.pdf
##   Supplementary_Fig3_micro_prot.pdf
##   Supplementary_Fig4_meta_prot.pdf
## plus the single-panel v1 PNGs in results/figures/cross_omics/.
## To restore it: uncomment everything below, comment out the active code
## above, and re-enable source("plot3in1_v3.R").
## ---------------------------------------------------------------------------
# ## 02_make_figures.R -- 3-in-1 figures for the blocks from 01_find_blocks.R.
# ##
# ##   Rscript code/cross_omics/02_make_figures.R
# ##
# ## Paper PDFs, written under results/figures/ with the names already used there:
# ##   Fig5_micro-metab_correlation.pdf.pdf
# ##   Supplementary_Fig3_micro_prot.pdf
# ##   Supplementary_Fig4_meta_prot.pdf
# ## Extra panels (png, and the single-panel v1 plots) go to results/figures/cross_omics/.
# ## Needs the ragg package for PNG. PDF uses cairo_pdf when X11 is available.
#
# cmd <- commandArgs(trailingOnly = FALSE)
# file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
# if (length(file_arg)) setwd(dirname(normalizePath(file_arg[[1]])))
#
# source("common.R")
# source("plot3in1_v1.R")
# source("plot3in1_v3.R")
#
# blocks_file <- file.path(tables_dir(), "cross_omics_blocks.rds")
# if (!file.exists(blocks_file)) {
#   stop("Missing ", blocks_file, ". Run 01_find_blocks.R first.")
# }
# blocks <- readRDS(blocks_file)
# extra <- file.path(figures_dir(), "cross_omics")
# dir.create(extra, recursive = TRUE, showWarnings = FALSE)
#
# v1_png <- function(file, X, Y, XY, res, step) {
#   ragg::agg_png(file.path(extra, file), width = 5.6, height = 4.2, units = "in", res = 300)
#   par(mar = c(0, 0, 0, 0))
#   plot3in1_v1(X, Y, XY, res, line = TRUE, step = step, bar = TRUE, cbar_pos = "right")
#   dev.off()
# }
#
# ## pair, x, y, v1 step, v3 step, v3 x / y label, v3 label offset, paper pdf name
# JOBS <- list(
#   list("micro_metab", "micro", "metab", c(5, 5),   c(5, 5),   "Metabolite", "Microbiome", 2.2,
#        "Fig5_micro-metab_correlation.pdf.pdf"),
#   list("micro_prot",  "micro", "prot",  c(20, 5),  c(20, 5),  "Proteomics", "Microbiome", 1,
#        "Supplementary_Fig3_micro_prot.pdf"),
#   list("metab_prot",  "metab", "prot",  c(20, 10), c(10, 10), "Proteomics", "Metabolite", 2.2,
#        "Supplementary_Fig4_meta_prot.pdf"))
#
# for (j in JOBS) {
#   nm <- j[[1]]
#   S   <- load_pair(j[[2]], j[[3]])
#   res <- lapply(blocks[[nm]], function(b) list(b$X, b$Y))
#   v1_png(paste0(nm, "_v1_Rmax.png"),     S$R_yy,  S$R_xx,  t(S$R_xy),  res, j[[4]])
#   v1_png(paste0(nm, "_v1_P.png"),        S$R_yy,  S$R_xx,  t(S$P_xy),  res, j[[4]])
#   v1_png(paste0(nm, "_v1_Active.png"),   S$R_yyA, S$R_xxA, t(S$R_xyA), res, j[[4]])
#   v1_png(paste0(nm, "_v1_Deceased.png"), S$R_yyD, S$R_xxD, t(S$R_xyD), res, j[[4]])
#   open_v3_device(file.path(figures_dir(), j[[9]]), xlabel_offset = j[[8]])
#   plot3in1_v3(S$R_yyA, S$R_xxA, t(S$R_xyA), S$R_yyD, S$R_xxD, t(S$R_xyD), res, line = TRUE,
#               step = j[[5]], bar = TRUE, xlabel_str = j[[6]], ylabel_str = j[[7]],
#               xlabel_offset = j[[8]])
#   dev.off()
#   open_v3_device(file.path(extra, paste0(nm, "_v3.png")), xlabel_offset = j[[8]])
#   plot3in1_v3(S$R_yyA, S$R_xxA, t(S$R_xyA), S$R_yyD, S$R_xxD, t(S$R_xyD), res, line = TRUE,
#               step = j[[5]], bar = TRUE, xlabel_str = j[[6]], ylabel_str = j[[7]],
#               xlabel_offset = j[[8]])
#   dev.off()
#   rm(S); gc()
# }
#
