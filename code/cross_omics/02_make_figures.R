## 02_make_figures.R -- 3-in-1 figures for the blocks from 01_find_blocks.R.
##
##   Rscript code/cross_omics/02_make_figures.R
##
## Paper PDFs, written under results/figures/ with the names already used there:
##   Fig5_micro-metab_correlation.pdf.pdf
##   Supplementary_Fig3_micro_prot.pdf
##   Supplementary_Fig4_meta_prot.pdf
## Extra panels (png, and the single-panel v1 plots) go to results/figures/cross_omics/.
## Needs the ragg package for PNG. PDF uses cairo_pdf when X11 is available.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
if (length(file_arg)) setwd(dirname(normalizePath(file_arg[[1]])))

source("common.R")
source("plot3in1_v1.R")
source("plot3in1_v3.R")

blocks_file <- file.path(tables_dir(), "cross_omics_blocks.rds")
if (!file.exists(blocks_file)) {
  stop("Missing ", blocks_file, ". Run 01_find_blocks.R first.")
}
blocks <- readRDS(blocks_file)
extra <- file.path(figures_dir(), "cross_omics")
dir.create(extra, recursive = TRUE, showWarnings = FALSE)

v1_png <- function(file, X, Y, XY, res, step) {
  ragg::agg_png(file.path(extra, file), width = 5.6, height = 4.2, units = "in", res = 300)
  par(mar = c(0, 0, 0, 0))
  plot3in1_v1(X, Y, XY, res, line = TRUE, step = step, bar = TRUE, cbar_pos = "right")
  dev.off()
}

## pair, x, y, v1 step, v3 step, v3 x / y label, v3 label offset, paper pdf name
JOBS <- list(
  list("micro_metab", "micro", "metab", c(5, 5),   c(5, 5),   "Metabolite", "Microbiome", 2.2,
       "Fig5_micro-metab_correlation.pdf.pdf"),
  list("micro_prot",  "micro", "prot",  c(20, 5),  c(20, 5),  "Proteomics", "Microbiome", 1,
       "Supplementary_Fig3_micro_prot.pdf"),
  list("metab_prot",  "metab", "prot",  c(20, 10), c(10, 10), "Proteomics", "Metabolite", 2.2,
       "Supplementary_Fig4_meta_prot.pdf"))

for (j in JOBS) {
  nm <- j[[1]]
  S   <- load_pair(j[[2]], j[[3]])
  res <- lapply(blocks[[nm]], function(b) list(b$X, b$Y))
  v1_png(paste0(nm, "_v1_Rmax.png"),     S$R_yy,  S$R_xx,  t(S$R_xy),  res, j[[4]])
  v1_png(paste0(nm, "_v1_P.png"),        S$R_yy,  S$R_xx,  t(S$P_xy),  res, j[[4]])
  v1_png(paste0(nm, "_v1_Active.png"),   S$R_yyA, S$R_xxA, t(S$R_xyA), res, j[[4]])
  v1_png(paste0(nm, "_v1_Deceased.png"), S$R_yyD, S$R_xxD, t(S$R_xyD), res, j[[4]])
  open_v3_device(file.path(figures_dir(), j[[9]]), xlabel_offset = j[[8]])
  plot3in1_v3(S$R_yyA, S$R_xxA, t(S$R_xyA), S$R_yyD, S$R_xxD, t(S$R_xyD), res, line = TRUE,
              step = j[[5]], bar = TRUE, xlabel_str = j[[6]], ylabel_str = j[[7]],
              xlabel_offset = j[[8]])
  dev.off()
  open_v3_device(file.path(extra, paste0(nm, "_v3.png")), xlabel_offset = j[[8]])
  plot3in1_v3(S$R_yyA, S$R_xxA, t(S$R_xyA), S$R_yyD, S$R_xxD, t(S$R_xyD), res, line = TRUE,
              step = j[[5]], bar = TRUE, xlabel_str = j[[6]], ylabel_str = j[[7]],
              xlabel_offset = j[[8]])
  dev.off()
  rm(S); gc()
}
