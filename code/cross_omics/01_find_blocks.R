## 01_find_blocks.R -- dense cross-omics blocks on the observed Active / Deceased split.
##
##   Rscript code/cross_omics/01_find_blocks.R
##
## About 20 seconds. Needs perm_workspace.RData in this directory (not committed).
## Writes cluster membership tables and results/tables/cross_omics_blocks.rds.
## 02_make_figures.R reads that rds file. The permute_*.R scripts do not.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
if (length(file_arg)) setwd(dirname(normalizePath(file_arg[[1]])))

source("common.R")
blocks <- list(); tabs <- list(); summ <- list()

## micro x metab, lambda 1.4
S  <- load_pair("micro", "metab", within = FALSE)
bl <- sort_result_matlab(S$R_xy, greedy_peeling_XY_all(S$P_xy, 1.4))
blocks$micro_metab <- bl
tabs$cross_omics_micro_metab_microbiome  <- member_table(bl, S$names_x, "X")
tabs$cross_omics_micro_metab_metabolites <- member_table(bl, S$names_y, "Y")
summ$micro_metab <- cbind(pair = "micro_metab", block_table(S, bl))

## micro x prot, lambda 1.3
S  <- load_pair("micro", "prot", within = FALSE)
bl <- sort_result_matlab(S$R_xy, greedy_peeling_XY_all(S$P_xy, 1.3))
blocks$micro_prot <- bl
tabs$cross_omics_micro_prot_microbiome <- member_table(bl, S$names_x, "X")
tabs$cross_omics_micro_prot_proteomics <- member_table(bl, S$names_y, "Y")
summ$micro_prot <- cbind(pair = "micro_prot", block_table(S, bl))

## metab x prot, lambda 1.5; block 4 split with lambda 1.4
S  <- load_pair("metab", "prot")
bl <- sort_result_matlab(S$R_xy, greedy_peeling_XY_all(S$P_xy, 1.5))
rx <- refine_signed(bl[[4]]$X, S$R_xx, 1.4)      # 40 -> 19 + 10 + 8
bl[[4]]$X <- rx$idx
bl[[4]]$Y <- refine_signed(bl[[4]]$Y, S$R_yy, 1.4)$idx   # 526 -> 358
blocks$metab_prot_all <- bl
blocks$metab_prot     <- bl[1:2]                 # figures use the first two blocks
tm <- member_table(bl, S$names_x, "X")
tm$ResultRow[tm$ResultRow == "4"] <- rep(paste0("4.", seq_along(rx$sizes)), rx$sizes)
tabs$cross_omics_metab_prot_metabolites <- tm
tabs$cross_omics_metab_prot_proteomics  <- member_table(bl, S$names_y, "Y")
summ$metab_prot <- cbind(pair = "metab_prot", block_table(S, bl))
rm(S); gc()

tab <- tables_dir()
summ <- do.call(rbind, summ)
saveRDS(blocks, file.path(tab, "cross_omics_blocks.rds"))
for (k in names(tabs)) write.csv(tabs[[k]], file.path(tab, paste0(k, ".csv")), row.names = FALSE)
write.csv(summ, file.path(tab, "cross_omics_block_summary.csv"), row.names = FALSE)
print(summ, digits = 3, row.names = FALSE)
