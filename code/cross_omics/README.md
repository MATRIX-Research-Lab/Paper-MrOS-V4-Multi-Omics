# Cross-omics blocks

Finds dense microbiome, metabolite, and protein correlation blocks and draws the
3-in-1 figures. Permutation scripts estimate whether a block this strong appears
after the Active / Deceased labels are shuffled.

Everything starts from `perm_workspace.RData` in this directory. That file is
not in the repository. Build it locally from the three modeling frames written
by `code/single_omics/preprocess_data.R`:

```r
source("code/single_omics/00_setup.R")
df_status_micro <- readRDS(file.path(paths$processed, "df_status_micro_ancombc.rds"))
df_status_metab <- readRDS(file.path(paths$processed, "df_status_metab.rds"))
df_status_prot  <- readRDS(file.path(paths$processed, "df_status_prot.rds"))
drop <- unique(c(analysis_covariates, analysis_non_feature_columns))
micro_cols <- setdiff(names(df_status_micro), drop)
metab_cols <- setdiff(names(df_status_metab), drop)
prot_cols  <- setdiff(names(df_status_prot), drop)
stopifnot(length(micro_cols) == 147L,
          length(metab_cols) == 1574L,
          length(prot_cols) == 7596L)
save(df_status_micro, df_status_metab, df_status_prot,
     micro_cols, metab_cols, prot_cols,
     file = "code/cross_omics/perm_workspace.RData")
```

The microbiome frame is stored as `df_status_micro`. Its rds file is named
`df_status_micro_ancombc.rds`.

## Run

From the repository root:

```bash
Rscript code/cross_omics/01_find_blocks.R    # cluster tables, about 20 seconds
Rscript code/cross_omics/02_make_figures.R   # Fig 5 and Supplementary Figs 3–4
```

`02` reads `results/tables/cross_omics_blocks.rds` from step 01. That rds is
gitignored. The three `permute_*.R` scripts do not read it. Each one shuffles
the labels 500 times for one pair:

```bash
Rscript code/cross_omics/permute_micro_metab.R
Rscript code/cross_omics/permute_micro_prot.R
Rscript code/cross_omics/permute_metab_prot.R
```

On a current laptop the three permutation jobs together take about three hours.
Split one job with `NODE=0 NPERM=50`. `perm_id = 0` is the observed split and
is written only by node 0. A finished `perm_id` is skipped on a rerun.

PNG output needs the R package `ragg`.

## Files

| File | What |
|---|---|
| `greedy.R` | peeling algorithm; sourced, not run |
| `common.R` | loads the workspace and builds the correlation matrices |
| `plot3in1_v1.R`, `plot3in1_v3.R` | plotting functions |
| `01_find_blocks.R` | observed blocks and cluster-name tables |
| `02_make_figures.R` | paper figures; run after 01 |
| `permute_micro_metab.R` | 500 label shuffles, microbiome × metabolomics |
| `permute_micro_prot.R` | 500 label shuffles, microbiome × proteomics |
| `permute_metab_prot.R` | 500 label shuffles, metabolomics × proteomics |

`01` and the permutation scripts both call `greedy.R`, but they do not call
each other. The permutation scripts re-find blocks on every shuffle and record
a Fisher combined statistic. The p-value is the permutation tail. The Chernoff
column in the output is a scale, not the p-value.
