## common.R -- data and helpers for the cross-omics block search.
## Entry scripts set the working directory to this folder before sourcing this file.
## Input: perm_workspace.RData in this directory (not committed).
##   Objects: df_status_{micro,metab,prot}, {micro,metab,prot}_cols

source("greedy.R")
if (!file.exists("perm_workspace.RData")) {
  stop("perm_workspace.RData is not in code/cross_omics/. ",
       "Build it locally from the three df_status_*.rds files written by ",
       "code/single_omics/preprocess_data.R. Do not commit it.")
}
load("perm_workspace.RData")

DF   <- list(micro = df_status_micro, metab = df_status_metab, prot = df_status_prot)
COLS <- list(micro = micro_cols, metab = metab_cols, prot = prot_cols)
ids  <- sort(Reduce(intersect, lapply(DF, `[[`, "ID")))          # 332 samples in all 3 omics
dec  <- DF$micro$status[match(ids, DF$micro$ID)] == "Deceased"    # 169 Deceased / 163 Active
feat <- function(d, k) as.matrix(DF[[d]][match(ids, DF[[d]]$ID), COLS[[d]]])[k, ]

## Pearson r (NaN -> 0) and -log10 p of cor.test (t, df = n - 2)
cor0 <- function(A, B = A) { r <- suppressWarnings(cor(A, B)); r[is.na(r)] <- 0; r }
nlp  <- function(r, n) { a <- pmin(abs(r), 1); -log10(2 * pt(-a * sqrt((n - 2) / (1 - a^2)), n - 2)) }

## |D| >= |A| -> Deceased value (ties go to Deceased)
cmb <- function(A, D) { k <- abs(D) >= abs(A); A[k] <- D[k]; A }

## All matrices for one pairing (P capped at 15)
load_pair <- function(x, y, within = TRUE) {
  S <- list(names_x = COLS[[x]], names_y = COLS[[y]])
  for (g in c("A", "D")) {
    k <- if (g == "D") dec else !dec
    X <- feat(x, k); Y <- feat(y, k)
    S[[paste0("R_xy", g)]] <- cor0(X, Y)
    S[[paste0("P_xy", g)]] <- nlp(S[[paste0("R_xy", g)]], sum(k))
    if (within) { S[[paste0("R_xx", g)]] <- cor0(X); S[[paste0("R_yy", g)]] <- cor0(Y) }
  }
  S$P_xy <- pmin(pmax(S$P_xyA, S$P_xyD), 15)
  S$R_xy <- cmb(S$R_xyA, S$R_xyD)
  if (within) { S$R_xx <- cmb(S$R_xxA, S$R_xxD); S$R_yy <- cmb(S$R_yyA, S$R_yyD) }
  S
}

## Sort features inside each block by |marginal|, then order blocks by mean R
sort_result_matlab <- function(Wp, bl) {
  bl <- sort_result(Wp, bl, symmetric = FALSE)
  bl[order(sapply(bl, function(b) mean(Wp[b$X, b$Y])), decreasing = TRUE)]
}

## metab x prot block-4 split: peel the signed within-omic matrix
refine_signed <- function(idx, R, lam) {
  sub  <- R[idx, idx]; sub <- (sub + t(sub)) / 2
  comm <- sort_result(sub, greedy_peeling_X_all(sub, lam), symmetric = TRUE)
  list(idx = idx[unlist(lapply(comm, `[[`, "X"))], sizes = sapply(comm, function(b) length(b$X)))
}

block_table <- function(S, bl) do.call(rbind, lapply(seq_along(bl), function(i) {
  I <- bl[[i]]$X; J <- bl[[i]]$Y
  data.frame(block = i, len1 = length(I), len2 = length(J), max_mean_r = mean(S$R_xy[I, J]),
             active_mean_r = mean(S$R_xyA[I, J]), deceased_mean_r = mean(S$R_xyD[I, J]),
             mean_negLog10P = mean(S$P_xy[I, J]))
}))

member_table <- function(bl, nm, w) do.call(rbind, lapply(seq_along(bl), function(i)
  data.frame(Element = nm[bl[[i]][[w]]], ResultRow = as.character(i))))

repo_root <- function() normalizePath(file.path(getwd(), "..", ".."))
tables_dir <- function() {
  d <- file.path(repo_root(), "results", "tables")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
figures_dir <- function() {
  d <- file.path(repo_root(), "results", "figures")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}
