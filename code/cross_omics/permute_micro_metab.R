## permute_micro_metab.R -- micro x metab
## Observed data + 500 permutations, one row per block.
## Does not read the blocks written by 01_find_blocks.R.
##
## Usage:
##   Rscript code/cross_omics/permute_micro_metab.R
##   NODE=0 NPERM=50 Rscript code/cross_omics/permute_micro_metab.R
##
## Output: results/tables/cross_omics_perm_micro_metab_node<NN>.csv
## perm_id = 0 is the observed run, produced by NODE 0 only.
## Place perm_workspace.RData in this directory first. It is not committed.

cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", cmd[grepl("^--file=", cmd)])
if (length(file_arg)) setwd(dirname(normalizePath(file_arg[[1]])))

## Method:
##   1) Per-cell Fisher z two-sample test,  H0: rho_D = rho_A
##   2) Fisher combining:  Psi_diff = -2 * sum(ln p)
##   3) Chernoff bound:  ln P(Psi >= t) <= -N (a - 1 - ln a),  a = Psi / (2N)
##   Inference comes from the permutation. The chi-square / Chernoff step is a
##   scale normalisation only, not evidence of significance.

source("greedy.R")
load("perm_workspace.RData")

NODE  <- as.integer(Sys.getenv("NODE",  "0"))
NPERM <- as.integer(Sys.getenv("NPERM", "500"))

CFG <- list(l = c(1.4, 1.5, 1.5), xs = c(0, 15), ys = c(0, 15))
CAP <- 15

DF   <- list(micro = df_status_micro, metab = df_status_metab, prot = df_status_prot)
COLS <- list(micro = micro_cols,      metab = metab_cols,      prot = prot_cols)

ids  <- sort(Reduce(intersect, lapply(DF, `[[`, "ID")))
feat <- function(d) as.matrix(DF[[d]][match(ids, DF[[d]]$ID), COLS[[d]], drop = FALSE])

X    <- feat("micro")
Y    <- feat("metab")
dec0 <- DF[["micro"]]$status[match(ids, DF[["micro"]]$ID)] == "Deceased"

m <- ncol(X); n <- ncol(Y)

stopifnot(nrow(X) == nrow(Y), length(dec0) == nrow(X),
          sum(dec0) == 169, sum(!dec0) == 163)

cat(sprintf("micro x metab | %d x %d | R_full %.0f MB | node %d | %d perms\n\n",
            m, n, (m + n)^2 * 8 / 1e6, NODE, NPERM))


## --- Correlations and p-values ----------------------------------------

cp <- function(A, Bm = NULL) {
  C <- if (is.null(Bm)) cor(A) else cor(A, Bm)
  C[is.na(C)] <- 0
  list(cor = pmin(pmax(C, -0.9999), 0.9999), k = nrow(A))
}

logp <- function(C, k) log(2) + pt(-abs(C) * sqrt((k - 2) / (1 - C^2)),
                                   df = k - 2, log.p = TRUE)

cmb <- function(A, D) { k <- abs(D) >= abs(A); A[k] <- D[k]; A }   # ties -> Deceased


## --- One full search --------------------------------------------------

search_blocks <- function(dec) {
  xyA <- cp(X[!dec, ], Y[!dec, ]); xyD <- cp(X[dec, ], Y[dec, ])
  R_xy <- cmb(xyA$cor, xyD$cor)
  R_xx <- cmb(cp(X[!dec, ])$cor, cp(X[dec, ])$cor)
  R_yy <- cmb(cp(Y[!dec, ])$cor, cp(Y[dec, ])$cor)

  P_xy <- pmax(-logp(xyA$cor, xyA$k), -logp(xyD$cor, xyD$k)) / log(10)
  P_xy[!is.finite(P_xy) | P_xy > CAP] <- CAP

  R_full <- matrix(0, m + n, m + n)
  R_full[1:m, 1:m]                 <- R_xx
  R_full[1:m, (m+1):(m+n)]         <- R_xy
  R_full[(m+1):(m+n), 1:m]         <- t(R_xy)
  R_full[(m+1):(m+n), (m+1):(m+n)] <- R_yy

  P_full <- matrix(0, m + n, m + n)
  P_full[1:m, (m+1):(m+n)] <- P_xy
  P_full[(m+1):(m+n), 1:m] <- t(P_xy)

  bl <- greedy_3step(P_full, R_full, m, n, CFG$l[1], CFG$l[2], CFG$l[3],
                     CFG$xs, CFG$ys, sort_xy = R_xy)

  rm(R_full, P_full, R_xx, R_yy, P_xy); gc(verbose = FALSE)

  list(blocks = bl, rD = xyD$cor, rA = xyA$cor, kD = xyD$k, kA = xyA$k)
}


## --- Chernoff bound ---------------------------------------------------
## For a <= 1 the bound degenerates to P <= 1, so record 0

chern <- function(Psi, N, df_per_cell = 2) {
  a <- Psi / (df_per_cell * N)
  ifelse(a > 1, -(df_per_cell / 2) * N * (a - 1 - log(a)), 0) / log(10)
}


## --- Block-level statistics -------------------------------------------
## Observed and permuted runs must go through this same function; any
## divergence silently corrupts the p-value

summarize <- function(S, perm_id) {
  if (length(S$blocks) == 0L) return(NULL)

  SE <- sqrt(1 / (S$kD - 3) + 1 / (S$kA - 3))

  do.call(rbind, lapply(seq_along(S$blocks), function(b) {
    I <- S$blocks[[b]]$X; J <- S$blocks[[b]]$Y
    stopifnot(all(I >= 1 & I <= m), all(J >= 1 & J <= n))

    rD <- S$rD[I, J, drop = FALSE]; rA <- S$rA[I, J, drop = FALSE]
    zD <- atanh(rD);                zA <- atanh(rA)

    Zc <- (zD - zA) / SE
    lP <- log(2) + pnorm(-abs(Zc), log.p = TRUE)      # per-cell difference p, natural log

    N        <- length(I) * length(J)
    Psi_diff <- -2 * sum(lP)
    sumZ2    <- sum(Zc^2)

    lD <- logp(rD, S$kD); lA <- logp(rA, S$kA)        # per-group Fisher, for comparison

    data.frame(
      perm_id = perm_id, block_id = b,
      p = length(I), q = length(J), N = N,
      n_valid = sum(is.finite(Zc)),

      Psi_diff = Psi_diff,
      a_diff   = Psi_diff / (2 * N),
      bound_diff = chern(Psi_diff, N, 2),

      sumZ2    = sumZ2,
      a_Z2     = sumZ2 / N,
      bound_Z2 = chern(sumZ2, N, 1),

      max_absZ = max(abs(Zc)), mean_absZ = mean(abs(Zc)),
      n_pos = sum(Zc > 0), n_neg = sum(Zc < 0),

      Psi_D = -2 * sum(lD), Psi_A = -2 * sum(lA),

      mean_z_D = mean(zD), mean_z_A = mean(zA),
      r_D = tanh(mean(zD)), r_A = tanh(mean(zA)),
      delta = mean(zD) - mean(zA),
      delta_absr = mean(abs(zD - zA)),
      min_absr = min(abs(tanh(mean(zD))), abs(tanh(mean(zA)))),
      max_absr = max(abs(tanh(mean(zD))), abs(tanh(mean(zA))))
    )
  }))
}


## --- Driver -----------------------------------------------------------

fn <- file.path(normalizePath(file.path(getwd(), "..", "..")), "results", "tables",
                sprintf("cross_omics_perm_micro_metab_node%02d.csv", NODE))
dir.create(dirname(fn), recursive = TRUE, showWarnings = FALSE)

wr <- function(x) {
  if (is.null(x)) return(invisible())
  write.table(x, fn, sep = ",", append = file.exists(fn),
              col.names = !file.exists(fn), row.names = FALSE)
}

done <- if (file.exists(fn)) unique(read.csv(fn)$perm_id) else integer(0)

## Observed run is done once, by NODE 0 only
if (NODE == 0L && !(0L %in% done)) {
  t0 <- Sys.time()
  o  <- summarize(search_blocks(dec0), 0L)
  wr(o)
  cat(sprintf("observed | %.1f s | %d blocks\n",
              difftime(Sys.time(), t0, units = "secs"),
              if (is.null(o)) 0L else nrow(o)))
  if (!is.null(o))
    print(o[, c("block_id","p","q","N","a_diff","bound_diff","r_D","r_A","delta")],
          digits = 4, row.names = FALSE)
  flush.console()
}

for (i in seq_len(NPERM)) {
  pid <- NODE * NPERM + i                 # globally unique; sets the seed
  if (pid %in% done) next
  set.seed(pid)
  t0 <- Sys.time()
  wr(summarize(search_blocks(sample(dec0)), pid))
  cat(sprintf("perm %4d | %.1f s\n", pid,
              difftime(Sys.time(), t0, units = "secs")))
  flush.console()
}

cat(sprintf("\ndone: %s\n", fn))
