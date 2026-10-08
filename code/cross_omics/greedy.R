## greedy.R -- 3-step greedy peeling, R port of the MATLAB reference.
##
## All indices are 1-based. The algorithm is fully deterministic: no RNG.
##
## Differences from the Python port, both matching MATLAB:
##   - greedy_peeling_XY_all tests "i < max_round", so it returns at most
##     max_round - 1 blocks (Python allowed max_round).
##   - greedy_3step's Y refinement returns Y-axis indices, not combined ones.
##     (MATLAB line 51 has an unfixed +m offset bug here.)
##
## Note: MATLAB runs these in single precision via single(Wp); this uses
## double. Results can differ where two marginals are within ~1e-3 of each
## other, since the peeling order is set by which.min.


## ---------------------------------------------------------------------
## Step 1: cross-omics (X x Y) peeling
## ---------------------------------------------------------------------

#' Peel one dense submatrix from a non-negative m x n matrix.
#' @return list(density, X, Y)
greedy_peeling_XY_one <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  m <- nrow(Wp); n <- ncol(Wp)

  C <- colSums(Wp)   # column marginals over remaining rows
  R <- rowSums(Wp)   # row marginals over remaining columns

  rem_X <- logical(m); rem_Y <- logical(n)
  best_rem_X <- logical(m); best_rem_Y <- logical(n)
  len_X <- m; len_Y <- n
  W_density <- 0

  for (it in seq_len(m + n)) {
    Cc <- C; Cc[rem_Y] <- Inf
    Rr <- R; Rr[rem_X] <- Inf
    IndT <- which.min(Cc); dT <- Cc[IndT]
    IndS <- which.min(Rr); dS <- Rr[IndS]
    if (!is.finite(dT) || !is.finite(dS)) break

    cratio <- if (len_Y > 0) len_X / len_Y else Inf

    if (cratio * dS <= dT) {          # remove a row (X feature)
      rem_X[IndS] <- TRUE
      C <- C - Wp[IndS, ]; C[rem_Y] <- 0
      Wp_sum <- sum(C); len_X <- len_X - 1L
    } else {                          # remove a column (Y feature)
      rem_Y[IndT] <- TRUE
      R <- R - Wp[, IndT]; R[rem_X] <- 0
      Wp_sum <- sum(R); len_Y <- len_Y - 1L
    }
    if (len_X <= 0 || len_Y <= 0) break

    score <- Wp_sum / (len_X * len_Y)^(lam / 2)
    if (is.finite(score) && score > W_density) {
      W_density <- score; best_rem_X <- rem_X; best_rem_Y <- rem_Y
    }
  }

  # W_density starts at 0, so if no candidate beats 0 the full matrix is
  # returned unpruned. Reachable whenever the input has negative entries.
  list(density = W_density, X = which(!best_rem_X), Y = which(!best_rem_Y))
}


#' Repeatedly peel dense X x Y blocks.
#' @return list of list(X, Y, density), indices into the original XY
greedy_peeling_XY_all <- function(XY, lam = 1.5, max_round = 5) {
  XY <- as.matrix(XY); storage.mode(XY) <- "double"
  m <- nrow(XY); n <- ncol(XY)
  org_density <- sum(XY) / (m * n)^(lam / 2)

  remain_X <- seq_len(m); remain_Y <- seq_len(n)
  Wp <- XY; result <- list()

  for (i in seq_len(max_round)) {
    if (length(remain_X) < 1L || length(remain_Y) < 1L) break
    pk <- greedy_peeling_XY_one(Wp, lam)

    if (pk$density > org_density && i < max_round &&
        length(pk$X) > 0 && length(pk$Y) > 0) {
      result[[length(result) + 1L]] <- list(
        X = sort(remain_X[pk$X]),
        Y = sort(remain_Y[pk$Y]),
        density = pk$density
      )
      remain_X <- setdiff(seq_len(m), unlist(lapply(result, `[[`, "X")))
      remain_Y <- setdiff(seq_len(n), unlist(lapply(result, `[[`, "Y")))
      if (length(remain_X) == 0L || length(remain_Y) == 0L) break
      Wp <- XY[remain_X, remain_Y, drop = FALSE]
    } else break
  }
  result
}


## ---------------------------------------------------------------------
## Step 2/3: within-omic (symmetric) peeling
## ---------------------------------------------------------------------

#' Peel one dense community from a symmetric non-negative matrix.
greedy_peeling_X_one <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  m <- nrow(Wp)
  C <- colSums(Wp)
  rem <- logical(m); best_rem <- logical(m)
  W_density <- 0

  for (i in seq_len(m)) {
    Cc <- C; Cc[rem] <- Inf
    Ind <- which.min(Cc)
    rem[Ind] <- TRUE
    C <- C - Wp[Ind, ]; C[rem] <- 0
    Wp_sum <- sum(C)

    denom <- (m - i)^lam
    # At i == m the denominator is 0; MATLAB yields NaN/Inf here and the
    # guard rejects it, so the all-removed state is never selected.
    score <- if (denom > 0) Wp_sum / denom else Inf
    if (is.finite(score) && score > W_density) {
      W_density <- score; best_rem <- rem
    }
  }
  list(density = W_density, X = which(!best_rem))
}


#' Repeatedly peel symmetric communities.
#'
#' WARNING: assumes NON-NEGATIVE input. On a signed correlation submatrix,
#' which.min picks the most anti-correlated node rather than the most weakly
#' connected one, and negative densities never beat the W_density = 0 start,
#' so every node comes back unpruned. Pass abs(R) unless the MATLAB reference
#' says otherwise.
greedy_peeling_X_all <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  stopifnot(isTRUE(all.equal(Wp, t(Wp), tolerance = 1e-8)))
  m <- nrow(Wp)
  org_density <- sum(Wp) / m^lam

  remain <- seq_len(m); cur <- Wp; result <- list()
  repeat {
    pk <- greedy_peeling_X_one(cur, lam)
    if (pk$density > org_density && length(pk$X) > 0) {
      result[[length(result) + 1L]] <- list(X = remain[pk$X], density = pk$density)
      remain <- setdiff(seq_len(m), unlist(lapply(result, `[[`, "X")))
      if (length(remain) == 0L) break
      cur <- Wp[remain, remain, drop = FALSE]
    } else break
  }
  result
}


## ---------------------------------------------------------------------
## sort_result
## ---------------------------------------------------------------------

#' Order features within each block by descending within-block marginal.
#'
#' Ties can order differently from MATLAB/numpy, whose sorts are not stable
#' while R's order() is. Membership is unaffected.
sort_result <- function(Wp, blocks, symmetric = NULL) {
  Wp <- as.matrix(Wp)
  if (is.null(symmetric)) symmetric <- nrow(Wp) == ncol(Wp)

  lapply(blocks, function(blk) {
    if (symmetric) {
      idx <- blk$X
      sub <- Wp[idx, idx, drop = FALSE]
      blk$X <- idx[order(rowSums(sub), decreasing = TRUE)]
    } else {
      xi <- blk$X; yi <- blk$Y
      sub <- Wp[xi, yi, drop = FALSE]
      blk$X <- xi[order(abs(rowSums(sub)), decreasing = TRUE)]
      blk$Y <- yi[order(abs(colSums(sub)), decreasing = TRUE)]
    }
    blk
  })
}


## ---------------------------------------------------------------------
## Driver
## ---------------------------------------------------------------------

#' Full 3-step greedy peeling.
#'
#' @param P either the m x n cross matrix (cross_only = TRUE) or the
#'   (m+n) x (m+n) combined matrix, of which only P[1:m, (m+1):(m+n)] is read
#' @param R (m+n) x (m+n) combined correlation matrix. Read only when a block
#'   exceeds X_size[2] or Y_size[2]; may be NULL when those are Inf.
#' @param X_size,Y_size c(min, max). The min filter is a strict >, so a block
#'   of exactly min is dropped. MATLAB defaults are c(10, 200) -- too large
#'   for small blocks, which get filtered away entirely.
#' @param sort_xy m x n matrix to sort blocks against, or NULL. MATLAB's
#'   greedy_3step has no XY-level sort; scripts call sort_result separately.
#' @return list of list(X, Y, density); X in 1..m, Y in 1..n
greedy_3step <- function(P, R = NULL, m, n,
                         lambda1 = 1.5, lambda2 = 1.5, lambda3 = 1.5,
                         X_size = c(10, 200), Y_size = c(10, 200),
                         max_round = 5, sort_xy = NULL, cross_only = FALSE) {
  P <- as.matrix(P); storage.mode(P) <- "double"

  Wp <- if (cross_only) {
    stopifnot(nrow(P) == m, ncol(P) == n); P
  } else {
    stopifnot(nrow(P) == m + n, ncol(P) == m + n)
    P[1:m, (m + 1):(m + n), drop = FALSE]
  }

  blocks <- greedy_peeling_XY_all(Wp, lambda1, max_round)

  if (!is.null(sort_xy)) {
    blocks <- sort_result(as.matrix(sort_xy), blocks, symmetric = FALSE)
  }

  blocks <- Filter(function(b)
    length(b$X) > X_size[1] && length(b$Y) > Y_size[1], blocks)
  if (length(blocks) == 0L) {
    warning("all blocks removed by the minimum-size filter (X_size[1] = ",
            X_size[1], ", Y_size[1] = ", Y_size[1], ")")
    return(list())
  }

  needs_R <- any(vapply(blocks, function(b)
    length(b$X) > X_size[2] || length(b$Y) > Y_size[2], logical(1)))
  if (needs_R) {
    if (is.null(R)) stop("a block exceeds its size maximum but R is NULL")
    R <- as.matrix(R); storage.mode(R) <- "double"
    stopifnot(nrow(R) == m + n, ncol(R) == m + n)
  }

  refine <- function(idx, sub, lam) {
    sub <- (sub + t(sub)) / 2
    comm <- greedy_peeling_X_all(sub, lam)
    if (length(comm) == 0L) return(idx)   # leave unrefined rather than empty
    comm <- sort_result(sub, comm, symmetric = TRUE)
    idx[unlist(lapply(comm, `[[`, "X"))]
  }

  for (k in seq_along(blocks)) {
    if (length(blocks[[k]]$X) > X_size[2]) {
      xi <- blocks[[k]]$X
      blocks[[k]]$X <- refine(xi, R[xi, xi, drop = FALSE], lambda2)
    }
    if (length(blocks[[k]]$Y) > Y_size[2]) {
      yi <- blocks[[k]]$Y
      iy <- m + yi                        # offset for indexing R only
      blocks[[k]]$Y <- refine(yi, R[iy, iy, drop = FALSE], lambda3)
    }
  }
  blocks
}
