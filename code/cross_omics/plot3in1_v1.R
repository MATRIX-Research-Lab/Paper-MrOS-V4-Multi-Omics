# R port of plot3in1_v1.m (base graphics; mirrors plot3in1_v1.py)
#
# Middle : rectangular heatmap of XY (selected X rows x selected Y cols)
# Top    : X(select_X, select_X) as a 45-degree rotated lower triangle (diamonds)
# Right  : Y(select_Y, select_Y) as a 45-degree rotated lower triangle (diamonds)
#
# MATLAB look: jet(256), caxis-style scaled colour mapping, inward ticks, no axis lines,
# Helvetica 10 pt, [0.15 0.15 0.15] text, default axes Position [0.13 0.11 0.775 0.815].
# Draws into the CURRENT device; open it with figsize 5.6 x 4.2 in (see plot_pattern_3in1.R).

MATLAB_GREY <- rgb(0.15, 0.15, 0.15)

matlab_jet <- function(m = 256) {            # exact port of MATLAB jet(m)
  n <- ceiling(m / 4)
  u <- c((1:n) / n, rep(1, n - 1), (n:1) / n)
  g <- ceiling(n / 2) - (m %% 4 == 1) + seq_along(u)
  r <- g + n; b <- g - n
  r <- r[r <= m]; g <- g[g <= m]; b <- b[b >= 1]
  J <- matrix(0, m, 3)
  J[r, 1] <- u[seq_along(r)]
  J[g, 2] <- u[seq_along(g)]
  J[b, 3] <- u[(length(u) - length(b) + 1):length(u)]
  rgb(J[, 1], J[, 2], J[, 3])
}

# MATLAB CDataMapping 'scaled': index = fix((c - cmin)/(cmax - cmin) * m) + 1, clamped
map_col <- function(v, clim, cmap) {
  m <- length(cmap)
  idx <- floor((v - clim[1]) / (clim[2] - clim[1]) * m) + 1
  cmap[pmin(pmax(idx, 1), m)]
}

unique_keep_last <- function(v) rev(rev(v)[!duplicated(rev(v))])

# Draw many 4-vertex polygons in one call (NA-separated), border = fill to hide seams
draw_quads <- function(px, py, cols) {       # px, py: 4 x k matrices
  k <- ncol(px)
  polygon(as.vector(rbind(px, NA))[-(5 * k)], as.vector(rbind(py, NA))[-(5 * k)],
          col = cols, border = cols, lwd = 0.25)
}

plot3in1_v1 <- function(X, Y, XY, result, line = FALSE, step = NULL, bar = FALSE,
                        one_based = TRUE, figsize = c(5.6, 4.2),
                        cbar_pos = c(0.15, 0.13, 0.013, 0.45),
                        clim_tri = NULL, clim_mid = NULL, box_lw = 2) {
  # result: list of length-2 lists, each = list(idx_col1, idx_col2) as in the MATLAB cell
  # array; columns are swapped internally exactly like result(:, [2 1]).
  # cbar_pos: c(left, bottom, width, height) in figure units, or "right" = just right of the plot,
  # bottom/height aligned with the middle square.
  # box_lw is in points (as MATLAB/matplotlib); R lwd = points / 0.75.
  result <- lapply(result, function(r) list(as.integer(unlist(r[[2]])), as.integer(unlist(r[[1]]))))
  off <- if (one_based) 0L else 1L
  select_X <- unique_keep_last(unlist(lapply(result, `[[`, 1)) + off)
  select_Y <- unique_keep_last(unlist(lapply(result, `[[`, 2)) + off)

  X <- as.matrix(X); Y <- as.matrix(Y); XY <- as.matrix(XY)
  C <- Y[select_Y, select_Y, drop = FALSE]
  B <- X[select_X, select_X, drop = FALSE]
  A <- XY[select_X, select_Y, drop = FALSE]
  m <- nrow(A); n <- ncol(A)
  ratio <- m / n / 1.2
  cmap <- matlab_jet(256)
  clim1 <- if (is.null(clim_tri)) c(-1, 1) else clim_tri
  clim2 <- if (!is.null(clim_mid)) clim_mid else if (max(A, na.rm = TRUE) > 1) c(0, 10) else c(-1, 1)

  # ---- axis equal tight: shrink the plot box inside the MATLAB axes Position ----
  xlim <- c(0, max(2 * m / ratio, n + 2 * m / ratio))
  ylim <- c(0, 2 * n + m / ratio)
  pos <- c(0.13, 0.11, 0.775, 0.815)
  box_w <- pos[3] * figsize[1]; box_h <- pos[4] * figsize[2]
  s <- min(box_w / diff(xlim), box_h / diff(ylim))         # inches per data unit
  w <- s * diff(xlim) / figsize[1]; h <- s * diff(ylim) / figsize[2]
  cx <- pos[1] + pos[3] / 2; cy <- pos[2] + pos[4] / 2
  op <- par(plt = c(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2), xaxs = "i", yaxs = "i",
            ps = 10, family = "Helvetica", col.axis = MATLAB_GREY, fg = MATLAB_GREY)
  on.exit(par(op), add = TRUE)
  plot.new(); plot.window(xlim, ylim)

  # ---- Clean tick labels (inward ticks, no axis line, MATLAB TickLength 0.01) ----
  # drawn BEFORE the patches: MATLAB Axes Layer = 'bottom' puts tick marks under the data
  xs <- if (is.null(step)) 100 else step[1]; ys <- if (is.null(step)) 10 else step[2]
  xTickIdx <- seq(xs, floor(m / xs) * xs, by = xs)
  yTickIdx <- rev(seq(ys, floor(n / ys) * ys, by = ys))
  tick_in <- 0.01 * max(w * figsize[1], h * figsize[2])     # inches
  tcl <- tick_in / par("csi")
  axis(1, at = (xTickIdx - 1) * 2 / ratio, labels = xTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
       tcl = tcl, mgp = c(0, -0.11, 0), col.ticks = MATLAB_GREY)
  axis(2, at = (n - yTickIdx) * 2 + 1, labels = yTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
       tcl = tcl, mgp = c(3, 0.23, 0), las = 1, col.ticks = MATLAB_GREY)

  # ---- Right: diamonds for Y ----
  ij <- do.call(rbind, lapply(1:n, function(i) cbind(i, i:n)))
  x0 <- ij[, 1] - 1 + m * 2 / ratio; y0 <- 2 * n - ij[, 2] * 2 + ij[, 1]
  vx <- c(0, 1, 0, -1); vy <- c(1, 0, -1, 0)
  draw_quads(outer(vx, x0, `+`), outer(vy, y0, `+`),
             map_col(C[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim1, cmap))

  # ---- Top: diamonds for X ----
  ij <- do.call(rbind, lapply(1:m, function(i) cbind(i, i:m)))
  x0 <- (ij[, 2] * 2 - ij[, 1]) / ratio; y0 <- (ij[, 1] - 1) / ratio + n * 2
  draw_quads(outer(vx / ratio, x0, `+`), outer(vy / ratio, y0, `+`),
             map_col(B[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim1, cmap))

  # ---- Middle square (drawn last = on top, as ax2 in MATLAB) ----
  ij <- expand.grid(j = 1:n, i = 1:m)
  x0 <- (ij$i - 1) * 2 / ratio; y0 <- (ij$j - 1) * 2
  cols <- map_col(A[cbind(ij$i, n + 1 - ij$j)], clim2, cmap)
  rect(x0, y0, x0 + 2 / ratio, y0 + 2, col = cols, border = cols, lwd = 0.25)

  # ---- Lines separating blocks ----
  if (line) {
    L <- t(sapply(result, function(r) c(length(r[[1]]), length(r[[2]]))))
    L[, 1] <- L[, 1] * 2 / ratio; L[, 2] <- L[, 2] * 2
    x <- 0; y <- sum(L[, 2]); i <- 1
    while (i <= nrow(L) && L[i, 1] != 0) {
      rect(x, y - L[i, 2], x + L[i, 1], y, border = "red", lwd = box_lw / 0.75)
      y <- y - L[i, 2]; x <- x + L[i, 1]; i <- i + 1
    }
  }

  # ---- Colorbar (middle panel), fixed figure position like MATLAB cb.Position ----
  if (bar) {
    if (is.character(cbar_pos) && cbar_pos == "right") {   # right of the plot, aligned with square
      sq_h <- h * (2 * n) / (2 * n + m / ratio)
      cbar_pos <- c(cx + w / 2 + 0.02, cy - h / 2, 0.013, sq_h)
    }
    par(plt = c(cbar_pos[1], cbar_pos[1] + cbar_pos[3], cbar_pos[2], cbar_pos[2] + cbar_pos[4]),
        new = TRUE)
    plot.new(); plot.window(c(0, 1), clim2)
    yb <- seq(clim2[1], clim2[2], length.out = length(cmap) + 1)
    rect(0, yb[-length(yb)], 1, yb[-1], col = cmap, border = NA)
    axis(4, at = pretty(clim2), labels = format(pretty(clim2), drop0trailing = TRUE, trim = TRUE),
         lwd = 0, lwd.ticks = 0.5 / 0.75, tcl = tcl, las = 1, mgp = c(3, 0.28, 0),
         col.ticks = MATLAB_GREY)
    box(lwd = 0.5 / 0.75, col = MATLAB_GREY)
  }
  invisible(list(select_X = select_X, select_Y = select_Y))
}
