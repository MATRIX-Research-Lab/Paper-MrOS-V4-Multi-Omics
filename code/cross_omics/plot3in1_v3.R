# R port of result/plot3in1_v3.m -- side-by-side (a) Active | (b) Deceased 3-in-1 plot
#
# Reuses the helpers of the illustration port plot3in1_v1.R (matlab_jet, map_col,
# unique_keep_last, draw_quads), which must be sourced first.
#
# Same call as MATLAB:
#   plot3in1_v3(R_yyA, R_xxA, t(R_xyA), R_yyD, R_xxD, t(R_xyD), result, TRUE, c(5, 5), TRUE,
#               "Metabolite", "Microbiome", 2.2)
# result: list of list(idx_col1, idx_col2) as the MATLAB p x 2 cell array; col1 indexes the
# right-triangle matrix (Y argument), col2 the top-triangle matrix (X argument).
#
# v3 vs v1 (as in the MATLAB code):
#   - two panels, letters a / b, x / y axis labels, colorbar only on panel b
#   - block boxes alternate blue [0 .447 .741] / orange [1 .5 0], 4 pt; empty rows are skipped
#     and the LAST box is extended down to the bottom edge of the rectangle
# Draws into the current device; open it at 20 x (9.5 + 0.5 * xlabel_offset) in, like the
# MATLAB PaperSize (see open_v3_device()).

open_v3_device <- function(file, xlabel_offset = 0, res = 300) {
  w <- 20; h <- 9.5 + 0.5 * xlabel_offset
  if (grepl("\\.pdf$", file)) {
    # cairo_pdf warns and opens no file when the cairo library is missing.
    # Fall back to pdf() so the named manuscript PDF is still written.
    opened <- FALSE
    suppressWarnings(tryCatch({
      grDevices::cairo_pdf(file, width = w, height = h, family = "Helvetica")
      opened <- file.exists(file) && grDevices::dev.cur() > 1L
    }, error = function(e) NULL))
    if (!isTRUE(opened)) {
      if (grDevices::dev.cur() > 1L) grDevices::dev.off()
      grDevices::pdf(file, width = w, height = h, family = "Helvetica")
    }
  } else {
    ragg::agg_png(file, width = w, height = h, units = "in", res = res, background = "white")
  }
  par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0))
  c(w, h)
}

# One panel. region = c(x0, x1, y0, y1) in figure units; margins (inches) leave room for the
# tick labels / axis labels. Returns the plot box in figure units (for the colorbar).
draw3in1_panel <- function(X, Y, XY, result, region, figsize, line = TRUE, step = c(100, 10),
                           clim_tri = c(-1, 1), clim_mid = NULL,
                           xlabel_str = NULL, ylabel_str = NULL, xlabel_offset = 0,
                           panel_label = NULL, mar_in = c(0.9, 1.1, 0.1, 0.1),
                           tick_ps = 18, label_ps = 16, letter_ps = 20, box_lw = 4) {
  result <- lapply(result, function(r) list(as.integer(unlist(r[[2]])), as.integer(unlist(r[[1]]))))
  select_X <- unique_keep_last(unlist(lapply(result, `[[`, 1)))
  select_Y <- unique_keep_last(unlist(lapply(result, `[[`, 2)))
  X <- as.matrix(X); Y <- as.matrix(Y); XY <- as.matrix(XY)
  C <- Y[select_Y, select_Y, drop = FALSE]
  B <- X[select_X, select_X, drop = FALSE]
  A <- XY[select_X, select_Y, drop = FALSE]
  m <- nrow(A); n <- ncol(A); ratio <- m / n / 1.2
  cmap  <- matlab_jet(256)
  clim2 <- if (!is.null(clim_mid)) clim_mid else if (max(A, na.rm = TRUE) > 1) c(0, 10) else c(-1, 1)

  # ---- axis equal tight inside the region (minus margins), centred ----
  xlim <- c(0, max(2 * m / ratio, n + 2 * m / ratio)); ylim <- c(0, 2 * n + m / ratio)
  bw <- (region[2] - region[1]) * figsize[1] - mar_in[2] - mar_in[4]
  bh <- (region[4] - region[3]) * figsize[2] - mar_in[1] - mar_in[3]
  s  <- min(bw / diff(xlim), bh / diff(ylim))
  w  <- s * diff(xlim) / figsize[1]; h <- s * diff(ylim) / figsize[2]
  cx <- region[1] + (mar_in[2] + bw / 2) / figsize[1]
  cy <- region[3] + (mar_in[1] + bh / 2) / figsize[2]
  plt <- c(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)
  par(plt = plt, xaxs = "i", yaxs = "i", ps = tick_ps, family = "Helvetica",
      col.axis = MATLAB_GREY, fg = MATLAB_GREY, new = TRUE, xpd = FALSE)
  plot.new(); plot.window(xlim, ylim)

  # ---- ticks first (MATLAB Layer = 'bottom') ----
  xTickIdx <- seq(step[1], floor(m / step[1]) * step[1], by = step[1])
  yTickIdx <- rev(seq(step[2], floor(n / step[2]) * step[2], by = step[2]))
  tcl <- (0.01 * max(w * figsize[1], h * figsize[2])) / par("csi")
  if (length(xTickIdx) && xTickIdx[1] <= m)
    axis(1, at = (xTickIdx - 1) * 2 / ratio, labels = xTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
         tcl = tcl, mgp = c(0, 0.1, 0), col.ticks = MATLAB_GREY)
  if (length(yTickIdx) && yTickIdx[1] <= n)
    axis(2, at = (n - yTickIdx) * 2 + 1, labels = yTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
         tcl = tcl, mgp = c(3, 0.35, 0), las = 1, col.ticks = MATLAB_GREY)

  vx <- c(0, 1, 0, -1); vy <- c(1, 0, -1, 0)
  # ---- 1) top diamonds (X) ----
  ij <- do.call(rbind, lapply(1:m, function(i) cbind(i, i:m)))
  x0 <- (ij[, 2] * 2 - ij[, 1]) / ratio; y0 <- (ij[, 1] - 1) / ratio + n * 2
  draw_quads(outer(vx / ratio, x0, `+`), outer(vy / ratio, y0, `+`),
             map_col(B[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim_tri, cmap))
  # ---- 2) right diamonds (Y) ----
  ij <- do.call(rbind, lapply(1:n, function(i) cbind(i, i:n)))
  x0 <- ij[, 1] - 1 + m * 2 / ratio; y0 <- 2 * n - ij[, 2] * 2 + ij[, 1]
  draw_quads(outer(vx, x0, `+`), outer(vy, y0, `+`),
             map_col(C[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim_tri, cmap))
  # ---- 3) centre rectangle last (on top) ----
  ij <- expand.grid(j = 1:n, i = 1:m)
  x0 <- (ij$i - 1) * 2 / ratio; y0 <- (ij$j - 1) * 2
  cols <- map_col(A[cbind(ij$i, n + 1 - ij$j)], clim2, cmap)
  rect(x0, y0, x0 + 2 / ratio, y0 + 2, col = cols, border = cols, lwd = 0.25)

  # ---- block boxes: blue / orange, last box extended to the bottom ----
  if (line) {
    L <- t(sapply(result, function(r) c(length(r[[1]]), length(r[[2]]))))
    L <- L[L[, 2] > 0 & L[, 1] > 0, , drop = FALSE]      # MATLAB: drop empty rows
    L[, 1] <- L[, 1] * 2 / ratio; L[, 2] <- L[, 2] * 2
    total_w <- m * 2 / ratio; total_h <- n * 2
    box_cols <- c(rgb(0, 0.447, 0.741), rgb(1, 0.5, 0))
    x <- 0; y <- total_h
    for (k in seq_len(nrow(L))) {
      rw <- min(L[k, 1], total_w - x)
      if (k == nrow(L)) { ry <- 0; rh <- min(y, total_h) }
      else { ry <- max(y - L[k, 2], 0); rh <- min(L[k, 2], y) }
      if (rw > 0 && rh > 0)
        rect(x, ry, x + rw, ry + rh, border = box_cols[(k - 1) %% 2 + 1], lwd = box_lw / 0.75)
      y <- y - L[k, 2]; x <- x + L[k, 1]
    }
  }

  # ---- axis labels and panel letter (MATLAB text() positions, data units) ----
  par(xpd = NA)
  if (!is.null(xlabel_str))
    text(m / ratio, -n * 0.1 - xlabel_offset, xlabel_str, adj = c(0.5, 1), cex = label_ps / tick_ps,
         col = MATLAB_GREY)
  if (!is.null(ylabel_str))
    text(-m * 0.15 / ratio, n, ylabel_str, srt = 90, adj = c(0.5, 0), cex = label_ps / tick_ps,
         col = MATLAB_GREY)
  if (!is.null(panel_label))
    text(0, 2 * n + m / ratio, panel_label, adj = c(0, 1), font = 2, cex = letter_ps / tick_ps,
         col = "black")
  par(xpd = FALSE)
  invisible(list(plt = plt, clim = clim2, sq_frac = (2 * n) / (2 * n + m / ratio),
                 select_X = select_X, select_Y = select_Y))
}

draw_colorbar <- function(pos, clim, cmap = matlab_jet(256), tick_ps = 18) {
  par(plt = pos, new = TRUE, xpd = FALSE, ps = tick_ps)
  plot.new(); plot.window(c(0, 1), clim, xaxs = "i", yaxs = "i")
  yb <- seq(clim[1], clim[2], length.out = length(cmap) + 1)
  rect(0, yb[-length(yb)], 1, yb[-1], col = cmap, border = NA)
  at <- pretty(clim, n = 10)
  axis(4, at = at, labels = format(at, drop0trailing = TRUE, trim = TRUE), lwd = 0,
       lwd.ticks = 0.5 / 0.75, tcl = -0.25, las = 1, mgp = c(3, 0.4, 0), col.ticks = MATLAB_GREY)
  box(lwd = 0.5 / 0.75, col = MATLAB_GREY)
}

plot3in1_v3 <- function(X_A, Y_A, XY_A, X_D, Y_D, XY_D, result, line = FALSE, step = c(100, 10),
                        bar = FALSE, xlabel_str = "Metabolite", ylabel_str = "Microbiome",
                        xlabel_offset = 0, figsize = c(20, 9.5 + 0.5 * xlabel_offset), ...) {
  regions <- list(a = c(0.00, 0.47, 0.00, 1.00), b = c(0.47, 0.94, 0.00, 1.00))
  mats <- list(a = list(X_A, Y_A, XY_A), b = list(X_D, Y_D, XY_D))
  info <- list()
  plot.new()                                   # open the page; panels then draw with new = TRUE
  for (p in names(regions)) {
    info[[p]] <- draw3in1_panel(mats[[p]][[1]], mats[[p]][[2]], mats[[p]][[3]], result,
                                region = regions[[p]], figsize = figsize, line = line, step = step,
                                xlabel_str = xlabel_str, ylabel_str = ylabel_str,
                                xlabel_offset = xlabel_offset, panel_label = p, ...)
  }
  if (bar) {
    pb <- info$b$plt
    draw_colorbar(c(0.955, 0.965, pb[3], pb[4]), info$b$clim)
  }
  invisible(info)
}
