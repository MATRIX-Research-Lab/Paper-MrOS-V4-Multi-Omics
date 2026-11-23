# plot_square_v4.R -- (a) Active | (b) Deceased cross-omics heatmaps without the
# within-omics triangles.
#
# Compared with plot3in1_v3.R:
#   - only the cross-correlation rectangle (the "middle square") is drawn
#   - every row and every column gets a tick mark; tick labels every `step`
#     cells (step = NULL picks the label step automatically)
#   - block boxes alternate blue / orange as before and carry a circled block
#     number inside the box: (1) in the first (blue) box, (2) in the second
#     (orange) box, ...
#   - axis titles "Number of <omics>"; legend = TRUE adds a strip under the two
#     panels giving the size of every block
#
# Needs matlab_jet(), map_col(), unique_keep_last() and MATLAB_GREY from
# plot3in1_v1.R: source that file first.
#
# plot_square_v4(R_A, R_D, blocks, ...)
#   R_A, R_D : Active / Deceased correlation matrices, rows = the omics on the
#              vertical axis, columns = the omics on the horizontal axis
#              (S$R_xyA / S$R_xyD from load_pair(); micro x metab, micro x prot,
#              metab x prot)
#   blocks   : list of list(X = row indices, Y = column indices), i.e. one
#              element of cross_omics_blocks.rds written by 01_find_blocks.R

BLOCK_COLS <- c(rgb(0, 0.447, 0.741), rgb(1, 0.5, 0))   # blue, orange, blue, ...

open_square_device <- function(file, figsize = c(20, 9.3), res = 300) {
  w <- figsize[1]; h <- figsize[2]
  if (grepl("\\.pdf$", file)) {
    # cairo_pdf warns and opens no file when the cairo library is missing;
    # fall back to pdf() so the named PDF is still written.
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
  invisible(figsize)
}

# Tick-label step: every cell when there are few cells, otherwise the smallest
# "nice" step that leaves at most `max_labels` labels on the axis.
auto_step <- function(count, max_labels = 20) {
  for (s in c(1, 2, 5, 10, 20, 25, 50, 100)) if (count / s <= max_labels) return(s)
  100
}

# Circled number: white disc with a coloured ring and the number inside.
# r_in = radius in inches, so the circle has the same size on every device and
# no special Unicode glyph is needed.
circled_number <- function(x, y, k, col, r_in = 0.13, cex = 0.85, lwd = 2) {
  ipx <- par("pin")[1] / diff(par("usr")[1:2])      # inches per user unit
  ipy <- par("pin")[2] / diff(par("usr")[3:4])
  t <- seq(0, 2 * pi, length.out = 72)
  polygon(x + r_in / ipx * cos(t), y + r_in / ipy * sin(t), col = "white", border = col,
          lwd = lwd, xpd = NA)
  text(x, y, labels = k, col = col, font = 2, cex = cex, xpd = NA)
}

# Width in inches of the widest label at `ps` points (Helvetica). Drivers use it to
# choose a page height that leaves room for feature names under the x axis.
label_extent_in <- function(labels, ps = 14, font = 1) {
  if (is.null(labels) || !length(labels)) return(0)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  par(ps = ps, family = "Helvetica"); plot.new()
  max(strwidth(labels, units = "inches", font = font))
}

# Vertical space in inches taken by column names drawn at `srt` degrees
# (90 = vertical, 45 / 30 = tilted up to the right, hanging down-left from the tick).
col_label_depth_in <- function(labels, ps = 14, srt = 90, font = 1) {
  if (is.null(labels) || !length(labels)) return(0)
  rad <- srt * pi / 180
  label_extent_in(labels, ps, font) * sin(rad) + (ps / 72) * cos(rad)
}

# One panel. region = c(x0, x1, y0, y1) in figure units; mar_in = margins in
# inches (bottom, left, top, right) for the tick labels and axis labels.
# The heatmap keeps the width / height ratio `panel_asp` (1.2 as in plot3in1_v3).
draw_square_panel <- function(R, blocks, region, figsize, step = NULL, clim = c(-1, 1),
                              xlabel_str = NULL, ylabel_str = NULL, panel_label = NULL,
                              row_labels = NULL, col_labels = NULL, name_ps = 14,
                              row_font = 1, col_font = 1, col_srt = 90,
                              panel_asp = 1.2, mar_in = c(1.3, 1.5, 0.6, 0.3),
                              tick_ps = 18, label_ps = 16, letter_ps = 20, box_lw = 4,
                              circle_r_in = 0.14, circle_inset_in = 0.32) {
  # row_labels / col_labels: feature names for the selected rows / columns (in the
  # order of `rows` / `cols` below). When given they replace the numeric tick labels
  # on that axis: row names horizontal on the left, column names below at `col_srt`
  # degrees (90 = vertical; 45 or 30 = tilted, which needs less vertical space).
  rows <- unique_keep_last(unlist(lapply(blocks, function(b) as.integer(b$X))))
  cols <- unique_keep_last(unlist(lapply(blocks, function(b) as.integer(b$Y))))
  M  <- as.matrix(R)[rows, cols, drop = FALSE]
  nr <- nrow(M); nc <- ncol(M)
  cmap <- matlab_jet(256)

  # ---- margins grow to fit the feature names (names + 0.25 in gap + axis title) ----
  gap_in <- 0.25; title_in <- 0.55
  ext_x <- ext_y <- depth_x <- 0
  rad <- col_srt * pi / 180
  if (!is.null(col_labels) || !is.null(row_labels)) {
    op <- par(ps = name_ps, family = "Helvetica")
    if (!is.null(col_labels)) ext_x <- max(strwidth(col_labels, units = "inches", font = col_font))
    if (!is.null(row_labels)) ext_y <- max(strwidth(row_labels, units = "inches", font = row_font))
    par(op)
    depth_x <- ext_x * sin(rad) + (name_ps / 72) * cos(rad)    # vertical extent of tilted names
    if (!is.null(col_labels)) mar_in[1] <- depth_x + gap_in + if (is.null(xlabel_str)) 0.15 else title_in
    if (!is.null(row_labels)) mar_in[2] <- ext_y + gap_in + if (is.null(ylabel_str)) 0.15 else title_in
  }

  # ---- plot box with aspect panel_asp, centred inside the region minus margins ----
  bw <- (region[2] - region[1]) * figsize[1] - mar_in[2] - mar_in[4]
  bh <- (region[4] - region[3]) * figsize[2] - mar_in[1] - mar_in[3]
  w_in <- min(bw, bh * panel_asp); h_in <- w_in / panel_asp
  cx <- region[1] + (mar_in[2] + bw / 2) / figsize[1]
  cy <- region[3] + (mar_in[1] + bh / 2) / figsize[2]
  w <- w_in / figsize[1]; h <- h_in / figsize[2]
  plt <- c(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2)
  par(plt = plt, xaxs = "i", yaxs = "i", ps = tick_ps, family = "Helvetica",
      col.axis = MATLAB_GREY, fg = MATLAB_GREY, new = TRUE, xpd = FALSE)
  plot.new(); plot.window(c(0, nc), c(0, nr))

  # ---- cells: column i spans [i - 1, i]; row j (counted from the top) spans [nr - j, nr - j + 1] ----
  ij <- expand.grid(j = 1:nr, i = 1:nc)
  cc <- map_col(M[cbind(ij$j, ij$i)], clim, cmap)
  rect(ij$i - 1, nr - ij$j, ij$i, nr - ij$j + 1, col = cc, border = cc, lwd = 0.25)

  # ---- ticks: a short tick on every cell, a longer labelled tick every `step` cells ----
  sx <- if (is.null(step)) auto_step(nc) else step[1]
  sy <- if (is.null(step)) auto_step(nr) else step[2]
  lx <- seq(sx, nc, by = sx); ly <- seq(sy, nr, by = sy)
  tk <- list(lwd = 0, lwd.ticks = 1 / 0.75, col.ticks = MATLAB_GREY)
  do.call(axis, c(list(1, at = (1:nc) - 0.5, labels = FALSE, tcl = -0.2), tk))
  do.call(axis, c(list(1, at = lx - 0.5, labels = FALSE, tcl = -0.45), tk))
  do.call(axis, c(list(2, at = nr - (1:nr) + 0.5, labels = FALSE, tcl = -0.2), tk))
  do.call(axis, c(list(2, at = nr - ly + 0.5, labels = FALSE, tcl = -0.45), tk))
  # labels via mtext: axis() silently drops every other label when 3-digit labels get close
  ipx <- par("pin")[1] / nc; ipy <- par("pin")[2] / nr                # inches per cell
  if (is.null(col_labels)) {
    mtext(lx, side = 1, at = lx - 0.5, line = 0.75, col = MATLAB_GREY)
  } else {
    text((1:nc) - 0.5, -gap_in / ipy, col_labels, srt = col_srt,
         adj = if (col_srt >= 89) c(1, 0.5) else c(1, 1), xpd = NA,
         cex = name_ps / tick_ps, font = col_font, col = MATLAB_GREY)
  }
  if (is.null(row_labels)) {
    mtext(ly, side = 2, at = nr - ly + 0.5, line = 0.75, las = 1, adj = 1, col = MATLAB_GREY)
  } else {
    text(-gap_in / ipx, nr - (1:nr) + 0.5, row_labels, adj = c(1, 0.5), xpd = NA,
         cex = name_ps / tick_ps, font = row_font, col = MATLAB_GREY)
  }

  # ---- block boxes (blue / orange, last box extended to the bottom edge as in v3)
  #      with the circled block number inside each box ----
  L <- t(sapply(blocks, function(b) c(length(b$Y), length(b$X))))   # width (cols), height (rows)
  x <- 0; y <- nr
  for (k in seq_len(nrow(L))) {
    rw <- min(L[k, 1], nc - x)
    if (k == nrow(L)) { ry <- 0; rh <- min(y, nr) }
    else { ry <- max(y - L[k, 2], 0); rh <- min(L[k, 2], y) }
    if (rw > 0 && rh > 0) {
      ck <- BLOCK_COLS[(k - 1) %% 2 + 1]
      rect(x, ry, x + rw, ry + rh, border = ck, lwd = box_lw / 0.75)
      circled_number(x + circle_inset_in / ipx, ry + rh - circle_inset_in / ipy, k, ck,
                     r_in = circle_r_in, cex = 0.9)
    }
    y <- y - L[k, 2]; x <- x + L[k, 1]
  }

  # ---- axis labels and panel letter ----
  if (!is.null(xlabel_str)) {
    if (is.null(col_labels))
      mtext(xlabel_str, side = 1, line = 3.0, cex = label_ps / tick_ps, col = MATLAB_GREY)
    else
      text(nc / 2, -(depth_x + gap_in + 0.3) / ipy, xlabel_str, adj = c(0.5, 1), xpd = NA,
           cex = label_ps / tick_ps, col = MATLAB_GREY)
  }
  if (!is.null(ylabel_str)) {
    if (is.null(row_labels))
      mtext(ylabel_str, side = 2, line = 3.6, cex = label_ps / tick_ps, col = MATLAB_GREY)
    else
      text(-(ext_y + gap_in + 0.3) / ipx, nr / 2, ylabel_str, srt = 90, adj = c(0.5, 0), xpd = NA,
           cex = label_ps / tick_ps, col = MATLAB_GREY)
  }
  if (!is.null(panel_label))
    mtext(panel_label, side = 3, line = 0.6, adj = 0, at = -0.08 * nc, font = 2,
          cex = letter_ps / tick_ps, col = "black")
  invisible(list(plt = plt, clim = clim, rows = rows, cols = cols))
}

draw_colorbar_v4 <- function(pos, clim, cmap = matlab_jet(256), tick_ps = 18) {
  par(plt = pos, new = TRUE, xpd = FALSE, ps = tick_ps, family = "Helvetica",
      col.axis = MATLAB_GREY, fg = MATLAB_GREY)
  plot.new(); plot.window(c(0, 1), clim, xaxs = "i", yaxs = "i")
  yb <- seq(clim[1], clim[2], length.out = length(cmap) + 1)
  rect(0, yb[-length(yb)], 1, yb[-1], col = cmap, border = NA)
  at <- pretty(clim, n = 10)
  axis(4, at = at, labels = format(at, drop0trailing = TRUE, trim = TRUE), lwd = 0,
       lwd.ticks = 0.5 / 0.75, tcl = -0.25, las = 1, mgp = c(3, 0.4, 0), col.ticks = MATLAB_GREY)
  box(lwd = 0.5 / 0.75, col = MATLAB_GREY)
}

# Legend strip: "(k) <entry k>" items, packed on one line when they fit, centred.
draw_block_legend <- function(entries, region, figsize, ps = 16, r_in = 0.13,
                              gap_in = 0.12, sep_in = 0.7) {
  par(plt = region, new = TRUE, xpd = NA, ps = ps, family = "Helvetica")
  plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
  W <- par("pin")[1]                                   # strip width in inches
  item <- 2 * r_in + gap_in + strwidth(entries, units = "inches")
  lines <- list(); cur <- integer(0); used <- 0
  for (k in seq_along(entries)) {
    add <- item[k] + if (length(cur)) sep_in else 0
    if (length(cur) && used + add > 0.98 * W) {
      lines[[length(lines) + 1]] <- cur; cur <- integer(0); used <- 0; add <- item[k]
    }
    cur <- c(cur, k); used <- used + add
  }
  lines[[length(lines) + 1]] <- cur
  nl <- length(lines); ys <- rev(seq_len(nl)) / (nl + 1)
  for (li in seq_len(nl)) {
    ks <- lines[[li]]
    x <- (W - sum(item[ks]) - sep_in * (length(ks) - 1)) / 2 / W
    for (k in ks) {
      circled_number(x + r_in / W, ys[li], k, BLOCK_COLS[(k - 1) %% 2 + 1], r_in = r_in, cex = 0.85)
      text(x + (2 * r_in + gap_in) / W, ys[li], entries[k], adj = c(0, 0.5), col = MATLAB_GREY)
      x <- x + (item[k] + sep_in) / W
    }
  }
}

plot_square_v4 <- function(R_A, R_D, blocks, step = NULL, bar = TRUE,
                           xlabel_str = "Number of metabolites", ylabel_str = "Number of microbes",
                           x_unit = "metabolites", y_unit = "microbes",
                           legend = FALSE,
                           legend_fmt = "Number of %s = %d, number of %s = %d",
                           figsize = c(20, 9.3), legend_h = 0.7, clim = c(-1, 1), ...) {
  blocks <- Filter(function(b) length(b$X) > 0 && length(b$Y) > 0, blocks)
  lb <- if (legend) legend_h / figsize[2] else 0
  regions <- list(a = c(0.00, 0.47, lb, 1.00), b = c(0.47, 0.94, lb, 1.00))
  mats <- list(a = R_A, b = R_D)
  info <- list()
  plot.new()                                   # open the page; panels then draw with new = TRUE
  for (p in names(regions)) {
    info[[p]] <- draw_square_panel(mats[[p]], blocks, regions[[p]], figsize, step = step,
                                   clim = clim, xlabel_str = xlabel_str, ylabel_str = ylabel_str,
                                   panel_label = p, ...)
  }
  if (bar) {
    pb <- info$b$plt
    draw_colorbar_v4(c(0.955, 0.965, pb[3], pb[4]), clim)
  }
  if (legend) {
    entries <- vapply(blocks, function(b) sprintf(legend_fmt, y_unit, length(b$X), x_unit, length(b$Y)), "")
    draw_block_legend(entries, c(0.03, 0.94, 0.0, lb), figsize)
  }
  invisible(info)
}
