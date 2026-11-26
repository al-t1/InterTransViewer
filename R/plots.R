#' Plot DEG counts with specific DEGs or R
#'
#' Two panels: the number of DEGs per experiment, and either the percent of
#' experiment-specific DEGs or the R metric. Draws on the current device.
#'
#' @param x An [itv] object, or a data frame from [itv_summary()].
#' @param value `"specific_percent"` or `"R"`.
#' @param main Title shown above both panels.
#' @param labels Optional character vector named by experiment, e.g. built
#'   from a table describing the experiments; see the example. Each bar is
#'   labelled with the experiment number and this label (by default, the
#'   experiment name). Entries for other experiments are ignored.
#' @return The summary data frame, invisibly.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' itv_plot_totals(x, "R", main = "Example")
#'
#' ## Labels from a description of the experiments
#' m <- read.delim(system.file("extdata", "example_description.csv",
#'                             package = "InterTransViewer"), sep = ";")
#' labels <- setNames(paste(m$Organ, m$Time, m$Method), m$experiment)
#' itv_plot_totals(x, "R", main = "Example", labels = labels)
#' @export
itv_plot_totals <- function(x, value = c("specific_percent", "R"), main = "",
                            labels = NULL) {
  value <- match.arg(value)
  tab <- if (is.data.frame(x)) x else itv_summary(x)
  value_title <- if (value == "R") "R metric" else "Experiment-specific DEGs, %"
  labels <- paste(tab$N, resolve_labels(labels, tab$experiment))
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  ## The label margin is sized to the longest label; the two bar panels share
  ## the rest of the device width equally.
  label_in <- max(graphics::strwidth(labels, units = "inches", cex = 0.8)) + 0.3
  plot_in <- max(graphics::par("din")[1] - label_in, 1)
  graphics::layout(matrix(1:2, 1), widths = c(plot_in / 2, plot_in / 2 + label_in))
  graphics::par(mar = c(3, 1, 4, 0.5))
  graphics::barplot(rev(tab$total), horiz = TRUE, col = "lightgreen",
                    xlim = c(max(tab$total, 1), 0),
                    main = paste0("Number of DEGs\n", main))
  graphics::par(mai = c(graphics::par("mai")[1], label_in, graphics::par("mai")[3:4]))
  b <- graphics::barplot(rev(tab[[value]]), horiz = TRUE, col = "grey",
                         main = paste0(value_title, "\n", main))
  graphics::mtext(rev(labels), side = 2, at = b, las = 2, line = 0.3, cex = 0.8)
  invisible(tab)
}

#' Heatmap of DEG list similarity
#'
#' Rows and columns are ordered by hierarchical clustering of the distances
#' 1 - S between experiments, where S is the similarity shown (or the one
#' in `order_by`, so that several heatmaps can share one order): neighbours
#' are experiments whose DEG lists overlap. Experiments without DEGs are left
#' out with a warning. Each row and column is labelled with the experiment
#' number (its position in `x`) and its label. Raw similarities are drawn
#' from 0 (light) to 1 (dark); chance-corrected ones from -1 (blue) through
#' 0 (white, chance level) to 1 (red), values below -1 drawn as -1. Draws on
#' the current device.
#'
#' @param x An [itv] object, or a matrix from [itv_similarity()].
#' @param metric,chance_corrected,universe,signed Used when `x` is an `itv`
#'   object, see [itv_similarity()]. A matrix from [itv_similarity()] carries
#'   its own metric and correction; for another matrix, `chance_corrected`
#'   sets the colour scale and `metric` the default title.
#' @param hclust_method Agglomeration method for [stats::hclust()].
#' @param main Plot title.
#' @param labels Optional character vector named by experiment, see
#'   [itv_plot_totals()]. By default, the experiment names.
#' @param order_by Optional similarity matrix (with the same experiment
#'   names) whose clustering (by 1 - similarity) sets the order of rows and
#'   columns.
#' @return The value of [stats::heatmap()] invisibly, or `NULL` if fewer than
#'   2 experiments have DEGs.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' itv_plot_heatmap(x)                                # corrected overlap
#' oc <- itv_similarity(x)
#' itv_plot_heatmap(x, chance_corrected = FALSE, order_by = oc)   # same order
#' @export
itv_plot_heatmap <- function(x, metric = c("overlap", "jaccard"), chance_corrected = TRUE,
                             hclust_method = "complete", main = NULL, labels = NULL,
                             universe = c("pair", "all", "none"), signed = FALSE,
                             order_by = NULL) {
  metric <- match.arg(metric)
  universe <- match.arg(universe)
  check_flag(chance_corrected, "chance_corrected")
  check_hclust_method(hclust_method, "hclust_method")
  S <- if (is.matrix(x) && is.numeric(x)) x else
    itv_similarity(x, metric, chance_corrected, universe, signed)
  if (!is.null(attr(S, "metric"))) metric <- attr(S, "metric")
  if (!is.null(attr(S, "chance_corrected"))) chance_corrected <- attr(S, "chance_corrected")
  if (is.null(rownames(S))) dimnames(S) <- rep(list(paste0("exp", seq_len(nrow(S)))), 2)
  exps <- rownames(S)
  lab <- paste(seq_len(nrow(S)), resolve_labels(labels, exps))
  if (is.null(main)) main <- paste0(metric, if (chance_corrected) ", chance-corrected")
  empty <- is.na(diag(S))
  if (any(empty)) {
    warning("Heatmap: experiments without DEGs were left out: ",
            paste(exps[empty], collapse = ", "), call. = FALSE)
    S <- S[!empty, !empty, drop = FALSE]
    lab <- lab[!empty]
    exps <- exps[!empty]
  }
  if (nrow(S) < 2) {
    warning("Heatmap needs at least 2 experiments with DEGs; not drawn.", call. = FALSE)
    return(invisible(NULL))
  }
  O <- S
  if (!is.null(order_by)) {
    if (!is.matrix(order_by) || !all(exps %in% rownames(order_by)) ||
        !all(exps %in% colnames(order_by)))
      stop("'order_by' must be a similarity matrix with the same experiments.",
           call. = FALSE)
    O <- order_by[exps, exps, drop = FALSE]
  }
  ## Cluster the distances 1 - S (a pair without a value gets the largest
  ## distance), then order the branches by mean similarity, as
  ## stats::heatmap() does.
  D <- 1 - unclass(O)[, ]
  D[is.na(D)] <- max(D, na.rm = TRUE)
  diag(D) <- 0
  dend <- stats::reorder(
    stats::as.dendrogram(stats::hclust(stats::as.dist(D), method = hclust_method)),
    rowMeans(O, na.rm = TRUE))
  name <- if (metric == "overlap") "I" else "Jaccard"
  if (chance_corrected) {
    S[!is.na(S) & S < -1] <- -1
    zlim <- c(-1, 1)
    col <- grDevices::hcl.colors(51, "Blue-Red 3")
    key_title <- paste0(name, ", corrected\n(0 = chance)")
  } else {
    zlim <- c(0, 1)
    col <- grDevices::hcl.colors(50, "YlOrRd", rev = TRUE)
    key_title <- name
  }
  dimnames(S) <- list(lab, lab)
  margin <- max(11, 0.6 * max(nchar(lab)))
  ## The title goes in an outer margin sized to its lines, so that a
  ## two-line title does not run into the dendrogram.
  n_lines <- length(strsplit(main, "\n", fixed = TRUE)[[1]])
  op <- graphics::par(oma = c(0, 0, 1.4 * n_lines + 0.6, 0))
  on.exit(graphics::par(op))
  h <- stats::heatmap(S, Rowv = dend, Colv = "Rowv", symm = TRUE, scale = "none",
                      zlim = zlim, col = col, margins = c(margin, margin))
  graphics::mtext(main, side = 3, outer = TRUE, line = 0.4, font = 2, cex = 1.2)
  draw_colour_key(col, zlim, key_title)
  invisible(h)
}

## A horizontal colour key in the empty top-left cell of stats::heatmap(),
## whose layout is 1 + 4 units wide and high (with both dendrograms), kept
## square and centred in the inner region.
draw_colour_key <- function(col, zlim, title) {
  omi <- graphics::par("omi")
  din <- graphics::par("din")
  W <- din[1] - omi[2] - omi[4]
  H <- din[2] - omi[1] - omi[3]
  side <- min(W, H)
  x0 <- (W - side) / 2
  y0 <- (H - side) / 2
  cell <- c(x0, x0 + side / 5, y0 + side * 4 / 5, y0 + side)
  fig <- c(cell[1:2] / W, cell[3:4] / H)
  old <- graphics::par(fig = fig, new = TRUE, mar = c(2.2, 1.2, 2.4, 1.2), mgp = c(3, 0.3, 0),
                       tcl = -0.2, cex = 0.7)
  on.exit(graphics::par(old))
  breaks <- seq(zlim[1], zlim[2], length.out = length(col) + 1)
  graphics::image(x = breaks, y = c(0, 1), z = matrix(utils::head(breaks, -1) + diff(breaks) / 2),
                  col = col, axes = FALSE, xlab = "", ylab = "", zlim = zlim)
  graphics::box()
  graphics::axis(1, at = c(zlim[1], mean(zlim), zlim[2]))
  graphics::mtext(title, side = 3, line = 0.3, cex = 0.7)
}

#' Plot bootstrap intervals
#'
#' One panel per set of experiments, side by side with a shared y axis; k
#' decreases from left to right. Draws on the current device.
#'
#' @param ci A data frame from [itv_ci()], or a named list of them.
#' @param main Title above all panels.
#' @return `ci`, invisibly.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' set.seed(1)
#' all <- itv_ci(itv_bootstrap(x, iter = 200))
#' root <- itv_ci(itv_bootstrap(itv_subset(x, 1:4), iter = 200))
#' itv_plot_ci(list(all = all, root = root), main = "Example")
#' @export
itv_plot_ci <- function(ci, main = "") {
  cis <- if (is.data.frame(ci)) list(ci) else ci
  if (is.null(names(cis))) names(cis) <- rep("", length(cis))
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(matrix(seq_along(cis), 1), widths = vapply(cis, nrow, integer(1)) + 3)
  graphics::par(oma = c(0, 0, 2.5, 0), mar = c(4, 4.5, 3, 1))
  ylim <- range(unlist(lapply(cis, function(d) c(d$lower, d$upper))), 0)
  for (i in seq_along(cis)) {
    d <- cis[[i]]
    nm <- names(cis)[i]
    plot(NA, xlim = c(max(d$k) + 0.5, min(d$k) - 0.5), ylim = ylim, xaxt = "n",
         xlab = "k", ylab = "d = N(m) - N(k)",
         main = sprintf("%s%sm = %s", nm, if (nzchar(nm)) ": " else "",
                        paste(unique(d$m), collapse = ", ")))
    graphics::axis(1, at = d$k)
    graphics::abline(h = 0, col = "red")
    graphics::segments(d$k, d$lower, d$k, d$upper, lwd = 6, col = "orange", lend = 1)
  }
  graphics::mtext(main, outer = TRUE, line = 0.7, font = 2)
  invisible(ci)
}

#' Histogram of one bootstrap comparison
#'
#' @param d Numeric vector, e.g. one column of [itv_bootstrap()] output.
#' @param main Plot title.
#' @param breaks Passed to [graphics::hist()].
#' @return The 2.5% and 97.5% quantiles, invisibly.
#' @export
itv_plot_hist <- function(d, main = "", breaks = 50) {
  q <- stats::quantile(d, c(0.025, 0.975), names = FALSE)
  h <- graphics::hist(d, breaks = breaks, plot = FALSE)
  h$counts <- h$counts / sum(h$counts)
  plot(h, col = "#279AF9", border = "#279AF1", xlab = "d", ylab = "Fraction",
       xlim = range(d, 0),
       main = sprintf("%s\n95%% percentile interval = [%g, %g]", main, q[1], q[2]),
       sub = paste(length(d), "iterations"))
  graphics::abline(v = q, col = "red", lty = "dotted", lwd = 1.5)
  graphics::abline(v = 0, col = "red", lty = "longdash", lwd = 1.5)
  invisible(q)
}
