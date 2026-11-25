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
#' Rows and columns are ordered by hierarchical clustering of the rows of the
#' similarity matrix (Euclidean distance), as in InterTransViewer 0.1.
#' Experiments without DEGs are left out with a warning. Each row and column
#' is labelled with the experiment number (its position in `x`) and its
#' label. Draws on the current device.
#'
#' @param x An [itv] object, or a matrix from [itv_similarity()].
#' @param metric Used when `x` is an `itv` object: `"overlap"` or `"jaccard"`.
#' @param hclust_method Agglomeration method for [stats::hclust()].
#' @param main Plot title.
#' @param labels Optional character vector named by experiment, see
#'   [itv_plot_totals()]. By default, the experiment names.
#' @return The value of [stats::heatmap()] invisibly, or `NULL` if fewer than
#'   2 experiments have DEGs.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' itv_plot_heatmap(itv_read(f))
#' @export
itv_plot_heatmap <- function(x, metric = c("overlap", "jaccard"),
                             hclust_method = "complete", main = NULL, labels = NULL) {
  metric <- match.arg(metric)
  check_hclust_method(hclust_method, "hclust_method")
  S <- if (is.matrix(x) && is.numeric(x)) x else itv_similarity(x, metric)
  if (is.null(rownames(S))) dimnames(S) <- rep(list(paste0("exp", seq_len(nrow(S)))), 2)
  lab <- paste(seq_len(nrow(S)), resolve_labels(labels, rownames(S)))
  if (is.null(main)) main <- metric
  empty <- is.na(diag(S))
  if (any(empty)) {
    warning("Heatmap: experiments without DEGs were left out: ",
            paste(rownames(S)[empty], collapse = ", "), call. = FALSE)
    S <- S[!empty, !empty, drop = FALSE]
    lab <- lab[!empty]
  }
  if (nrow(S) < 2) {
    warning("Heatmap needs at least 2 experiments with DEGs; not drawn.", call. = FALSE)
    return(invisible(NULL))
  }
  dimnames(S) <- list(lab, lab)
  margin <- max(11, 0.6 * max(nchar(lab)))
  invisible(stats::heatmap(S, symm = TRUE, scale = "none", zlim = c(0, 1),
                           hclustfun = function(d) stats::hclust(d, method = hclust_method),
                           col = grDevices::hcl.colors(50, "YlOrRd", rev = TRUE),
                           margins = c(margin, margin), main = main,
                           xlab = "colour: 0 (light) to 1 (dark)"))
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
