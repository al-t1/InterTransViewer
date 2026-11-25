#' Plot DEG counts with specific DEGs or R
#'
#' Two panels: the number of DEGs per experiment, and either the percent of
#' experiment-specific DEGs or the R metric. Draws on the current device.
#'
#' @param x An [itv] object, or a data frame from [itv_summary()].
#' @param value `"specific_percent"` or `"R"`.
#' @param main Title shown above both panels.
#' @return The summary data frame, invisibly.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' itv_plot_totals(itv_read(f), "R", main = "Example")
#' @export
itv_plot_totals <- function(x, value = c("specific_percent", "R"), main = "") {
  value <- match.arg(value)
  tab <- if (is.data.frame(x)) x else itv_summary(x)
  value_title <- if (value == "R") "R metric" else "Experiment-specific DEGs, %"
  op <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(op))
  graphics::layout(matrix(1:2, 1), widths = c(1, 1.4))
  graphics::par(mar = c(3, 1, 4, 0.5))
  graphics::barplot(rev(tab$total), horiz = TRUE, col = "lightgreen",
                    xlim = c(max(tab$total, 1), 0),
                    main = paste0("Number of DEGs\n", main))
  graphics::par(mar = c(3, 13, 4, 1))
  b <- graphics::barplot(rev(tab[[value]]), horiz = TRUE, col = "grey",
                         main = paste0(value_title, "\n", main))
  graphics::mtext(rev(paste(tab$N, tab$experiment)), side = 2, at = b,
                  las = 2, line = 0.3, cex = 0.8)
  invisible(tab)
}

#' Heatmap of DEG list similarity
#'
#' Rows and columns are ordered by hierarchical clustering of the rows of the
#' similarity matrix (Euclidean distance), as in InterTransViewer 0.1.
#' Experiments without DEGs are left out with a warning. Draws on the
#' current device.
#'
#' @param x An [itv] object, or a matrix from [itv_similarity()].
#' @param metric Used when `x` is an `itv` object: `"overlap"` or `"jaccard"`.
#' @param hclust_method Agglomeration method for [stats::hclust()].
#' @param main Plot title.
#' @return The value of [stats::heatmap()] invisibly, or `NULL` if fewer than
#'   2 experiments have DEGs.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' itv_plot_heatmap(itv_read(f))
#' @export
itv_plot_heatmap <- function(x, metric = c("overlap", "jaccard"),
                             hclust_method = "complete", main = NULL) {
  metric <- match.arg(metric)
  S <- if (is.matrix(x) && is.numeric(x)) x else itv_similarity(x, metric)
  if (is.null(main)) main <- metric
  empty <- is.na(diag(S))
  if (any(empty)) {
    warning("Heatmap: experiments without DEGs were left out: ",
            paste(rownames(S)[empty], collapse = ", "), call. = FALSE)
    S <- S[!empty, !empty, drop = FALSE]
  }
  if (nrow(S) < 2) {
    warning("Heatmap needs at least 2 experiments with DEGs; not drawn.", call. = FALSE)
    return(invisible(NULL))
  }
  invisible(stats::heatmap(S, symm = TRUE, scale = "none", zlim = c(0, 1),
                           hclustfun = function(d) stats::hclust(d, method = hclust_method),
                           col = grDevices::hcl.colors(50, "YlOrRd", rev = TRUE),
                           margins = c(11, 11), main = main,
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
