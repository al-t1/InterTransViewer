#' Cluster experiments by logFC
#'
#' Hierarchical clustering of experiments by the logFC of genes measured in
#' all of them, with bootstrap support for each clade.
#'
#' Only genes with a finite logFC and a known FDR in every experiment are
#' used, optionally only those that are DE in at least `min_exps`
#' experiments. The logFC values are then, in this order and each if asked:
#' \enumerate{
#'   \item divided by the range of each experiment's logFC (`scale_exps_by_range`),
#'     so that experiments with stronger responses overall do not dominate;
#'   \item centred: each gene's mean across experiments is subtracted
#'     (`center_genes`);
#'   \item divided by each gene's standard deviation across experiments
#'     (`scale_genes_by_SD`); genes with zero variance are dropped.
#' }
#' With `center_genes` and `scale_genes_by_SD` each gene is z-scored. The defaults
#' (`scale_exps_by_range` and `scale_genes_by_SD`, no centring) are the method of
#' InterTransViewer 0.1 and the 2023 paper.
#'
#' An experiment whose logFC does not vary (for example, one without DEGs
#' when `insignif_to_zero = TRUE`) is left out with a warning. Distances
#' between experiments are Euclidean, or 1 - r for Pearson or Spearman
#' correlation across genes; correlation ignores `scale_exps_by_range`. The tree is
#' built with [stats::hclust()].
#'
#' Clade support is the fraction of `boot_iter` trees, built from genes
#' resampled with replacement, that contain the clade. Genes are not
#' independent (they respond in groups), so this support overstates how sure
#' the clustering is; read it as relative, not as a probability.
#'
#' @param x An [itv] object.
#' @param min_exps Use genes that are DE in at least this many experiments.
#' @param insignif_to_zero If `TRUE`, logFC values that are not significant
#'   are set to 0 first.
#' @param scale_exps_by_range Divide each experiment's logFC by its range.
#' @param center_genes Subtract each gene's mean across experiments.
#' @param scale_genes_by_SD Divide each gene by its standard deviation across
#'   experiments.
#' @param distance `"euclidean"`, `"pearson"` (1 - Pearson r) or
#'   `"spearman"` (1 - Spearman rho).
#' @param method Agglomeration method for [stats::hclust()].
#' @param boot_iter Number of bootstrap replicates for clade support.
#' @return An object of class `itv_cluster`: a list with `tree` (an
#'   `hclust` object), `experiments` (names of the experiments used; leaf i
#'   of the tree is `experiments[i]`), `numbers` (their positions in `x`),
#'   `support` (one value per merge, in the order of
#'   `tree$merge`), `n_genes` (genes used), `n_complete` (genes measured in
#'   all experiments), `n_total` (genes in `x`) and the settings.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' set.seed(1)
#' cl <- itv_cluster(x, boot_iter = 20)
#' cl
#' plot(cl)
#' ## Genes z-scored, 1 - Pearson r
#' plot(itv_cluster(x, center_genes = TRUE, distance = "pearson", boot_iter = 20))
#' @export
itv_cluster <- function(x, min_exps = 0, insignif_to_zero = FALSE,
                        scale_exps_by_range = TRUE, center_genes = FALSE, scale_genes_by_SD = TRUE,
                        distance = c("euclidean", "pearson", "spearman"),
                        method = "ward.D2", boot_iter = 100) {
  check_itv(x)
  check_number(min_exps, "min_exps", 0, whole = TRUE)
  check_flag(insignif_to_zero, "insignif_to_zero")
  check_flag(scale_exps_by_range, "scale_exps_by_range")
  check_flag(center_genes, "center_genes")
  check_flag(scale_genes_by_SD, "scale_genes_by_SD")
  distance <- match.arg(distance)
  check_hclust_method(method, "method")
  check_number(boot_iter, "boot_iter", 1, whole = TRUE)
  complete <- rowSums(!is.finite(x$logFC) | is.na(x$FDR)) == 0
  keep <- complete & rowSums(x$deg) >= min_exps
  X <- x$logFC[keep, , drop = FALSE]
  if (insignif_to_zero) X[!x$deg[keep, , drop = FALSE]] <- 0
  if (nrow(X) < 2) stop("Too few genes left for clustering (", nrow(X), ").", call. = FALSE)

  rng <- apply(X, 2, function(v) diff(range(v)))
  flat <- !(is.finite(rng) & rng > 0)
  if (any(flat)) {
    warning("Clustering: experiments with no logFC variation were left out: ",
            paste(colnames(X)[flat], collapse = ", "), call. = FALSE)
    X <- X[, !flat, drop = FALSE]
  }
  if (ncol(X) < 2) stop("Clustering needs at least 2 experiments.", call. = FALSE)

  S <- t(X)                                     # experiments x genes
  if (scale_exps_by_range) S <- S / rng[!flat]
  if (center_genes) S <- sweep(S, 2, colMeans(S), "-")
  if (scale_genes_by_SD) {
    sds <- apply(S, 2, stats::sd)
    ok <- is.finite(sds) & sds > 0
    S <- sweep(S[, ok, drop = FALSE], 2, sds[ok], "/")
  }
  if (ncol(S) < 2) stop("Too few genes left for clustering (", ncol(S), ").", call. = FALSE)

  dist_fun <- cluster_distance(distance)
  tree <- stats::hclust(dist_fun(S), method = method)
  structure(list(tree = tree, experiments = rownames(S),
                 numbers = match(rownames(S), colnames(x$deg)),
                 support = clade_support(S, tree, method, boot_iter, dist_fun),
                 n_genes = ncol(S), n_complete = sum(complete), n_total = nrow(x$deg),
                 min_exps = min_exps, boot_iter = boot_iter, distance = distance,
                 scale_exps_by_range = scale_exps_by_range, center_genes = center_genes,
                 scale_genes_by_SD = scale_genes_by_SD),
            class = "itv_cluster")
}

## Distances between the rows of an experiments x genes matrix.
cluster_distance <- function(distance) {
  if (distance == "euclidean") return(stats::dist)
  function(S) {
    r <- suppressWarnings(stats::cor(t(S), method = distance))
    r[is.na(r)] <- 0          # an experiment constant on these genes
    stats::as.dist(1 - r)
  }
}

#' @export
print.itv_cluster <- function(x, ...) {
  cat("<itv_cluster> ", length(x$experiments), " experiments, ", x$tree$method, " on ",
      x$n_genes, " genes, ", cluster_distance_name(x), " distance, ", x$boot_iter,
      " bootstrap replicates\n", sep = "")
  if (!is.null(x$n_total))
    cat("Genes measured in all experiments: ", x$n_complete, " of ", x$n_total, "\n",
        sep = "")
  if (!is.null(x$scale_exps_by_range))
    cat("Scaling: ", cluster_scaling_text(x), "\n", sep = "")
  cat("Clade support, %: ", paste(round(100 * x$support), collapse = " "), "\n", sep = "")
  invisible(x)
}

#' Plot a clustering of experiments
#'
#' Draws the dendrogram with the bootstrap support of each clade, in
#' percent, at its node.
#'
#' @param x An `itv_cluster` object from [itv_cluster()].
#' @param labels Optional character vector named by experiment, see
#'   [itv_plot_totals()]. Each leaf is labelled with the experiment number
#'   (its position in the data that was clustered) and this label (by
#'   default, the experiment name).
#' @param main Plot title.
#' @param support_col Colour of the support values; `NA` hides them.
#' @param ... Passed to [plot()] for the dendrogram.
#' @return `x`, invisibly.
#' @export
plot.itv_cluster <- function(x, labels = NULL, main = NULL,
                             support_col = "firebrick", ...) {
  numbers <- if (is.null(x$numbers)) seq_along(x$experiments) else x$numbers
  labels <- paste(numbers, resolve_labels(labels, x$experiments))
  if (is.null(main))
    main <- sprintf("%s, %s distance, on %d genes%s\nsupport, %%, from %d replicates",
                    x$tree$method, cluster_distance_name(x), x$n_genes,
                    if (is.null(x$n_total)) "" else
                      sprintf(" (%d of %d measured in all)", x$n_complete, x$n_total),
                    x$boot_iter)
  op <- graphics::par(mar = c(max(6, 0.6 * max(nchar(labels))), 4, 5, 1), cex = 0.85)
  on.exit(graphics::par(op))
  tree <- x$tree
  tree$labels <- labels
  plot(stats::as.dendrogram(tree), main = main, ylab = "Height", ...)
  if (!is.na(support_col))
    graphics::text(node_x(x$tree), x$tree$height, round(100 * x$support),
                   pos = 3, cex = 0.75, col = support_col, xpd = NA)
  invisible(x)
}

cluster_distance_name <- function(x) if (is.null(x$distance)) "euclidean" else x$distance

cluster_scaling_text <- function(x) {
  steps <- c(if (x$scale_exps_by_range) "experiments / range", if (x$center_genes) "genes centred",
             if (x$scale_genes_by_SD) "genes / SD")
  if (length(steps)) paste(steps, collapse = ", ") else "none"
}

## Each clade as a string of sorted leaf indices, one per merge step.
clade_keys <- function(tree) {
  m <- tree$merge
  members <- vector("list", nrow(m))
  side <- function(i, k) if (m[i, k] < 0) -m[i, k] else members[[m[i, k]]]
  for (i in seq_len(nrow(m))) members[[i]] <- sort(c(side(i, 1), side(i, 2)))
  vapply(members, paste, character(1), collapse = " ")
}

## Fraction of bootstrap trees (genes resampled) that contain each clade.
clade_support <- function(S, tree, method, iter, dist_fun = stats::dist) {
  ref <- clade_keys(tree)
  hits <- numeric(length(ref))
  for (b in seq_len(iter)) {
    Sb <- S[, sample.int(ncol(S), replace = TRUE), drop = FALSE]
    hits <- hits + (ref %in% clade_keys(stats::hclust(dist_fun(Sb), method = method)))
  }
  hits / iter
}

## x position of each merge node, as drawn by plot.hclust and plot.dendrogram.
node_x <- function(tree) {
  leaf_x <- match(seq_along(tree$order), tree$order)
  m <- tree$merge
  x <- numeric(nrow(m))
  pos <- function(v) if (v < 0) leaf_x[-v] else x[v]
  for (i in seq_len(nrow(m))) x[i] <- (pos(m[i, 1]) + pos(m[i, 2])) / 2
  x
}
