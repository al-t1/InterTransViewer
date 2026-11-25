#' Cluster experiments by logFC
#'
#' Hierarchical clustering of experiments by the logFC of genes measured in
#' all of them, with bootstrap support for each clade.
#'
#' Only genes with a finite logFC and a known FDR in every experiment are
#' used, optionally only those that are DE in at least `min_exps`
#' experiments. Each experiment is divided by the range of its logFC, then
#' each gene by its standard deviation across experiments; genes with zero
#' variance are dropped. An experiment whose logFC does not vary (for
#' example, one without DEGs when `insignif_to_zero = TRUE`) is left out with
#' a warning. Euclidean distances are clustered with [stats::hclust()].
#' Clade support is the fraction of `boot_iter` trees, built from genes
#' resampled with replacement, that contain the clade.
#'
#' @param x An [itv] object.
#' @param min_exps Use genes that are DE in at least this many experiments.
#' @param insignif_to_zero If `TRUE`, logFC values that are not significant
#'   are set to 0 first.
#' @param method Agglomeration method for [stats::hclust()].
#' @param boot_iter Number of bootstrap replicates for clade support.
#' @return An object of class `itv_cluster`: a list with `tree` (an
#'   `hclust` object), `experiments` (names of the experiments used; leaf i
#'   of the tree is `experiments[i]`), `support` (one value per merge, in the order of
#'   `tree$merge`), `n_genes` (genes used) and `n_complete` (genes measured in
#'   all experiments).
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' set.seed(1)
#' cl <- itv_cluster(itv_read(f), boot_iter = 20)
#' cl
#' plot(cl)
#' @export
itv_cluster <- function(x, min_exps = 0, insignif_to_zero = FALSE,
                        method = "ward.D2", boot_iter = 100) {
  check_itv(x)
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

  S <- t(X)
  S <- S / rng[!flat]
  sds <- apply(S, 2, stats::sd)
  ok <- is.finite(sds) & sds > 0
  S <- sweep(S[, ok, drop = FALSE], 2, sds[ok], "/")
  if (ncol(S) < 2) stop("Too few genes left for clustering (", ncol(S), ").", call. = FALSE)

  tree <- stats::hclust(stats::dist(S), method = method)
  structure(list(tree = tree, experiments = rownames(S),
                 support = clade_support(S, tree, method, boot_iter),
                 n_genes = ncol(S), n_complete = sum(complete),
                 min_exps = min_exps, boot_iter = boot_iter),
            class = "itv_cluster")
}

#' @export
print.itv_cluster <- function(x, ...) {
  cat("<itv_cluster> ", length(x$experiments), " experiments, ", x$tree$method,
      " on ", x$n_genes, " genes, ", x$boot_iter, " bootstrap replicates\n", sep = "")
  cat("Clade support, %: ", paste(round(100 * x$support), collapse = " "), "\n", sep = "")
  invisible(x)
}

#' Plot a clustering of experiments
#'
#' Draws the dendrogram with the bootstrap support of each clade, in
#' percent, at its node.
#'
#' @param x An `itv_cluster` object from [itv_cluster()].
#' @param labels Leaf labels, one per experiment in `x$experiments`.
#' @param main Plot title.
#' @param support_col Colour of the support values; `NA` hides them.
#' @param ... Passed to [plot()] for the dendrogram.
#' @return `x`, invisibly.
#' @export
plot.itv_cluster <- function(x, labels = x$experiments, main = NULL,
                             support_col = "firebrick", ...) {
  if (length(labels) != length(x$experiments))
    stop("labels must have one entry per experiment in x$experiments.", call. = FALSE)
  if (is.null(main))
    main <- sprintf("%s on %d genes; support, %%, from %d replicates",
                    x$tree$method, x$n_genes, x$boot_iter)
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

## Each clade as a string of sorted leaf indices, one per merge step.
clade_keys <- function(tree) {
  m <- tree$merge
  members <- vector("list", nrow(m))
  side <- function(i, k) if (m[i, k] < 0) -m[i, k] else members[[m[i, k]]]
  for (i in seq_len(nrow(m))) members[[i]] <- sort(c(side(i, 1), side(i, 2)))
  vapply(members, paste, character(1), collapse = " ")
}

## Fraction of bootstrap trees (genes resampled) that contain each clade.
clade_support <- function(S, tree, method, iter) {
  ref <- clade_keys(tree)
  hits <- numeric(length(ref))
  for (b in seq_len(iter)) {
    Sb <- S[, sample.int(ncol(S), replace = TRUE), drop = FALSE]
    hits <- hits + (ref %in% clade_keys(stats::hclust(stats::dist(Sb), method = method)))
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
