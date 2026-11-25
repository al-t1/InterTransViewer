#' Bootstrap test of homogeneity
#'
#' Compares the number of genes that are DE in at least one of `size1`
#' experiments drawn from `x` with the same number for `size2` experiments
#' drawn from `y`. Experiments are drawn with replacement. The default
#' compares all m experiments of `x` with k = m-1, ..., 1 experiments of `x`
#' (Tyapkin et al. 2023). If the interval of d = N(m) - N(k) from [itv_ci()]
#' includes 0, adding experiments does not reliably add DEGs.
#'
#' @param x An [itv] object, or a logical matrix of DEG calls.
#' @param y The second set; `x` by default.
#' @param size1 Number of experiments drawn from `x`; all of them by default.
#' @param sizes2 Numbers of experiments drawn from `y`; by default
#'   `size1 - 1` down to 1.
#' @param iter Number of bootstrap iterations for each size.
#' @param progress Show a progress bar.
#' @return A numeric matrix of d, `iter` rows and one column per size in
#'   `sizes2`, named like `"16vs15"`, with attributes `size1` and `sizes2`.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' set.seed(1)
#' b <- itv_bootstrap(x, iter = 200)
#' itv_ci(b)
#' @export
itv_bootstrap <- function(x, y = x, size1 = NULL, sizes2 = NULL, iter = 1000,
                          progress = interactive()) {
  deg1 <- deg_matrix(x)
  deg2 <- deg_matrix(y)
  check_number(iter, "iter", 1, whole = TRUE)
  check_flag(progress, "progress")
  if (!missing(y)) warn_gene_sets(deg1, deg2)
  if (is.null(size1)) size1 <- ncol(deg1)
  if (is.null(sizes2)) {
    if (size1 < 2) stop("The default sizes need at least 2 experiments.", call. = FALSE)
    sizes2 <- (size1 - 1):1
  }
  check_number(size1, "size1", 1, whole = TRUE)
  if (!is.numeric(sizes2) || length(sizes2) == 0 || anyNA(sizes2) ||
      any(sizes2 != round(sizes2)))
    stop("'sizes2' must be whole numbers.", call. = FALSE)
  if (any(sizes2 < 1)) stop("Sample sizes must be >= 1.", call. = FALSE)

  ## Each experiment as the row numbers of its DEGs; N is the size of their union.
  idx1 <- lapply(seq_len(ncol(deg1)), function(j) which(deg1[, j]))
  idx2 <- lapply(seq_len(ncol(deg2)), function(j) which(deg2[, j]))
  N <- function(idx, v, n) {
    hit <- logical(n)
    hit[unlist(idx[unique(v)], use.names = FALSE)] <- TRUE
    sum(hit)
  }
  out <- matrix(NA_real_, iter, length(sizes2),
                dimnames = list(NULL, paste0(size1, "vs", sizes2)))
  if (progress) pb <- utils::txtProgressBar(max = length(out), style = 3)
  for (s in seq_along(sizes2)) {
    for (b in seq_len(iter)) {
      v1 <- sample.int(ncol(deg1), size1, replace = TRUE)
      v2 <- sample.int(ncol(deg2), sizes2[s], replace = TRUE)
      out[b, s] <- N(idx1, v1, nrow(deg1)) - N(idx2, v2, nrow(deg2))
      if (progress) utils::setTxtProgressBar(pb, (s - 1) * iter + b)
    }
  }
  if (progress) close(pb)
  attr(out, "size1") <- size1
  attr(out, "sizes2") <- sizes2
  out
}

## The size of a union of DEG lists grows with the number of genes tested,
## so comparing sets measured on different genes is confounded.
warn_gene_sets <- function(deg1, deg2) {
  g1 <- rownames(deg1)
  g2 <- rownames(deg2)
  same <- if (is.null(g1) || is.null(g2)) nrow(deg1) == nrow(deg2) else setequal(g1, g2)
  if (!same)
    warning("x and y do not have the same genes (", nrow(deg1), " and ", nrow(deg2),
            if (!is.null(g1) && !is.null(g2))
              paste0(", ", length(intersect(g1, g2)), " shared"),
            "). The number of DEGs in a union depends on how many genes were ",
            "tested, so d also reflects this difference. Consider keeping only the ",
            "shared genes in both.", call. = FALSE)
}

#' Percentile intervals of bootstrap results
#'
#' @param boot A matrix from [itv_bootstrap()].
#' @param probs Lower and upper quantiles.
#' @return A data frame with one row per column of `boot`: `comparison`,
#'   `m` and `k` (sample sizes), `lower`, `upper`, and `includes_zero`.
#' @export
itv_ci <- function(boot, probs = c(0.025, 0.975)) {
  boot <- as.matrix(boot)
  q <- t(apply(boot, 2, stats::quantile, probs = probs, names = FALSE))
  m <- attr(boot, "size1")
  k <- attr(boot, "sizes2")
  if (is.null(m) || is.null(k)) {
    parts <- strsplit(colnames(boot), "vs", fixed = TRUE)
    m <- suppressWarnings(as.integer(vapply(parts, `[`, "", 1)))
    k <- suppressWarnings(as.integer(vapply(parts, function(p) p[2], "")))
  }
  data.frame(comparison = colnames(boot), m = m, k = k,
             lower = q[, 1], upper = q[, 2],
             includes_zero = q[, 1] <= 0 & q[, 2] >= 0, row.names = NULL)
}
