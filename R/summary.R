#' DEG counts, specific DEGs and the R metric
#'
#' For each experiment: `total`, the number of DEGs (n); `specific`, the
#' number of DEGs found in no other experiment; `specific_percent`, their
#' share (100 delta); and `R = delta * N / n`, where N is the number of genes
#' that are DE in at least one experiment. `specific_percent` and `R` are
#' `NA` for an experiment without DEGs.
#'
#' @param x An [itv] object, or a logical genes x experiments matrix of DEG calls.
#' @return A data frame with columns `N` (experiment number), `experiment`,
#'   `total`, `specific`, `specific_percent` and `R`.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' itv_summary(itv_read(f))
#' @export
itv_summary <- function(x) {
  deg <- deg_matrix(x)
  n_exp_per_gene <- rowSums(deg)
  N_all    <- sum(n_exp_per_gene > 0)
  total    <- colSums(deg)
  specific <- colSums(deg & n_exp_per_gene == 1)
  frac     <- ifelse(total > 0, specific / total, NA_real_)
  data.frame(N = seq_along(total), experiment = colnames(deg),
             total = unname(total), specific = unname(specific),
             specific_percent = unname(round(100 * frac, 2)),
             R = unname(frac / (total / N_all)),
             row.names = NULL)
}

#' Genes that are DE in one experiment only
#'
#' @param x An [itv] object.
#' @return A data frame with columns `ID`, `experiment` and `logFC`, ordered
#'   by experiment.
#' @export
itv_specific <- function(x) {
  check_itv(x)
  idx <- which(x$deg & rowSums(x$deg) == 1, arr.ind = TRUE)
  data.frame(ID = rownames(x$deg)[idx[, 1]], experiment = colnames(x$deg)[idx[, 2]],
             logFC = x$logFC[idx], row.names = NULL)
}

#' Significant logFC values
#'
#' @param x An [itv] object.
#' @return A matrix of the genes that are DE in at least one experiment,
#'   holding the logFC where the gene is DE and `NA` elsewhere. A DEG without
#'   a logFC is also `NA`; use `x$deg` for the calls themselves.
#' @export
itv_signif <- function(x) {
  check_itv(x)
  s <- x$logFC
  s[!x$deg] <- NA
  s[rowSums(x$deg) > 0, , drop = FALSE]
}

#' Number of experiments in which each gene is DE
#'
#' @param x An [itv] object, or a logical matrix of DEG calls with gene IDs
#'   as row names.
#' @return A data frame with columns `ID` and `n_experiments`, for the genes
#'   that are DE in at least one experiment.
#' @export
itv_joint_degs <- function(x) {
  deg <- deg_matrix(x)
  n <- rowSums(deg)
  data.frame(ID = rownames(deg)[n > 0], n_experiments = unname(n[n > 0]),
             row.names = NULL)
}

#' Pairwise similarity of DEG lists
#'
#' With c the number of shared DEGs, the overlap index (the similarity
#' index I of Tyapkin et al. 2023) is c / min(|A|, |B|), and the Jaccard
#' index is c / |A or B|. Pairs involving an experiment without DEGs are `NA`.
#'
#' @param x An [itv] object, or a logical matrix of DEG calls.
#' @param metric `"overlap"` or `"jaccard"`.
#' @return A symmetric experiments x experiments matrix.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' round(itv_similarity(itv_read(f)), 2)
#' @export
itv_similarity <- function(x, metric = c("overlap", "jaccard")) {
  metric <- match.arg(metric)
  d <- deg_matrix(x) * 1
  inter <- crossprod(d)
  n <- diag(inter)
  den <- if (metric == "overlap") outer(n, n, pmin) else outer(n, n, "+") - inter
  S <- inter / den
  S[den == 0] <- NA
  S
}
