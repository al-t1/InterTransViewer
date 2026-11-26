#' DEG counts, specific DEGs and the R metric
#'
#' For each experiment:
#' \describe{
#'   \item{`total`}{n, the number of DEGs;}
#'   \item{`specific`}{the DEGs that are DE in no other experiment;}
#'   \item{`specific_percent`}{their share, 100 delta;}
#'   \item{`specific_tested`, `specific_untested`}{`specific` split by
#'     whether every other experiment tested the gene (it has an FDR there)
#'     and found it not DE, or at least one other experiment never tested it.
#'     An untested gene, e.g. one that is not on a microarray, may only look
#'     specific;}
#'   \item{`R`}{delta N / n, where N is the number of genes that are DE in
#'     at least one experiment;}
#'   \item{`specific_oe`}{observed / expected specific DEGs: `specific_tested`
#'     divided by the number expected if the DEG lists were independent,
#'     \eqn{m \prod_{j \ne A} (1 - n_j / G_j)}, with m the experiment's DEGs
#'     tested in every other experiment, and \eqn{n_j} and \eqn{G_j} the
#'     DEGs and genes tested of each other experiment. Above 1: more
#'     specific DEGs than chance.}
#' }
#' `specific_percent`, `R` and `specific_oe` are `NA` for an experiment
#' without DEGs.
#'
#' @param x An [itv] object, or a logical genes x experiments matrix of DEG
#'   calls (then every gene counts as tested in every experiment).
#' @return A data frame with columns `N` (experiment number), `experiment`,
#'   `total`, `specific`, `specific_percent`, `specific_tested`,
#'   `specific_untested`, `R` and `specific_oe`.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' itv_summary(itv_read(f))
#' @export
itv_summary <- function(x) {
  deg <- deg_matrix(x)
  meas <- measured_matrix(x, deg)
  k <- ncol(deg)
  n_exp_per_gene <- rowSums(deg)
  everywhere <- rowSums(meas) == k
  N_all    <- sum(n_exp_per_gene > 0)
  total    <- colSums(deg)
  specific <- colSums(deg & n_exp_per_gene == 1)
  tested   <- colSums(deg & n_exp_per_gene == 1 & everywhere)
  frac     <- ifelse(total > 0, specific / total, NA_real_)

  ## Expected specific DEGs if the lists were independent.
  p_de <- total / (colSums(meas) + extra_genes(x))
  m <- colSums(deg & everywhere)
  expected <- m * vapply(seq_len(k), function(a) prod(1 - p_de[-a]), numeric(1))
  oe <- ifelse(total > 0 & expected > 0, tested / expected, NA_real_)

  data.frame(N = seq_along(total), experiment = colnames(deg),
             total = unname(total), specific = unname(specific),
             specific_percent = unname(round(100 * frac, 2)),
             specific_tested = unname(tested),
             specific_untested = unname(specific - tested),
             R = unname(frac / (total / N_all)),
             specific_oe = unname(oe),
             row.names = NULL)
}

#' Genes that are DE in one experiment only
#'
#' @param x An [itv] object.
#' @return A data frame with columns `ID`, `experiment`, `logFC` and
#'   `tested` (`TRUE` if every other experiment tested the gene; see
#'   [itv_summary()]), ordered by experiment.
#' @export
itv_specific <- function(x) {
  check_itv(x)
  idx <- which(x$deg & rowSums(x$deg) == 1, arr.ind = TRUE)
  everywhere <- rowSums(x$measured) == ncol(x$deg)
  data.frame(ID = rownames(x$deg)[idx[, 1]], experiment = colnames(x$deg)[idx[, 2]],
             logFC = x$logFC[idx], tested = unname(everywhere[idx[, 1]]),
             row.names = NULL)
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

#' Similarity of DEG lists
#'
#' Compares the DEG lists of each pair of experiments A and B, on the genes
#' tested in both (by default). With c the number of shared DEGs:
#' \describe{
#'   \item{overlap}{the similarity index I of Tyapkin et al. (2023),
#'     c / min(|A|, |B|): the Szymkiewicz-Simpson overlap coefficient;}
#'   \item{jaccard}{c / |A or B|.}
#' }
#' With `chance_corrected = TRUE` (the default), each value S becomes
#' (S - E) / (1 - E), where E is the value expected if the two lists were
#' drawn at random from the G genes tested in both: 0 means chance level, 1
#' the largest similarity possible, negative less overlap than chance.
#'
#' By chance, c is hypergeometric with mean \eqn{E[c] = |A| |B| / G}, so the
#' expected overlap is max(|A|, |B|) / G: two long lists overlap a lot even
#' without any shared biology. The expected Jaccard is approximated by
#' \eqn{E[c] / (|A| + |B| - E[c])}. These expectations treat all genes as
#' equally likely to be DE, so they are a lower bar: genes that respond in
#' many conditions raise the real overlap.
#'
#' Chance correction matters when lists are long compared with the genes
#' tested; for short lists E is small and the corrected values are close to
#' the raw ones.
#'
#' `itv_similarity()` gives one metric as a matrix; `itv_pairs()` gives all
#' of them, raw and corrected, with the expected values and a
#' hypergeometric p-value, as one table.
#'
#' @param x An [itv] object, or a logical matrix of DEG calls (then every
#'   gene counts as tested everywhere, and `signed` is not available).
#' @param metric `"overlap"` (index I) or `"jaccard"`.
#' @param chance_corrected Correct the metric for chance, see Details.
#' @param universe The genes each pair is compared on: `"pair"`, those
#'   tested (with an FDR) in both experiments; `"all"`, those tested in
#'   every experiment; `"none"`, all genes, an untested gene counting as not
#'   DE. Genes outside the table (see `n_genes` in [itv()]) count as tested
#'   everywhere.
#' @param signed If `TRUE`, a gene is shared only if it changes in the same
#'   direction in both experiments (a DEG without a logFC is then never
#'   shared), and the expectation is \eqn{(up_A up_B + down_A down_B) / G}.
#' @return `itv_similarity()`: a symmetric experiments x experiments matrix,
#'   with attributes `metric` and `chance_corrected`; pairs with an
#'   experiment without DEGs are `NA`. `itv_pairs()`: a data
#'   frame with one row per pair: `experiment_a`, `experiment_b`, `genes`
#'   (G), `n_a`, `n_b` (DEGs among the G genes), `shared`, `shared_same`,
#'   `shared_opposite` (shared DEGs, and those changing in the same and in
#'   opposite directions), `overlap`, `overlap_expected`,
#'   `overlap_corrected`, `jaccard`, `jaccard_expected`,
#'   `jaccard_corrected`, `enrichment` (shared / expected shared) and
#'   `p_value` (hypergeometric, one-sided, for `shared`; `NA` if
#'   `signed = TRUE`).
#' @references Simpson G.G. (1943) Mammals and the nature of continents.
#'   American Journal of Science 241:1-31.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' round(itv_similarity(x), 2)                            # corrected overlap
#' round(itv_similarity(x, chance_corrected = FALSE), 2)  # index I as published
#' round(itv_similarity(x, "jaccard"), 2)
#' head(itv_pairs(x))
#' head(itv_pairs(x, signed = TRUE))
#' @export
itv_similarity <- function(x, metric = c("overlap", "jaccard"), chance_corrected = TRUE,
                           universe = c("pair", "all", "none"), signed = FALSE) {
  metric <- match.arg(metric)
  universe <- match.arg(universe)
  check_flag(chance_corrected, "chance_corrected")
  S <- pair_stats(x, universe, signed)[[sim_key(metric, chance_corrected)]]
  attr(S, "metric") <- metric
  attr(S, "chance_corrected") <- chance_corrected
  S
}

sim_key <- function(metric, chance_corrected) {
  if (chance_corrected) paste0(metric, "_corrected") else metric
}

#' @rdname itv_similarity
#' @export
itv_pairs <- function(x, universe = c("pair", "all", "none"), signed = FALSE) {
  universe <- match.arg(universe)
  s <- pair_stats(x, universe, signed)
  exps <- colnames(s$overlap)
  idx <- which(upper.tri(s$overlap), arr.ind = TRUE)
  idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
  get <- function(m) unname(m[idx])
  data.frame(experiment_a = exps[idx[, 1]], experiment_b = exps[idx[, 2]],
             genes = get(s$genes), n_a = get(s$n_a), n_b = get(s$n_b),
             shared = get(s$shared), shared_same = get(s$shared_same),
             shared_opposite = get(s$shared_opposite),
             overlap = get(s$overlap), overlap_expected = get(s$overlap_expected),
             overlap_corrected = get(s$overlap_corrected),
             jaccard = get(s$jaccard), jaccard_expected = get(s$jaccard_expected),
             jaccard_corrected = get(s$jaccard_corrected),
             enrichment = get(s$enrichment), p_value = get(s$p_value),
             row.names = NULL)
}

## All pairwise counts and similarities, as experiments x experiments
## matrices. In a matrix m, m[a, b] refers to the genes tested in both a and
## b; n_a[a, b] is the number of DEGs of a among them, n_b[a, b] those of b.
pair_stats <- function(x, universe, signed) {
  check_flag(signed, "signed")
  deg <- deg_matrix(x)
  meas <- measured_matrix(x, deg)
  lfc <- if (inherits(x, "itv")) x$logFC else NULL
  if (signed && is.null(lfc))
    stop("signed = TRUE needs an itv object with logFC values.", call. = FALSE)
  if (universe == "all") {
    rows <- rowSums(meas) == ncol(deg)
    deg <- deg[rows, , drop = FALSE]; meas <- meas[rows, , drop = FALSE]
    if (!is.null(lfc)) lfc <- lfc[rows, , drop = FALSE]
  } else if (universe == "none") {
    meas[] <- TRUE
  }
  D <- deg * 1
  M <- meas * 1
  G <- crossprod(M) + extra_genes(x)
  n_a <- crossprod(D, M)
  n_b <- t(n_a)
  shared <- crossprod(D)
  if (!is.null(lfc)) {
    up <- (deg & !is.na(lfc) & lfc > 0) * 1
    down <- (deg & !is.na(lfc) & lfc < 0) * 1
    same <- crossprod(up) + crossprod(down)
    opposite <- crossprod(up, down) + crossprod(down, up)
  } else {
    same <- opposite <- shared * NA
  }
  if (signed) {
    up_a <- crossprod(up, M)
    down_a <- crossprod(down, M)
    c_obs <- same
    c_exp <- (up_a * t(up_a) + down_a * t(down_a)) / G
  } else {
    c_obs <- shared
    c_exp <- n_a * n_b / G
  }
  c_exp[G == 0] <- NA

  lo <- pmin(n_a, n_b)
  overlap <- c_obs / lo
  overlap_expected <- c_exp / lo
  overlap[lo == 0] <- overlap_expected[lo == 0] <- NA
  union <- n_a + n_b - c_obs
  jaccard <- c_obs / union
  jaccard_expected <- c_exp / (n_a + n_b - c_exp)
  jaccard[union == 0] <- jaccard_expected[union == 0] <- NA
  correct <- function(s, e) {
    out <- (s - e) / (1 - e)
    out[!is.na(e) & e >= 1] <- NA
    out
  }
  p <- if (signed) shared * NA else
    stats::phyper(shared - 1, n_a, G - n_a, n_b, lower.tail = FALSE)
  p[lo == 0] <- NA
  enrichment <- c_obs / c_exp
  enrichment[lo == 0] <- NA
  keep_names <- function(m) { dimnames(m) <- list(colnames(deg), colnames(deg)); m }
  lapply(list(genes = G, n_a = n_a, n_b = n_b, shared = shared, shared_same = same,
              shared_opposite = opposite, overlap = overlap,
              overlap_expected = overlap_expected,
              overlap_corrected = correct(overlap, overlap_expected),
              jaccard = jaccard, jaccard_expected = jaccard_expected,
              jaccard_corrected = correct(jaccard, jaccard_expected),
              enrichment = enrichment, p_value = p), keep_names)
}

## The "tested" matrix of an itv object; for a bare DEG matrix, everything.
measured_matrix <- function(x, deg) {
  if (inherits(x, "itv")) x$measured else array(TRUE, dim(deg), dimnames(deg))
}
