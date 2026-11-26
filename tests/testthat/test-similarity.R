## Brute-force pair statistics on the genes in `U`.
brute_pair <- function(x, a, b, U, signed = FALSE) {
  A <- x$deg[U, a]; B <- x$deg[U, b]
  fa <- x$logFC[U, a]; fb <- x$logFC[U, b]
  G <- sum(U)
  if (signed) {
    ua <- A & fa > 0; da <- A & fa < 0; ub <- B & fb > 0; db <- B & fb < 0
    shared <- sum(ua & ub, na.rm = TRUE) + sum(da & db, na.rm = TRUE)
    e <- (sum(ua, na.rm = TRUE) * sum(ub, na.rm = TRUE) +
            sum(da, na.rm = TRUE) * sum(db, na.rm = TRUE)) / G
  } else {
    shared <- sum(A & B)
    e <- sum(A) * sum(B) / G
  }
  lo <- min(sum(A), sum(B))
  I <- shared / lo; EI <- e / lo
  J <- shared / (sum(A) + sum(B) - shared); EJ <- e / (sum(A) + sum(B) - e)
  list(genes = G, n_a = sum(A), n_b = sum(B), overlap = I, overlap_expected = EI,
       overlap_corrected = (I - EI) / (1 - EI), jaccard = J, jaccard_expected = EJ,
       jaccard_corrected = (J - EJ) / (1 - EJ), enrichment = shared / e)
}

check_against_brute <- function(p, x, universe_fun, signed = FALSE) {
  for (i in seq_len(nrow(p))) {
    a <- p$experiment_a[i]; b <- p$experiment_b[i]
    ref <- brute_pair(x, a, b, universe_fun(a, b), signed)
    for (nm in names(ref)) expect_equal(p[[nm]][i], ref[[nm]], info = paste(a, b, nm))
  }
}

test_that("itv_pairs() matches a direct computation for each universe", {
  x <- itv_read(example_file())
  m <- x$measured
  p <- itv_pairs(x)
  expect_equal(nrow(p), choose(8, 2))
  check_against_brute(p, x, function(a, b) m[, a] & m[, b])
  check_against_brute(itv_pairs(x, "all"), x, function(a, b) rowSums(m) == ncol(m))
  check_against_brute(itv_pairs(x, "none"), x, function(a, b) rep(TRUE, nrow(m)))
  check_against_brute(itv_pairs(x, signed = TRUE), x, function(a, b) m[, a] & m[, b],
                      signed = TRUE)
  ## Pairs with the microarray study B are compared on its 1050 genes.
  expect_equal(unique(p$genes[p$experiment_b == "root_4h_B"]), 1050)
  expect_equal(unique(p$genes[p$experiment_a == "root_1h_A" &
                                p$experiment_b == "root_4h_A"]), 1500)
})

test_that("the hypergeometric p-value and the shared counts are right", {
  x <- itv_read(example_file())
  p <- itv_pairs(x)
  i <- which(p$experiment_a == "root_1h_A" & p$experiment_b == "shoot_1h_C")
  expect_equal(p$p_value[i], stats::phyper(p$shared[i] - 1, p$n_a[i],
                                           p$genes[i] - p$n_a[i], p$n_b[i],
                                           lower.tail = FALSE))
  expect_equal(p$shared, p$shared_same + p$shared_opposite)
  expect_true(all(is.na(itv_pairs(x, signed = TRUE)$p_value)))
})

test_that("signed similarity does not count genes changing in opposite directions", {
  lfc <- cbind(a = c(2, 2, -2, 1, 0), b = c(2, -2, -2, 1, 0))
  fdr <- cbind(a = c(0.01, 0.01, 0.01, 0.5, 0.5), b = c(0.01, 0.01, 0.01, 0.5, 0.5))
  x <- itv(lfc, fdr)
  p <- itv_pairs(x)
  expect_equal(c(p$shared, p$shared_same, p$shared_opposite), c(3, 2, 1))
  expect_equal(p$overlap, 1)
  ps <- itv_pairs(x, signed = TRUE)
  expect_equal(ps$overlap, 2 / 3)
  expect_equal(ps$jaccard, 2 / 4)      # (gene, direction) sets: {1+,2+,3-} vs {1+,2-,3-}
  expect_equal(ps$overlap_expected, (2 * 1 + 1 * 2) / 5 / 3)
  expect_error(itv_pairs(x$deg, signed = TRUE), "needs an itv object")
  expect_error(itv_pairs(x, signed = NA), "'signed'")
})

test_that("chance correction: 0 at the expected value, 1 at full overlap", {
  ## a's DEGs are inside b's, so I = 1 and the corrected I is 1.
  deg <- cbind(a = c(TRUE, TRUE, FALSE, FALSE), b = c(TRUE, TRUE, TRUE, FALSE))
  p <- itv_pairs(deg)
  expect_equal(p$overlap, 1)
  expect_equal(p$overlap_corrected, 1)
  expect_equal(p$overlap_expected, 3 / 4)
  ## Exactly the expected overlap gives 0.
  deg <- cbind(a = rep(c(TRUE, FALSE), 5), b = rep(c(TRUE, TRUE, FALSE, FALSE, FALSE), 2))
  p <- itv_pairs(deg)                 # shared 2, E = 5 * 4 / 10 = 2
  expect_equal(p$enrichment, 1)
  expect_equal(p$overlap_corrected, 0)
  expect_equal(p$jaccard_corrected, 0)
})

test_that("itv_similarity() returns each metric as a matrix", {
  x <- itv_read(example_file())
  p <- itv_pairs(x)
  for (m in c("overlap", "jaccard")) for (cc in c(TRUE, FALSE)) {
    S <- itv_similarity(x, m, chance_corrected = cc)
    expect_equal(dim(S), c(8L, 8L))
    expect_equal(attr(S, "metric"), m)
    expect_equal(attr(S, "chance_corrected"), cc)
    expect_equal(unname(S[cbind(match(p$experiment_a, rownames(S)),
                                match(p$experiment_b, rownames(S)))]),
                 p[[if (cc) paste0(m, "_corrected") else m]])
  }
  ## The default is the chance-corrected overlap.
  oc <- itv_similarity(x)
  expect_equal(oc, itv_similarity(x, "overlap", chance_corrected = TRUE))
  expect_equal(unname(diag(oc)), rep(1, 8))
  expect_error(itv_similarity(x, "overlap_corrected"))
  expect_error(itv_similarity(x, chance_corrected = "yes"), "'chance_corrected'")
  ## Within-organ pairs are well above chance, between-organ pairs near it.
  expect_gt(min(oc[1:4, 1:4]), 0.5)
  expect_lt(max(abs(oc[1:4, 5:7])), 0.1)
})

test_that("n_genes adds tested, non-DE genes to the chance levels only", {
  x <- itv_read(example_file())
  y <- itv_read(example_file(), n_genes = 3000)
  expect_equal(y$n_genes, 3000)
  py <- itv_pairs(y); px <- itv_pairs(x)
  expect_equal(py$genes, px$genes + 1500)
  expect_equal(py$shared, px$shared)
  expect_equal(py$overlap, px$overlap)
  expect_true(all(py$overlap_expected < px$overlap_expected))
  expect_equal(itv_subset(y, 1:2)$n_genes, 3000)
  expect_output(print(y), "Genes tested: 3000")
  expect_error(itv_read(example_file(), n_genes = 10), "'n_genes'")

  ## A table of DEGs only is noticed.
  only <- x$deg[rowSums(x$deg) > 0, ]
  expect_warning(itv(x$logFC[rownames(only), ], x$FDR[rownames(only), ]),
                 "looks filtered to DEGs")
  expect_silent(itv(x$logFC[rownames(only), ], x$FDR[rownames(only), ], n_genes = 1500))
})

test_that("itv_summary() splits specific DEGs into tested and untested", {
  x <- itv_read(example_file())
  s <- itv_summary(x)
  expect_equal(s$specific, s$specific_tested + s$specific_untested)
  ## Study B (microarray) measures every gene it reports, and the others
  ## measure all genes, so B's specific DEGs are all tested; root_1h_A's
  ## specific DEGs off the array are not.
  everywhere <- rowSums(x$measured) == 8
  spec <- x$deg & rowSums(x$deg) == 1
  expect_equal(s$specific_tested, unname(colSums(spec & everywhere)))
  expect_equal(s$specific_untested[3:4], c(0, 0))
  expect_gt(s$specific_untested[1], 0)
  sp <- itv_specific(x)
  expect_equal(as.vector(tapply(sp$tested, factor(sp$experiment, colnames(x$deg)), sum)),
               s$specific_tested)

  ## specific_oe by hand for experiment 1.
  n <- colSums(x$deg); G <- colSums(x$measured)
  e1 <- sum(x$deg[, 1] & everywhere) * prod(1 - n[-1] / G[-1])
  expect_equal(s$specific_oe[1], s$specific_tested[1] / e1)

  ## A bare DEG matrix counts everything as tested.
  s2 <- itv_summary(x$deg)
  expect_equal(s2$specific_untested, rep(0, 8))
})

test_that("the heatmap takes the metric, the correction and a shared order", {
  x <- itv_read(example_file())
  pdf(NULL); on.exit(grDevices::dev.off())
  oc <- itv_similarity(x)
  h1 <- itv_plot_heatmap(oc)
  h2 <- itv_plot_heatmap(x, "jaccard", chance_corrected = FALSE, order_by = oc)
  expect_equal(h1$rowInd, h2$rowInd)
  expect_equal(h1, itv_plot_heatmap(x))
  expect_silent(itv_plot_heatmap(x, "jaccard", universe = "all", signed = TRUE))
  ## A plain matrix uses chance_corrected for its colour scale.
  expect_silent(itv_plot_heatmap(unclass(oc)[, ], chance_corrected = FALSE))
  expect_error(itv_plot_heatmap(x, order_by = oc[1:3, 1:3]), "'order_by'")
})

test_that("itv_run() writes one similarity table and one heatmap", {
  out <- tempfile("itv_sim_")
  res <- itv_run(example_file(), out_dir = out, condition = "ex", cluster_boot_iter = 3,
                 bootstrap = FALSE, progress = FALSE)
  expect_equal(res$similarity, itv_similarity(res$itv))
  expect_equal(res$pairs, itv_pairs(res$itv))
  sim_files <- grep("Sim|Heatmap", list.files(out), value = TRUE)
  expect_setequal(sim_files, c("ex_SimMatrix_overlap_corrected.tsv",
                               "ex_Heatmap_overlap_corrected.pdf"))
  tab <- utils::read.delim(file.path(out, "ex_SimMatrix_overlap_corrected.tsv"),
                           check.names = FALSE)
  expect_equal(unname(as.matrix(tab[, -1])), unname(round(unclass(res$similarity)[, ], 4)))
  f <- file.path(out, "ex_Heatmap_overlap_corrected.pdf")
  expect_length(grepRaw("/Type /Page /", readBin(f, "raw", file.size(f)),
                        fixed = TRUE, all = TRUE), 1)
  info <- readLines(file.path(out, "ex_run_info.txt"))
  expect_true(paste("similarity: overlap, chance_corrected = TRUE, universe = pair,",
                    "signed = FALSE") %in% info)

  one <- itv_run(example_file(), out_dir = tempfile(), condition = "one",
                 similarity_metric = "jaccard", chance_corrected = FALSE,
                 universe = "none", signed = TRUE,
                 cluster_boot_iter = 3, bootstrap = FALSE, progress = FALSE)
  expect_equal(one$similarity,
               itv_similarity(one$itv, "jaccard", FALSE, universe = "none", signed = TRUE))
  expect_true(file.exists(file.path(one$out_dir, "one_Heatmap_jaccard.pdf")))
  expect_true(file.exists(file.path(one$out_dir, "one_SimMatrix_jaccard.tsv")))
  expect_equal(one$pairs, itv_pairs(one$itv, "none", signed = TRUE))
  expect_error(itv_run(example_file(), out_dir = tempfile(), universe = "both"))
  expect_error(itv_run(example_file(), out_dir = tempfile(), chance_corrected = NA),
               "'chance_corrected'")
})
