four_points <- function() {
  stats::hclust(stats::dist(c(a = 0, b = 1, c = 10, d = 11)), method = "complete")
}

test_that("clade_keys lists the members of each merge", {
  expect_equal(clade_keys(four_points()), c("1 2", "3 4", "1 2 3 4"))
})

test_that("node_x places each merge midway between its children", {
  tree <- four_points()
  expect_equal(tree$order, 1:4)
  expect_equal(node_x(tree), c(1.5, 3.5, 2.5))
})

random_itv <- function(n = 60, k = 5, seed = 2) {
  set.seed(seed)
  lfc <- matrix(stats::rnorm(n * k), n, k,
                dimnames = list(paste0("G", seq_len(n)), paste0("e", seq_len(k))))
  itv(lfc, matrix(stats::runif(n * k), n, k), FDR_cutoff = 0.2)
}

test_that("itv_cluster returns a tree, support values and gene counts", {
  x <- random_itv()
  set.seed(1)
  cl <- itv_cluster(x, boot_iter = 20)
  expect_s3_class(cl, "itv_cluster")
  expect_s3_class(cl$tree, "hclust")
  expect_equal(cl$experiments, itv_experiments(x))
  expect_length(cl$support, 4)
  expect_true(all(cl$support >= 0 & cl$support <= 1))
  expect_equal(cl$support[4], 1)   # the root is in every tree
  expect_equal(cl$n_genes, 60)
  expect_equal(cl$n_complete, 60)
  expect_match(capture.output(print(cl))[1], "5 experiments, ward.D2 on 60 genes")
})

test_that("clade support is reproducible with a seed", {
  x <- random_itv()
  set.seed(5); a <- itv_cluster(x, boot_iter = 15)
  set.seed(5); b <- itv_cluster(x, boot_iter = 15)
  expect_equal(a$support, b$support)
})

test_that("genes with NA, NaN, Inf or -Inf, or missing FDR, are left out", {
  x <- random_itv()
  x$logFC[1, 1] <- Inf
  x$logFC[2, 2] <- -Inf
  x$logFC[3, 3] <- NA
  x$FDR[4, 4]   <- NA
  cl <- itv_cluster(x, boot_iter = 2)
  expect_equal(cl$n_complete, 56)
  expect_equal(cl$n_genes, 56)
})

test_that("min_exps keeps only genes DE in enough experiments", {
  x <- random_itv()
  cl <- itv_cluster(x, min_exps = 2, boot_iter = 2)
  expect_equal(cl$n_genes, sum(rowSums(x$deg) >= 2))
})

test_that("an experiment without logFC variation is left out with a warning", {
  x <- random_itv()
  x$logFC[, "e3"] <- 0.5
  expect_warning(cl <- itv_cluster(x, boot_iter = 2), "no logFC variation.*e3")
  expect_equal(cl$experiments, c("e1", "e2", "e4", "e5"))
  expect_equal(sort(cl$tree$order), 1:4)
})

test_that("insignif_to_zero with a DEG-free experiment does not fail", {
  set.seed(2)
  lfc <- matrix(stats::rnorm(300), 60, 5, dimnames = list(NULL, paste0("e", 1:5)))
  fdr <- matrix(stats::runif(300), 60, 5)
  fdr[, 2] <- 0.9
  x <- itv(lfc, fdr, FDR_cutoff = 0.2)
  expect_warning(cl <- itv_cluster(x, insignif_to_zero = TRUE, boot_iter = 2),
                 "no logFC variation.*e2")
  expect_false("e2" %in% cl$experiments)
})

test_that("clustering stops clearly when too little is left", {
  x <- random_itv(k = 2)
  x$logFC[, 2] <- 1
  expect_error(suppressWarnings(itv_cluster(x, boot_iter = 2)), "at least 2 experiments")
  x <- random_itv()
  x$logFC[-1, 1] <- NA
  expect_error(itv_cluster(x, boot_iter = 2), "Too few genes")
})

test_that("itv_cluster uses the method argument", {
  cl <- itv_cluster(random_itv(), method = "average", boot_iter = 2)
  expect_equal(cl$tree$method, "average")
})

test_that("the example data separates root and shoot experiments", {
  set.seed(1)
  cl <- itv_cluster(itv_read(example_file()), boot_iter = 20)
  groups <- stats::cutree(cl$tree, k = 2)
  names(groups) <- cl$experiments
  expect_equal(length(unique(groups[1:4])), 1)
  expect_equal(length(unique(groups[5:7])), 1)
  expect_false(groups[1] == groups[5])
})
