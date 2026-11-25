test_that("itv calls DEGs from FDR, and logFC only when a cutoff is set", {
  lfc <- matrix(c(2, 0.1, NA, -3, 1), ncol = 1)
  fdr <- matrix(c(0.01, 0.01, 0.01, NA, 0.05), ncol = 1)
  expect_equal(unname(itv(lfc, fdr)$deg[, 1]), c(TRUE, TRUE, TRUE, FALSE, FALSE))
  expect_equal(unname(itv(lfc, fdr, logFC_cutoff = log2(1.5))$deg[, 1]),
               c(TRUE, FALSE, FALSE, FALSE, FALSE))
  expect_true(itv(matrix(1), matrix(0.01), logFC_cutoff = 1)$deg[1, 1])
  expect_true(itv(matrix(-Inf), matrix(0.01), logFC_cutoff = 1)$deg[1, 1])
  expect_false(itv(matrix(1), matrix(0.01), FDR_cutoff = 0.01)$deg[1, 1])
})

test_that("itv fills in names and stores the DEG definition", {
  x <- itv(matrix(1:4, 2), matrix(c(0.01, 0.5, 0.01, 0.01), 2))
  expect_equal(itv_experiments(x), c("exp1", "exp2"))
  expect_equal(rownames(x$deg), c("gene1", "gene2"))
  expect_equal(dimnames(x$FDR), dimnames(x$logFC))
  expect_equal(x$FDR_cutoff, 0.05)
  expect_equal(x$logFC_cutoff, 0)
  expect_null(x$metadata)
  expect_true(all(x$measured))
})

test_that("itv rejects bad input", {
  expect_error(itv(matrix("a"), matrix(0.1)), "numeric")
  expect_error(itv(matrix(1:4, 2), matrix(0.1, 2, 3)), "same dimensions")
  expect_error(itv(matrix(1:2, 1, dimnames = list("g", c("a", "a"))), matrix(0.1, 1, 2)),
               "Experiment names must be unique")
  expect_error(itv(matrix(1:2, 2, dimnames = list(c("g", "g"), "a")), matrix(0.1, 2, 1)),
               "Gene IDs must be unique")
  expect_error(itv(matrix(1), matrix(0.1), FDR_cutoff = 0), "FDR_cutoff")
  expect_error(itv(matrix(1), matrix(0.1), FDR_cutoff = c(0.1, 0.2)), "FDR_cutoff")
  expect_error(itv(matrix(1), matrix(0.1), logFC_cutoff = -1), "logFC_cutoff")
})

test_that("print shows the size, the DEG definition and metadata", {
  x <- itv_read(example_file())
  out <- capture.output(print(x))
  expect_match(out[1], "1500 genes x 8 experiments")
  expect_match(out[2], "FDR < 0.05")
  expect_true(any(grepl("Metadata: none", out)))
  x2 <- itv_read(example_file(), logFC_cutoff = 1,
                 metadata = example_file("example_description.csv"))
  out2 <- capture.output(print(x2))
  expect_match(out2[2], "\\|logFC\\| >= 1")
  expect_true(any(grepl("Metadata: 5 columns", out2)))
})

test_that("itv_subset keeps experiments by name or number, with their metadata", {
  x <- itv_read(example_file(), metadata = example_file("example_description.csv"))
  s <- itv_subset(x, c("shoot_1h_C", "root_1h_A"))
  expect_equal(itv_experiments(s), c("shoot_1h_C", "root_1h_A"))
  expect_equal(s$logFC, x$logFC[, c(5, 1)])
  expect_equal(s$deg, x$deg[, c(5, 1)])
  expect_equal(s$measured, x$measured[, c(5, 1)])
  expect_equal(s$metadata$Organ, c("Shoot", "Root"))
  expect_equal(itv_subset(x, c(5, 1)), s)
})

test_that("itv_subset rejects unknown or repeated experiments", {
  x <- toy_itv()
  expect_error(itv_subset(x, "nope"), "not found: nope")
  expect_error(itv_subset(x, 9), "no experiment number 9")
  expect_error(itv_subset(x, c("e1", "e1")), "more than once")
  expect_error(itv_subset(x, character(0)), "empty")
  expect_error(itv_subset(list(), 1), "itv object")
})

test_that("select_experiments refuses numbers unless allowed", {
  av <- c("a", "b", "c")
  expect_equal(select_experiments(c("c", "a"), av), c(3L, 1L))
  expect_error(select_experiments(c(1, 3), av, "boot_subsets$x"), "must list experiment names")
  expect_equal(select_experiments(c(1, 3), av, allow_numbers = TRUE), c(1L, 3L))
  expect_error(select_experiments(list("a"), av), "character vector")
})
