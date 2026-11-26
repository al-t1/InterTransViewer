test_that("itv_summary gives the hand-computed totals, specific DEGs and R", {
  s <- itv_summary(toy_itv())
  expect_equal(s$experiment, c("e1", "e2", "e3"))
  expect_equal(s$total, c(3, 3, 1))
  expect_equal(s$specific, c(1, 1, 0))
  expect_equal(s$specific_percent, c(33.33, 33.33, 0))
  ## N_all = 4 genes DE somewhere; R = (specific / total) / (total / N_all).
  expect_equal(s$R, c((1/3) / (3/4), (1/3) / (3/4), 0))
  expect_equal(itv_summary(toy_deg()), s)
})

test_that("itv_summary handles an experiment without DEGs", {
  s <- itv_summary(toy_itv(extra_empty = TRUE))
  expect_equal(s$total[4], 0)
  expect_true(is.na(s$specific_percent[4]))
  expect_true(is.na(s$R[4]))
  expect_equal(s$R[1:3], itv_summary(toy_itv())$R)
})

test_that("itv_specific lists genes DE in one experiment only", {
  sp <- itv_specific(toy_itv())
  expect_equal(sp$ID, c("G2", "G3"))
  expect_equal(sp$experiment, c("e1", "e2"))
  expect_equal(sp$logFC, c(2, 8))
})

test_that("itv_signif and itv_joint_degs cover genes DE somewhere", {
  x <- toy_itv()
  s <- itv_signif(x)
  expect_equal(rownames(s), c("G1", "G2", "G3", "G5"))
  expect_equal(unname(s["G1", ]), c(1, 6, NA))
  expect_equal(unname(s["G5", ]), c(5, 10, 15))
  j <- itv_joint_degs(x)
  expect_equal(j$ID, c("G1", "G2", "G3", "G5"))
  expect_equal(j$n_experiments, c(2, 1, 1, 3))
})

test_that("a DEG without logFC is counted but shown as NA", {
  x <- itv(matrix(c(NA, 1), 2, dimnames = list(c("g1", "g2"), "a")),
           matrix(c(0.01, 0.01), 2))
  expect_equal(itv_summary(x)$total, 2)
  expect_true(is.na(itv_signif(x)["g1", "a"]))
  expect_match(capture.output(print(x))[4], "DEG calls without a logFC: 1")
})

test_that("itv_similarity gives the overlap index I and Jaccard", {
  I <- itv_similarity(toy_itv(), "overlap", chance_corrected = FALSE)
  expect_equal(unname(diag(I)), c(1, 1, 1))
  expect_equal(I["e1", "e2"], 2 / 3)   # c = 2, min(3, 3) = 3
  expect_equal(I["e1", "e3"], 1)       # c = 1, min(3, 1) = 1
  expect_true(isSymmetric(I))
  J <- itv_similarity(toy_itv(), "jaccard", chance_corrected = FALSE)
  expect_equal(J["e1", "e2"], 2 / 4)   # 2 / (3 + 3 - 2)
  expect_equal(J["e1", "e3"], 1 / 3)
  expect_error(itv_similarity(toy_itv(), "cosine"))
})

test_that("itv_similarity gives NA for an experiment without DEGs", {
  I <- itv_similarity(toy_itv(extra_empty = TRUE))
  expect_true(all(is.na(I["e4", ])))
  expect_equal(I[1:3, 1:3], unclass(itv_similarity(toy_itv()))[, ])
})

test_that("itv_similarity matches a direct set computation", {
  x <- itv_read(example_file())
  I <- itv_similarity(x, chance_corrected = FALSE)
  a <- which(x$deg[, "root_1h_A"]); b <- which(x$deg[, "shoot_4h_D"])
  expect_equal(I["root_1h_A", "shoot_4h_D"], length(intersect(a, b)) / min(length(a), length(b)))
})

test_that("functions that need logFC refuse a bare DEG matrix", {
  expect_error(itv_specific(toy_deg()), "itv object")
  expect_error(itv_summary(matrix(1:4, 2)), "logical matrix")
})
