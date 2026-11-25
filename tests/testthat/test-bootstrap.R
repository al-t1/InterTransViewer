random_deg <- function(n = 200, k = 6, p = 0.2, seed = 4) {
  set.seed(seed)
  matrix(stats::runif(n * k) < p, n, k, dimnames = list(NULL, paste0("e", seq_len(k))))
}

## Straightforward version of the same sampling, for comparison.
naive_bootstrap <- function(deg1, deg2, size1, sizes2, iter) {
  N <- function(deg, v) sum(rowSums(deg[, v, drop = FALSE]) > 0)
  sapply(sizes2, function(s) vapply(seq_len(iter), function(b) {
    v1 <- sample.int(ncol(deg1), size1, replace = TRUE)
    v2 <- sample.int(ncol(deg2), s, replace = TRUE)
    N(deg1, v1) - N(deg2, v2)
  }, numeric(1)))
}

test_that("itv_bootstrap has one column per k, named m vs k", {
  b <- itv_bootstrap(random_deg(), iter = 10, progress = FALSE)
  expect_equal(dim(b), c(10L, 5L))
  expect_equal(colnames(b), paste0("6vs", 5:1))
  expect_equal(attr(b, "size1"), 6)
  expect_equal(attr(b, "sizes2"), 5:1)
})

test_that("itv_bootstrap matches the straightforward computation", {
  deg <- random_deg()
  set.seed(11); fast <- itv_bootstrap(deg, iter = 30, progress = FALSE)
  set.seed(11); ref <- naive_bootstrap(deg, deg, 6, 5:1, 30)
  expect_equal(unname(unclass(fast)[, ]), unname(ref), ignore_attr = TRUE)
})

test_that("itv objects, two different sets and custom sizes work", {
  x <- itv_read(example_file())
  d2 <- random_deg(k = 4, n = 1500, seed = 2)
  set.seed(3)
  fast <- itv_bootstrap(x, d2, size1 = 4, sizes2 = 4, iter = 25, progress = FALSE)
  set.seed(3)
  ref <- naive_bootstrap(x$deg, d2, 4, 4, 25)
  expect_equal(colnames(fast), "4vs4")
  expect_equal(as.vector(fast), as.vector(ref))
})

test_that("identical experiments give d = 0", {
  deg <- random_deg(k = 1)[, rep(1, 5)]
  colnames(deg) <- paste0("e", 1:5)
  expect_true(all(itv_bootstrap(deg, iter = 20, progress = FALSE) == 0))
})

test_that("experiments without DEGs are allowed", {
  deg <- cbind(random_deg(k = 3), e4 = FALSE)
  expect_silent(b <- itv_bootstrap(deg, iter = 5, progress = FALSE))
  expect_false(anyNA(b))
})

test_that("itv_bootstrap is reproducible with a seed", {
  deg <- random_deg()
  set.seed(9); a <- itv_bootstrap(deg, iter = 10, progress = FALSE)
  set.seed(9); b <- itv_bootstrap(deg, iter = 10, progress = FALSE)
  expect_identical(a, b)
})

test_that("bad sample sizes are rejected", {
  expect_error(itv_bootstrap(random_deg(k = 1), iter = 2, progress = FALSE), "at least 2")
  expect_error(itv_bootstrap(random_deg(), sizes2 = 0, iter = 2, progress = FALSE), ">= 1")
})

test_that("itv_ci gives percentile intervals, sizes and zero flags", {
  b <- cbind("3vs2" = 1:100, "3vs1" = -50:49)
  ci <- itv_ci(b)
  expect_equal(ci$comparison, c("3vs2", "3vs1"))
  expect_equal(ci$m, c(3L, 3L))
  expect_equal(ci$k, c(2L, 1L))
  expect_equal(ci$lower, unname(c(quantile(1:100, 0.025), quantile(-50:49, 0.025))))
  expect_equal(ci$upper, unname(c(quantile(1:100, 0.975), quantile(-50:49, 0.975))))
  expect_equal(ci$includes_zero, c(FALSE, TRUE))
  b2 <- itv_bootstrap(random_deg(), iter = 10, progress = FALSE)
  expect_equal(itv_ci(b2)$k, 5:1)
})
