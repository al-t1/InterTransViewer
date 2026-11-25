## Draw into a temporary PDF; return its path after closing it.
draw <- function(expr) {
  f <- tempfile(fileext = ".pdf")
  grDevices::pdf(f)
  dev <- grDevices::dev.cur()
  on.exit(grDevices::dev.off(dev))
  force(expr)
  f
}
nonempty_pdf <- function(f) file.exists(f) && file.size(f) > 1000

test_that("plotting functions draw without errors and restore par", {
  x <- itv_read(example_file())
  draw({
    mar <- graphics::par("mar")
    itv_plot_totals(x, "specific_percent", "example")
    itv_plot_totals(itv_summary(x), "R")
    expect_equal(graphics::par("mar"), mar)
    itv_plot_heatmap(x)
    itv_plot_heatmap(itv_similarity(x, "jaccard"), main = "jaccard")
    set.seed(1)
    cl <- itv_cluster(x, boot_iter = 3)
    plot(cl)
    plot(cl, labels = paste("label", 1:8), support_col = NA)
    expect_equal(graphics::par("mar"), mar)
    ci <- itv_ci(itv_bootstrap(x, iter = 20, progress = FALSE))
    itv_plot_ci(ci)
    itv_plot_ci(list(all = ci, sub = ci[6:7, ]), main = "two panels")
    itv_plot_hist(stats::rnorm(100), "hist")
    expect_equal(graphics::par("mar"), mar)
  })
  expect_error(plot(itv_cluster(x, boot_iter = 1), labels = "too few"), "one entry")
})

test_that("heatmap leaves out experiments without DEGs instead of failing", {
  x <- toy_itv(extra_empty = TRUE)
  f <- draw(expect_warning(itv_plot_heatmap(x), "without DEGs were left out: e4"))
  expect_true(nonempty_pdf(f))
})

test_that("heatmap with fewer than 2 usable experiments is skipped", {
  x <- itv(matrix(1, 2, 3), matrix(c(0.01, 0.5, rep(0.5, 4)), 2, 3))
  r <- collect_warnings(itv_plot_heatmap(x))
  expect_null(r$value)
  expect_true(any(grepl("at least 2 experiments", r$warnings)))
})

test_that("totals plot works with an experiment without DEGs", {
  expect_silent(draw(itv_plot_totals(toy_itv(extra_empty = TRUE), "R")))
})
