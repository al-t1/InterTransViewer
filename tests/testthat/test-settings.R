quiet_run <- function(...) {
  itv_run(example_file(), out_dir = tempfile("itv_set_"), cluster_boot_iter = 3,
          boot_iter = 20, progress = FALSE, ...)
}

test_that("itv_cluster() scaling and distance options", {
  x <- itv_read(example_file())
  set.seed(1)
  ref <- itv_cluster(x, boot_iter = 2)
  set.seed(1)
  expect_equal(itv_cluster(x, scale_exps_by_range = TRUE, center_genes = FALSE, scale_genes_by_SD = TRUE,
                           distance = "euclidean", boot_iter = 2)$tree$merge, ref$tree$merge)

  ## z-scored genes and 1 - Pearson r, as by hand.
  cl <- itv_cluster(x, center_genes = TRUE, distance = "pearson", boot_iter = 2)
  keep <- rowSums(!is.finite(x$logFC) | is.na(x$FDR)) == 0
  S <- t(x$logFC[keep, ])
  S <- S / apply(S, 1, function(v) diff(range(v)))
  S <- scale(S)
  ref_tree <- stats::hclust(stats::as.dist(1 - stats::cor(t(S))), method = "ward.D2")
  expect_equal(cl$tree$merge, ref_tree$merge)
  expect_equal(cl$tree$height, ref_tree$height)
  expect_equal(cl$distance, "pearson")
  expect_equal(cl$n_total, 1500)
  expect_output(print(cl), "pearson distance")
  expect_output(print(cl), "measured in all experiments: 1050 of 1500")
  expect_output(print(cl), "genes centred")

  ## No scaling at all; Spearman.
  raw <- itv_cluster(x, scale_exps_by_range = FALSE, scale_genes_by_SD = FALSE, boot_iter = 2)
  expect_equal(raw$tree$merge,
               stats::hclust(stats::dist(t(x$logFC[keep, ])), "ward.D2")$merge)
  expect_s3_class(itv_cluster(x, distance = "spearman", boot_iter = 2), "itv_cluster")
  expect_error(itv_cluster(x, distance = "manhattan"))
  expect_error(itv_cluster(x, center_genes = NA), "'center_genes'")
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(plot(cl))
})

test_that("itv_settings() keeps named itv_run() arguments", {
  s <- itv_settings(FDR_cutoff = 0.01, seed = NULL, cluster_distance = "pearson")
  expect_s3_class(s, "itv_settings")
  expect_equal(names(s), c("FDR_cutoff", "seed", "cluster_distance"))
  expect_null(s$seed)
  s2 <- itv_settings(s, boot_iter = 10, FDR_cutoff = 0.1)
  expect_equal(s2$FDR_cutoff, 0.1)
  expect_equal(s2$boot_iter, 10)
  expect_true("seed" %in% names(s2))
  expect_output(print(s), "cluster_distance = \"pearson\"")
  expect_output(print(itv_settings()), "all defaults")
  expect_error(itv_settings(FDR = 0.1), "Not settings of itv_run\\(\\): FDR")
  expect_error(itv_settings(file = "a.csv"), "Not settings")
  expect_error(itv_settings(0.1), "must be named")
  expect_error(itv_settings(boot_iter = 1, boot_iter = 2), "Repeated")
})

test_that("itv_run() takes settings, and arguments in the call override them", {
  s <- itv_settings(FDR_cutoff = 0.01, cluster_distance = "pearson", bootstrap = FALSE,
                    logFC_cutoff = 1)
  res <- quiet_run(settings = s, logFC_cutoff = 0)
  expect_equal(res$itv$FDR_cutoff, 0.01)          # from the settings
  expect_equal(res$itv$logFC_cutoff, 0)           # the call wins
  expect_equal(res$cluster$distance, "pearson")
  expect_length(res$ci, 0)
  expect_equal(res$settings$FDR_cutoff, 0.01)
  expect_equal(res$settings$boot_iter, 20)
  expect_false("out_dir" %in% names(res$settings))
  ## The settings of a run reproduce it.
  again <- quiet_run(settings = res$settings)
  expect_equal(again$summary, res$summary)
  expect_equal(again$similarity, res$similarity)

  expect_error(quiet_run(settings = list(FDR_cutoff = 0.1)), "itv_settings")
  expect_error(quiet_run(settings = itv_settings(cluster_distance = "manhattan")),
               "'cluster_distance' must be one of")
  expect_error(quiet_run(settings = itv_settings(FDR_cutoff = 2)), "'FDR_cutoff'")
})

test_that("printing a run shows the summary, results and warnings", {
  res <- suppressWarnings(quiet_run(boot_subsets = list(bad = c("root_1h_A", "nope"))))
  expect_s3_class(res, "itv_run")
  out <- paste(utils::capture.output(print(res)), collapse = "\n")
  expect_match(out, "<itv_run> example: 8 experiments, 1500 genes", fixed = TRUE)
  expect_match(out, "seedling_2h_E", fixed = TRUE)
  expect_match(out, "Similarity: overlap, chance-corrected, genes tested in both", fixed = TRUE)
  expect_match(out, "Clustering: ward.D2, euclidean distance", fixed = TRUE)
  expect_match(out, "all (8 experiments)", fixed = TRUE)
  expect_match(out, "Bootstrap subset skipped", fixed = TRUE)
  expect_true(any(grepl("experiments not found: nope", res$notes)))
  ## Warnings still reach the user.
  expect_warning(quiet_run(boot_subsets = list(bad = "nope")), "skipped")
})

test_that("the heatmap draws a colour key without disturbing the device", {
  x <- itv_read(example_file())
  pdf(NULL); on.exit(grDevices::dev.off())
  op <- graphics::par(no.readonly = TRUE)
  itv_plot_heatmap(x)
  itv_plot_heatmap(x, "jaccard", chance_corrected = FALSE)
  expect_equal(graphics::par("fig"), op$fig)
  expect_equal(graphics::par("oma"), op$oma)
})

test_that("the heatmap is ordered by clustering 1 - similarity", {
  x <- itv_read(example_file())
  pdf(NULL); on.exit(grDevices::dev.off())
  oc <- unclass(itv_similarity(x))[, ]
  h <- itv_plot_heatmap(x)
  ref <- stats::reorder(stats::as.dendrogram(stats::hclust(stats::as.dist(1 - oc), "complete")),
                        rowMeans(oc))
  expect_equal(h$rowInd, stats::order.dendrogram(ref))
  ## Root (1-4) and shoot (5-7) experiments each sit together.
  pos <- match(1:8, h$rowInd)
  expect_equal(diff(range(pos[1:4])), 3)
  expect_equal(diff(range(pos[5:7])), 2)
})
