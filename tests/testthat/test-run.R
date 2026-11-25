expected_files <- function(prefix) {
  paste0(prefix, "_", c("SignifTable.csv", "JointDEGs.csv", "TotalSpecificR.csv",
                        "SpecificDEGs.csv", "TotalSpecificDEG.pdf", "Rmetric.pdf",
                        "SimMatrix_overlap.csv", "Heatmap_overlap.pdf",
                        "Clustering_min0.pdf", "Bootstrap.rds", "Bootstrap_CI.csv",
                        "Bootstrap_CI.pdf", "run_info.txt"))
}
quick_run <- function(file, ...) {
  collect_warnings(itv_run(file, out_dir = tempfile("itv_run_"), cluster_boot_iter = 3,
                           boot_iter = 20, progress = FALSE, ...))
}

test_that("itv_run writes all outputs for the example data", {
  r <- quick_run(example_file(), metadata = example_file("example_description.csv"),
                 boot_subsets = list(root = c("root_1h_A", "root_4h_A", "root_4h_B",
                                              "root_8h_B")))
  expect_length(r$warnings, 0)
  out <- r$value$out_dir
  expect_true(all(file.exists(file.path(out, expected_files("example")))))
  expect_true(file.exists(file.path(out, "example_Bootstrap_CI_root.csv")))
  expect_true(file.exists(file.path(out, "example_Description_Updated.csv")))

  tab <- utils::read.csv(file.path(out, "example_TotalSpecificR.csv"), sep = ";")
  expect_equal(tab$total, r$value$summary$total)
  meta <- utils::read.csv(file.path(out, "example_Description_Updated.csv"), sep = ";",
                          colClasses = "character")
  expect_equal(as.numeric(meta$DEGs), tab$total)
  expect_equal(meta$Organ[1], "Root")
  ci <- utils::read.csv(file.path(out, "example_Bootstrap_CI_root.csv"), sep = ";")
  expect_equal(ci$k, 3:1)
  expect_equal(names(r$value$ci), c("all", "root"))
})

test_that("itv_run gives the same results for the same seed", {
  a <- quick_run(example_file())$value
  b <- quick_run(example_file())$value
  expect_equal(a$cluster$support, b$cluster$support)
  expect_equal(a$bootstrap$all, b$bootstrap$all)
})

test_that("itv_run restores the random state", {
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  quick_run(example_file(), bootstrap = FALSE)
  expect_equal(stats::runif(1), before)
})

test_that("itv_run survives an experiment without DEGs", {
  f <- file.path(tempdir(), "Toy_Profiles.csv")
  write_toy(f, k = 4, n = 60, empty = 4)
  r <- quick_run(f, cluster_insignif_to_zero = TRUE,
                 boot_subsets = list(gone = c("exp1", "nope")))
  w <- r$warnings
  expect_true(any(grepl("Experiments without DEGs: exp4", w)))
  expect_true(any(grepl("Heatmap: experiments without DEGs were left out: exp4", w)))
  expect_true(any(grepl("no logFC variation were left out: exp4", w)))
  expect_true(any(grepl("Bootstrap subset skipped.*not found: nope", w)))
  expect_true(all(file.exists(file.path(r$value$out_dir, expected_files("Toy")))))
})

test_that("itv_run: numeric boot_subsets are refused, experiments take numbers", {
  r <- quick_run(example_file(), boot_subsets = list(nums = c(1, 2)))
  expect_true(any(grepl("must list experiment names", r$warnings)))
  expect_false("nums" %in% names(r$value$ci))

  r <- quick_run(example_file(), experiments = c(5, 1, 2),
                 metadata = example_file("example_description.csv"))
  expect_equal(r$value$summary$experiment, c("shoot_1h_C", "root_1h_A", "root_4h_A"))
  expect_equal(r$value$itv$metadata$Organ, c("Shoot", "Root", "Root"))
})

test_that("itv_run continues when the metadata does not fit", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Order;Organ", "1;Root", "2;Root"), f)
  r <- quick_run(example_file(), metadata = f, bootstrap = FALSE)
  expect_true(any(grepl("Metadata not used: The metadata has 2 rows", r$warnings)))
  expect_false(file.exists(file.path(r$value$out_dir, "example_Description_Updated.csv")))
  expect_true(file.exists(file.path(r$value$out_dir, "example_TotalSpecificR.csv")))
})

test_that("itv_run uses the condition and heatmap metric in file names", {
  r <- quick_run(example_file(), condition = "demo", heatmap_metric = "jaccard",
                 bootstrap = FALSE)
  out <- r$value$out_dir
  expect_true(file.exists(file.path(out, "demo_SimMatrix_jaccard.csv")))
  expect_true(file.exists(file.path(out, "demo_Heatmap_jaccard.pdf")))
  expect_false(file.exists(file.path(out, "demo_Bootstrap.rds")))
})
