expected_files <- function(prefix) {
  paste0(prefix, "_", c("SignifTable.tsv", "JointDEGs.tsv", "TotalSpecificR.tsv",
                        "SpecificDEGs.tsv", "TotalSpecificDEG.pdf", "Rmetric.pdf",
                        "SimMatrix_overlap_corrected.tsv", "Heatmap_overlap_corrected.pdf",
                        "Clustering_min0.pdf", "Bootstrap.rds", "Bootstrap_CI.tsv",
                        "Bootstrap_CI.pdf", "run_info.txt"))
}
quick_run <- function(file, ...) {
  collect_warnings(itv_run(file, out_dir = tempfile("itv_run_"), cluster_boot_iter = 3,
                           boot_iter = 20, progress = FALSE, ...))
}

test_that("itv_run writes all outputs for the example data", {
  m <- utils::read.delim(example_file("example_description.csv"), sep = ";")
  lab <- stats::setNames(paste(m$Organ, m$Time, m$Method), m$experiment)
  r <- quick_run(example_file(), labels = lab,
                 boot_subsets = list(root = c("root_1h_A", "root_4h_A", "root_4h_B",
                                              "root_8h_B")))
  expect_length(r$warnings, 0)
  out <- r$value$out_dir
  expect_setequal(list.files(out), c(expected_files("example"),
                                     "example_Bootstrap.rds", "example_Bootstrap_CI.tsv",
                                     "example_Bootstrap_root.rds",
                                     "example_Bootstrap_CI_root.tsv"))

  tab <- utils::read.csv(file.path(out, "example_TotalSpecificR.tsv"), sep = "\t", quote = "")
  expect_equal(tab$total, r$value$summary$total)
  joined <- merge(m, r$value$summary, by = "experiment")
  expect_equal(nrow(joined), 8)
  info <- readLines(file.path(out, "example_run_info.txt"))
  expect_true(any(grepl("root_1h_A = Root 1h RNA-seq", info, fixed = TRUE)))
  ci <-utils::read.csv(file.path(out, "example_Bootstrap_CI_root.tsv"), sep = "\t", quote = "")
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

  r <- quick_run(example_file(), experiments = c(5, 1, 2))
  expect_equal(r$value$summary$experiment, c("shoot_1h_C", "root_1h_A", "root_4h_A"))
})

test_that("itv_run: labels need an entry for every experiment analysed", {
  lab <- c(shoot_1h_C = "Shoot", root_1h_A = "Root")
  r <- quick_run(example_file(), experiments = c("shoot_1h_C", "root_1h_A"),
                 labels = lab, bootstrap = FALSE)
  expect_length(r$warnings, 0)
  out <- tempfile("itv_run_")
  expect_error(itv_run(example_file(), out_dir = out, labels = lab, progress = FALSE),
               "no entry for: root_4h_A")
  expect_false(dir.exists(out))
  expect_error(itv_run(example_file(), out_dir = out, labels = c("a", "b")),
               "named by experiment")
})

test_that("itv_run uses the condition and similarity metric in file names", {
  r <- quick_run(example_file(), condition = "demo", similarity_metric = "jaccard",
                 chance_corrected = FALSE, bootstrap = FALSE)
  out <- r$value$out_dir
  expect_true(file.exists(file.path(out, "demo_SimMatrix_jaccard.tsv")))
  expect_true(file.exists(file.path(out, "demo_Heatmap_jaccard.pdf")))
  expect_false(file.exists(file.path(out, "demo_Bootstrap.rds")))
})

test_that("itv_run writes tab-separated tables that keep semicolons in names", {
  f <- tempfile(fileext = ".tsv")
  write_toy(f, sep = "\t", exps = c("a;1", "b,2", "c"))
  r <- quick_run(f, sep = "\t", condition = "toy", bootstrap = FALSE)
  out <- file.path(r$value$out_dir, "toy_TotalSpecificR.tsv")
  expect_equal(unique(utils::count.fields(out, sep = "\t", quote = "")), 9L)
  tab <- utils::read.delim(out, quote = "")
  expect_equal(tab$experiment, c("a;1", "b,2", "c"))

  g <- tempfile()
  write_tsv(data.frame(x = c("tab\there", "line\nbreak"), y = 1:2), g)
  expect_equal(readLines(g), c("x\ty", "tab here\t1", "line break\t2"))
})

test_that("itv_run warns when the output folder is not empty", {
  out <- tempfile("itv_run_")
  dir.create(out)
  writeLines("old", file.path(out, "old_file.txt"))
  r <- collect_warnings(itv_run(example_file(), out_dir = out, bootstrap = FALSE,
                                cluster_boot_iter = 3, progress = FALSE))
  expect_true(any(grepl("is not empty", r$warnings)))
  expect_true(file.exists(file.path(out, "old_file.txt")))
  r2 <- quick_run(example_file(), bootstrap = FALSE)
  expect_false(any(grepl("is not empty", r2$warnings)))
})

test_that("itv_run checks its settings before writing anything", {
  bad <- list(list(cluster_min_exps = 1.5, "cluster_min_exps"),
              list(cluster_boot_iter = 0, "cluster_boot_iter"),
              list(boot_iter = "100", "boot_iter"),
              list(FDR_cutoff = 5, "FDR_cutoff"),
              list(logFC_cutoff = NA, "logFC_cutoff"),
              list(seed = 1.5, "seed"),
              list(cluster_method = "ward", "cluster_method"),
              list(heatmap_hclust = "wardD2", "heatmap_hclust"),
              list(bootstrap = "yes", "bootstrap"),
              list(cluster_insignif_to_zero = NA, "cluster_insignif_to_zero"),
              list(condition = "a/b", "condition"),
              list(boot_subsets = list(c("root_1h_A", "root_4h_A")), "boot_subsets"),
              list(boot_subsets = list(all = c("root_1h_A", "root_4h_A")), "boot_subsets"),
              list(boot_subsets = c(a = "root_1h_A"), "boot_subsets"))
  for (b in bad) {
    out <- tempfile("itv_run_")
    args <- c(list(example_file(), out_dir = out, progress = FALSE), b[1])
    expect_error(do.call(itv_run, args), b[[2]], fixed = TRUE)
    expect_false(dir.exists(out))
  }
  expect_error(itv_run(tempfile()), "File not found")
})
