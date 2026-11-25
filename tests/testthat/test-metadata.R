toy_meta <- function() {
  data.frame(Order = c("1", "2", "3"), Experiment = c("GSE1", "GSE2", "SRP3"),
             Time = c("0.5h", "4h", "NaN"), Age = c("?d", "7d", "7d"),
             DEGs = c("10", "20", "30"), Specific_DEGs = c("1%", "2%", "3%"),
             stringsAsFactors = FALSE)
}

test_that("itv_read_metadata keeps every value as text", {
  f <- tempfile(fileext = ".csv")
  writeLines(c("Order;Experiment;Time;Age;Conc", "1;GSE1;0.5h;?d;NaN", "2;GSE2;4h;7d;"), f)
  m <- itv_read_metadata(f)
  expect_equal(m$Order, c("1", "2"))
  expect_equal(m$Time, c("0.5h", "4h"))
  expect_equal(m$Age, c("?d", "7d"))
  expect_equal(m$Conc, c("NaN", ""))
})

test_that("itv_read_metadata detects comma and tab separators", {
  f1 <- tempfile(); writeLines(c("Order,Experiment", "1,GSE1"), f1)
  f2 <- tempfile(); writeLines(c("Order\tExperiment", "1\tGSE1"), f2)
  expect_equal(itv_read_metadata(f1)$Experiment, "GSE1")
  expect_equal(itv_read_metadata(f2)$Experiment, "GSE1")
  f3 <- tempfile(); writeLines(c("Order Experiment", "1 GSE1"), f3)
  expect_error(itv_read_metadata(f3), "Cannot tell")
  expect_error(itv_read_metadata(tempfile()), "File not found")
})

test_that("metadata rows are matched by order and get this run's counts", {
  x <- itv_add_metadata(toy_itv(), toy_meta())
  m <- itv_update_metadata(x)
  expect_equal(m$Experiment, c("GSE1", "GSE2", "SRP3"))
  expect_equal(m$experiment, c("e1", "e2", "e3"))
  expect_equal(m$DEGs, c(3, 3, 1))
  expect_equal(m$Specific_DEGs, c("33.33%", "33.33%", "0.00%"))
  expect_equal(m$Time, c("0.5h", "4h", "NaN"))
  expect_equal(names(m), c(names(toy_meta()), "experiment"))
})

test_that("metadata follows itv_subset", {
  x <- itv_subset(itv_add_metadata(toy_itv(), toy_meta()), c("e3", "e1"))
  m <- itv_update_metadata(x)
  expect_equal(m$Experiment, c("SRP3", "GSE1"))
  expect_equal(m$DEGs, c(1, 3))
})

test_that("metadata is matched by a column of experiment names", {
  meta <- toy_meta()
  meta$name <- c("e3", "e1", "e2")   # rows in a different order
  x <- itv_add_metadata(toy_itv(), meta)
  expect_equal(x$metadata$Experiment, c("GSE2", "SRP3", "GSE1"))
  expect_equal(x$metadata$name, c("e1", "e2", "e3"))
})

test_that("metadata with the wrong number of rows is refused", {
  expect_error(itv_add_metadata(toy_itv(), toy_meta()[1:2, ]),
               "2 rows but there are 3 experiments")
  expect_error(itv_add_metadata(toy_itv(), 1:3), "data frame")
  expect_error(itv_update_metadata(toy_itv()), "no metadata")
})

test_that("an experiment without DEGs gets NA for Specific_DEGs", {
  meta <- rbind(toy_meta(), toy_meta()[1, ])
  m <- itv_update_metadata(itv_add_metadata(toy_itv(extra_empty = TRUE), meta))
  expect_equal(m$DEGs[4], 0)
  expect_true(is.na(m$Specific_DEGs[4]))
})

test_that("itv_read attaches metadata from a file", {
  x <- itv_read(example_file(), metadata = example_file("example_description.csv"))
  expect_equal(nrow(x$metadata), 8)
  expect_equal(x$metadata$Organ[c(1, 5, 8)], c("Root", "Shoot", "Seedling"))
})
