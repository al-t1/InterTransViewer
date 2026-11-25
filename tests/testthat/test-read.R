test_that("semicolon file is read into an itv object", {
  f <- tempfile(fileext = ".csv")
  ref <- write_toy(f, k = 3)
  x <- itv_read(f)
  expect_s3_class(x, "itv")
  expect_equal(dim(x$logFC), c(40L, 3L))
  expect_equal(itv_experiments(x), c("exp1", "exp2", "exp3"))
  expect_equal(rownames(x$logFC), ref$ids)
  expect_equal(unname(x$logFC), ref$logFC)
  expect_equal(unname(x$FDR), ref$FDR)
  expect_equal(x$deg, !is.na(x$FDR) & x$FDR < 0.05)
})

test_that("comma and tab separators are detected with a message", {
  f0 <- tempfile(fileext = ".csv"); write_toy(f0)
  x0 <- itv_read(f0)
  f1 <- tempfile(fileext = ".csv"); write_toy(f1, sep = ",")
  f2 <- tempfile(fileext = ".tsv"); write_toy(f2, sep = "\t")
  f3 <- tempfile(fileext = ".txt"); write_toy(f3, sep = "\t")
  expect_message(x1 <- itv_read(f1), "using comma")
  expect_message(x2 <- itv_read(f2), "using tab")
  expect_message(x3 <- itv_read(f3), "using tab")
  expect_equal(x1, x0); expect_equal(x2, x0); expect_equal(x3, x0)
  expect_silent(itv_read(f2, sep = "\t"))
})

test_that("decimal commas are detected", {
  f <- tempfile(fileext = ".csv"); ref <- write_toy(f, sep = ";", dec = ",")
  expect_message(x <- itv_read(f), "Decimal commas")
  expect_equal(unname(x$logFC), ref$logFC)
  expect_equal(unname(x$FDR), ref$FDR)
})

test_that("a UTF-8 byte order mark (Excel CSV) is handled", {
  f <- tempfile(fileext = ".csv"); write_toy(f, bom = TRUE)
  x <- itv_read(f)
  expect_equal(rownames(x$logFC)[1], "AT1G00001")
  expect_equal(itv_experiments(x), c("exp1", "exp2", "exp3"))
})

test_that("bad files give clear errors", {
  expect_error(itv_read(tempfile()), "File not found")
  f <- tempfile(); writeLines(character(0), f)
  expect_error(itv_read(f), "empty")
  f <- tempfile(); writeLines("ID;logFC_a;FDR_a", f)
  expect_error(itv_read(f), "no data rows")
  f <- tempfile(); writeLines(c("ID;logFC_a;FDR_a;logFC_b", "G1;1;0.01;2"), f)
  expect_error(itv_read(f), "Cannot tell")
  f <- tempfile(); writeLines(c("ID logFC_a FDR_a", "G1 1 0.01"), f)
  expect_error(itv_read(f), "Cannot tell")
  f <- tempfile(); writeLines(c("ID;a,b;c,d", "G1;1,2;3,4"), f)
  expect_error(itv_read(f, sep = "\t"), "Cannot tell")
})

test_that("missing and duplicated IDs are removed with warnings", {
  f <- tempfile()
  writeLines(c("ID;logFC_a;FDR_a", "G1;1;0.01", "G1;2;0.02", ";3;0.03", "G2;4;0.04"), f)
  r <- collect_warnings(itv_read(f))
  expect_true(any(grepl("without a gene ID", r$warnings)))
  expect_true(any(grepl("duplicated gene IDs", r$warnings)))
  expect_equal(r$value$logFC[, 1], c(G1 = 1, G2 = 4))
})

test_that("empty cells and NA become NA; text in numeric columns warns", {
  f <- tempfile()
  writeLines(c("ID;logFC_a;FDR_a", "G1;;0.01", "G2;NA;NA", "G3;oops;0.5"), f)
  expect_warning(x <- itv_read(f), "non-numeric")
  expect_true(all(is.na(x$logFC[, 1])))
  expect_equal(unname(x$FDR[, 1]), c(0.01, NA, 0.5))
  expect_equal(unname(x$measured[, 1]), c(TRUE, FALSE, TRUE))
})

test_that("Inf and -Inf logFC are kept as numbers", {
  f <- tempfile()
  writeLines(c("ID;logFC_a;FDR_a", "G1;Inf;0.01", "G2;-Inf;0.01", "G3;1;0.5"), f)
  expect_equal(unname(itv_read(f)$logFC[, 1]), c(Inf, -Inf, 1))
})

test_that("an all-NA experiment is read as numeric NA", {
  f <- tempfile()
  writeLines(c("ID;logFC_a;FDR_a;logFC_b;FDR_b", "G1;1;0.01;;", "G2;2;0.5;;"), f)
  x <- itv_read(f)
  expect_type(x$logFC, "double")
  expect_false(any(x$measured[, "b"]))
})

test_that("suspicious columns give warnings", {
  f <- tempfile()
  writeLines(c("ID;FDR_a;logFC_a", "G1;0.01;2.5", "G2;0.2;-1"), f)
  r <- collect_warnings(itv_read(f))
  expect_true(any(grepl("do not all look like logFC/FDR", r$warnings)))
  expect_true(any(grepl("outside \\[0, 1\\]", r$warnings)))
})

test_that("experiment names are taken from the logFC headers", {
  f <- tempfile()
  writeLines(c("ID;logFC_4h_X;FDR_4h_X;2h_Y.logFC;2h_Y.FDR;logFC;FDR;logFC_4h_X;FDR_4h_X",
               "G1;1;0.01;1;0.01;1;0.01;1;0.01"), f)
  expect_equal(itv_experiments(itv_read(f)), c("4h_X", "2h_Y", "exp3", "4h_X.1"))
})

test_that("the example data is read without warnings", {
  expect_silent(x <- itv_read(example_file()))
  expect_equal(dim(x$deg), c(1500L, 8L))
  expect_equal(sum(!x$measured[, "root_4h_B"]), 450)
})
