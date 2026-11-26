de_folder <- function() {
  system.file("extdata", "de_results", package = "InterTransViewer", mustWork = TRUE)
}
de_file <- function(name) file.path(de_folder(), name)

test_that("limma, DESeq2 and edgeR files are detected and fully joined", {
  expect_message(x <- itv_weave(de_folder()), "\\(DESeq2\\)")
  expect_s3_class(x, "itv")
  expect_equal(itv_experiments(x), c("root_1h_limma", "root_4h_DESeq2", "shoot_4h_edgeR"))
  expect_equal(nrow(x$deg), 400)
  expect_equal(unname(colSums(x$measured)), c(300, 400, 377))

  limma <- utils::read.csv(de_file("root_1h_limma.csv"), row.names = 1)
  deseq <- utils::read.csv(de_file("root_4h_DESeq2.csv"), row.names = 1)
  edger <- utils::read.delim(de_file("shoot_4h_edgeR.tsv"))
  expect_equal(unname(x$logFC[rownames(limma), 1]), limma$logFC)
  expect_equal(unname(x$FDR[rownames(limma), 1]), limma$adj.P.Val)
  expect_equal(unname(x$logFC[rownames(deseq), 2]), deseq$log2FoldChange)
  expect_equal(unname(x$FDR[rownames(deseq), 2]), deseq$padj)
  expect_equal(unname(x$FDR[rownames(edger), 3]), edger$FDR)
  ## Genes a study did not report are NA there, and not DE.
  off_array <- setdiff(rownames(x$deg), rownames(limma))
  expect_length(off_array, 100)
  expect_true(all(is.na(x$FDR[off_array, 1])))
  expect_false(any(x$deg[off_array, 1]))
})

test_that("the message lists the columns used", {
  msg <- testthat::capture_messages(itv_weave(de_folder()))
  expect_match(msg, "root_1h_limma: 300 genes; ID = row names, logFC = 'logFC', FDR = 'adj.P.Val' (limma)",
               fixed = TRUE, all = FALSE)
  expect_match(msg, "shoot_4h_edgeR: 377 genes; ID = row names, logFC = 'logFC', FDR = 'FDR' (edgeR)",
               fixed = TRUE, all = FALSE)
  expect_match(msg, "Joined: 400 genes x 3 experiments; 118 genes", fixed = TRUE, all = FALSE)
  expect_silent(itv_weave(de_folder(), verbose = FALSE))
})

test_that("files, names and data frames give the same result", {
  x <- itv_weave(de_folder(), verbose = FALSE)
  files <- file.path(de_folder(), c("root_1h_limma.csv", "root_4h_DESeq2.csv",
                                    "shoot_4h_edgeR.tsv"))
  expect_equal(itv_weave(files, verbose = FALSE), x)
  named <- itv_weave(stats::setNames(files, c("a", "b", "c")), verbose = FALSE)
  expect_equal(itv_experiments(named), c("a", "b", "c"))
  dfs <- list(root_1h_limma = utils::read.csv(files[1], row.names = 1),
              root_4h_DESeq2 = utils::read.csv(files[2], row.names = 1),
              shoot_4h_edgeR = utils::read.delim(files[3]))
  expect_equal(itv_weave(dfs, verbose = FALSE), x)
  ## A TopTags-like list holding the table.
  z <- itv_weave(list(e = list(table = dfs[[3]])), verbose = FALSE)
  expect_equal(unname(z$FDR[, 1]), dfs[[3]]$FDR)
  expect_equal(rownames(z$deg), rownames(dfs[[3]]))
})

test_that("ID columns, separators and decimal commas are handled", {
  f <- tempfile(fileext = ".txt")
  writeLines(c("gene_id;Symbol;log2FC;pval;qvalue", "G1;A;1,5;0,001;0,01",
               "G2;B;-2;0,2;0,4", ";C;1;0,1;0,2", "G1;D;3;0,01;0,03"), f)
  expect_message(suppressWarnings(itv_weave(f, verbose = FALSE)), "Decimal commas")
  r <- collect_warnings(suppressMessages(itv_weave(f, verbose = FALSE)))
  expect_true(any(grepl("without a gene ID", r$warnings)))
  expect_true(any(grepl("duplicated gene IDs", r$warnings)))
  x <- r$value
  expect_equal(rownames(x$deg), c("G1", "G2"))
  expect_equal(unname(x$logFC[, 1]), c(1.5, -2))
  expect_equal(unname(x$FDR[, 1]), c(0.01, 0.4))

  ## Numeric IDs in a column called "gene"; tab-separated; gzipped.
  g <- tempfile(fileext = ".tsv.gz")
  con <- gzfile(g, "w")
  writeLines(c("gene\tlogFC\tFDR", "101\t1\t0.01", "102\t2\t0.5"), con)
  close(con)
  y <- itv_weave(g, verbose = FALSE)
  expect_equal(rownames(y$deg), c("101", "102"))
  expect_equal(itv_experiments(y), sub("\\.tsv\\.gz$", "", basename(g)))
})

test_that("upper_ids joins IDs written in different case", {
  a <- data.frame(ID = c("At1g01010", "AT1G01020"), logFC = 1:2, FDR = c(0.01, 0.02))
  b <- data.frame(ID = c("AT1G01010", "at1g01020"), logFC = 3:4, FDR = c(0.01, 0.02))
  expect_equal(nrow(itv_weave(list(a = a, b = b), verbose = FALSE)$deg), 4)
  x <- itv_weave(list(a = a, b = b), upper_ids = TRUE, verbose = FALSE)
  expect_equal(rownames(x$deg), c("AT1G01010", "AT1G01020"))
  expect_true(all(x$deg))
})

test_that("column choices override detection, for all or per experiment", {
  files <- file.path(de_folder(), c("root_1h_limma.csv", "root_4h_DESeq2.csv"))
  raw <- itv_weave(files, FDR = c(root_4h_DESeq2 = "pvalue"), verbose = FALSE)
  deseq <- utils::read.csv(files[2], row.names = 1)
  expect_equal(unname(raw$FDR[rownames(deseq), 2]), deseq$pvalue)
  expect_equal(raw$FDR[, 1], itv_weave(files, verbose = FALSE)$FDR[, 1])
  ## By number (the row-name column is number 1) and by name, case-insensitive.
  num <- itv_weave(files[2], id = 1, logFC = 3, FDR = 7, verbose = FALSE)
  nam <- itv_weave(files[2], id = "", logFC = "LOG2FOLDCHANGE", FDR = "padj",
                   verbose = FALSE)
  expect_equal(num, nam)
  expect_equal(num, itv_weave(files[2], verbose = FALSE))
  ## One per experiment, in order.
  both <- itv_weave(files, FDR = c("adj.P.Val", "pvalue"), verbose = FALSE)
  expect_equal(both, raw)
})

test_that("unclear inputs give clear errors", {
  d <- data.frame(ID = c("g1", "g2"), logFC = 1:2, P.Value = c(0.01, 0.5))
  expect_error(itv_weave(list(a = d)), "Only raw p-values were found")
  expect_error(itv_weave(list(a = d)), "Columns: 1 ID, 2 logFC, 3 P.Value", fixed = TRUE)
  expect_error(itv_weave(list(a = d), FDR = "qval"), "no column 'qval'")
  expect_error(itv_weave(list(a = d), FDR = 9), "no column number 9")
  expect_error(itv_weave(list(a = d), FDR = c(b = "P.Value")), "not experiment names")
  expect_error(itv_weave(list(a = d), FDR = c("x", "y")), "one per experiment")
  expect_error(itv_weave(list(a = d), FDR = "logFC"), "must be different")
  expect_error(itv_weave(list(a = d, a = d), FDR = "P.Value"), "must be unique")
  expect_error(itv_weave(list(a = d), names = c("x", "y"), FDR = "P.Value"), "'names'")
  expect_error(itv_weave(d), "list of data frames")
  expect_error(itv_weave(tempfile()), "File not found")
  empty <- tempfile(); dir.create(empty)
  expect_error(itv_weave(empty), "No .csv")
  expect_warning(itv_weave(list(a = transform(d, P.Value = c(2, 0.1))), FDR = "P.Value",
                           verbose = FALSE), "outside \\[0, 1\\]")
})

test_that("itv_write() writes a table that itv_read() reads back", {
  x <- itv_weave(de_folder(), logFC_cutoff = 1, verbose = FALSE)
  f <- tempfile(fileext = ".csv")
  expect_equal(itv_write(x, f), f)
  expect_equal(readLines(f, 1),
               paste("ID", "logFC_root_1h_limma", "FDR_root_1h_limma", "logFC_root_4h_DESeq2",
                     "FDR_root_4h_DESeq2", "logFC_shoot_4h_edgeR", "FDR_shoot_4h_edgeR",
                     sep = ";"))
  expect_equal(itv_read(f, logFC_cutoff = 1), x)

  y <- itv(matrix(1, 1, 1, dimnames = list("g;1", "a")), matrix(0.01, 1, 1))
  expect_error(itv_write(y, f), "sep = \"\\t\"", fixed = TRUE)
  itv_write(y, f, sep = "\t")
  expect_equal(suppressMessages(itv_read(f, sep = "\t")), y)
})

test_that("itv_run() takes an itv object", {
  x <- itv_weave(de_folder(), verbose = FALSE)
  out <- file.path(tempfile(), "res")
  expect_error(itv_run(x, out_dir = out), "Give 'condition'")
  res <- itv_run(x, out_dir = out, condition = "DE", logFC_cutoff = 1,
                 cluster_boot_iter = 5, boot_iter = 20, progress = FALSE)
  expect_equal(res$itv$logFC_cutoff, 1)
  expect_equal(res$summary$total, unname(colSums(itv(x$logFC, x$FDR, logFC_cutoff = 1)$deg)))
  expect_true(file.exists(file.path(out, "DE_TotalSpecificR.tsv")))
  expect_match(readLines(file.path(out, "DE_run_info.txt"))[2],
               "file = an itv object (400 genes x 3 experiments)", fixed = TRUE)
})
