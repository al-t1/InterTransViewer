example_description <- function() {
  utils::read.delim(example_file("example_description.csv"), sep = ";")
}

test_that("itv_filter() keeps experiments by DEG count, bounds inclusive", {
  x <- itv_read(example_file())
  n <- colSums(x$deg)    # 311 331 175 174 265 271 265 122
  expect_message(y <- itv_filter(x, min_degs = 175), "Removed 2 of 8")
  expect_equal(itv_experiments(y), names(n)[n >= 175])
  expect_message(itv_filter(x, min_degs = 175), "root_8h_B (174 DEGs < min_degs)",
                 fixed = TRUE)
  expect_message(y <- itv_filter(x, max_degs = 265), "> max_degs")
  expect_equal(itv_experiments(y), names(n)[n <= 265])
  expect_silent(y <- itv_filter(x))
  expect_equal(y, x)
  expect_error(itv_filter(x, min_degs = 1000), "No experiments left")
  expect_error(itv_filter(x, min_degs = -1), "'min_degs'")
  expect_error(itv_filter(x, max_degs = NA), "'max_degs'")
})

test_that("itv_filter() evaluates conditions in the description", {
  x <- itv_read(example_file())
  d <- example_description()
  y <- suppressMessages(itv_filter(x, Organ == "Root", Method == "RNA-seq", description = d))
  expect_equal(itv_experiments(y), c("root_1h_A", "root_4h_A"))
  expect_message(itv_filter(x, Organ == "Root", description = d),
                 "shoot_1h_C (not Organ == \"Root\")", fixed = TRUE)
  ## Variables from the calling environment; order of the description does not matter.
  organs <- c("Shoot", "Seedling")
  y <- suppressMessages(itv_filter(x, Organ %in% organs, description = d[8:1, ]))
  expect_equal(itv_experiments(y), c("shoot_1h_C", "shoot_4h_C", "shoot_4h_D", "seedling_2h_E"))
  ## NA counts as not met.
  d$Time[1] <- NA
  y <- suppressMessages(itv_filter(x, Time == "1h", description = d))
  expect_equal(itv_experiments(y), "shoot_1h_C")

  expect_error(itv_filter(x, Organ == "Root"), "need a 'description'")
  expect_error(itv_filter(x, Organ == "Root", description = d[-1, ]), "no row for: root_1h_A")
  expect_error(itv_filter(x, Organ, description = d), "TRUE or FALSE")
  expect_error(itv_filter(x, Organ == "Root", description = d[, -1]), "'experiment' column")
})

test_that("itv_labels() pastes columns or fills a template", {
  d <- example_description()
  l <- itv_labels(d, c("Organ", "Time"))
  expect_equal(names(l), d$experiment)
  expect_equal(unname(l[1:2]), c("Root 1h", "Root 4h"))
  expect_equal(unname(itv_labels(d, c("Organ", "Time"), sep = ", ")[1]), "Root, 1h")
  t <- itv_labels(d, "{Organ}, {Time} ({Method}) {Organ}")
  expect_equal(unname(t[3]), "Root, 4h (Microarray) Root")
  ## Missing values are left out.
  d$Time[1] <- NA
  expect_equal(unname(itv_labels(d, c("Organ", "Time", "Method"))[1]), "Root RNA-seq")
  ## Usable as plot labels.
  x <- itv_read(example_file())
  pdf(NULL); on.exit(grDevices::dev.off())
  expect_silent(itv_plot_totals(itv_summary(x), labels = itv_labels(d, "Organ")))

  expect_error(itv_labels(d, "Tissue"), "No column 'Tissue'")
  expect_error(itv_labels(d, "{Tissue}"), "No column 'Tissue'")
  expect_error(itv_labels(d, 1), "'columns'")
  expect_error(itv_labels(d[, -1], "Organ"), "'experiment' column")
  expect_error(itv_labels(rbind(d, d), "Organ"), "unique")
})

test_that("itv_combine() joins genes and experiments", {
  x <- itv_read(example_file())
  a <- itv_subset(x, 1:3)
  b <- itv_subset(x, 4:8)
  expect_equal(itv_combine(a, b), x)
  expect_equal(itv_combine(list(a, b)), x)

  ## Different genes: all genes, NA where not measured; or only shared ones.
  b2 <- b
  for (m in c("logFC", "FDR", "deg", "measured")) b2[[m]] <- b[[m]][1:1000, , drop = FALSE]
  a2 <- a
  for (m in c("logFC", "FDR", "deg", "measured")) a2[[m]] <- a[[m]][501:1500, , drop = FALSE]
  all <- itv_combine(a2, b2)
  expect_equal(nrow(all$deg), 1500)
  expect_true(all(is.na(all$FDR[rownames(a$deg)[1:500], 1:3])))
  expect_equal(all$deg[rownames(a2$deg), 1:3], a2$deg)
  shared <- itv_combine(a2, b2, genes = "shared")
  expect_equal(rownames(shared$deg), rownames(x$deg)[501:1000])

  ## Repeated names need prefixes.
  expect_error(itv_combine(a, a), "Name the objects")
  ab <- itv_combine(first = a, second = a)
  expect_equal(itv_experiments(ab)[c(1, 4)], c("first_root_1h_A", "second_root_1h_A"))

  ## Cutoffs must agree unless given.
  a15 <- itv(a$logFC, a$FDR, logFC_cutoff = 1)
  expect_error(itv_combine(a15, b), "different logFC_cutoff")
  expect_equal(itv_combine(a15, b, logFC_cutoff = 1)$logFC_cutoff, 1)
  expect_error(itv_combine(a, list(1)), "Expected an itv object")
  expect_error(itv_combine(), "at least one")
})
