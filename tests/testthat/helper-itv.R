## Helpers for the tests.

example_file <- function(name = "example_profiles.csv") {
  system.file("extdata", name, package = "InterTransViewer", mustWork = TRUE)
}

## Write a small profiles table: ID, then logFC_<exp>/FDR_<exp> pairs.
## `empty` = experiments (by number) with no significant FDR.
write_toy <- function(path, sep = ";", dec = ".", n = 40, k = 3, seed = 1,
                      empty = integer(0), exps = paste0("exp", seq_len(k)),
                      bom = FALSE) {
  set.seed(seed)
  fc  <- matrix(round(stats::rnorm(n * k, sd = 2), 3), n, k)
  fdr <- matrix(round(stats::runif(n * k)^2, 4), n, k)   # about 22% below 0.05
  fdr[, empty] <- 0.9
  df <- data.frame(ID = sprintf("AT1G%05d", seq_len(n)))
  for (j in seq_len(k)) {
    df[[paste0("logFC_", exps[j])]] <- fc[, j]
    df[[paste0("FDR_", exps[j])]]   <- fdr[, j]
  }
  utils::write.table(df, path, sep = sep, dec = dec, quote = FALSE, row.names = FALSE)
  if (bom) {
    body <- readBin(path, "raw", file.size(path))
    writeBin(c(as.raw(c(0xef, 0xbb, 0xbf)), body), path)
  }
  invisible(list(logFC = fc, FDR = fdr, ids = df$ID))
}

## The small DEG matrix used for hand-checked numbers.
##        e1 e2 e3
##   G1    T  T  F
##   G2    T  F  F
##   G3    F  T  F
##   G4    F  F  F
##   G5    T  T  T
toy_deg <- function() {
  matrix(c(TRUE, TRUE, FALSE, FALSE, TRUE,
           TRUE, FALSE, TRUE, FALSE, TRUE,
           FALSE, FALSE, FALSE, FALSE, TRUE), 5, 3,
         dimnames = list(paste0("G", 1:5), paste0("e", 1:3)))
}

## An itv object whose DEG calls are exactly toy_deg().
toy_itv <- function(extra_empty = FALSE) {
  deg <- toy_deg()
  if (extra_empty) deg <- cbind(deg, e4 = FALSE)
  lfc <- matrix(seq_along(deg), nrow(deg), dimnames = dimnames(deg))
  itv(lfc, ifelse(deg, 0.01, 0.5))
}

## Collect warnings instead of failing on them.
collect_warnings <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}
