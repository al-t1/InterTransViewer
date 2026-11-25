#' Create an itv object
#'
#' Holds the logFC and FDR of each gene in each experiment and the DEG calls
#' derived from them. Usually created by [itv_read()].
#'
#' A gene is a DEG in an experiment when its FDR is below `FDR_cutoff` and,
#' if `logFC_cutoff > 0`, its absolute logFC is at least `logFC_cutoff`.
#' A missing FDR means "not DE". A missing logFC with a significant FDR counts
#' as DE only when there is no logFC cutoff.
#'
#' @param logFC,FDR Numeric matrices (or data frames) of the same size,
#'   genes in rows and experiments in columns. Row names are gene IDs and
#'   column names experiment names; both must be unique. Missing names are
#'   filled in as `gene1, ...` and `exp1, ...`.
#' @param FDR_cutoff A number in (0, 1].
#' @param logFC_cutoff A number >= 0, on the scale of `logFC` (usually log2),
#'   e.g. `log2(1.5)`.
#' @return An object of class `itv`: a list with
#'   \item{logFC, FDR}{the input matrices;}
#'   \item{deg}{logical matrix of DEG calls;}
#'   \item{measured}{logical matrix, `TRUE` where the FDR is not missing.
#'     Kept for reference; the indicators treat unmeasured genes as not DE;}
#'   \item{metadata}{experiment descriptions or `NULL`, see [itv_add_metadata()];}
#'   \item{FDR_cutoff, logFC_cutoff}{the DEG definition.}
#' @examples
#' lfc <- matrix(c(2, -1, 0.1, 3, 0.5, -2), 3,
#'               dimnames = list(c("g1", "g2", "g3"), c("a", "b")))
#' fdr <- matrix(c(0.01, 0.2, 0.04, 0.001, 0.5, 0.03), 3, dimnames = dimnames(lfc))
#' x <- itv(lfc, fdr)
#' x$deg
#' itv(lfc, fdr, logFC_cutoff = 1)$deg
#' @export
itv <- function(logFC, FDR, FDR_cutoff = 0.05, logFC_cutoff = 0) {
  logFC <- as.matrix(logFC)
  FDR <- as.matrix(FDR)
  if (!is.numeric(logFC) || !is.numeric(FDR))
    stop("logFC and FDR must be numeric.", call. = FALSE)
  if (!identical(dim(logFC), dim(FDR)))
    stop("logFC and FDR must have the same dimensions.", call. = FALSE)
  if (ncol(logFC) < 1 || nrow(logFC) < 1)
    stop("logFC and FDR need at least one gene and one experiment.", call. = FALSE)
  if (is.null(colnames(logFC))) colnames(logFC) <- paste0("exp", seq_len(ncol(logFC)))
  if (is.null(rownames(logFC))) rownames(logFC) <- paste0("gene", seq_len(nrow(logFC)))
  if (anyDuplicated(colnames(logFC))) stop("Experiment names must be unique.", call. = FALSE)
  if (anyDuplicated(rownames(logFC))) stop("Gene IDs must be unique.", call. = FALSE)
  dimnames(FDR) <- dimnames(logFC)
  if (!is.numeric(FDR_cutoff) || length(FDR_cutoff) != 1 || is.na(FDR_cutoff) ||
      FDR_cutoff <= 0 || FDR_cutoff > 1)
    stop("FDR_cutoff must be a number in (0, 1].", call. = FALSE)
  if (!is.numeric(logFC_cutoff) || length(logFC_cutoff) != 1 || is.na(logFC_cutoff) ||
      logFC_cutoff < 0)
    stop("logFC_cutoff must be a number >= 0.", call. = FALSE)

  structure(list(logFC = logFC, FDR = FDR,
                 deg = call_degs(logFC, FDR, FDR_cutoff, logFC_cutoff),
                 measured = !is.na(FDR), metadata = NULL,
                 FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff),
            class = "itv")
}

call_degs <- function(logFC, FDR, FDR_cutoff = 0.05, logFC_cutoff = 0) {
  deg <- !is.na(FDR) & FDR < FDR_cutoff
  if (logFC_cutoff > 0)
    deg <- deg & !is.na(logFC) & abs(logFC) >= logFC_cutoff
  deg
}

#' @export
print.itv <- function(x, ...) {
  exps <- colnames(x$deg)
  cat("<itv> ", nrow(x$deg), " genes x ", length(exps), " experiments\n", sep = "")
  cat("DEG: FDR < ", x$FDR_cutoff,
      if (x$logFC_cutoff > 0) paste0(" and |logFC| >= ", signif(x$logFC_cutoff, 3)),
      "\n", sep = "")
  cat("Genes DE in at least one experiment: ", sum(rowSums(x$deg) > 0), "\n", sep = "")
  no_fc <- sum(x$deg & is.na(x$logFC))
  if (no_fc > 0) cat("DEG calls without a logFC: ", no_fc, "\n", sep = "")
  cat("Experiments: ", paste(utils::head(exps, 6), collapse = ", "),
      if (length(exps) > 6) ", ...", "\n", sep = "")
  cat("Metadata: ", if (is.null(x$metadata)) "none" else
        paste(ncol(x$metadata), "columns"), "\n", sep = "")
  invisible(x)
}

#' Experiment names of an itv object
#'
#' @param x An [itv] object.
#' @return A character vector.
#' @export
itv_experiments <- function(x) {
  check_itv(x)
  colnames(x$deg)
}

#' Keep some experiments
#'
#' @param x An [itv] object.
#' @param experiments Experiment names, or numbers in the order of
#'   [itv_experiments()]. The result follows this order. Unknown names or
#'   numbers and repeated experiments are errors.
#' @return An [itv] object with the selected experiments (and metadata rows).
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' itv_subset(x, c("root_1h_A", "root_4h_A"))
#' itv_subset(x, 1:3)
#' @export
itv_subset <- function(x, experiments) {
  check_itv(x)
  idx <- select_experiments(experiments, colnames(x$deg), "experiments",
                            allow_numbers = TRUE)
  for (m in c("logFC", "FDR", "deg", "measured")) x[[m]] <- x[[m]][, idx, drop = FALSE]
  if (!is.null(x$metadata)) {
    x$metadata <- x$metadata[idx, , drop = FALSE]
    rownames(x$metadata) <- NULL
  }
  x
}

## Turn a selection of experiments into column numbers. Names are matched
## against `available`; numbers are allowed only if allow_numbers = TRUE.
select_experiments <- function(sel, available, what = "selection",
                               allow_numbers = FALSE) {
  if (is.numeric(sel)) {
    if (!allow_numbers)
      stop(what, " must list experiment names, e.g. c(\"",
           paste(utils::head(available, 2), collapse = "\", \""), "\").", call. = FALSE)
    bad <- sel[!sel %in% seq_along(available)]
    if (length(bad))
      stop(what, ": no experiment number ", paste(bad, collapse = ", "),
           " (there are ", length(available), ").", call. = FALSE)
    idx <- as.integer(sel)
  } else if (is.character(sel)) {
    missing <- setdiff(sel, available)
    if (length(missing))
      stop(what, ": experiments not found: ", paste(missing, collapse = ", "), call. = FALSE)
    idx <- match(sel, available)
  } else {
    stop(what, " must be a character vector of experiment names.", call. = FALSE)
  }
  if (length(idx) == 0) stop(what, " is empty.", call. = FALSE)
  if (anyDuplicated(idx)) stop(what, " lists an experiment more than once.", call. = FALSE)
  idx
}

check_itv <- function(x) {
  if (!inherits(x, "itv")) stop("Expected an itv object; see itv_read().", call. = FALSE)
  invisible(x)
}

## The DEG matrix of an itv object, or a logical matrix given directly.
deg_matrix <- function(x) {
  if (inherits(x, "itv")) return(x$deg)
  if (is.matrix(x) && is.logical(x)) {
    if (is.null(colnames(x))) colnames(x) <- paste0("exp", seq_len(ncol(x)))
    return(x)
  }
  stop("Expected an itv object or a logical matrix of DEG calls.", call. = FALSE)
}
