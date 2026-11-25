#' Read experiment descriptions
#'
#' Reads a table with one row per experiment (e.g. method, organ, time).
#' All values are kept as text, so entries such as "0.5h", "?d" or "NaN"
#' stay as they are.
#'
#' @param file Path to the table.
#' @param sep Column separator, a semicolon by default. If it does not split
#'   the header, the one of semicolon, comma and tab that gives the most
#'   columns is used.
#' @return A data frame of character columns.
#' @examples
#' f <- system.file("extdata", "example_description.csv", package = "InterTransViewer")
#' itv_read_metadata(f)
#' @export
itv_read_metadata <- function(file, sep = ";") {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  header <- readLines(file, n = 1, warn = FALSE, encoding = "UTF-8")
  if (length(header) == 0 || !nzchar(trimws(header)))
    stop("'", basename(file), "' is empty or has no header line.", call. = FALSE)
  n_fields <- function(s) length(strsplit(header, s, fixed = TRUE)[[1]])
  if (n_fields(sep) < 2) {
    counts <- vapply(c(";", ",", "\t"), n_fields, integer(1))
    if (max(counts) < 2)
      stop("Cannot tell how the columns of '", basename(file), "' are separated.",
           call. = FALSE)
    sep <- names(counts)[which.max(counts)]
  }
  utils::read.table(file, header = TRUE, sep = sep, quote = "\"",
                    colClasses = "character", na.strings = character(0),
                    check.names = FALSE, comment.char = "",
                    fileEncoding = "UTF-8-BOM")
}

#' Attach experiment descriptions
#'
#' Rows are matched to experiments by a column that contains every
#' experiment name, if there is one. Otherwise the rows must follow the order
#' of the experiments in `x`, and their number must equal the number of
#' experiments. Attach metadata before [itv_subset()], which keeps the
#' matching rows.
#'
#' @param x An [itv] object.
#' @param metadata A data frame, or the path of a file read with
#'   [itv_read_metadata()].
#' @return `x` with `x$metadata` set, one row per experiment.
#' @export
itv_add_metadata <- function(x, metadata) {
  check_itv(x)
  if (is.character(metadata) && length(metadata) == 1)
    metadata <- itv_read_metadata(metadata)
  if (!is.data.frame(metadata))
    stop("metadata must be a data frame or a file path.", call. = FALSE)
  exps <- colnames(x$deg)
  is_key <- vapply(metadata, function(v) all(exps %in% v) && !anyDuplicated(v),
                   logical(1))
  if (any(is_key)) {
    rows <- match(exps, metadata[[which(is_key)[1]]])
  } else {
    if (nrow(metadata) != length(exps))
      stop("The metadata has ", nrow(metadata), " rows but there are ",
           length(exps), " experiments. Rows must follow the order of the ",
           "experiments, or one column must hold the experiment names.", call. = FALSE)
    rows <- seq_along(exps)
  }
  metadata <- metadata[rows, , drop = FALSE]
  rownames(metadata) <- NULL
  x$metadata <- metadata
  x
}

#' Experiment descriptions with DEG counts
#'
#' Adds the columns `experiment`, `DEGs` and `Specific_DEGs` (percent, e.g.
#' "1.44%") to the attached metadata. Existing `DEGs` and `Specific_DEGs`
#' columns are overwritten.
#'
#' @param x An [itv] object with metadata, see [itv_add_metadata()].
#' @return A data frame.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' m <- system.file("extdata", "example_description.csv", package = "InterTransViewer")
#' itv_update_metadata(itv_read(f, metadata = m))
#' @export
itv_update_metadata <- function(x) {
  check_itv(x)
  if (is.null(x$metadata)) stop("x has no metadata; see itv_add_metadata().", call. = FALSE)
  s <- itv_summary(x)
  m <- x$metadata
  m$experiment    <- s$experiment
  m$DEGs          <- s$total
  m$Specific_DEGs <- ifelse(is.na(s$specific_percent), NA,
                            sprintf("%.2f%%", s$specific_percent))
  m
}
