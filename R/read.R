#' Read differential expression profiles
#'
#' Reads a table with one row per gene: an ID column followed by one
#' logFC/FDR column pair per experiment, and returns an [itv] object.
#'
#' Experiment names come from the logFC column headers, without a "logFC_"
#' prefix (or a "_logFC"/".logFC" suffix). Empty cells and "NA" are missing
#' values. Decimal commas and the byte order mark of Excel CSV files are
#' handled. Rows without an ID are removed, and duplicated IDs keep their
#' first row, both with a warning. Warnings are also given when the column
#' names do not look like logFC/FDR pairs or FDR values fall outside \[0, 1\].
#'
#' @param file Path to the table.
#' @param sep Column separator. A semicolon by default. If it does not fit
#'   the header, semicolon, comma and tab are tried, and the only one that
#'   fits is used, with a message.
#' @param FDR_cutoff,logFC_cutoff,n_genes DEG definition and, for a table
#'   of DEGs only, the number of genes tested; see [itv()].
#' @return An object of class `itv`.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' x
#' @export
itv_read <- function(file, sep = ";", FDR_cutoff = 0.05, logFC_cutoff = 0,
                     n_genes = NULL) {
  p <- read_profiles(file, sep)
  itv(p$logFC, p$FDR, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff,
      n_genes = n_genes)
}

## Pick the column separator for a profiles table: `sep` if it splits the
## header into an ID column plus logFC/FDR pairs, otherwise the only one of
## ";", "," and tab that does.
detect_sep <- function(file, sep = ";") {
  header <- readLines(file, n = 1, warn = FALSE, encoding = "UTF-8")
  if (length(header) == 0 || !nzchar(trimws(header)))
    stop("'", basename(file), "' is empty or has no header line.", call. = FALSE)
  n_fields <- function(s) length(strsplit(header, s, fixed = TRUE)[[1]])
  fits <- function(s) { n <- n_fields(s); n >= 3 && (n - 1) %% 2 == 0 }
  if (fits(sep)) return(sep)

  candidates <- c(semicolon = ";", comma = ",", tab = "\t")
  ok <- candidates[vapply(candidates, fits, logical(1))]
  if (length(ok) == 1) {
    message("The header of '", basename(file), "' is not separated by '", sep,
            "'; using ", names(ok), " instead.")
    return(unname(ok))
  }
  stop("Cannot tell how the columns of '", basename(file), "' are separated. ",
       "Expected an ID column followed by logFC/FDR pairs, separated by ",
       "semicolons, commas or tabs. Header starts with: ", substr(header, 1, 80),
       call. = FALSE)
}

## Read the profiles table into two gene x experiment matrices, logFC and FDR.
read_profiles <- function(file, sep = ";") {
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  sep <- detect_sep(file, sep)
  read <- function(dec) {
    utils::read.table(file, header = TRUE, sep = sep, dec = dec, quote = "\"",
                      na.strings = c("NA", ""), check.names = FALSE,
                      stringsAsFactors = FALSE, comment.char = "",
                      fileEncoding = "UTF-8-BOM")
  }
  x <- read(".")
  if (nrow(x) == 0)
    stop("'", basename(file), "' has a header but no data rows.", call. = FALSE)
  if (ncol(x) < 3 || (ncol(x) - 1) %% 2 != 0)
    stop("Expected an ID column followed by logFC/FDR column pairs; found ",
         ncol(x), " columns.", call. = FALSE)

  ## Decimal commas (e.g. "0,05"), common in files saved by Excel.
  chr <- vapply(x[-1], is.character, logical(1))
  if (sep != "," && any(chr) &&
      any(grepl("^-?[0-9]*,[0-9]+([eE][-+]?[0-9]+)?$", unlist(x[-1][chr])))) {
    message("Decimal commas detected in '", basename(file), "'; reading them as decimals.")
    x <- read(",")
  }

  ids <- trimws(as.character(x[[1]]))
  bad <- is.na(ids) | ids == ""
  if (any(bad)) {
    warning(sum(bad), " rows without a gene ID were removed.", call. = FALSE)
    x <- x[!bad, , drop = FALSE]; ids <- ids[!bad]
  }
  dup <- duplicated(ids)
  if (any(dup)) {
    warning(sum(dup), " duplicated gene IDs; keeping the first occurrence.", call. = FALSE)
    x <- x[!dup, , drop = FALSE]; ids <- ids[!dup]
  }

  as_num <- function(v, name) {
    if (is.numeric(v)) return(v)
    if (is.logical(v)) return(as.numeric(v))   # a column that is entirely NA
    out <- suppressWarnings(as.numeric(v))
    if (any(is.na(out) & !is.na(v)))
      warning("Column '", name, "': non-numeric values were set to NA.", call. = FALSE)
    out
  }
  to_matrix <- function(cols) {
    m <- do.call(cbind, Map(as_num, x[cols], names(x)[cols]))
    dimnames(m) <- list(ids, NULL)
    m
  }
  fc_cols <- seq(2, ncol(x), by = 2)
  logFC <- to_matrix(fc_cols)
  FDR   <- to_matrix(fc_cols + 1)

  fc_names  <- names(x)[fc_cols]
  fdr_names <- names(x)[fc_cols + 1]
  if (!all(grepl("fc|fold", fc_names, ignore.case = TRUE)) ||
      !all(grepl("fdr|padj|adj|p.?val|q.?val", fdr_names, ignore.case = TRUE)))
    warning("Column names do not all look like logFC/FDR pairs ",
            "(expected: ID, logFC_<exp1>, FDR_<exp1>, logFC_<exp2>, ...). ",
            "Check the column order.", call. = FALSE)

  exps <- sub("^logFC[._ ]?|[._ ]?logFC$", "", fc_names, ignore.case = TRUE)
  check_pairs(exps, fc_names, fdr_names, fc_cols)
  exps[exps == ""] <- paste0("exp", which(exps == ""))
  exps <- make.unique(exps)
  colnames(logFC) <- colnames(FDR) <- exps

  if (any(FDR < 0 | FDR > 1, na.rm = TRUE))
    warning("Some FDR values are outside [0, 1]; check the column order.", call. = FALSE)
  list(logFC = logFC, FDR = FDR)
}

## Stop if an FDR column names a different experiment than the logFC column
## before it, e.g. logFC_a followed by FDR_b. Only checked where the FDR
## and logFC headers both have a recognised prefix or suffix and both names
## are left non-empty.
check_pairs <- function(fc_exps, fc_names, fdr_names, fc_cols) {
  tag <- "(FDR|padj|p[._]?adj|adj[._]?P[._]?Val|q[._]?val(ue)?|p[._]?val(ue)?)"
  fdr_exps <- sub(paste0("^", tag, "[._ ]?|[._ ]?", tag, "$"), "", fdr_names,
                  ignore.case = TRUE)
  comparable <- fdr_exps != fdr_names & fc_exps != fc_names &
    nzchar(fdr_exps) & nzchar(fc_exps)
  bad <- comparable & fc_exps != fdr_exps
  if (any(bad))
    stop("Name mismatch between ",
         paste0("'", fc_names[bad], "' (column ", fc_cols[bad], ") and '",
                fdr_names[bad], "' (column ", fc_cols[bad] + 1, ")", collapse = "; "),
         ". Each logFC column must be followed by the FDR column of the same ",
         "experiment.", call. = FALSE)
}
