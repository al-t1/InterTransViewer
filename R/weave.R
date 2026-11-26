#' Combine differential expression results into one itv object
#'
#' Takes one differential expression result per experiment, e.g. the output
#' of limma, DESeq2 or edgeR, finds the gene ID, logFC and FDR columns, and
#' joins the experiments by gene ID. Genes missing from some experiments are
#' kept, with `NA` there (a full join), so they count as "not DE" in those
#' experiments. Use [itv_write()] to save the result as a profiles table.
#'
#' The columns are found by name, ignoring case and punctuation:
#' \describe{
#'   \item{logFC}{`logFC` (limma, edgeR), `log2FoldChange` (DESeq2),
#'     `log2FC`, `lfc`, `logFoldChange`;}
#'   \item{FDR}{`adj.P.Val` (limma), `padj` (DESeq2), `FDR` (edgeR),
#'     `qvalue`, `qval`, `p.adjust`;}
#'   \item{ID}{the row names (the unnamed first column written by
#'     `write.csv()`, or the row names of a data frame), or a column called
#'     `ID`, `GeneID`, `gene_id`, `gene`, `feature`, `locus`, `AGI`,
#'     `Symbol.AGI`, `ProbeID`, `Symbol`, ...; failing that, the first text
#'     column.}
#' }
#' Raw p-values are never used as FDR unless you name them in `FDR`. A
#' message lists the columns used for each experiment; give `id`, `logFC`
#' or `FDR` to override them.
#'
#' Files may be separated by commas, semicolons or tabs (detected from the
#' header) and may be gzipped. A header with one field fewer than the data
#' rows, as written by `write.table()`, means the first column holds row
#' names. Rows without an ID are removed, and duplicated IDs keep their
#' first row, both with a warning.
#'
#' @param sources The results, one per experiment: a character vector of
#'   file paths, a folder (every `.csv`, `.tsv`, `.txt` or `.tab` file in
#'   it, gzipped or not), or a list of data frames, such as
#'   `limma::topTable(fit, number = Inf)`,
#'   `as.data.frame(DESeq2::results(dds))` or
#'   `edgeR::topTags(qlf, n = Inf)` (the `table` is used).
#' @param names Experiment names. By default the names of the list, or the
#'   file names without the extension.
#' @param id,logFC,FDR Optional column choices, to use instead of the
#'   detected ones: a column name or number, the same for all experiments,
#'   or a vector with one per experiment (in the order of `sources`, or
#'   named by experiment; experiments left out are detected). Column numbers
#'   count the row-name column of a file, if there is one. `id = ""` means
#'   the row names.
#' @param sep Column separator of the files, or `NULL` to detect it.
#' @param upper_ids Convert gene IDs to upper case before joining, e.g. when
#'   one study writes `At1g01010` and another `AT1G01010`.
#' @param FDR_cutoff,logFC_cutoff,n_genes DEG definition and, for results
#'   of DEGs only, the number of genes tested; see [itv()].
#' @param verbose Print the columns used and the size of the result.
#' @return An object of class `itv`, with the genes in order of first
#'   appearance.
#' @seealso [itv_write()] to save the joined table, [itv_combine()] to join
#'   itv objects.
#' @examples
#' ## limma, DESeq2 and edgeR results for three experiments
#' folder <- system.file("extdata", "de_results", package = "InterTransViewer")
#' list.files(folder)
#' x <- itv_weave(folder)
#' x
#'
#' ## Other names; choose a column by hand
#' x <- itv_weave(folder, names = c("Root_1h", "Root_4h", "Shoot_4h"),
#'                FDR = c(Root_4h = "padj"))
#'
#' ## Data frames, e.g. straight from limma::topTable()
#' d <- read.csv(file.path(folder, "root_1h_limma.csv"), row.names = 1)
#' itv_weave(list(root_1h = d))
#'
#' ## Save it as one profiles table, for itv_read() or itv_run()
#' f <- tempfile(fileext = ".csv")
#' itv_write(x, f)
#' @export
itv_weave <- function(sources, names = NULL, id = NULL, logFC = NULL, FDR = NULL,
                      sep = NULL, upper_ids = FALSE, FDR_cutoff = 0.05,
                      logFC_cutoff = 0, n_genes = NULL, verbose = TRUE) {
  check_flag(upper_ids, "upper_ids")
  check_flag(verbose, "verbose")
  check_number(FDR_cutoff, "FDR_cutoff", 0, 1, min_open = TRUE)
  check_number(logFC_cutoff, "logFC_cutoff", 0)
  if (!is.null(sep)) check_string(sep, "sep")

  tables <- weave_sources(sources, sep)
  exps <- if (is.null(names)) base::names(tables) else names
  if (!is.character(exps) || length(exps) != length(tables) || anyNA(exps) ||
      any(!nzchar(exps)))
    stop("'names' must be ", length(tables), " non-empty experiment names, one per ",
         "result.", call. = FALSE)
  if (anyDuplicated(exps))
    stop("Experiment names must be unique; repeated: ",
         paste(unique(exps[duplicated(exps)]), collapse = ", "),
         ". Give them with 'names'.", call. = FALSE)
  base::names(tables) <- exps
  choice <- list(id = per_experiment(id, exps, "id"),
                 logFC = per_experiment(logFC, exps, "logFC"),
                 FDR = per_experiment(FDR, exps, "FDR"))

  parts <- vector("list", length(exps))
  base::names(parts) <- exps
  used <- character(length(exps))
  for (i in seq_along(exps)) {
    p <- weave_one(tables[[i]], exps[i], choice$id[[i]], choice$logFC[[i]],
                   choice$FDR[[i]], upper_ids)
    parts[[i]] <- p
    used[i] <- p$used
  }

  ids <- unique(unlist(lapply(parts, `[[`, "ids"), use.names = FALSE))
  fill <- function(what) {
    m <- matrix(NA_real_, length(ids), length(exps), dimnames = list(ids, exps))
    for (e in exps) m[match(parts[[e]]$ids, ids), e] <- parts[[e]][[what]]
    m
  }
  x <- itv(fill("logFC"), fill("FDR"), FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff,
           n_genes = n_genes)
  if (verbose) {
    partial <- sum(rowSums(x$measured) < length(exps))
    message(paste(c("Columns used:", paste0("  ", used),
                    sprintf("Joined: %d genes x %d experiment%s; %d genes %s.",
                            length(ids), length(exps), if (length(exps) > 1) "s" else "",
                            partial, "have no FDR in at least one experiment")),
                  collapse = "\n"))
  }
  x
}

#' Write an itv object as a profiles table
#'
#' Writes the table read by [itv_read()]: an `ID` column, then `logFC_<exp>`
#' and `FDR_<exp>` for each experiment. Missing values are written as `NA`.
#'
#' @param x An [itv] object, e.g. from [itv_weave()].
#' @param file Path of the file to write.
#' @param sep Column separator. A semicolon, as [itv_read()] expects by
#'   default; use a tab (`"\\t"`) if gene IDs or experiment names contain
#'   semicolons.
#' @return `file`, invisibly.
#' @examples
#' folder <- system.file("extdata", "de_results", package = "InterTransViewer")
#' x <- itv_weave(folder, verbose = FALSE)
#' f <- tempfile(fileext = ".csv")
#' itv_write(x, f)
#' all.equal(itv_read(f), x)
#' @export
itv_write <- function(x, file, sep = ";") {
  check_itv(x)
  check_string(file, "file")
  check_string(sep, "sep")
  exps <- colnames(x$deg)
  text <- c(rownames(x$deg), exps)
  bad <- grepl(sep, text, fixed = TRUE) | grepl("[\"\r\n]", text)
  if (any(bad))
    stop("Gene IDs or experiment names contain the separator, quotes or line breaks: ",
         paste(utils::head(text[bad], 3), collapse = ", "),
         if (sep != "\t") "; try sep = \"\\t\"", ".", call. = FALSE)
  df <- data.frame(ID = rownames(x$deg), check.names = FALSE, stringsAsFactors = FALSE)
  for (e in exps) {
    df[[paste0("logFC_", e)]] <- unname(x$logFC[, e])
    df[[paste0("FDR_", e)]] <- unname(x$FDR[, e])
  }
  utils::write.table(df, file, sep = sep, quote = FALSE, row.names = FALSE)
  invisible(file)
}

## ---- helpers -------------------------------------------------------------

## Column name candidates, compared after weave_key(); in order of preference.
weave_candidates <- list(
  logFC = c("logfc", "log2foldchange", "log2fc", "lfc", "logfoldchange", "log2fold",
            "log2ratio", "logratio"),
  FDR = c("adjpval", "padj", "fdr", "qvalue", "qval", "padjust", "padjusted",
          "adjpvalue", "adjp", "fdrbh"),
  id = c("", "rownames", "rowname", "id", "geneid", "geneids", "gene", "genes",
         "feature", "featureid", "locus", "locusid", "agi", "symbolagi", "ensembl",
         "ensemblid", "probeid", "probe", "probesetid", "x", "v1", "symbol",
         "genename", "genesymbol", "name"),
  pvalue = c("pvalue", "pval", "p")
)
weave_key <- function(s) tolower(gsub("[^[:alnum:]]", "", s))

## Turn `sources` into a named list of data frames. Row names of data frames
## become a first column with an empty name, as in files written by
## write.csv().
weave_sources <- function(sources, sep) {
  if (is.character(sources)) {
    if (length(sources) == 0) stop("'sources' is empty.", call. = FALSE)
    if (length(sources) == 1 && dir.exists(sources)) {
      files <- sort(list.files(sources, "\\.(csv|tsv|txt|tab)(\\.gz)?$",
                               full.names = TRUE, ignore.case = TRUE))
      if (length(files) == 0)
        stop("No .csv, .tsv, .txt or .tab files in '", sources, "'.", call. = FALSE)
    } else {
      files <- sources
    }
    missing <- files[!file.exists(files) | dir.exists(files)]
    if (length(missing))
      stop("File not found: ", paste(missing, collapse = ", "), call. = FALSE)
    out <- lapply(files, read_de_table, sep = sep)
    nms <- base::names(sources)
    if (is.null(nms) || length(files) != length(sources) || any(!nzchar(nms)))
      nms <- sub("\\.[^.]*$", "", sub("\\.gz$", "", basename(files), ignore.case = TRUE))
    base::names(out) <- nms
    return(out)
  }
  if (is.data.frame(sources) || !is.list(sources) || length(sources) == 0)
    stop("'sources' must be file paths, a folder, or a list of data frames, ",
         "e.g. list(exp1 = topTable(fit1, number = Inf), exp2 = ...).", call. = FALSE)
  nms <- base::names(sources)
  if (is.null(nms)) nms <- paste0("exp", seq_along(sources))
  nms[is.na(nms) | !nzchar(nms)] <- paste0("exp", which(is.na(nms) | !nzchar(nms)))
  out <- lapply(seq_along(sources), function(i) {
    d <- sources[[i]]
    if (is.list(d) && !is.data.frame(d) && is.data.frame(d$table)) d <- d$table  # TopTags
    d <- tryCatch(as.data.frame(d, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(d) || ncol(d) == 0)
      stop("Result ", i, " ('", nms[i], "') is not a table.", call. = FALSE)
    if (.row_names_info(d) > 0) {        # real row names, not 1..n
      d <- cbind(data.frame(id = rownames(d), stringsAsFactors = FALSE), d)
      names(d)[1] <- ""
    }
    rownames(d) <- NULL
    d
  })
  base::names(out) <- nms
  out
}

## Read one result file. A header with one field fewer than the rows means
## the first column holds row names (write.table()); it gets an empty name.
read_de_table <- function(file, sep = NULL) {
  first <- readLines(file, n = 2, warn = FALSE, encoding = "UTF-8")
  if (length(first) == 0 || !nzchar(trimws(first[1])))
    stop("'", basename(file), "' is empty.", call. = FALSE)
  if (is.null(sep)) {
    counts <- vapply(c("\t", ";", ","), function(s)
      length(strsplit(first[1], s, fixed = TRUE)[[1]]), integer(1))
    if (max(counts) < 2)
      stop("Cannot tell how the columns of '", basename(file), "' are separated; ",
           "give 'sep'.", call. = FALSE)
    sep <- names(counts)[which.max(counts)]
  }
  read <- function(skip, header, dec = ".") {
    utils::read.table(file, header = FALSE, sep = sep, dec = dec, quote = "\"",
                      skip = skip, nrows = if (header) 1 else -1,
                      na.strings = if (header) character(0) else c("NA", ""),
                      colClasses = if (header) "character" else NA,
                      stringsAsFactors = FALSE, comment.char = "",
                      fileEncoding = "UTF-8-BOM", check.names = FALSE)
  }
  hdr <- unlist(read(0, TRUE), use.names = FALSE)
  if (length(first) < 2)
    stop("'", basename(file), "' has a header but no data rows.", call. = FALSE)
  body <- read(1, FALSE)
  chr <- vapply(body, is.character, logical(1))
  if (sep != "," && any(chr[-1]) &&
      any(grepl("^-?[0-9]*,[0-9]+([eE][-+]?[0-9]+)?$", unlist(body[-1][chr[-1]])))) {
    message("Decimal commas detected in '", basename(file), "'; reading them as decimals.")
    body <- read(1, FALSE, dec = ",")
  }
  if (ncol(body) == length(hdr) + 1) hdr <- c("", hdr)
  if (ncol(body) != length(hdr))
    stop("'", basename(file), "': the header has ", length(hdr), " columns but the ",
         "rows have ", ncol(body), ".", call. = FALSE)
  names(body) <- hdr
  body
}

## A column choice per experiment: NULL (detect), one value for all, a value
## per experiment in order, or values named by experiment (others detected).
per_experiment <- function(v, exps, what) {
  if (is.null(v)) return(vector("list", length(exps)))
  ok <- (is.character(v) || is.numeric(v)) && !anyNA(v)
  nms <- names(v)
  if (ok && !is.null(nms)) {
    unknown <- setdiff(nms, exps)
    if (length(unknown) || anyDuplicated(nms))
      stop("'", what, "' is named by experiment, but these are not experiment names ",
           "(or repeat): ", paste(unique(c(unknown, nms[duplicated(nms)])), collapse = ", "),
           ". The experiments are: ", paste(exps, collapse = ", "), ".", call. = FALSE)
    out <- vector("list", length(exps))
    out[match(nms, exps)] <- as.list(unname(v))
    return(out)
  }
  if (ok && length(v) == 1) return(rep(list(v), length(exps)))
  if (ok && length(v) == length(exps)) return(as.list(v))
  stop("'", what, "' must be one column name or number for all experiments, one ",
       "per experiment (", length(exps), "), or a vector named by experiment.",
       call. = FALSE)
}

## Pick the ID, logFC and FDR columns of one table and return its values.
weave_one <- function(d, exp, id, logFC, FDR, upper_ids) {
  cols <- names(d)
  keys <- weave_key(cols)
  numeric_col <- vapply(d, function(v) is.numeric(v) || all(is.na(v)), logical(1))
  pick <- function(choice, what) {
    if (!is.null(choice)) {
      if (is.numeric(choice)) {
        if (length(choice) != 1 || !choice %in% seq_along(cols))
          stop(exp, ": no column number ", choice, " for ", what, " (there are ",
               length(cols), ").", call. = FALSE)
        return(as.integer(choice))
      }
      j <- match(choice, cols)
      if (is.na(j)) j <- match(weave_key(choice), keys)
      if (is.na(j))
        stop(exp, ": no column '", choice, "' for ", what, ". Columns: ",
             col_list(cols), ".", call. = FALSE)
      return(j)
    }
    cand <- weave_candidates[[what]]
    allowed <- if (what == "id") !numeric_col | keys %in% c("", "id", "geneid", "gene")
               else numeric_col
    hits <- match(cand, ifelse(allowed, keys, NA))
    hits <- hits[!is.na(hits)]
    if (length(hits)) return(hits[1])
    if (what == "id" && any(!numeric_col)) return(which(!numeric_col)[1])
    hint <- if (what == "FDR" && any(keys %in% weave_candidates$pvalue))
      " Only raw p-values were found; they are not used as FDR unless you ask for it."
    else ""
    stop(exp, ": cannot find the ", what, " column.", hint, " Columns: ",
         col_list(cols), ". Name it with ", what, " = c(", exp, " = \"...\").",
         call. = FALSE)
  }
  j_id <- pick(id, "id"); j_fc <- pick(logFC, "logFC"); j_fdr <- pick(FDR, "FDR")
  if (anyDuplicated(c(j_id, j_fc, j_fdr)))
    stop(exp, ": the ID, logFC and FDR columns must be different (columns ",
         j_id, ", ", j_fc, ", ", j_fdr, ").", call. = FALSE)

  ids <- trimws(as.character(d[[j_id]]))
  if (upper_ids) ids <- toupper(ids)
  fc <- weave_numeric(d[[j_fc]], exp, cols[j_fc])
  fdr <- weave_numeric(d[[j_fdr]], exp, cols[j_fdr])
  bad <- is.na(ids) | !nzchar(ids)
  if (any(bad))
    warning(exp, ": ", sum(bad), " rows without a gene ID were removed.", call. = FALSE)
  dup <- !bad & duplicated(ids)
  if (any(dup))
    warning(exp, ": ", sum(dup), " duplicated gene IDs; keeping the first occurrence.",
            call. = FALSE)
  keep <- !bad & !dup
  if (!any(keep)) stop(exp, ": no rows with a gene ID.", call. = FALSE)
  if (any(fdr < 0 | fdr > 1, na.rm = TRUE))
    warning(exp, ": some FDR values in '", cols[j_fdr], "' are outside [0, 1].",
            call. = FALSE)

  label <- function(j) if (nzchar(cols[j])) paste0("'", cols[j], "'") else "row names"
  tool <- if (all(c("log2foldchange", "padj") %in% keys)) " (DESeq2)"
    else if (all(c("logfc", "adjpval") %in% keys)) " (limma)"
    else if (all(c("logfc", "logcpm", "fdr") %in% keys)) " (edgeR)"
    else ""
  list(ids = ids[keep], logFC = fc[keep], FDR = fdr[keep],
       used = sprintf("%s: %d genes; ID = %s, logFC = %s, FDR = %s%s", exp, sum(keep),
                      label(j_id), label(j_fc), label(j_fdr), tool))
}

weave_numeric <- function(v, exp, name) {
  if (is.numeric(v)) return(as.numeric(v))
  if (is.logical(v) && all(is.na(v))) return(rep(NA_real_, length(v)))
  out <- suppressWarnings(as.numeric(as.character(v)))
  if (any(is.na(out) & !is.na(v)))
    warning(exp, ": column '", name, "' has non-numeric values; they were set to NA.",
            call. = FALSE)
  out
}

col_list <- function(cols) {
  shown <- ifelse(nzchar(cols), cols, "<row names>")
  paste0(seq_along(cols), " ", shown, collapse = ", ")
}
