#' Keep experiments by DEG count or by their description
#'
#' Keeps the experiments with at least `min_degs` and at most `max_degs`
#' DEGs and, if conditions are given, those whose row of `description`
#' meets all of them. A message lists the experiments removed and why.
#' To pick experiments by name or number, use [itv_subset()].
#'
#' @param x An [itv] object.
#' @param ... Conditions on the columns of `description`, e.g.
#'   `Organ == "Root"`, `Method %in% c("RNA-seq", "Microarray")`. An
#'   experiment is kept if it meets all of them; `NA` counts as not met.
#' @param min_degs,max_degs Bounds (inclusive) on the number of DEGs.
#' @param description A data frame with one row per experiment and an
#'   `experiment` column of experiment names. Needed for the conditions.
#' @return An [itv] object with the experiments kept, in their order.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' d <- read.delim(system.file("extdata", "example_description.csv",
#'                             package = "InterTransViewer"), sep = ";")
#' x <- itv_read(f)
#' itv_filter(x, min_degs = 200)
#' itv_filter(x, Organ == "Root", Method == "RNA-seq", description = d)
#' @export
itv_filter <- function(x, ..., min_degs = 0, max_degs = Inf, description = NULL) {
  check_itv(x)
  check_number(min_degs, "min_degs", 0)
  if (!identical(max_degs, Inf)) check_number(max_degs, "max_degs", 0)
  exps <- colnames(x$deg)
  n_deg <- colSums(x$deg)
  why <- ifelse(n_deg < min_degs, sprintf("%d DEGs < min_degs", n_deg),
                ifelse(n_deg > max_degs, sprintf("%d DEGs > max_degs", n_deg), ""))

  conds <- eval(substitute(alist(...)))
  if (length(conds)) {
    if (is.null(description))
      stop("Conditions need a 'description' table to be evaluated in.", call. = FALSE)
    d <- description_rows(description, exps)
    for (cond in conds) {
      ok <- eval(cond, d, parent.frame())
      if (!is.logical(ok) || !length(ok) %in% c(1, nrow(d)))
        stop("The condition ", paste(deparse(cond), collapse = " "),
             " must give TRUE or FALSE for each experiment.", call. = FALSE)
      ok <- rep_len(ok, nrow(d))
      fail <- why == "" & (is.na(ok) | !ok)
      why[fail] <- paste("not", paste(deparse(cond), collapse = " "))
    }
  }

  keep <- why == ""
  if (!any(keep)) stop("No experiments left after filtering.", call. = FALSE)
  if (any(!keep))
    message("Removed ", sum(!keep), " of ", length(exps), " experiments: ",
            paste0(exps[!keep], " (", why[!keep], ")", collapse = ", "), ".")
  itv_subset(x, exps[keep])
}

#' Plot labels from a description table
#'
#' Builds the `labels` used by the plots and [itv_run()] from a table that
#' describes the experiments, one row each.
#'
#' @param description A data frame with an `experiment` column of
#'   experiment names and the columns to use.
#' @param columns Either the names of the columns to paste together, e.g.
#'   `c("Organ", "Time")`, or one template with column names in curly
#'   braces, as in the examples.
#' @param sep Separator between columns, when `columns` lists names.
#' @return A character vector of labels named by experiment. Missing values
#'   are left out of the labels.
#' @examples
#' d <- read.delim(system.file("extdata", "example_description.csv",
#'                             package = "InterTransViewer"), sep = ";")
#' itv_labels(d, c("Organ", "Time"))
#' itv_labels(d, "{Organ}, {Time} ({Method})")
#'
#' x <- itv_read(system.file("extdata", "example_profiles.csv",
#'                           package = "InterTransViewer"))
#' itv_plot_totals(itv_summary(x), labels = itv_labels(d, c("Organ", "Time")))
#' @export
itv_labels <- function(description, columns, sep = " ") {
  d <- check_description(description)
  if (!is.character(columns) || length(columns) == 0 || anyNA(columns))
    stop("'columns' must be column names or a template such as \"{Organ}, {Time}\".",
         call. = FALSE)
  check_string(sep, "sep")
  value <- function(col) {
    if (!col %in% names(d))
      stop("No column '", col, "' in the description. Columns: ",
           paste(names(d), collapse = ", "), ".", call. = FALSE)
    v <- as.character(d[[col]])
    v[is.na(v)] <- ""
    v
  }
  if (length(columns) == 1 && grepl("{", columns, fixed = TRUE)) {
    fields <- regmatches(columns, gregexpr("\\{[^{}]+\\}", columns))[[1]]
    if (length(fields) == 0)
      stop("The template has no {column} fields.", call. = FALSE)
    out <- rep(columns, nrow(d))
    for (f in unique(fields)) {
      v <- value(substr(f, 2, nchar(f) - 1))
      out <- vapply(seq_along(out), function(i) gsub(f, v[i], out[i], fixed = TRUE),
                    character(1))
    }
  } else {
    parts <- lapply(columns, value)
    out <- vapply(seq_len(nrow(d)), function(i) {
      p <- vapply(parts, `[`, character(1), i)
      paste(p[nzchar(p)], collapse = sep)
    }, character(1))
  }
  out <- trimws(gsub("[ ]{2,}", " ", out))
  stats::setNames(out, d$experiment)
}

#' Join itv objects
#'
#' Puts the experiments of several [itv] objects side by side, e.g. auxin and
#' ethylene experiments, to analyse them together. Genes are joined by ID.
#'
#' @param ... [itv] objects, or one list of them. If they are named, the
#'   names are put in front of the experiment names (`auxin_4h_X`), which
#'   is needed when two objects have experiments of the same name.
#' @param genes `"all"` keeps every gene, with `NA` (not DE) in the
#'   experiments that did not measure it; `"shared"` keeps only the genes
#'   present in every object.
#' @param FDR_cutoff,logFC_cutoff DEG definition of the result. By default
#'   that of the objects, which must then agree.
#' @param n_genes Number of genes tested, see [itv()]. By default the
#'   largest `n_genes` of the objects (or the number of genes joined, if
#'   larger), or `NULL` if none of them has one.
#' @return An [itv] object.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' x <- itv_read(f)
#' roots <- itv_subset(x, 1:4)
#' shoots <- itv_subset(x, 5:7)
#' itv_combine(roots, shoots)
#' itv_combine(a = roots, b = roots)   # same experiments, told apart by prefix
#' @export
itv_combine <- function(..., genes = c("all", "shared"), FDR_cutoff = NULL,
                        logFC_cutoff = NULL, n_genes = NULL) {
  genes <- match.arg(genes)
  objs <- list(...)
  if (length(objs) == 1 && is.list(objs[[1]]) && !inherits(objs[[1]], "itv"))
    objs <- objs[[1]]
  if (length(objs) == 0) stop("Give at least one itv object.", call. = FALSE)
  for (o in objs) check_itv(o)
  cutoff <- function(given, field) {
    if (!is.null(given)) return(given)
    v <- unique(vapply(objs, `[[`, numeric(1), field))
    if (length(v) > 1)
      stop("The objects use different ", field, " values (",
           paste(v, collapse = ", "), "); set ", field, ".", call. = FALSE)
    v
  }
  FDR_cutoff <- cutoff(FDR_cutoff, "FDR_cutoff")
  logFC_cutoff <- cutoff(logFC_cutoff, "logFC_cutoff")

  prefix <- names(objs)
  exps <- unlist(lapply(seq_along(objs), function(i) {
    e <- colnames(objs[[i]]$deg)
    if (!is.null(prefix) && nzchar(prefix[i])) paste0(prefix[i], "_", e) else e
  }))
  if (anyDuplicated(exps))
    stop("Experiment names repeat: ", paste(unique(exps[duplicated(exps)]), collapse = ", "),
         ". Name the objects, e.g. itv_combine(auxin = x, ethylene = y), to put the ",
         "names in front.", call. = FALSE)
  all_ids <- lapply(objs, function(o) rownames(o$deg))
  ids <- if (genes == "all") unique(unlist(all_ids)) else Reduce(intersect, all_ids)
  if (length(ids) == 0) stop("The objects share no genes.", call. = FALSE)
  join <- function(field) {
    m <- do.call(cbind, lapply(objs, function(o) o[[field]][match(ids, rownames(o$deg)), ,
                                                            drop = FALSE]))
    dimnames(m) <- list(ids, exps)
    m
  }
  if (is.null(n_genes)) {
    given <- unlist(lapply(objs, `[[`, "n_genes"))
    if (length(given)) n_genes <- max(given, length(ids))
  }
  itv(join("logFC"), join("FDR"), FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff,
      n_genes = n_genes)
}

## The rows of a description table for `exps`, in their order.
description_rows <- function(description, exps) {
  d <- check_description(description)
  missing <- setdiff(exps, d$experiment)
  if (length(missing))
    stop("The description has no row for: ", paste(missing, collapse = ", "), ".",
         call. = FALSE)
  d[match(exps, d$experiment), , drop = FALSE]
}

check_description <- function(description) {
  if (!is.data.frame(description) || !"experiment" %in% names(description))
    stop("'description' must be a data frame with an 'experiment' column of ",
         "experiment names.", call. = FALSE)
  e <- as.character(description$experiment)
  if (anyNA(e) || any(!nzchar(e)) || anyDuplicated(e))
    stop("The 'experiment' column of the description must hold unique, ",
         "non-empty names.", call. = FALSE)
  description$experiment <- e
  description
}
