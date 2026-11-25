## Argument checks shared by the exported functions. Each stops with a
## message naming the argument, or returns the value invisibly.

check_number <- function(x, name, min = -Inf, max = Inf, whole = FALSE,
                         min_open = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x) &&
    (if (min_open) x > min else x >= min) && x <= max &&
    (!whole || x == round(x))
  if (!ok) {
    what <- if (whole) "a whole number" else "a number"
    range <- if (is.finite(min) && is.finite(max)) {
      sprintf(" in %s%g, %g]", if (min_open) "(" else "[", min, max)
    } else if (is.finite(min)) {
      sprintf(" %s %g", if (min_open) ">" else ">=", min)
    } else ""
    stop("'", name, "' must be ", what, range, ".", call. = FALSE)
  }
  invisible(x)
}

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x))
    stop("'", name, "' must be TRUE or FALSE.", call. = FALSE)
  invisible(x)
}

check_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x))
    stop("'", name, "' must be a single non-empty character string.", call. = FALSE)
  invisible(x)
}

hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                    "mcquitty", "median", "centroid")

check_hclust_method <- function(x, name) {
  check_string(x, name)
  if (!x %in% hclust_methods)
    stop("'", name, "' must be one of: ", paste(hclust_methods, collapse = ", "), ".",
         call. = FALSE)
  invisible(x)
}
