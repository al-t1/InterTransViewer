#' Run the whole analysis on one file
#'
#' Reads a profiles table, computes all indicators, and writes tables
#' (tab-separated, unquoted) and PDFs to `out_dir`, with names starting with
#' `condition`:
#' `SignifTable.tsv`, `JointDEGs.tsv`, `TotalSpecificR.tsv`,
#' `SpecificDEGs.tsv`, `TotalSpecificDEG.pdf`, `Rmetric.pdf`,
#' `SimMatrix_<metric>.tsv` and `Heatmap_<metric>.pdf` (e.g.
#' `overlap_corrected`, or `overlap` without chance correction),
#' `Clustering_min<min_exps>.pdf`,
#' `Bootstrap[_<subset>].rds`, `Bootstrap_CI[_<subset>].tsv`,
#' `Bootstrap_CI.pdf` and `run_info.txt` (settings and session info).
#' Tabs and line breaks inside text values are written as spaces.
#'
#' Tables identify experiments by name; the plots label them with their
#' number and `labels`. To add a description of the experiments to the
#' results, join it to `summary` by experiment name, e.g.
#' `merge(description, res$summary, by = "experiment")`.
#'
#' All settings are checked before any output is written. If `out_dir`
#' already holds files, a warning says so; files with the same names are
#' overwritten and other files are left as they are.
#' Problems in one step (clustering, a bootstrap subset) give a warning and
#' the other steps still run.
#'
#' @param file Path to the profiles table, see [itv_read()], or an [itv]
#'   object, e.g. from [itv_weave()] or [itv_combine()]; its DEGs are called
#'   again with `FDR_cutoff` and `logFC_cutoff`.
#' @param settings Optional settings from [itv_settings()]; arguments given
#'   in the call override them.
#' @param out_dir Output folder; created if needed.
#' @param condition Prefix for output files. By default the file name
#'   without "_Profiles" and the extension; needed when `file` is an itv
#'   object.
#' @param sep Column separator of `file`, see [itv_read()].
#' @param labels Optional plot labels: a character vector named by
#'   experiment, see [itv_plot_totals()]. It must have an entry for every
#'   experiment analysed; other entries are ignored.
#' @param experiments Optional experiments to analyse, by name or number,
#'   see [itv_subset()].
#' @param FDR_cutoff,logFC_cutoff,n_genes DEG definition and, for a table
#'   of DEGs only, the number of genes tested; see [itv()]. For an itv
#'   object, `n_genes = NULL` keeps its own.
#' @param seed Random seed for the bootstraps, or `NULL` to use the current
#'   random state. The previous random state is restored afterwards.
#' @param similarity_metric,chance_corrected,universe,signed How DEG lists
#'   are compared, see [itv_similarity()]. By default the overlap index I,
#'   corrected for chance, on the genes tested in both experiments of a
#'   pair. All metrics, raw and corrected, are in `res$pairs`.
#' @param heatmap_hclust Agglomeration method used to order the heatmap.
#' @param cluster_method,cluster_min_exps,cluster_insignif_to_zero,cluster_scale_exps_by_range,cluster_center_genes,cluster_scale_genes_by_SD,cluster_distance,cluster_boot_iter
#'   Passed to [itv_cluster()] as `method`, `min_exps`, `insignif_to_zero`,
#'   `scale_exps_by_range`, `center_genes`, `scale_genes_by_SD`, `distance` and
#'   `boot_iter`.
#' @param bootstrap Run the homogeneity bootstrap.
#' @param boot_iter Iterations for each k.
#' @param boot_subsets Named list of extra sets of experiments to bootstrap,
#'   each a character vector of experiment names. A set with unknown names is
#'   skipped with a warning.
#' @param progress Show bootstrap progress bars.
#' @return Invisibly, an object of class `itv_run`: a list with `condition`,
#'   `itv` (the data), `summary`, `similarity` (the matrix written), `pairs`
#'   (all metrics for every pair, see [itv_pairs()]; not written to a file),
#'   `cluster`, `bootstrap` and `ci` (lists by set), `out_dir`, `settings`
#'   (the settings used, see [itv_settings()]) and `notes` (the warnings
#'   given). Printing it shows the summary table, the main results and the
#'   warnings.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' m <- read.delim(system.file("extdata", "example_description.csv",
#'                             package = "InterTransViewer"), sep = ";")
#' res <- itv_run(f, out_dir = file.path(tempdir(), "itv_example"),
#'                labels = itv_labels(m, c("Organ", "Time", "Method")),
#'                cluster_boot_iter = 20, boot_iter = 100,
#'                boot_subsets = list(root = c("root_1h_A", "root_4h_A",
#'                                             "root_4h_B", "root_8h_B")))
#' res$summary
#' list.files(res$out_dir)
#'
#' ## The description with this run's results
#' merge(m, res$summary, by = "experiment")
#' @export
itv_run <- function(file, settings = NULL,
                    out_dir = file.path("Results", condition),
                    condition = default_condition(file),
                    sep = ";", labels = NULL, experiments = NULL,
                    FDR_cutoff = 0.05, logFC_cutoff = 0, n_genes = NULL, seed = 42,
                    similarity_metric = c("overlap", "jaccard"), chance_corrected = TRUE,
                    universe = c("pair", "all", "none"), signed = FALSE,
                    heatmap_hclust = "complete",
                    cluster_method = "ward.D2", cluster_min_exps = 0,
                    cluster_insignif_to_zero = FALSE, cluster_scale_exps_by_range = TRUE,
                    cluster_center_genes = FALSE, cluster_scale_genes_by_SD = TRUE,
                    cluster_distance = c("euclidean", "pearson", "spearman"),
                    cluster_boot_iter = 100,
                    bootstrap = TRUE, boot_iter = 1000, boot_subsets = list(),
                    progress = interactive()) {
  ## Settings fill in the arguments not given in this call.
  if (!is.null(settings)) {
    if (!inherits(settings, "itv_settings"))
      stop("'settings' must come from itv_settings().", call. = FALSE)
    for (nm in names(settings))
      if (eval(call("missing", as.name(nm)))) assign(nm, settings[[nm]])
  }
  similarity_metric <- check_choice(similarity_metric, c("overlap", "jaccard"),
                                    "similarity_metric")
  universe <- check_choice(universe, c("pair", "all", "none"), "universe")
  cluster_distance <- check_choice(cluster_distance, c("euclidean", "pearson", "spearman"),
                                   "cluster_distance")
  args <- mget(setdiff(names(formals(itv_run)), "settings"))

  ## Run, keeping a copy of every warning for print().
  notes <- character(0)
  res <- withCallingHandlers(
    do.call(run_analysis, args),
    warning = function(w) notes <<- c(notes, conditionMessage(w)))
  res$settings <- structure(args[setting_names], class = "itv_settings")
  res$notes <- notes
  invisible(structure(res, class = "itv_run"))
}

## The settings kept in a run's result and reusable for other data: all but
## the data-specific ones.
setting_names <- c("FDR_cutoff", "logFC_cutoff", "seed", "similarity_metric",
                   "chance_corrected", "universe", "signed", "heatmap_hclust",
                   "cluster_method", "cluster_min_exps", "cluster_insignif_to_zero",
                   "cluster_scale_exps_by_range", "cluster_center_genes", "cluster_scale_genes_by_SD",
                   "cluster_distance", "cluster_boot_iter", "bootstrap", "boot_iter",
                   "progress")

## The analysis of itv_run(), with every setting resolved.
run_analysis <- function(file, out_dir, condition, sep, labels, experiments,
                         FDR_cutoff, logFC_cutoff, n_genes, seed, similarity_metric,
                         chance_corrected, universe, signed, heatmap_hclust,
                         cluster_method, cluster_min_exps, cluster_insignif_to_zero,
                         cluster_scale_exps_by_range, cluster_center_genes, cluster_scale_genes_by_SD,
                         cluster_distance, cluster_boot_iter, bootstrap, boot_iter,
                         boot_subsets, progress) {
  check_run_args(file = file, out_dir = out_dir, condition = condition, sep = sep,
                 FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff, seed = seed,
                 chance_corrected = chance_corrected, signed = signed,
                 heatmap_hclust = heatmap_hclust, cluster_method = cluster_method,
                 cluster_min_exps = cluster_min_exps,
                 cluster_insignif_to_zero = cluster_insignif_to_zero,
                 cluster_scale_exps_by_range = cluster_scale_exps_by_range,
                 cluster_center_genes = cluster_center_genes,
                 cluster_scale_genes_by_SD = cluster_scale_genes_by_SD,
                 cluster_boot_iter = cluster_boot_iter, bootstrap = bootstrap,
                 boot_iter = boot_iter, boot_subsets = boot_subsets, progress = progress)
  if (!is.null(labels)) check_labels(labels)

  ## 1. Data and DEG calls; labels are checked against the experiments before
  ## anything is written.
  x <- if (inherits(file, "itv")) {
    itv(file$logFC, file$FDR, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff,
        n_genes = if (is.null(n_genes)) file$n_genes else n_genes)
  } else {
    itv_read(file, sep = sep, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff,
             n_genes = n_genes)
  }
  if (!is.null(experiments)) x <- itv_subset(x, experiments)
  deg <- x$deg
  exp_labels <- resolve_labels(labels, colnames(deg))
  names(exp_labels) <- colnames(deg)
  no_fc <- sum(deg & is.na(x$logFC))
  if (no_fc > 0)
    warning(no_fc, " DEG calls have a significant FDR but no logFC; ",
            "they are counted as DEGs, with logFC shown as NA.", call. = FALSE)
  no_deg <- colnames(deg)[colSums(deg) == 0]
  if (length(no_deg) > 0)
    warning("Experiments without DEGs: ", paste(no_deg, collapse = ", "), call. = FALSE)

  if (!is.null(seed)) {
    restore_rng <- save_rng()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE)))
    warning("Output folder '", out_dir, "' is not empty. Files with the same names ",
            "are overwritten; other files are left as they are.", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(name) file.path(out_dir, paste0(condition, "_", name))
  pdf_to <- function(name, width, height) {
    grDevices::pdf(out(name), width = width, height = height)
    dev <- grDevices::dev.cur()
    function() grDevices::dev.off(dev)
  }

  ## 2. Gene tables.
  signif <- itv_signif(x)
  write_tsv(data.frame(ID = rownames(signif), signif, check.names = FALSE),
             out("SignifTable.tsv"))
  write_tsv(itv_joint_degs(x), out("JointDEGs.tsv"))

  ## 3. Totals, specific DEGs and R.
  tsr <- itv_summary(x)
  write_tsv(tsr, out("TotalSpecificR.tsv"))
  write_tsv(itv_specific(x), out("SpecificDEGs.tsv"))
  height <- max(4, 0.28 * nrow(tsr) + 2)
  width <- max(10, 6 + 0.08 * max(nchar(paste(tsr$N, exp_labels))))
  close_dev <- pdf_to("TotalSpecificDEG.pdf", width, height)
  itv_plot_totals(tsr, "specific_percent", condition, labels = exp_labels)
  close_dev()
  close_dev <- pdf_to("Rmetric.pdf", width, height)
  itv_plot_totals(tsr, "R", condition, labels = exp_labels)
  close_dev()

  ## 4. Similarity of DEG lists: one metric, one table, one heatmap.
  sim_name <- sim_key(similarity_metric, chance_corrected)
  sim <- itv_similarity(x, similarity_metric, chance_corrected, universe, signed)
  write_tsv(data.frame(experiment = rownames(sim), round(sim, 4), check.names = FALSE),
            out(paste0("SimMatrix_", sim_name, ".tsv")))
  if (sum(!is.na(diag(sim))) >= 2) {
    side <- max(7.5, 5 + 0.1 * max(nchar(exp_labels)))
    close_dev <- pdf_to(paste0("Heatmap_", sim_name, ".pdf"), side + 0.5, side)
    itv_plot_heatmap(sim, hclust_method = heatmap_hclust,
                     main = paste0(condition, ": ", sim_titles[[similarity_metric]],
                                   if (chance_corrected) ", chance-corrected", "\n",
                                   universe_title[[universe]],
                                   if (signed) ", same direction only"),
                     labels = exp_labels)
    close_dev()
  } else {
    itv_plot_heatmap(sim)   # gives the warnings, draws nothing
  }

  ## 5. Clustering.
  cl <- tryCatch(
    itv_cluster(x, min_exps = cluster_min_exps, insignif_to_zero = cluster_insignif_to_zero,
                scale_exps_by_range = cluster_scale_exps_by_range, center_genes = cluster_center_genes,
                scale_genes_by_SD = cluster_scale_genes_by_SD, distance = cluster_distance,
                method = cluster_method, boot_iter = cluster_boot_iter),
    error = function(e) {
      warning("Clustering skipped: ", conditionMessage(e), call. = FALSE)
      NULL
    })
  if (!is.null(cl)) {
    leaf_labels <- stats::setNames(sprintf("%s (%d DEGs, %.1f%% spec.)", exp_labels,
                                    tsr$total, tsr$specific_percent), tsr$experiment)
    close_dev <- pdf_to(paste0("Clustering_min", cluster_min_exps, ".pdf"),
                        max(8, 0.5 * length(cl$experiments) + 3),
                        max(9, 6 + 0.1 * max(nchar(leaf_labels))))
    plot(cl, labels = leaf_labels,
         main = sprintf(paste0("%s: %s, %s distance, on %d genes (DE in >= %d ",
                               "experiments; %d of %d genes measured in all experiments)\n",
                               "scaling: %s; red numbers: bootstrap support, %%, ",
                               "%d replicates"),
                        condition, cluster_method, cluster_distance, cl$n_genes,
                        cluster_min_exps, cl$n_complete, nrow(deg),
                        cluster_scaling_text(cl), cluster_boot_iter))
    close_dev()
  }

  ## 6. Homogeneity bootstrap for all experiments and each subset.
  boots <- list()
  cis <- list()
  if (bootstrap) {
    sets <- list(all = seq_len(ncol(deg)))
    for (nm in names(boot_subsets)) {
      sets[[nm]] <- tryCatch(
        select_experiments(boot_subsets[[nm]], colnames(deg), paste0("boot_subsets$", nm)),
        error = function(e) {
          warning("Bootstrap subset skipped. ", conditionMessage(e), call. = FALSE)
          NULL
        })
    }
    for (nm in names(sets)) {
      cols <- sets[[nm]]
      if (is.null(cols)) next
      if (length(cols) < 2) {
        warning("Bootstrap set '", nm, "' skipped; it needs at least 2 experiments.",
                call. = FALSE)
        next
      }
      b <- itv_bootstrap(deg[, cols, drop = FALSE], iter = boot_iter, progress = progress)
      tag <- if (nm == "all") "" else paste0("_", nm)
      saveRDS(b, out(paste0("Bootstrap", tag, ".rds")))
      ci <- itv_ci(b)
      write_tsv(ci, out(paste0("Bootstrap_CI", tag, ".tsv")))
      boots[[nm]] <- b
      cis[[nm]] <- ci
    }
    if (length(cis) > 0) {
      close_dev <- pdf_to("Bootstrap_CI.pdf",
                          2 + 0.3 * sum(vapply(cis, nrow, integer(1)) + 3), 5)
      itv_plot_ci(cis, main = sprintf("%s: 95%% percentile intervals of d, %d iterations",
                                      condition, boot_iter))
      close_dev()
    }
  }

  ## 7. Settings and session.
  writeLines(c(sprintf("InterTransViewer %s", utils::packageVersion("InterTransViewer")),
               sprintf("file = %s", if (inherits(file, "itv"))
                 sprintf("an itv object (%d genes x %d experiments)",
                         nrow(file$deg), ncol(file$deg)) else file),
               sprintf("labels = %s", if (is.null(labels)) "NULL" else
                         paste(paste0(names(exp_labels), " = ", exp_labels), collapse = "; ")),
               sprintf("experiments = %s", paste(deparse(experiments), collapse = "")),
               sprintf("FDR_cutoff = %g; logFC_cutoff = %g; n_genes = %s; seed = %s",
                       FDR_cutoff, logFC_cutoff,
                       if (is.null(x$n_genes)) "NULL" else x$n_genes,
                       if (is.null(seed)) "NULL" else seed),
               sprintf("similarity: %s, chance_corrected = %s, universe = %s, signed = %s",
                       similarity_metric, chance_corrected, universe, signed),
               sprintf("heatmap: %s", heatmap_hclust),
               sprintf(paste0("clustering: %s, distance = %s, min_exps = %s, ",
                              "insignif_to_zero = %s, scale_exps_by_range = %s, center_genes = %s, ",
                              "scale_genes_by_SD = %s, boot_iter = %s"),
                       cluster_method, cluster_distance, cluster_min_exps,
                       cluster_insignif_to_zero, cluster_scale_exps_by_range, cluster_center_genes,
                       cluster_scale_genes_by_SD, cluster_boot_iter),
               sprintf("bootstrap: %s, iter = %s", bootstrap, boot_iter),
               sprintf("boot_subsets = %s", paste(deparse(boot_subsets), collapse = "")),
               "", utils::capture.output(utils::sessionInfo())),
             out("run_info.txt"))

  list(condition = condition, itv = x, summary = tsr, similarity = sim,
       pairs = itv_pairs(x, universe = universe, signed = signed), cluster = cl,
       bootstrap = boots, ci = cis, out_dir = out_dir)
}

#' @export
print.itv_run <- function(x, ...) {
  s <- x$settings
  d <- x$itv
  cat("<itv_run> ", x$condition, ": ", ncol(d$deg), " experiments, ", nrow(d$deg),
      " genes", if (!is.null(d$n_genes)) paste0(" (", d$n_genes, " tested)"), "\n", sep = "")
  cat("DEG: FDR < ", d$FDR_cutoff,
      if (d$logFC_cutoff > 0) paste0(" and |logFC| >= ", signif(d$logFC_cutoff, 3)), "\n",
      sep = "")
  cat("Output: ", x$out_dir, "\n\n", sep = "")

  tab <- x$summary
  tab$R <- round(tab$R, 2)
  tab$specific_oe <- round(tab$specific_oe, 2)
  print(tab, row.names = FALSE)

  sim <- x$similarity
  off <- sim[upper.tri(sim)]
  cat("\nSimilarity: ", attr(sim, "metric"),
      if (isTRUE(attr(sim, "chance_corrected"))) ", chance-corrected", ", ",
      universe_title[[s$universe]], if (isTRUE(s$signed)) ", same direction only", "\n",
      sep = "")
  if (any(!is.na(off))) {
    best <- which(sim == max(off, na.rm = TRUE) & upper.tri(sim), arr.ind = TRUE)[1, ]
    cat(sprintf("  pairs: median %.2f, range %.2f to %.2f; most similar: %s and %s\n",
                stats::median(off, na.rm = TRUE), min(off, na.rm = TRUE),
                max(off, na.rm = TRUE), rownames(sim)[best[1]], colnames(sim)[best[2]]))
  }

  cl <- x$cluster
  if (is.null(cl)) {
    cat("Clustering: skipped\n")
  } else {
    cat(sprintf("Clustering: %s, %s distance, %d genes (%s); clade support %d-%d%%\n",
                cl$tree$method, cl$distance, cl$n_genes, cluster_scaling_text(cl),
                round(100 * min(cl$support)), round(100 * max(cl$support))))
  }

  if (length(x$ci)) {
    cat("Homogeneity bootstrap (95% interval of d excludes 0 at k =):\n")
    for (nm in names(x$ci)) {
      ci <- x$ci[[nm]]
      ks <- ci$k[!ci$includes_zero]
      cat(sprintf("  %s (%d experiments): %s\n", nm, ci$m[1],
                  if (length(ks)) paste(ks, collapse = ", ") else "none"))
    }
  } else {
    cat("Homogeneity bootstrap: not run\n")
  }

  if (length(x$notes)) {
    cat("\nWarnings:\n")
    cat(paste0("  - ", x$notes, "\n"), sep = "")
  }
  invisible(x)
}

#' Settings for itv_run()
#'
#' Bundles settings of [itv_run()] to reuse them across datasets. Give only
#' the settings you want to fix; the others keep their defaults. Arguments
#' typed in an [itv_run()] call override the settings, so one set can be
#' varied run by run:
#' `itv_run(f, settings = s, boot_iter = 100)`.
#'
#' The values are checked when [itv_run()] starts, before anything is
#' written. The result of [itv_run()] holds the settings it used, except the
#' data-specific ones (file, output folder, condition, labels, experiments,
#' subsets, `n_genes`, `sep`), as `res$settings`, ready to apply to the next
#' dataset.
#'
#' @param ... Named arguments of [itv_run()] (all but `file` and
#'   `settings`). An unnamed first argument may be an `itv_settings` object
#'   to start from; the named ones then replace or add to it.
#' @return An object of class `itv_settings`, a named list.
#' @examples
#' s <- itv_settings(FDR_cutoff = 0.01, logFC_cutoff = log2(1.5),
#'                   cluster_distance = "pearson", cluster_center_genes = TRUE)
#' s
#' itv_settings(s, boot_iter = 200)         # the same, with fewer iterations
#'
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' res <- itv_run(f, settings = s, out_dir = file.path(tempdir(), "itv_settings"),
#'                cluster_boot_iter = 10, boot_iter = 50)   # these two override s
#' res
#' res$settings
#' @export
itv_settings <- function(...) {
  args <- list(...)
  nms <- names(args)
  if (is.null(nms)) nms <- rep("", length(args))
  base <- list()
  if (length(args) && !nzchar(nms[1]) && inherits(args[[1]], "itv_settings")) {
    base <- unclass(args[[1]])
    args <- args[-1]
    nms <- nms[-1]
  }
  if (any(!nzchar(nms)))
    stop("Settings must be named, e.g. itv_settings(FDR_cutoff = 0.01).", call. = FALSE)
  valid <- setdiff(names(formals(itv_run)), c("file", "settings"))
  bad <- setdiff(nms, valid)
  if (length(bad))
    stop("Not settings of itv_run(): ", paste(bad, collapse = ", "), ". See ?itv_run.",
         call. = FALSE)
  if (anyDuplicated(nms))
    stop("Repeated settings: ", paste(unique(nms[duplicated(nms)]), collapse = ", "), ".",
         call. = FALSE)
  base[nms] <- args
  structure(base, class = "itv_settings")
}

#' @export
print.itv_settings <- function(x, ...) {
  cat("<itv_settings>", if (!length(x)) "none (all defaults)", "\n")
  for (nm in names(x))
    cat("  ", nm, " = ", paste(deparse(x[[nm]]), collapse = " "), "\n", sep = "")
  invisible(x)
}

sim_titles <- c(overlap = "overlap index I", jaccard = "Jaccard index")
universe_title <- c(pair = "genes tested in both", all = "genes tested in all",
                    none = "all genes")

## The default condition: the file name without "_Profiles" and the extension.
default_condition <- function(file) {
  if (inherits(file, "itv"))
    stop("Give 'condition' (the prefix of the output files) when 'file' is an ",
         "itv object.", call. = FALSE)
  if (!is.character(file) || length(file) != 1 || is.na(file))
    stop("'file' must be a file path or an itv object.", call. = FALSE)
  sub("(_profiles)?\\.[a-z]+$", "", basename(file), ignore.case = TRUE)
}

## Tab-separated and unquoted. A tab or line break inside a value would shift
## the columns, so it is written as a space.
write_tsv <- function(x, file) {
  clean <- function(v) gsub("[\t\r\n]+", " ", v)
  chr <- vapply(x, function(v) is.character(v) || is.factor(v), logical(1))
  x[chr] <- lapply(x[chr], function(v) clean(as.character(v)))
  names(x) <- clean(names(x))
  utils::write.table(x, file, sep = "\t", quote = FALSE, row.names = FALSE)
}

## Check every itv_run() setting before any work is done.
check_run_args <- function(file, out_dir, condition, sep, FDR_cutoff, logFC_cutoff,
                           seed, chance_corrected, signed, heatmap_hclust,
                           cluster_method, cluster_min_exps, cluster_insignif_to_zero,
                           cluster_scale_exps_by_range, cluster_center_genes, cluster_scale_genes_by_SD,
                           cluster_boot_iter, bootstrap, boot_iter, boot_subsets, progress) {
  check_flag(chance_corrected, "chance_corrected")
  check_flag(signed, "signed")
  check_flag(cluster_scale_exps_by_range, "cluster_scale_exps_by_range")
  check_flag(cluster_center_genes, "cluster_center_genes")
  check_flag(cluster_scale_genes_by_SD, "cluster_scale_genes_by_SD")
  if (!inherits(file, "itv")) {
    check_string(file, "file")
    if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
  }
  check_string(out_dir, "out_dir")
  check_string(condition, "condition")
  if (grepl("[/\\\\:*?\"<>|]", condition))
    stop("'condition' is used in file names and must not contain any of / \\ : * ? \" < > |",
         call. = FALSE)
  check_string(sep, "sep")
  check_number(FDR_cutoff, "FDR_cutoff", 0, 1, min_open = TRUE)
  check_number(logFC_cutoff, "logFC_cutoff", 0)
  if (!is.null(seed)) check_number(seed, "seed", whole = TRUE)
  check_hclust_method(heatmap_hclust, "heatmap_hclust")
  check_hclust_method(cluster_method, "cluster_method")
  check_number(cluster_min_exps, "cluster_min_exps", 0, whole = TRUE)
  check_flag(cluster_insignif_to_zero, "cluster_insignif_to_zero")
  check_number(cluster_boot_iter, "cluster_boot_iter", 1, whole = TRUE)
  check_flag(bootstrap, "bootstrap")
  check_number(boot_iter, "boot_iter", 1, whole = TRUE)
  check_flag(progress, "progress")
  nms <- names(boot_subsets)
  if (!is.list(boot_subsets) ||
      (length(boot_subsets) > 0 &&
       (is.null(nms) || any(!nzchar(nms)) || anyDuplicated(nms) || "all" %in% nms)))
    stop("'boot_subsets' must be a list with a unique name, other than \"all\", for ",
         "each set, e.g. list(root = c(\"exp1\", \"exp2\")).", call. = FALSE)
  invisible(TRUE)
}

## Save the global random state; the returned function restores it.
save_rng <- function() {
  env <- globalenv()
  if (exists(".Random.seed", envir = env, inherits = FALSE)) {
    seed <- get(".Random.seed", envir = env, inherits = FALSE)
    function() assign(".Random.seed", seed, envir = env)
  } else {
    function() if (exists(".Random.seed", envir = env, inherits = FALSE))
      rm(".Random.seed", envir = env)
  }
}
