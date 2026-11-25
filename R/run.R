#' Run the whole analysis on one file
#'
#' Reads a profiles table, computes all indicators, and writes tables
#' (tab-separated, unquoted) and PDFs to `out_dir`, with names starting with
#' `condition`:
#' `SignifTable.tsv`, `JointDEGs.tsv`, `TotalSpecificR.tsv`,
#' `SpecificDEGs.tsv`, `TotalSpecificDEG.pdf`, `Rmetric.pdf`,
#' `SimMatrix_<metric>.tsv`,
#' `Heatmap_<metric>.pdf`, `Clustering_min<min_exps>.pdf`,
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
#' @param file Path to the profiles table, see [itv_read()].
#' @param out_dir Output folder; created if needed.
#' @param condition Prefix for output files. By default the file name
#'   without "_Profiles" and the extension.
#' @param sep Column separator of `file`, see [itv_read()].
#' @param labels Optional plot labels: a character vector named by
#'   experiment, see [itv_plot_totals()]. It must have an entry for every
#'   experiment analysed; other entries are ignored.
#' @param experiments Optional experiments to analyse, by name or number,
#'   see [itv_subset()].
#' @param FDR_cutoff,logFC_cutoff DEG definition, see [itv()].
#' @param seed Random seed for the bootstraps, or `NULL` to use the current
#'   random state. The previous random state is restored afterwards.
#' @param heatmap_metric `"overlap"` or `"jaccard"`.
#' @param heatmap_hclust Agglomeration method used to order the heatmap.
#' @param cluster_method,cluster_min_exps,cluster_insignif_to_zero,cluster_boot_iter
#'   Passed to [itv_cluster()] as `method`, `min_exps`, `insignif_to_zero`
#'   and `boot_iter`.
#' @param bootstrap Run the homogeneity bootstrap.
#' @param boot_iter Iterations for each k.
#' @param boot_subsets Named list of extra sets of experiments to bootstrap,
#'   each a character vector of experiment names. A set with unknown names is
#'   skipped with a warning.
#' @param progress Show bootstrap progress bars.
#' @return Invisibly, a list with `itv` (the data), `summary`, `similarity`,
#'   `cluster`, `bootstrap` and `ci` (lists by set) and `out_dir`.
#' @examples
#' f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
#' m <- read.delim(system.file("extdata", "example_description.csv",
#'                             package = "InterTransViewer"), sep = ";")
#' res <- itv_run(f, out_dir = file.path(tempdir(), "itv_example"),
#'                labels = setNames(paste(m$Organ, m$Time, m$Method), m$experiment),
#'                cluster_boot_iter = 20, boot_iter = 100,
#'                boot_subsets = list(root = c("root_1h_A", "root_4h_A",
#'                                             "root_4h_B", "root_8h_B")))
#' res$summary
#' list.files(res$out_dir)
#'
#' ## The description with this run's results
#' merge(m, res$summary, by = "experiment")
#' @export
itv_run <- function(file,
                    out_dir = file.path("Results", condition),
                    condition = sub("(_profiles)?\\.[a-z]+$", "", basename(file),
                                    ignore.case = TRUE),
                    sep = ";", labels = NULL, experiments = NULL,
                    FDR_cutoff = 0.05, logFC_cutoff = 0, seed = 42,
                    heatmap_metric = c("overlap", "jaccard"), heatmap_hclust = "complete",
                    cluster_method = "ward.D2", cluster_min_exps = 0,
                    cluster_insignif_to_zero = FALSE, cluster_boot_iter = 100,
                    bootstrap = TRUE, boot_iter = 1000, boot_subsets = list(),
                    progress = interactive()) {
  heatmap_metric <- match.arg(heatmap_metric)
  force(condition)
  check_run_args(file = file, out_dir = out_dir, condition = condition, sep = sep,
                 FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff, seed = seed,
                 heatmap_hclust = heatmap_hclust, cluster_method = cluster_method,
                 cluster_min_exps = cluster_min_exps,
                 cluster_insignif_to_zero = cluster_insignif_to_zero,
                 cluster_boot_iter = cluster_boot_iter, bootstrap = bootstrap,
                 boot_iter = boot_iter, boot_subsets = boot_subsets, progress = progress)
  if (!is.null(labels)) check_labels(labels)

  ## 1. Data and DEG calls; labels are checked against the experiments before
  ## anything is written.
  x <- itv_read(file, sep = sep, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff)
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

  ## 4. Similarity of DEG lists.
  sim <- itv_similarity(x, heatmap_metric)
  write_tsv(data.frame(experiment = rownames(sim), round(sim, 4), check.names = FALSE),
             out(paste0("SimMatrix_", heatmap_metric, ".tsv")))
  if (sum(!is.na(diag(sim))) >= 2) {
    side <- max(7.5, 5 + 0.1 * max(nchar(exp_labels)))
    close_dev <- pdf_to(paste0("Heatmap_", heatmap_metric, ".pdf"), side + 0.5, side)
    itv_plot_heatmap(sim, hclust_method = heatmap_hclust,
                     main = paste0(condition, ": ", heatmap_metric), labels = exp_labels)
    close_dev()
  } else {
    itv_plot_heatmap(sim)   # gives the warnings, draws nothing
  }

  ## 5. Clustering.
  cl <- tryCatch(
    itv_cluster(x, min_exps = cluster_min_exps, insignif_to_zero = cluster_insignif_to_zero,
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
         main = sprintf(paste0("%s: %s on %d genes (DE in >= %d experiments; %d of %d ",
                               "genes measured in all experiments)\n",
                               "red numbers: bootstrap support, %%, %d replicates"),
                        condition, cluster_method, cl$n_genes, cluster_min_exps,
                        cl$n_complete, nrow(deg), cluster_boot_iter))
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
               sprintf("file = %s", file),
               sprintf("labels = %s", if (is.null(labels)) "NULL" else
                         paste(paste0(names(exp_labels), " = ", exp_labels), collapse = "; ")),
               sprintf("experiments = %s", paste(deparse(experiments), collapse = "")),
               sprintf("FDR_cutoff = %g; logFC_cutoff = %g; seed = %s",
                       FDR_cutoff, logFC_cutoff, if (is.null(seed)) "NULL" else seed),
               sprintf("heatmap: %s, %s", heatmap_metric, heatmap_hclust),
               sprintf("clustering: %s, min_exps = %s, insignif_to_zero = %s, boot_iter = %s",
                       cluster_method, cluster_min_exps, cluster_insignif_to_zero,
                       cluster_boot_iter),
               sprintf("bootstrap: %s, iter = %s", bootstrap, boot_iter),
               sprintf("boot_subsets = %s", paste(deparse(boot_subsets), collapse = "")),
               "", utils::capture.output(utils::sessionInfo())),
             out("run_info.txt"))

  invisible(list(itv = x, summary = tsr, similarity = sim, cluster = cl,
                 bootstrap = boots, ci = cis, out_dir = out_dir))
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
                           seed, heatmap_hclust, cluster_method, cluster_min_exps,
                           cluster_insignif_to_zero, cluster_boot_iter, bootstrap,
                           boot_iter, boot_subsets, progress) {
  check_string(file, "file")
  if (!file.exists(file)) stop("File not found: ", file, call. = FALSE)
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
