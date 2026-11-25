#' Run the whole analysis on one file
#'
#' Reads a profiles table, computes all indicators, and writes tables
#' (semicolon-separated) and PDFs to `out_dir`, with names starting with
#' `condition`:
#' `SignifTable.csv`, `JointDEGs.csv`, `TotalSpecificR.csv`,
#' `SpecificDEGs.csv`, `TotalSpecificDEG.pdf`, `Rmetric.pdf`,
#' `Description_Updated.csv` (with `metadata`), `SimMatrix_<metric>.csv`,
#' `Heatmap_<metric>.pdf`, `Clustering_min<min_exps>.pdf`,
#' `Bootstrap[_<subset>].rds`, `Bootstrap_CI[_<subset>].csv`,
#' `Bootstrap_CI.pdf` and `run_info.txt` (settings and session info).
#'
#' Problems in one step (metadata, clustering, a bootstrap subset) give a
#' warning and the other steps still run.
#'
#' @param file Path to the profiles table, see [itv_read()].
#' @param out_dir Output folder; created if needed.
#' @param condition Prefix for output files. By default the file name
#'   without "_Profiles" and the extension.
#' @param sep Column separator of `file`, see [itv_read()].
#' @param metadata Optional path of an experiment description table, see
#'   [itv_add_metadata()]. Rows follow the experiments in `file` (before
#'   `experiments` is applied) or are matched by a column of experiment names.
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
#' m <- system.file("extdata", "example_description.csv", package = "InterTransViewer")
#' res <- itv_run(f, out_dir = file.path(tempdir(), "itv_example"), metadata = m,
#'                cluster_boot_iter = 20, boot_iter = 100,
#'                boot_subsets = list(root = c("root_1h_A", "root_4h_A",
#'                                             "root_4h_B", "root_8h_B")))
#' res$summary
#' list.files(res$out_dir)
#' @export
itv_run <- function(file,
                    out_dir = file.path("Results", condition),
                    condition = sub("(_profiles)?\\.[a-z]+$", "", basename(file),
                                    ignore.case = TRUE),
                    sep = ";", metadata = NULL, experiments = NULL,
                    FDR_cutoff = 0.05, logFC_cutoff = 0, seed = 42,
                    heatmap_metric = c("overlap", "jaccard"), heatmap_hclust = "complete",
                    cluster_method = "ward.D2", cluster_min_exps = 0,
                    cluster_insignif_to_zero = FALSE, cluster_boot_iter = 100,
                    bootstrap = TRUE, boot_iter = 1000, boot_subsets = list(),
                    progress = interactive()) {
  heatmap_metric <- match.arg(heatmap_metric)
  force(condition)
  if (!is.null(seed)) {
    restore_rng <- save_rng()
    on.exit(restore_rng(), add = TRUE)
    set.seed(seed)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out <- function(name) file.path(out_dir, paste0(condition, "_", name))
  pdf_to <- function(name, width, height) {
    grDevices::pdf(out(name), width = width, height = height)
    dev <- grDevices::dev.cur()
    function() grDevices::dev.off(dev)
  }

  ## 1. Data, metadata and DEG calls.
  x <- itv_read(file, sep = sep, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff)
  if (!is.null(metadata)) {
    x <- tryCatch(itv_add_metadata(x, metadata), error = function(e) {
      warning("Metadata not used: ", conditionMessage(e), call. = FALSE)
      x
    })
  }
  if (!is.null(experiments)) x <- itv_subset(x, experiments)
  deg <- x$deg
  no_fc <- sum(deg & is.na(x$logFC))
  if (no_fc > 0)
    warning(no_fc, " DEG calls have a significant FDR but no logFC; ",
            "they are counted as DEGs, with logFC shown as NA.", call. = FALSE)
  no_deg <- colnames(deg)[colSums(deg) == 0]
  if (length(no_deg) > 0)
    warning("Experiments without DEGs: ", paste(no_deg, collapse = ", "), call. = FALSE)

  ## 2. Gene tables.
  signif <- itv_signif(x)
  write_csv2(data.frame(ID = rownames(signif), signif, check.names = FALSE),
             out("SignifTable.csv"))
  write_csv2(itv_joint_degs(x), out("JointDEGs.csv"))

  ## 3. Totals, specific DEGs, R and metadata.
  tsr <- itv_summary(x)
  write_csv2(tsr, out("TotalSpecificR.csv"))
  write_csv2(itv_specific(x), out("SpecificDEGs.csv"))
  height <- max(4, 0.28 * nrow(tsr) + 2)
  close_dev <- pdf_to("TotalSpecificDEG.pdf", 10, height)
  itv_plot_totals(tsr, "specific_percent", condition)
  close_dev()
  close_dev <- pdf_to("Rmetric.pdf", 10, height)
  itv_plot_totals(tsr, "R", condition)
  close_dev()
  if (!is.null(x$metadata))
    write_csv2(itv_update_metadata(x), out("Description_Updated.csv"))

  ## 4. Similarity of DEG lists.
  sim <- itv_similarity(x, heatmap_metric)
  write_csv2(data.frame(experiment = rownames(sim), round(sim, 4), check.names = FALSE),
             out(paste0("SimMatrix_", heatmap_metric, ".csv")))
  if (sum(!is.na(diag(sim))) >= 2) {
    close_dev <- pdf_to(paste0("Heatmap_", heatmap_metric, ".pdf"), 8, 7.5)
    itv_plot_heatmap(sim, hclust_method = heatmap_hclust,
                     main = paste0(condition, ": ", heatmap_metric))
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
    labels <- sprintf("%d %s (%d DEGs, %.1f%% spec.)", tsr$N, tsr$experiment,
                      tsr$total, tsr$specific_percent)[match(cl$experiments, tsr$experiment)]
    close_dev <- pdf_to(paste0("Clustering_min", cluster_min_exps, ".pdf"),
                        max(8, 0.5 * length(labels) + 3), 9)
    plot(cl, labels = labels,
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
      write_csv2(ci, out(paste0("Bootstrap_CI", tag, ".csv")))
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
               sprintf("metadata = %s", if (is.null(metadata)) "NULL" else metadata),
               sprintf("experiments = %s", paste(deparse(experiments), collapse = "")),
               sprintf("FDR_cutoff = %g; logFC_cutoff = %g; seed = %s",
                       FDR_cutoff, logFC_cutoff, if (is.null(seed)) "NULL" else seed),
               sprintf("heatmap: %s, %s", heatmap_metric, heatmap_hclust),
               sprintf("clustering: %s, min_exps = %d, insignif_to_zero = %s, boot_iter = %d",
                       cluster_method, cluster_min_exps, cluster_insignif_to_zero,
                       cluster_boot_iter),
               sprintf("bootstrap: %s, iter = %d", bootstrap, boot_iter),
               sprintf("boot_subsets = %s", paste(deparse(boot_subsets), collapse = "")),
               "", utils::capture.output(utils::sessionInfo())),
             out("run_info.txt"))

  invisible(list(itv = x, summary = tsr, similarity = sim, cluster = cl,
                 bootstrap = boots, ci = cis, out_dir = out_dir))
}

write_csv2 <- function(x, file) {
  utils::write.table(x, file, sep = ";", quote = FALSE, row.names = FALSE)
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
