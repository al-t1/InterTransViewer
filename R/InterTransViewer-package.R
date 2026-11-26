#' InterTransViewer: compare differential expression profiles across experiments
#'
#' A typical analysis:
#' \enumerate{
#'   \item [itv_read()] reads a table of logFC/FDR pairs into an `itv` object
#'     and calls DEGs; [itv_weave()] builds one from separate limma, DESeq2
#'     or edgeR results, and [itv_write()] saves it. [itv_subset()] and
#'     [itv_filter()] keep some experiments; [itv_combine()] joins objects.
#'   \item [itv_summary()] gives the number of DEGs, the number and percent of
#'     experiment-specific DEGs, and the R metric per experiment;
#'     [itv_specific()], [itv_signif()] and [itv_joint_degs()] list the genes.
#'   \item [itv_pairs()] and [itv_similarity()] compare DEG lists (overlap
#'     index I and Jaccard, raw and chance-corrected);
#'     [itv_plot_heatmap()] draws them.
#'   \item [itv_cluster()] clusters experiments by logFC with bootstrap clade
#'     support.
#'   \item [itv_bootstrap()] and [itv_ci()] test the homogeneity of a set of
#'     experiments; [itv_plot_ci()] draws the intervals.
#' }
#' [itv_run()] does all of this for one file and writes tables and PDFs.
#'
#' Experiments are named by the logFC column headers of the profiles table.
#' The plots take an optional `labels` argument, a character vector named by
#' experiment, for longer labels built from your own description of the
#' experiments with [itv_labels()].
#'
#' @references Tyapkin A.V., Lavrekha V.V., Ubogoeva E.V., Oshchepkov D.Yu.,
#'   Omelyanchuk N.A., Zemlyanskaya E.V. InterTransViewer: a comparative
#'   description of differential gene expression profiles from different
#'   experiments. Vavilov Journal of Genetics and Breeding. 2023;27(8):1042-1052.
#'   \doi{10.18699/VJGB-23-119}
#' @keywords internal
"_PACKAGE"
