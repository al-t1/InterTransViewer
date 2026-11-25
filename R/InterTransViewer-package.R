#' InterTransViewer: compare differential expression profiles across experiments
#'
#' A typical analysis:
#' \enumerate{
#'   \item [itv_read()] reads a table of logFC/FDR pairs into an `itv` object
#'     and calls DEGs; [itv_subset()] keeps some experiments;
#'     [itv_add_metadata()] attaches experiment descriptions.
#'   \item [itv_summary()] gives the number of DEGs, the number and percent of
#'     experiment-specific DEGs, and the R metric per experiment;
#'     [itv_specific()], [itv_signif()] and [itv_joint_degs()] list the genes.
#'   \item [itv_similarity()] compares DEG lists (overlap index I or Jaccard);
#'     [itv_plot_heatmap()] draws them.
#'   \item [itv_cluster()] clusters experiments by logFC with bootstrap clade
#'     support.
#'   \item [itv_bootstrap()] and [itv_ci()] test the homogeneity of a set of
#'     experiments; [itv_plot_ci()] draws the intervals.
#' }
#' [itv_run()] does all of this for one file and writes tables and PDFs.
#'
#' @references Tyapkin A.V., Lavrekha V.V., Ubogoeva E.V., Oshchepkov D.Yu.,
#'   Omelyanchuk N.A., Zemlyanskaya E.V. InterTransViewer: a comparative
#'   description of differential gene expression profiles from different
#'   experiments. Vavilov Journal of Genetics and Breeding. 2023;27(8):1042-1052.
#'   \doi{10.18699/VJGB-23-119}
#' @keywords internal
"_PACKAGE"
