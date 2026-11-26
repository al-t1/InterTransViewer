# InterTransViewer 0.4.0

## Changes to results

* Similarity is computed on the genes tested (with an FDR) in both
  experiments of a pair, `universe = "pair"`. Before, a gene that one study
  did not test (e.g. not on a microarray) counted as "not DE" there, which
  lowered the similarity of that study to the others. `universe = "none"`
  gives the old values, `universe = "all"` uses the genes tested in every
  experiment. Pairs of experiments that tested the same genes are
  unchanged.
* Similarity is corrected for chance by default (`chance_corrected = TRUE`
  in `itv_similarity()`, `itv_plot_heatmap()` and `itv_run()`): each value
  S becomes (S - E) / (1 - E), where E is its expected value if the two DEG
  lists were drawn at random from the genes tested in both. 0 means chance
  level. Two long lists overlap a lot even at random, so this matters when
  lists are long compared with the genes tested; for short lists the
  corrected values are close to the raw ones. `chance_corrected = FALSE`
  gives the index I as published. `itv_run()` now writes
  `SimMatrix_overlap_corrected.tsv` and `Heatmap_overlap_corrected.pdf` by
  default; its `heatmap_metric` argument is replaced by `similarity_metric`
  (`"overlap"` or `"jaccard"`) and `chance_corrected`.

## New features

* `itv_pairs()` gives a table of all pairs of experiments: genes compared,
  DEGs, shared DEGs (same and opposite direction), the overlap index I (the
  Szymkiewicz-Simpson overlap coefficient) and the Jaccard index with
  their expected and chance-corrected values, the enrichment of shared DEGs
  over chance and a hypergeometric p-value. `itv_run()` returns it as
  `res$pairs`. Chance-corrected heatmaps use a -1 to 1 scale with white at
  chance level.
* `signed = TRUE` counts a gene as shared only if it changes in the same
  direction in both experiments.
* `itv_summary()` splits the specific DEGs into `specific_tested` (every
  other experiment tested the gene and found it not DE) and
  `specific_untested` (some experiment never tested it, so it may only look
  specific), and adds `specific_oe`, observed / expected specific DEGs if
  the lists were independent, next to R. `itv_specific()` gains a `tested`
  column.
* `n_genes` in `itv()`, `itv_read()`, `itv_weave()`, `itv_combine()` and
  `itv_run()`: the number of genes tested, for tables holding only DEGs,
  used by the chance levels. A table of 100 or more genes that are all DEGs
  gives a warning suggesting it.
* `itv_plot_heatmap()` orders experiments by clustering the distances
  1 - S, where S is the similarity shown, so that neighbours are
  experiments whose DEG lists overlap. Before, it clustered Euclidean
  distances between rows of the similarity matrix (`stats::heatmap()`'s
  default, as in the paper), which can put together two experiments that
  share few DEGs but are equally unlike the rest. It also draws a colour
  key, gains `order_by`, to draw several heatmaps in one order, and puts
  the title in the outer margin so that it can have two lines.
* `itv_cluster()` gains `scale_exps_by_range` (divide each experiment by
  its logFC range), `center_genes` and `scale_genes_by_SD` (subtract each
  gene's mean, divide by its SD; both together z-score the genes) and
  `distance`
  (`"euclidean"`, or 1 - r for `"pearson"` or `"spearman"`). The defaults
  are the method of the paper. Printing and plotting a clustering show how
  many genes were measured in all experiments, and the help notes that
  gene-resampling support overstates certainty. `itv_run()` passes them on
  as `cluster_scale_exps_by_range`, `cluster_center_genes`,
  `cluster_scale_genes_by_SD` and `cluster_distance`.
* `itv_settings()` bundles `itv_run()` settings to reuse across datasets:
  `itv_run(f, settings = s)`. Give only the settings you want to fix;
  arguments in the call override them. `res$settings` holds the settings a
  run used, ready for the next dataset.
* `itv_run()` results print as a report: the summary table, the
  similarity, clustering and bootstrap results, and any warnings, which are
  also kept in `res$notes`.

* `itv_weave()` builds an itv object from one differential expression
  result per experiment: files, a folder, or data frames such as
  `limma::topTable()`, `DESeq2::results()` or `edgeR::topTags()` output.
  It finds the ID, logFC and FDR columns (or takes them by name or number,
  for all experiments or per experiment) and joins the experiments by gene
  ID, keeping genes missing from some of them. Replaces `DE_weave()` of the
  0.1 script. Example results are in `inst/extdata/de_results/`.
* `itv_write()` saves an itv object as a profiles table for `itv_read()`.
* `itv_run()` also takes an itv object, e.g. from `itv_weave()`.
* `itv_filter()` keeps experiments by DEG count (`min_degs`, `max_degs`) or
  by conditions on a description table, e.g.
  `itv_filter(x, Organ == "Root", description = d)`. Replaces `KeepExps()`.
* `itv_labels()` builds plot labels from a description table, from column
  names or a template such as `"{Organ}, {Time}"`.
* `itv_combine()` joins itv objects, e.g. auxin and ethylene experiments,
  on all genes or the shared ones.

# InterTransViewer 0.3.0

## Breaking changes

* Experiment descriptions (metadata) are no longer handled by the package.
  `itv_read_metadata()`, `itv_add_metadata()` and `itv_update_metadata()`
  are removed, as are the `metadata` arguments of `itv_read()` and
  `itv_run()`, `x$metadata`, and the `Description_Updated` output.
  Experiments are named by the logFC column headers; to add a description
  to the results, join it by name, e.g.
  `merge(description, res$summary, by = "experiment")`.
* `plot.itv_cluster()` takes `labels` as a vector named by experiment (was
  unnamed, in leaf order).

## New features

* `labels` argument in `itv_plot_totals()`, `itv_plot_heatmap()`,
  `plot.itv_cluster()` and `itv_run()`: a character vector named by
  experiment, e.g. built from your own description table. Entries for other
  experiments are ignored, so one vector serves every subset. All three plots
  now label each experiment with its number and its label, and size their
  margins to the labels.
* `itv_cluster()` results record the experiment numbers (`numbers`).

## Other changes

* `itv_run()` writes its tables as tab-separated `.tsv` files without quotes
  (were semicolon-separated `.csv`), so semicolons in experiment names or
  metadata no longer shift columns. Tabs and line breaks inside values are
  written as spaces.
* `itv_run()` checks all its settings before reading or writing anything;
  before, a setting such as `cluster_min_exps = 1.5` failed at the very end,
  after most files were written. `itv()`, `itv_cluster()`, `itv_bootstrap()`
  and `itv_plot_heatmap()` check their arguments too.
* `itv_run()` warns when `out_dir` already holds files.
* `itv_read()` stops when a logFC column is followed by the FDR column of a
  different experiment, naming both columns and their positions, e.g.
  "Name mismatch between 'logFC_a' (column 2) and 'FDR_b' (column 3)".
* `itv_bootstrap(x, y)` warns when `x` and `y` were measured on different genes.
* `itv_plot_totals()` sizes the label margin to the longest experiment name.

# InterTransViewer 0.2.0

First release as an R package; a rewrite of the 2023 script (version 0.1)
using base R only.

* Same indicators as the 2023 paper: DEG counts, experiment-specific DEGs,
  R metric, similarity index I, Ward.D2 clustering with bootstrap clade
  support, and the homogeneity bootstrap.
* Reads one profiles table (ID, logFC/FDR pairs). Semicolon, comma and tab
  separators, decimal commas and Excel byte order marks are detected.
* Fixes relative to 0.1:
  - a missing logFC with a significant FDR no longer becomes logFC = 1
    (it passed any logFC cutoff below 1);
  - clustering removes `-Inf` as well as `Inf`, and uses its `method`;
  - experiments without DEGs no longer break the heatmap, and experiments
    without logFC variation are left out of clustering with a warning;
  - the bootstrap returns a matrix and is about 40 times faster, with the
    same results for the same seed;
  - results are reproducible (`seed`), and settings and session information
    are saved with the results.
* Results come from the profiles table, so DEG counts can differ slightly
  from the paper, whose supplementary tables used rounded FDR values.
