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
