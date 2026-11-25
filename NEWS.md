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
