# InterTransViewer

<!-- badges: start -->
[![R-CMD-check](https://github.com/al-t1/InterTransViewer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/al-t1/InterTransViewer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

InterTransViewer describes and compares differential gene expression profiles
from several experiments, as proposed in
[Tyapkin et al. (2023)](https://doi.org/10.18699/VJGB-23-119). It helps decide
which experiments can be analysed together. For each experiment it reports:

- **n**: the number of differentially expressed genes (DEGs);
- **specific DEGs**: DEGs found in no other experiment, as a count and a
  percent, split into those every other experiment tested and those some
  experiment never tested (e.g. genes missing from a microarray);
- **R**: the specific fraction scaled by the experiment's share of all DEGs,
  and the observed / expected number of specific DEGs;
- **similarity** of DEG lists between experiments: the overlap index I and
  Jaccard, each also corrected for chance (two long lists overlap a lot even
  at random), on the genes tested in both;
- **clustering** of experiments by logFC, with bootstrap support for each clade;
- a **homogeneity bootstrap**: whether adding experiments still adds DEGs.

> **Note:** InterTransViewer compares the *results* of differential expression
> analyses (logFC and FDR per gene, e.g. from limma, DESeq2 or edgeR), not the
> expression data itself. It does not take raw counts, intensities or
> individual microarray/RNA-seq samples, and it does not normalise data or
> correct batch effects. Run the differential expression analysis of each
> experiment first, then compare the results here.

The package uses base R only, with no other dependencies.

## Installation

```r
# install.packages("remotes")
remotes::install_github("al-t1/InterTransViewer")
```

## Input

One table with a row per gene: an ID column, then a logFC and an FDR column
for each experiment. Columns can be separated by semicolons (default), commas
or tabs.

```
ID;logFC_root_1h_A;FDR_root_1h_A;logFC_root_4h_A;FDR_root_4h_A;...
GENE0001;-4.657;1.849e-51;0.03106;0.9708;...
```

Experiment names come from the logFC headers, so name your columns the way
you want experiments to appear in the results. An empty cell or `NA` means
the gene is not DE in that experiment.

### From limma, DESeq2 or edgeR results

If you have one result per experiment instead, `itv_weave()` finds the ID,
logFC and FDR columns of limma, DESeq2 and edgeR output (and of most other
tables), and joins the experiments by gene ID. Genes missing from a study are
kept, as not DE there.

```r
x <- itv_weave("DE_results/")          # a folder, file paths, or a list of data frames
x <- itv_weave("DE_results/", FDR = c(study_3 = "q_value"))   # name a column by hand
itv_write(x, "Hormone_Profiles.csv")   # save it as one profiles table
res <- itv_run(x, condition = "Hormone")
```

`itv_combine(auxin = x, ethylene = y)` joins itv objects the same way.

## Usage

All results for one file, written to a folder as tab-separated tables and PDFs:

```r
library(InterTransViewer)
f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")

res <- itv_run(f, out_dir = "Results/example",
               boot_subsets = list(root = c("root_1h_A", "root_4h_A",
                                            "root_4h_B", "root_8h_B")))
res$summary
```

Printing `res` gives a short report. To use the same settings for several
datasets, bundle them; arguments in the call still override them:

```r
s <- itv_settings(FDR_cutoff = 0.01, cluster_distance = "pearson",
                  cluster_center_genes = TRUE)
res_a <- itv_run("Auxin_Profiles.csv", settings = s)
res_b <- itv_run("Ethylene_Profiles.csv", settings = s, boot_iter = 200)
```

Or step by step:

```r
x <- itv_read(f, FDR_cutoff = 0.05)
itv_summary(x)                      # n, specific DEGs, R
itv_plot_totals(x, "R")
itv_plot_heatmap(x)                 # similarity index I, corrected for chance
itv_pairs(x)                        # all pairs: I, Jaccard, chance levels, p-values
plot(itv_cluster(x))                # clustering with bootstrap support
ci <- itv_ci(itv_bootstrap(x))      # homogeneity bootstrap
itv_plot_ci(ci)
```

### Describing the experiments

The package works with experiment names. If you have a table that describes
the experiments (one row each, with an `experiment` column of their names),
use it to label the plots, to choose experiments, and to read the results:

```r
m <- read.delim(system.file("extdata", "example_description.csv",
                            package = "InterTransViewer"), sep = ";")
labels <- itv_labels(m, c("Organ", "Time", "Method"))   # or "{Organ}, {Time} ({Method})"

x <- itv_read(f)
roots <- itv_filter(x, Organ == "Root", min_degs = 100, description = m)

res <- itv_run(f, out_dir = "Results/example", labels = labels)
merge(m, res$summary, by = "experiment")   # descriptions with the results
```

A full analysis script to copy and edit:
`system.file("scripts", "itv_analysis.R", package = "InterTransViewer")`.

The example data are simulated: eight experiments, a root and a shoot group,
and one microarray study that measures 70% of the genes
(see `data-raw/make_example.R`).

## Citation

If you use InterTransViewer, please cite:

> Tyapkin A.V., Lavrekha V.V., Ubogoeva E.V., Oshchepkov D.Yu., Omelyanchuk N.A.,
> Zemlyanskaya E.V. InterTransViewer: a comparative description of differential
> gene expression profiles from different experiments. *Vavilov Journal of
> Genetics and Breeding*. 2023;27(8):1042-1052.
> [doi:10.18699/VJGB-23-119](https://doi.org/10.18699/VJGB-23-119)

## License

GPL (>= 3). Copyright (C) 2026 A. V. Tyapkin.

This package is a new implementation of the methods of InterTransViewer
(C) 2023 Institute of Cytology and Genetics SB RAS (Rospatent registration
No. 2023681227). It does not contain code or data from that program.
