# InterTransViewer

<!-- badges: start -->
[![R-CMD-check](https://github.com/al-t1/InterTransViewer/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/al-t1/InterTransViewer/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

InterTransViewer describes and compares differential gene expression profiles
from several experiments, as proposed in
[Tyapkin et al. (2023)](https://doi.org/10.18699/VJGB-23-119). It helps decide
which experiments can be analysed together. For each experiment it reports:

- **n**: the number of differentially expressed genes (DEGs);
- **specific DEGs**: DEGs found in no other experiment, as a count and a percent;
- **R**: the specific fraction scaled by the experiment's share of all DEGs;
- **similarity** of DEG lists between experiments (the overlap index I, or Jaccard);
- **clustering** of experiments by logFC, with bootstrap support for each clade;
- a **homogeneity bootstrap**: whether adding experiments still adds DEGs.

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

Experiment names come from the logFC headers. An empty cell or `NA` means the
gene is not DE in that experiment. Optionally, a second table describes the
experiments (method, organ, time, ...), one row per experiment in the same order.

## Usage

All results for one file, written to a folder:

```r
library(InterTransViewer)
f <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
m <- system.file("extdata", "example_description.csv", package = "InterTransViewer")

res <- itv_run(f, out_dir = "Results/example", metadata = m,
               boot_subsets = list(root = c("root_1h_A", "root_4h_A",
                                            "root_4h_B", "root_8h_B")))
res$summary
```

Or step by step:

```r
x <- itv_read(f, FDR_cutoff = 0.05, metadata = m)
itv_summary(x)                      # n, specific DEGs, R
itv_plot_totals(x, "R")
itv_plot_heatmap(x)                 # similarity index I
plot(itv_cluster(x))                # clustering with bootstrap support
ci <- itv_ci(itv_bootstrap(x))      # homogeneity bootstrap
itv_plot_ci(ci)
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
