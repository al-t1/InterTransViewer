## InterTransViewer: step-by-step analysis script.
## Copy this file next to your data and edit the settings. Every step can be
## run and inspected on its own; itv_run() does the same in one call.
## Find this file with: system.file("scripts", "itv_analysis.R", package = "InterTransViewer")

library(InterTransViewer)

## ===== Settings ==============================================================
input_file   <- system.file("extdata", "example_profiles.csv", package = "InterTransViewer")
FDR_cutoff   <- 0.05
logFC_cutoff <- 0        # e.g. log2(1.5)
set.seed(42)

## Optional plot labels, named by experiment, from your own description of
## the experiments. Use labels <- NULL to label plots with experiment names.
description <- read.delim(system.file("extdata", "example_description.csv",
                                      package = "InterTransViewer"), sep = ";")
labels <- itv_labels(description, c("Organ", "Time", "Method"))

## ===== 1. Data ===============================================================
x <- itv_read(input_file, FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff)
## Or, from one limma/DESeq2/edgeR result per experiment in a folder:
## x <- itv_weave("DE_results", FDR_cutoff = FDR_cutoff, logFC_cutoff = logFC_cutoff)
x
## Keep some experiments:
## x <- itv_subset(x, c("root_1h_A", "root_4h_A"))
## x <- itv_filter(x, Organ == "Root", min_degs = 100, description = description)

## ===== 2. DEG counts, specific DEGs, R =======================================
tsr <- itv_summary(x)
tsr
merge(description, tsr, by = "experiment")          # descriptions with the results
head(itv_specific(x))
itv_plot_totals(x, "specific_percent", labels = labels)
itv_plot_totals(x, "R", labels = labels)

## ===== 3. Similarity of DEG lists ============================================
## On the genes tested in both experiments of each pair (universe = "pair");
## signed = TRUE counts only genes changing in the same direction.
pairs <- itv_pairs(x)
head(pairs[order(-pairs$overlap_corrected), ])   # the pairs most above chance
## Overlap index I corrected for chance (0 = chance level); metric = "jaccard"
## for Jaccard, chance_corrected = FALSE for the raw values.
sim <- itv_similarity(x, "overlap", chance_corrected = TRUE)
round(sim, 2)
itv_plot_heatmap(sim, labels = labels)

## ===== 4. Clustering of experiments ==========================================
cl <- itv_cluster(x, min_exps = 0, boot_iter = 100)
plot(cl, labels = labels)

## ===== 5. Homogeneity bootstrap ==============================================
boot_all  <- itv_bootstrap(x, iter = 1000)
boot_root <- itv_bootstrap(itv_subset(x, c("root_1h_A", "root_4h_A", "root_4h_B",
                                           "root_8h_B")), iter = 1000)
ci <- list(all = itv_ci(boot_all), root = itv_ci(boot_root))
ci
itv_plot_ci(ci)

## Two sets of the same size, e.g. two conditions:
## d <- itv_bootstrap(x_auxin, x_ethylene, size1 = 16, sizes2 = 16)
## itv_plot_hist(d[, 1], "16 auxin - 16 ethylene")
