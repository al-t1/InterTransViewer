## Makes inst/extdata/de_results/: one differential expression result per
## experiment, as written by limma, DESeq2 and edgeR, for the itv_weave()
## examples and tests. Run from the package root:
##   Rscript data-raw/make_de_results.R
## Needs the Bioconductor packages limma, edgeR and DESeq2 (not needed to
## use InterTransViewer).

suppressPackageStartupMessages({
  library(limma); library(edgeR); library(DESeq2)
})
set.seed(2026)
out <- file.path("inst", "extdata", "de_results")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(out, full.names = TRUE))

n <- 400
ids <- sprintf("AT%dG%04d0", rep(1:5, each = n / 5), rep(seq_len(n / 5), 5))
## Genes 1-60 respond in every experiment; 61-100 only in the root ones.
effect <- c(rep(2, 30), rep(-2, 30), rep(1.5, 40), rep(0, n - 100))
group <- factor(rep(c("control", "treated"), each = 3), levels = c("control", "treated"))
design <- model.matrix(~group)

## Counts for one RNA-seq experiment: negative binomial, 3 vs 3.
sim_counts <- function(effect) {
  mu <- 2^stats::runif(n, 3, 11)
  m <- sapply(seq_along(group), function(j)
    stats::rnbinom(n, mu = mu * 2^(effect * (group[j] == "treated")), size = 10))
  dimnames(m) <- list(ids, paste0(group, rep(1:3, 2)))
  m
}

## limma: a microarray covering 300 of the genes, root, 1 h.
## topTable() keeps the gene IDs as row names; write.csv() writes them in a
## first column with an empty header.
expr <- matrix(stats::rnorm(n * 6, 8, 0.4), n, 6, dimnames = list(ids, NULL))
expr[, group == "treated"] <- expr[, group == "treated"] + effect
on_array <- sort(sample(n, 300))
fit <- eBayes(lmFit(expr[on_array, ], design))
tt <- topTable(fit, coef = 2, number = Inf)
write.csv(signif(tt, 5), file.path(out, "root_1h_limma.csv"))

## DESeq2: root, 4 h. as.data.frame(results()) written with write.csv();
## some padj are NA (independent filtering and outliers).
dds <- DESeqDataSetFromMatrix(sim_counts(effect), data.frame(group), ~group)
dds <- DESeq(dds, quiet = TRUE)
res <- as.data.frame(results(dds))
write.csv(signif(res, 5), file.path(out, "root_4h_DESeq2.csv"))

## edgeR: shoot, 4 h; only the shared response. Lowly expressed genes are
## removed by filterByExpr(); topTags()$table written with write.table(),
## tab-separated, so the header has no field for the row names.
y <- DGEList(sim_counts(ifelse(seq_len(n) <= 60, effect, 0)), group = group)
y <- y[filterByExpr(y), , keep.lib.sizes = FALSE]
y <- normLibSizes(y)
qlf <- glmQLFTest(glmQLFit(y, design), coef = 2)
tab <- topTags(qlf, n = Inf)$table
write.table(signif(tab, 5), file.path(out, "shoot_4h_edgeR.tsv"), sep = "\t")

for (f in list.files(out, full.names = TRUE))
  cat(basename(f), ":", length(readLines(f)) - 1, "genes,", file.size(f), "bytes\n")
