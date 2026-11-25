## Builds the simulated example files in inst/extdata/.
## Run from the package folder:  Rscript data-raw/make_example.R
##
## The data are simulated, not real measurements. Eight "experiments" on
## 1500 genes: a response shared by all experiments, a root program (four
## root experiments), a shoot program (three shoot experiments), and
## experiment-specific responses. One microarray experiment measures only
## 70% of the genes (the rest are empty), and one seedling experiment shares
## little with the others.

set.seed(2026)
G <- 1500
ids <- sprintf("GENE%04d", seq_len(G))

exps <- data.frame(
  name   = c("root_1h_A", "root_4h_A", "root_4h_B", "root_8h_B",
             "shoot_1h_C", "shoot_4h_C", "shoot_4h_D", "seedling_2h_E"),
  Study  = c("A", "A", "B", "B", "C", "C", "D", "E"),
  Method = c("RNA-seq", "RNA-seq", "Microarray", "Microarray",
             "RNA-seq", "RNA-seq", "RNA-seq", "RNA-seq"),
  Organ  = c("Root", "Root", "Root", "Root", "Shoot", "Shoot", "Shoot", "Seedling"),
  Time   = c("1h", "4h", "4h", "8h", "1h", "4h", "4h", "2h"),
  noise  = c(0.30, 0.25, 0.40, 0.40, 0.30, 0.25, 0.30, 0.35),
  stringsAsFactors = FALSE)

core  <- 1:120          # respond everywhere
root  <- 121:400        # respond in root experiments
shoot <- 401:620        # respond in shoot experiments
other <- 621:G          # background; source of experiment-specific responses
effect <- function(n) sample(c(-1, 1), n, replace = TRUE) * stats::rgamma(n, 4, 2.5)
base <- numeric(G)
base[c(core, root, shoot)] <- effect(length(c(core, root, shoot)))

logFC <- FDR <- matrix(NA_real_, G, nrow(exps))
for (j in seq_len(nrow(exps))) {
  mu <- numeric(G)
  on <- core[stats::runif(length(core)) < 0.8]
  if (exps$Organ[j] == "Root")  on <- c(on, root[stats::runif(length(root)) < 0.85])
  if (exps$Organ[j] == "Shoot") on <- c(on, shoot[stats::runif(length(shoot)) < 0.85])
  if (exps$Organ[j] == "Seedling") on <- on[stats::runif(length(on)) < 0.4]
  mu[on] <- base[on] * stats::runif(length(on), 0.6, 1.4)
  own <- sample(other, if (exps$Organ[j] == "Seedling") 120 else 25)
  mu[own] <- effect(length(own))
  se <- exps$noise[j]
  lfc <- mu + stats::rnorm(G, 0, se)
  p <- 2 * stats::pnorm(-abs(lfc) / se)
  logFC[, j] <- lfc
  FDR[, j] <- stats::p.adjust(p, "BH")
}
## The microarray study (B) measures only 70% of the genes.
not_on_array <- sample(G, round(0.3 * G))
logFC[not_on_array, exps$Study == "B"] <- NA
FDR[not_on_array, exps$Study == "B"] <- NA

tab <- data.frame(ID = ids)
for (j in seq_len(nrow(exps))) {
  tab[[paste0("logFC_", exps$name[j])]] <- signif(logFC[, j], 4)
  tab[[paste0("FDR_", exps$name[j])]]   <- signif(FDR[, j], 4)
}
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
utils::write.table(tab, "inst/extdata/example_profiles.csv", sep = ";",
                   quote = FALSE, row.names = FALSE, na = "")
utils::write.table(cbind(Order = seq_len(nrow(exps)), exps[, c("Study", "Method", "Organ", "Time")]),
                   "inst/extdata/example_description.csv", sep = ";",
                   quote = FALSE, row.names = FALSE)
