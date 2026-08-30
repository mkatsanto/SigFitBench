library(limSolve)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

mut_mat <- as.matrix(read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE))
cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep=",", row.names=1, header=TRUE, check.names=FALSE))
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

# Reimplementation of YAPSA::LCD(), which just loops lsei() over samples with
# G = identity / H = 0 (i.e. plain non-negative least squares) -- avoided
# pulling in bioconductor-yapsa itself since it has no osx-arm64 conda
# build path (r-nmf) and its other imports (e.g. BSgenome.Hsapiens.UCSC.hg19,
# KEGGREST) are unrelated to this call.
lcd <- function(mutation_catalogue, signatures) {
  n_sigs <- ncol(signatures)
  G <- diag(n_sigs)
  H <- rep(0, n_sigs)
  exposures <- matrix(0, nrow = n_sigs, ncol = ncol(mutation_catalogue),
                       dimnames = list(colnames(signatures), colnames(mutation_catalogue)))
  for (i in seq_len(ncol(mutation_catalogue))) {
    fit <- lsei(A = signatures, B = mutation_catalogue[, i], G = G, H = H, verbose = FALSE)
    exposures[, i] <- fit$X
  }
  exposures
}

res <- lcd(mut_mat, cosmic3)
res_to_save <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits=4)
write.csv(res_to_save, file=snakemake@output[["result"]])
