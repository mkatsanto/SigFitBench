#!/usr/bin/env Rscript
# mSigAct SparseAssignActivity: fit all signatures by maximum likelihood (negative binomial by
# default), then remove signatures while a likelihood-ratio test (p.thresh) finds the fit not
# significantly worse. Output: mutations per signature (signatures x samples).

sink(file(snakemake@log[["stdout"]], open = "wt"), type = "output")
sink(file(snakemake@log[["stderr"]], open = "wt"), type = "message")

suppressPackageStartupMessages({
  library(mSigAct)
  library(ICAMS)
})

counts <- read.table(snakemake@input[["data"]], sep = "\t", header = TRUE, row.names = 1,
                     check.names = FALSE)
sigs <- read.table(snakemake@input[["signature_db"]], sep = "\t", header = TRUE, row.names = 1,
                   check.names = FALSE)

# ICAMS names the 96 SBS channels as 5' base, reference, 3' base, alternative:
# "A[C>A]A" -> "ACAA"
to_icams <- function(x) paste0(substr(x, 1, 1), substr(x, 3, 3), substr(x, 7, 7), substr(x, 5, 5))
order96 <- ICAMS::catalog.row.order$SBS96
rownames(counts) <- to_icams(rownames(counts))
rownames(sigs) <- to_icams(rownames(sigs))
if (!setequal(rownames(counts), order96) || !setequal(rownames(sigs), order96)) {
  stop("Mutation channels do not match the 96 ICAMS SBS channels")
}

spectra <- ICAMS::as.catalog(as.matrix(counts[order96, , drop = FALSE]),
                             region = "genome", catalog.type = "counts")
sig_catalog <- ICAMS::as.catalog(as.matrix(sigs[order96, , drop = FALSE]),
                                 region = "genome", catalog.type = "counts.signature")

res <- SparseAssignActivity(
  spectra = spectra,
  sigs = sig_catalog,
  output.dir = "msigact_output",
  max.level = snakemake@params[["max_level"]],
  p.thresh = snakemake@params[["p_thresh"]],
  m.opts = DefaultManyOpts(likelihood.dist = snakemake@params[["likelihood"]]),
  num.parallel.samples = 1,
  mc.cores.per.sample = snakemake@threads,
  seed = 0,
  drop.low.mut.samples = FALSE
)
if (!is.null(res$error.messages)) message("mSigAct errors: ", paste(res$error.messages, collapse = "; "))

# proposed.assignment holds mutation counts for the selected signatures only
assignment <- as.matrix(res$proposed.assignment)
assignment[is.na(assignment)] <- 0
out <- matrix(0, nrow = ncol(sigs), ncol = ncol(counts),
              dimnames = list(colnames(sigs), colnames(counts)))
out[rownames(assignment), colnames(assignment)] <- assignment
write.csv(round(out, 4), file = snakemake@output[["result"]])
