#!/usr/bin/env Rscript
# mSigAct SparseAssignActivity: fit all signatures by maximum likelihood (negative binomial by
# default), then remove signatures while a likelihood-ratio test (p.thresh) finds the fit not
# significantly worse. Output: mutations per signature (signatures x samples).

# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

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

sig_catalog <- ICAMS::as.catalog(as.matrix(sigs[order96, , drop = FALSE]),
                                 region = "genome", catalog.type = "counts.signature")

# Fit one sample at a time
out <- matrix(0, nrow = ncol(sigs), ncol = ncol(counts),
              dimnames = list(colnames(sigs), colnames(counts)))
for (s in colnames(counts)) {
  arguments <- list(
    spectra = ICAMS::as.catalog(as.matrix(counts[order96, s, drop = FALSE]),
                                region = "genome", catalog.type = "counts"),
    sigs = sig_catalog,
    output.dir = file.path("msigact_output", s),
    max.level = snakemake@params[["max_level"]],
    p.thresh = snakemake@params[["p_thresh"]],
    m.opts = DefaultManyOpts(likelihood.dist = snakemake@params[["likelihood"]]),
    num.parallel.samples = 1,
    mc.cores.per.sample = snakemake@threads,
    seed = 0,
    drop.low.mut.samples = FALSE
  )
  # extra arguments from the config (tool_overrides: msigact: extra) are added or replace these
  extra <- snakemake@params[["extra"]]
  if (length(extra) > 0) arguments[names(extra)] <- extra
  res <- do.call(SparseAssignActivity, arguments)
  if (!is.null(res$error.messages)) message("mSigAct errors for ", s, ": ", paste(res$error.messages, collapse = "; "))

  # mSigAct gives no solution for a sample when it reaches its search limit (max.subsets),
  # and then returns missing values or all zeros; report that as a failure instead of
  # writing zeros
  assignment <- if (is.null(res$proposed.assignment)) NULL else as.matrix(res$proposed.assignment)
  if (is.null(assignment) || !(s %in% colnames(assignment)) ||
      (sum(assignment[, s], na.rm = TRUE) == 0 && sum(counts[, s]) > 0)) {
    stop("mSigAct found no solution for sample ", s, ": ",
         if (is.null(res$error.messages)) "no reason given" else res$error.messages[[1]])
  }

  # proposed.assignment holds mutation counts for the selected signatures only
  assignment[is.na(assignment)] <- 0
  out[rownames(assignment), s] <- assignment[, s]
}
write.csv(round(out, 4), file = snakemake@output[["result"]])
