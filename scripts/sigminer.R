# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(sigminer))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic <- read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=",")
context_file <- read.table(snakemake@input[["data"]], row.names=1, header=TRUE, check.names=FALSE)
counts <- as.matrix(context_file)
counts <- counts[rownames(cosmic), , drop = FALSE]

# Fit one sample at a time. sig_fit returns the exposures directly: mutations per signature
fit <- matrix(0, nrow=ncol(cosmic), ncol=ncol(counts),
              dimnames=list(colnames(cosmic), colnames(counts)))
for (s in colnames(counts)) {
  one <- sig_fit(counts[, s, drop=FALSE], cosmic, method=snakemake@params[["method"]])
  fit[rownames(one), s] <- one[, 1]
}
write.csv(fit, file=snakemake@output[["result"]])
