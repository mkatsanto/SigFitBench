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
  arguments <- list(catalogue_matrix=counts[, s, drop=FALSE], sig=cosmic, method=snakemake@params[["method"]])
  # extra arguments from the config (tool_overrides: sigminer: extra) are added or replace these
  extra <- snakemake@params[["extra"]]
  arguments[names(extra)] <- extra
  one <- do.call(sig_fit, arguments)
  fit[rownames(one), s] <- one[, 1]
}
write.csv(fit, file=snakemake@output[["result"]])
