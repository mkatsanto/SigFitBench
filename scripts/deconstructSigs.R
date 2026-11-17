# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(deconstructSigs))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], sep="\t", row.names=1, header=TRUE, check.names=FALSE)
cosmic3 <- as.data.frame(t(cosmic3))
mut_matrix <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE, check.names=FALSE)
mut_matrix <- as.data.frame(t(mut_matrix))

# Match the mutation types by name, so the order of the two files does not matter
stopifnot(setequal(colnames(cosmic3), colnames(mut_matrix)))
mut_matrix <- mut_matrix[, colnames(cosmic3), drop = FALSE]

df <- data.frame()
for (sample in rownames(mut_matrix)) {
  cat("Processing sample", sample, "\n")
  arguments <- list(tumor.ref=mut_matrix, sample.id=sample, signatures.ref=cosmic3,
                    contexts.needed=TRUE, signature.cutoff=0)
  # extra arguments from the config (tool_overrides: deconstructsigs: extra) are added or replace these
  extra <- snakemake@params[["extra"]]
  arguments[names(extra)] <- extra
  res <- do.call(whichSignatures, arguments)
  # deconstructSigs reports shares of the sample; convert to numbers of mutations
  vals <- as.data.frame(c(res$weights[1,] * sum(mut_matrix[sample, ])))
  rownames(vals) <- sample
  df <- rbind(df, vals)
}
df <- df %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(df), file=snakemake@output[["result"]])
