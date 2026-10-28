# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(MutationalPatterns))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep="\t", header=TRUE, row.names=1, check.names=FALSE))
mut_mat <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE)
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

# Fit one sample at a time
fit <- matrix(0, nrow=ncol(cosmic3), ncol=ncol(mut_mat),
              dimnames=list(colnames(cosmic3), colnames(mut_mat)))
for (s in colnames(mut_mat)) {
  one <- as.matrix(mut_mat[, s, drop=FALSE]) + 0.0001
  if (snakemake@params[["strict"]]) {
    contribution <- fit_to_signatures_strict(one, cosmic3)$fit_res$contribution
  } else {
    contribution <- fit_to_signatures(one, cosmic3)$contribution
  }
  fit[rownames(contribution), s] <- contribution[, 1]
}

res_to_save <- as.data.frame(fit) %>% mutate_if(is.numeric, round, digits=4)
write.csv(res_to_save, file=snakemake@output[["result"]])
