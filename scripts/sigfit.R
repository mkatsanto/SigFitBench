# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(sigfit))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=",")
mut_mat <- read.table(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE, check.names=FALSE)

# sigfit wants samples (counts) and signatures as rows, mutation types as columns, same order
counts <- t(as.matrix(mut_mat[rownames(cosmic3), , drop=FALSE]))
signatures <- t(as.matrix(cosmic3))

# Fit one sample at a time
mean_rows <- list()
lower_rows <- list()
for (s in rownames(counts)) {
  fit <- fit_signatures(counts=counts[s, , drop=FALSE], signatures=signatures, chains=1, seed=0)
  exposures <- retrieve_pars(fit, par="exposures", hpd_prob=0.90)
  mean_rows[[s]] <- as.data.frame(exposures$mean, row.names=s)
  lower_rows[[s]] <- as.data.frame(exposures$lower_90, row.names=s)
}

# sigfit reports shares of each sample; convert to numbers of mutations
totals <- rowSums(counts)
x <- as.data.frame(do.call(rbind, mean_rows) * totals) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(x), file=snakemake@output[["result"]])

y <- as.data.frame(do.call(rbind, lower_rows) * totals) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(y), file=snakemake@output[["lower90"]])
