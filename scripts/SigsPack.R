# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(SigsPack))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep=",", row.names=1, header=TRUE, check.names=FALSE))
mut_mat <- as.matrix(read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE))
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

# Fit one sample at a time
exposures <- matrix(0, nrow=ncol(cosmic3), ncol=ncol(mut_mat),
                    dimnames=list(colnames(cosmic3), colnames(mut_mat)))
for (s in colnames(mut_mat)) {
  one <- signature_exposure(mut_mat[, s, drop=FALSE], P=cosmic3)$exposures
  # SigsPack reports shares of the sample; convert to numbers of mutations
  exposures[rownames(one), s] <- one[, 1] * sum(mut_mat[, s])
}
res_to_save <- as.data.frame(exposures) %>% mutate_if(is.numeric, round, digits=4)
colnames(res_to_save) <- paste0("S", seq(0, ncol(res_to_save) - 1))
write.csv(res_to_save, file=snakemake@output[["result"]])
