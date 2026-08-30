library(sigminer)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic <- read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=",")
context_file <- read.table(snakemake@input[["data"]], row.names=1, header=TRUE, check.names=FALSE)
counts <- as.matrix(context_file)
counts <- counts[rownames(cosmic), , drop = FALSE]

fit <- sig_fit(counts, cosmic, method=snakemake@params[["method"]])
write.csv(fit$exposure, file=snakemake@output[["result"]])
