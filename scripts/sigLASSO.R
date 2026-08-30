library(siglasso)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic <- read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=",")
context_file <- read.table(snakemake@input[["data"]], row.names=1, header=TRUE, check.names=FALSE)
counts <- as.matrix(context_file)
counts <- counts[rownames(cosmic), , drop = FALSE]

fit <- siglasso(counts, cosmic)
res_round <- as.data.frame(fit$exposure) %>% mutate_if(is.numeric, round, digits=4)
write.csv(res_round, file=snakemake@output[["result"]])
