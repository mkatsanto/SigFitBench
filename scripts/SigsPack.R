library(SigsPack)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep=",", header=TRUE, check.names=FALSE))
mut_mat <- as.matrix(read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE))

res <- signature_exposure(mut_mat, P=cosmic3)
res_to_save <- as.data.frame(res$exposures) %>% mutate_if(is.numeric, round, digits=4)
colnames(res_to_save) <- paste0("S", seq(0, ncol(res_to_save) - 1))
write.csv(res_to_save, file=snakemake@output[["result"]])
