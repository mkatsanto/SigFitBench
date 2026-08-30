library(MutationalPatterns)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep="\t", header=TRUE, row.names=1, check.names=FALSE))
mut_mat <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE)
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

if (snakemake@params[["strict"]]) {
  res <- fit_to_signatures_strict(mut_mat + 0.0001, cosmic3)
  res_to_save <- as.data.frame(res$fit_res$contribution) %>% mutate_if(is.numeric, round, digits=4)
} else {
  res <- fit_to_signatures(mut_mat + 0.0001, cosmic3)
  res_to_save <- as.data.frame(res$contribution) %>% mutate_if(is.numeric, round, digits=4)
}

write.csv(res_to_save, file=snakemake@output[["result"]])
