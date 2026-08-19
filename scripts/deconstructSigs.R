library(deconstructSigs)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], sep="\t", row.names=1, header=TRUE, check.names=FALSE)
cosmic3 <- as.data.frame(t(cosmic3))
mut_matrix <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE, check.names=FALSE)
mut_matrix <- as.data.frame(t(mut_matrix))

stopifnot(identical(colnames(cosmic3), colnames(mut_matrix)))

df <- data.frame()
for (sample in rownames(mut_matrix)) {
  cat("Processing sample", sample, "\n")
  res <- whichSignatures(tumor.ref=mut_matrix, sample.id=sample, signatures.ref=cosmic3,
                         contexts.needed=TRUE, signature.cutoff=0)
  vals <- as.data.frame(c(res$weights[1,]))
  rownames(vals) <- sample
  df <- rbind(df, vals)
}
df <- df %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(df), file=snakemake@output[["result"]])
