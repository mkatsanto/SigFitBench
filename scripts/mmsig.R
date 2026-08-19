library(mmsig)
library(dplyr)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], sep=",", header=TRUE, check.names=FALSE)
mut_mat <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE)

set.seed(0)
sig_out <- mm_fit_signatures(muts.input=mut_mat,
                             sig.input=cosmic3,
                             input.format="classes",
                             sample.sigt.profs=NULL,
                             strandbias=FALSE,
                             bootstrap=FALSE,
                             refcheck=FALSE,
                             cos_sim_threshold=0.01,
                             dbg=FALSE)
df <- sig_out$estimate
df <- df[, !names(df) %in% c("mutations")]
res_to_save <- as.data.frame(df) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(res_to_save), file=snakemake@output[["result"]])
