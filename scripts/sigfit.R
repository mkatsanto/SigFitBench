library(dplyr)
library(sigfit)

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=",")
context_file <- read.table(snakemake@input[["data"]], header=TRUE, check.names=FALSE)
counts <- build_catalogues(context_file)

fit <- fit_signatures(counts=counts, signatures=cosmic3, chains=1, seed=0)
exposures <- retrieve_pars(fit, par="exposures", hpd_prob=0.90)

x <- as.data.frame(exposures$mean) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(x), file=snakemake@output[["result"]])

y <- as.data.frame(exposures$lower_90) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(y), file=snakemake@output[["lower90"]])
