# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(mmsig))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

cosmic3 <- read.table(snakemake@input[["cosmic"]], sep=",", header=TRUE, check.names=FALSE)
mut_mat <- read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE)

# Fit one sample at a time
set.seed(0)
df <- do.call(rbind, lapply(colnames(mut_mat), function(s) {
  sig_out <- mm_fit_signatures(muts.input=mut_mat[, s, drop=FALSE],
                               sig.input=cosmic3,
                               input.format="classes",
                               sample.sigt.profs=NULL,
                               strandbias=FALSE,
                               bootstrap=FALSE,
                               refcheck=FALSE,
                               cos_sim_threshold=0.01,
                               dbg=FALSE)
  sig_out$estimate
}))
# mmsig reports shares of the sample; convert to numbers of mutations
df <- df[, !names(df) %in% c("mutations")] * df$mutations
res_to_save <- as.data.frame(df) %>% mutate_if(is.numeric, round, digits=4)
write.csv(t(res_to_save), file=snakemake@output[["result"]])
