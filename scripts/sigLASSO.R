# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(siglasso))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

# siglasso expects the signatures as a numeric matrix
cosmic <- as.matrix(read.table(snakemake@input[["cosmic"]], row.names=1, header=TRUE, check.names=FALSE, sep=","))
context_file <- read.table(snakemake@input[["data"]], row.names=1, header=TRUE, check.names=FALSE)
counts <- as.matrix(context_file)
counts <- counts[rownames(cosmic), , drop = FALSE]

# Fit one sample at a time (equivalent to the batch call, since siglasso fits each sample
# independently), so a sample that fails gets all-zero contributions instead of stopping the run.
sig_names <- colnames(cosmic)
zero_col <- setNames(rep(0, length(sig_names)), sig_names)

fit_cols <- lapply(colnames(counts), function(s) {
  res <- tryCatch(
    siglasso(counts[, s, drop = FALSE], cosmic, plot = FALSE),
    error = function(e) {
      message(sprintf("siglasso failed on sample %s: %s -- writing all-zero contribution", s, conditionMessage(e)))
      NULL
    }
  )
  # siglasso returns proportions; convert them to mutation counts
  if (is.null(res)) zero_col else res[sig_names, 1] * sum(counts[, s])
})

fit <- do.call(cbind, fit_cols)
colnames(fit) <- colnames(counts)
rownames(fit) <- sig_names

res_round <- as.data.frame(fit) %>% mutate_if(is.numeric, round, digits=4)
write.csv(res_round, file=snakemake@output[["result"]])
