# If anything below fails, write "# TOOL FAILED: <reason>" to the result file(s) and finish,
# so the rest of the workflow and the HTML report still run. The error stays in the log.
options(error = function() {
  reason <- gsub("\\s+", " ", trimws(geterrmessage()))
  for (f in unlist(snakemake@output)) writeLines(paste("# TOOL FAILED:", reason), f)
  try(sink(type = "message"), silent = TRUE)
  message("TOOL FAILED (", snakemake@rule, "): ", reason)
  quit(save = "no", status = 0)
})

suppressPackageStartupMessages(library(limSolve))
suppressPackageStartupMessages(library(dplyr))

sink(file(snakemake@log[["stdout"]], open="wt"), type="output")
sink(file(snakemake@log[["stderr"]], open="wt"), type="message")

mut_mat <- as.matrix(read.csv(snakemake@input[["data"]], sep="\t", row.names=1, header=TRUE))
cosmic3 <- as.matrix(read.table(snakemake@input[["cosmic"]], sep=",", row.names=1, header=TRUE, check.names=FALSE))
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

# The fit of YAPSA's LCD(): for each sample, least squares with no negative contributions
# (lsei with G = identity, H = 0). It is written out here because the YAPSA package itself
# cannot be installed with conda on Apple Silicon Macs.
lcd <- function(mutation_catalogue, signatures) {
  n_sigs <- ncol(signatures)
  G <- diag(n_sigs)
  H <- rep(0, n_sigs)
  exposures <- matrix(0, nrow = n_sigs, ncol = ncol(mutation_catalogue),
                       dimnames = list(colnames(signatures), colnames(mutation_catalogue)))
  for (i in seq_len(ncol(mutation_catalogue))) {
    arguments <- list(A = signatures, B = mutation_catalogue[, i], G = G, H = H, verbose = FALSE)
    # extra arguments from the config (tool_overrides: yapsa: extra) are added or replace these
    extra <- snakemake@params[["extra"]]
    arguments[names(extra)] <- extra
    fit <- do.call(lsei, arguments)
    exposures[, i] <- fit$X
  }
  exposures
}

res <- lcd(mut_mat, cosmic3)
res_to_save <- as.data.frame(res) %>% mutate_if(is.numeric, round, digits=4)
write.csv(res_to_save, file=snakemake@output[["result"]])
