suppressMessages(library(optparse))

option_list <- list(
  make_option("--data", type = "character"),
  make_option("--cosmic", type = "character"),
  make_option("--output", type = "character"),
  make_option("--method", type = "character", default = "KLD"),
  make_option("--filter-type", type = "character", default = "giniScaledThreshold", dest = "filter_type"),
  make_option("--nboot", type = "integer", default = 100),
  make_option("--nparallel", type = "integer", default = 1)
)
opt <- parse_args(OptionParser(option_list = option_list))

# Safety net for signature.tools.lib (installed via GitHub, not conda): if
# for any reason it isn't already installed, install it here instead.
# Parallel job instances of this rule can hit this at the same time, so a
# lock directory guards against concurrent installs corrupting each other in
# the shared conda env. The lock is reclaimed if it's older than
# STALE_LOCK_SECS: a from-source install of these packages should never take
# that long, so an older lock means its owner was killed before it could
# clean up (e.g. an interrupted run) rather than one still legitimately
# working -- otherwise every future run would wait out the timeout and
# silently fall through to the library() call below, failing with a
# misleading "no package called" error instead of actually installing.
if (!requireNamespace("signature.tools.lib", quietly = TRUE)) {
  lock_dir <- file.path(Sys.getenv("CONDA_PREFIX"), ".signature_tools_lib_install.lock")
  STALE_LOCK_SECS <- 900
  acquired <- FALSE
  for (i in 1:600) {
    if (dir.create(lock_dir, showWarnings = FALSE)) {
      acquired <- TRUE
      break
    }
    if (requireNamespace("signature.tools.lib", quietly = TRUE)) break
    lock_age <- tryCatch(
      as.numeric(Sys.time() - file.info(lock_dir)$mtime, units = "secs"),
      error = function(e) 0
    )
    if (!is.na(lock_age) && lock_age > STALE_LOCK_SECS) {
      message("[signature_tools] reclaiming stale install lock (", round(lock_age), "s old)")
      unlink(lock_dir, recursive = TRUE)
      next
    }
    Sys.sleep(1)
  }
  if (acquired) {
    on.exit(unlink(lock_dir, recursive = TRUE), add = TRUE)
    if (!requireNamespace("signature.tools.lib", quietly = TRUE)) {
      remotes::install_cran("NMF", upgrade = "never")
      remotes::install_github("linxihui/NNLM", upgrade = "never")
      remotes::install_github("Nik-Zainal-Group/indelsig.tools.lib", upgrade = "never")
      remotes::install_github("Nik-Zainal-Group/signature.tools.lib", upgrade = "never")
    }
  }
}

suppressPackageStartupMessages(library(signature.tools.lib))
suppressPackageStartupMessages(library(dplyr))

mut_mat <- read.csv(opt$data, sep = "\t", row.names = 1, header = TRUE)
cosmic3 <- as.matrix(read.table(opt$cosmic, sep = "\t", header = TRUE, row.names = 1, check.names = FALSE))
mut_mat <- mut_mat[rownames(cosmic3), , drop = FALSE]

# Fit one sample at a time
exposures <- do.call(rbind, lapply(colnames(mut_mat), function(s) {
  res <- Fit(catalogues = mut_mat[, s, drop = FALSE], signatures = cosmic3,
             method = opt$method,
             exposureFilterType = opt$filter_type,
             useBootstrap = TRUE, nboot = opt$nboot, nparallel = opt$nparallel)
  one <- res$exposures[1, , drop = FALSE]
  rownames(one) <- s
  one
}))
res_to_save <- t(exposures)
res_to_save <- res_to_save[!(row.names(res_to_save) == "unassigned"), , drop = FALSE]
res_to_save <- as.data.frame(res_to_save) %>% mutate_if(is.numeric, round, digits = 4)
write.csv(res_to_save, file = opt$output)
