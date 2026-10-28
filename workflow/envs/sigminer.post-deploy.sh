#!/usr/bin/env bash
# Runs once, right after Snakemake creates the sigminer env (Snakemake runs
# "<env>.post-deploy.sh" automatically). Installs maftools (Bioconductor), NMF and sigminer.
# The last line stops with an error if anything failed, so a broken env is not used.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  BiocManager::install("maftools", update = FALSE, ask = FALSE)
  remotes::install_cran(c("NMF", "sigminer"), upgrade = "never")
  library(maftools); library(NMF); library(nnls); library(sigminer)
'

# Put this env's R first on PATH when the env is activated, so another R installed on the
# machine (e.g. from Homebrew) is never used instead
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
