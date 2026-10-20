#!/usr/bin/env bash
# Runs once after Snakemake builds msigact.yaml: installs mSigAct and the packages it needs that
# are not available on conda, at the versions mSigAct requires.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages(c("philentropy", "sets", "cosmicsig"), dependencies = FALSE)
  remotes::install_github(c("steverozen/ICAMS@v3.0.6-branch",
                            "steverozen/PCAWG7@v0.1.3-branch",
                            "Rozen-Lab/mSigTools@v1.0.8-branch"),
                          upgrade = "never", dependencies = FALSE)
  remotes::install_github("steverozen/mSigAct@v3.0.1-branch",
                          upgrade = "never", dependencies = FALSE)
  library(mSigAct)
'

# Use this env's R first on PATH whenever the env is activated
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
