#!/usr/bin/env bash
# Runs once, right after Snakemake creates the sigfit env (Snakemake runs
# "<env>.post-deploy.sh" automatically). sigfit is only on GitHub, so it is installed here.
# The last line stops with an error if the install failed, so a broken env is not used.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  remotes::install_github("kgori/sigfit", upgrade = "never")
  library(sigfit)
'

# Put this env's R first on PATH when the env is activated, so another R installed on the
# machine (e.g. from Homebrew) is never used instead
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
