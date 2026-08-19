#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds the
# sigminer.yaml conda env (Snakemake convention: a script named
# "<envfile>.post-deploy.sh" next to "<envfile>.yaml").
#
# sigminer hard-Imports NMF, which conda-forge only ships for
# osx-64/linux/win (no osx-arm64), so conda refuses to solve
# bioconductor/conda builds that depend on it at all on Apple Silicon.
# NMF is installed here from CRAN source instead (the "compilers" package
# in sigminer.yaml provides the toolchain for this), then sigminer itself
# is installed from CRAN -- once NMF already exists as a normal installed
# R package, install.packages() doesn't care that it didn't come from conda.
#
# Uses "$CONDA_PREFIX/bin/Rscript" by absolute path rather than the bare
# "Rscript" a plain call would use: Snakemake's own conda activation appends
# (rather than prepends) this env's bin/ to PATH, so a same-named binary
# earlier on PATH (e.g. a Homebrew R install) can otherwise silently shadow
# this env's own Rscript -- installing these packages into (or loading them
# from) the wrong R entirely.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  remotes::install_cran("NMF", upgrade = "never")
  remotes::install_cran("sigminer", upgrade = "never")
'

# Also fix Rscript resolution for future activations of this env (the
# sigminer rule invokes R via a plain `script:` directive, which relies on
# PATH -- see the note above).
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
