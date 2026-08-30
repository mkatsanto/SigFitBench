#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds the sigfit.yaml
# conda env (Snakemake convention: a script named "<envfile>.post-deploy.sh"
# next to "<envfile>.yaml").
#
# sigfit is GitHub-only (kgori/sigfit), not on CRAN/Bioconductor/conda at
# all, so it can't be listed in the yaml itself and has to be installed
# this way. Its own .stan model files are compiled at install time via
# rstan/StanHeaders (provided by conda in sigfit.yaml) -- this is
# inherently slow and CPU-heavy, not something this script can speed up.
# build_vignettes is left off (the package default) to avoid also needing
# pandoc just to install the package.
#
# Uses "$CONDA_PREFIX/bin/Rscript" by absolute path rather than the bare
# "Rscript" a plain call would use: Snakemake's own conda activation appends
# (rather than prepends) this env's bin/ to PATH, so a same-named binary
# earlier on PATH (e.g. a Homebrew R install) can otherwise silently shadow
# this env's own Rscript -- installing this package into (or loading it
# from) the wrong R entirely.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  remotes::install_github("kgori/sigfit", upgrade = "never")
'

# Also fix Rscript resolution for future activations of this env (the
# sigfit rule invokes R via a plain `script:` directive, which relies on
# PATH -- see the note above).
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
