#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds the mmsig.yaml
# conda env (Snakemake convention: a script named "<envfile>.post-deploy.sh"
# next to "<envfile>.yaml").
#
# mmsig is GitHub-only (evenrus/mmsig), not on CRAN/Bioconductor/conda at
# all, so it can't be listed in the yaml itself and has to be installed
# this way. It hard-Imports MutationalPatterns, which itself can't come
# from conda on osx-arm64 (no native build for its NMF dependency) -- so
# NMF is installed here from CRAN source first (the "compilers" package in
# mmsig.yaml provides the toolchain for this), then MutationalPatterns via
# BiocManager, then mmsig itself. mmsig's other Bioconductor-only Depends
# (BSgenome.Hsapiens.UCSC.hg19, TxDb.Hsapiens.UCSC.hg19.knownGene) are left
# for install_github to pull in automatically as unmet Imports, via R's own
# BiocManager-backed installer (r-biocmanager is in the yaml) -- same
# reasoning as signature_tools.post-deploy.sh.
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
  BiocManager::install("MutationalPatterns", update = FALSE, ask = FALSE)
  remotes::install_github("evenrus/mmsig", upgrade = "never")
'

# Also fix Rscript resolution for future activations of this env (the
# mmsig rule invokes R via a plain `script:` directive, which relies on
# PATH -- see the note above).
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
