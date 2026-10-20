#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds the
# signature_tools.yaml conda env (Snakemake convention: a script named
# "<envfile>.post-deploy.sh" next to "<envfile>.yaml").
#
# signature.tools.lib and two of its dependencies (NNLM, indelsig.tools.lib)
# are not on conda-forge/bioconda, only on GitHub, so they can't be listed
# in the yaml itself and have to be installed this way.
#
# NMF is also installed here rather than via conda (r-nmf): conda-forge only
# ships osx-64/linux/win builds of it, no osx-arm64, so on Apple Silicon it
# has to be compiled from CRAN source instead (the "compilers" package in
# signature_tools.yaml provides the toolchain for this).
#
# The BSgenome.Hsapiens.* genome data packages (hg38, hs37d5, hg19) aren't
# listed in signature_tools.yaml either: conda's post-link step for these
# downloads the genome data with a bare curl call and no retry logic, so
# any network hiccup fails the whole env build. install_github below pulls
# them in automatically as unmet Imports of signature.tools.lib instead,
# via R's own BiocManager-backed installer (r-biocmanager is in the yaml).
#
# Uses "$CONDA_PREFIX/bin/Rscript" by absolute path rather than the bare
# "Rscript" a plain call would use: Snakemake's own conda activation appends
# (rather than prepends) this env's bin/ to PATH, so a same-named binary
# earlier on PATH (e.g. a Homebrew R install) can otherwise silently shadow
# this env's own Rscript -- installing these packages into (or loading them
# from) the wrong R entirely.
set -euo pipefail

"$CONDA_PREFIX/bin/Rscript" --vanilla -e '
  # BSgenome.Hsapiens.* genome data packages are large (hundreds of MB) and
  # Rs default 60s download.file() timeout is too short for a slow mirror.
  # 600s gives a slow connection a real chance to finish.
  options(timeout = 600)
  remotes::install_cran("NMF", upgrade = "never")
  remotes::install_github("linxihui/NNLM", upgrade = "never")
  remotes::install_github("Nik-Zainal-Group/indelsig.tools.lib", upgrade = "never")
  remotes::install_github("Nik-Zainal-Group/signature.tools.lib", upgrade = "never")
'

# Also fix Rscript resolution for future activations of this env (e.g. if
# the signature_tools rule is ever changed back to a plain `script:`
# directive instead of its current explicit-path `shell:` invocation).
mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
