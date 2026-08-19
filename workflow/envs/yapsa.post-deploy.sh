#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds this conda env
# (Snakemake convention: a script named "<envfile>.post-deploy.sh" next to
# "<envfile>.yaml").
#
# Snakemake's own conda activation (a bare `source <conda_root>/bin/activate
# '<prefix>'`, not the shell-integrated `conda activate` from `conda init`)
# APPENDS this env's bin/ to PATH instead of prepending it. Any same-named
# binary earlier on PATH (e.g. a Homebrew R install's Rscript) then silently
# shadows this env's own Rscript for every rule using the `script:`
# directive -- causing "there is no package called X" errors for packages
# that ARE installed, just loaded under the wrong R (and potentially
# segfaulting, since compiled packages aren't ABI-compatible across R
# builds). This activation hook re-prepends the env's own bin/ on every
# future activation of this env, fixing Rscript resolution regardless of
# what else is on the user's PATH.
set -euo pipefail

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
