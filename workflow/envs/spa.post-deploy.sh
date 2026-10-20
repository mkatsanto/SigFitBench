#!/usr/bin/env bash
# Runs once, automatically, right after Snakemake builds the spa.yaml conda
# env (Snakemake convention: a script named "<envfile>.post-deploy.sh" next
# to "<envfile>.yaml"). Used by both the spa and spa_low rules, which share
# this env.
#
# Snakemake's own conda activation (a bare `source <conda_root>/bin/activate
# '<prefix>'`, not the shell-integrated `conda activate` from `conda init`)
# APPENDS this env's bin/ to PATH instead of prepending it. If anything
# earlier on PATH already provides a `python` (e.g. an active virtualenv in
# the launching shell), that shadows this env's own Python entirely --
# running SPA.py under a completely different interpreter, with unrelated
# (and likely mismatched) numpy/scipy/pandas builds. That mismatch is a
# common cause of a native segfault on import, with no Python traceback at
# all, rather than a clean ImportError. This activation hook re-prepends the
# env's own bin/ on every future activation of this env, fixing Python
# resolution regardless of what else is on the user's PATH.
set -euo pipefail

mkdir -p "$CONDA_PREFIX/etc/conda/activate.d"
cat > "$CONDA_PREFIX/etc/conda/activate.d/99-fix-path-priority.sh" <<'HOOK'
export PATH="$CONDA_PREFIX/bin:$PATH"
HOOK
