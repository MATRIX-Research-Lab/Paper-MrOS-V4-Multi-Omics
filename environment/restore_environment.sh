#!/usr/bin/env bash
#
# Rebuild the recorded cluster analysis environment on a linux-64 machine.
#
#   bash environment/restore_environment.sh "$HOME/scratch/conda/envs/r-4.4.3"
#
# Recreates the conda environment from the hash-pinned lock, then installs the
# CRAN-only packages at their recorded versions. Verify afterwards with
# environment/capture_environment.sh, which re-checks the contracts.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

TARGET_PREFIX="${1:-}"
if [[ -z "$TARGET_PREFIX" ]]; then
  echo "usage: bash environment/restore_environment.sh <conda-env-prefix>" >&2
  exit 2
fi

LOCK="$SCRIPT_DIR/conda-linux-64-r-4.4.3.lock"
[[ -r "$LOCK" ]] || { echo "ERROR: missing $LOCK" >&2; exit 2; }

if [[ "$(uname -s)-$(uname -m)" != "Linux-x86_64" ]]; then
  echo "WARNING: the explicit lock is linux-64 only." >&2
  echo "On another platform, solve environment/conda-r-4.4.3.yml instead." >&2
  echo "Package versions will differ and results may not match exactly." >&2
  exit 2
fi

command -v conda >/dev/null || { echo "ERROR: conda is not on PATH." >&2; exit 2; }

# Refuse to touch an existing environment. conda create would fail here anyway,
# but with an error that reads like a bug rather than a decision. Rebuilding is
# for a fresh machine or a deliberate reset -- an environment that already
# works needs nothing from this script.
if [[ -e "$TARGET_PREFIX/conda-meta" ]]; then
  echo "ERROR: a conda environment already exists at:" >&2
  echo "         $TARGET_PREFIX" >&2
  echo >&2
  echo "This script builds an environment; it does not update one in place." >&2
  echo "If the environment already works, you do not need to run this at all." >&2
  echo "To verify it instead:" >&2
  echo "         bash environment/capture_environment.sh" >&2
  echo "To rebuild from scratch, remove it first, or restore to a new prefix:" >&2
  echo "         conda env remove --prefix \"$TARGET_PREFIX\"" >&2
  exit 2
fi

echo "==> creating conda environment at $TARGET_PREFIX"
conda create --yes --prefix "$TARGET_PREFIX" --file "$LOCK"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$TARGET_PREFIX"

echo "==> installing CRAN-only packages at recorded versions"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript "$SCRIPT_DIR/restore_cran_packages.R"

echo
echo "Environment restored. Verify with:"
echo "  bash environment/capture_environment.sh"
