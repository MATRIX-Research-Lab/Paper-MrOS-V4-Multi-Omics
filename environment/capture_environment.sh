#!/usr/bin/env bash
#
# Verify and record the cluster analysis environment.
#
# Run from an activated conda environment on the cluster:
#   conda activate "$HOME/scratch/conda/envs/r-4.4.3"
#   bash environment/capture_environment.sh
#
# Refuses to record an environment that does not satisfy the package and API
# contracts in code/single_omics/00_setup.R, so a recorded environment is
# always one the pipeline can actually run.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_ROOT"

OUT="$REPO_ROOT/environment"
mkdir -p "$OUT"

if [[ -z "${CONDA_PREFIX:-}" ]]; then
  echo "ERROR: no conda environment is active." >&2
  echo "  conda activate \"\$HOME/scratch/conda/envs/r-4.4.3\"" >&2
  exit 2
fi

R_BIN="$(command -v R || true)"
if [[ -z "$R_BIN" ]]; then
  echo "ERROR: R is not on PATH." >&2
  exit 2
fi
case "$R_BIN" in
  "$CONDA_PREFIX"/*) ;;
  *)
    echo "ERROR: R at $R_BIN is not the conda R in $CONDA_PREFIX." >&2
    echo "A cluster R module is probably loaded. Start a clean shell." >&2
    exit 2
    ;;
esac

echo "==> conda environment: $CONDA_PREFIX"
echo "==> R: $R_BIN"

# 1. Exact, hash-pinned environment. linux-64 only; use this to rebuild.
echo "==> writing conda-linux-64-r-4.4.3.lock"
conda list --prefix "$CONDA_PREFIX" --explicit > "$OUT/conda-linux-64-r-4.4.3.lock"

# 2. Human-readable, cross-platform spec of what was explicitly requested.
echo "==> writing conda-r-4.4.3.yml"
conda env export --prefix "$CONDA_PREFIX" --from-history > "$OUT/conda-r-4.4.3.yml"

# 3. Verify the contracts, then record the R side. Must run after step 1:
#    the R script reads the lock to work out which packages came from conda
#    and which came from CRAN.
echo "==> verifying contracts and recording R packages"
Rscript "$SCRIPT_DIR/record_r_packages.R"

# 4. Provenance of the code and platform this environment goes with.
echo "==> writing MANIFEST.txt"
{
  echo "MrOS v4 mortality — cluster environment record"
  echo
  echo "captured           : $(date -u '+%Y-%m-%dT%H:%M:%SZ') (UTC)"
  echo "host               : $(hostname)"
  echo "kernel             : $(uname -sr)"
  if [[ -r /etc/os-release ]]; then
    echo "distribution       : $(. /etc/os-release && echo "$PRETTY_NAME")"
  fi
  echo "conda prefix       : $CONDA_PREFIX"
  echo "conda version      : $(conda --version 2>/dev/null || echo unknown)"
  echo "R                  : $("$R_BIN" --version | head -1)"
  echo "git commit         : $(git rev-parse HEAD 2>/dev/null || echo unavailable)"
  echo "git branch         : $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo unavailable)"
  if ! git diff-index --quiet HEAD -- 2>/dev/null; then
    echo "git state          : DIRTY — uncommitted changes present"
  else
    echo "git state          : clean"
  fi
} > "$OUT/MANIFEST.txt"

echo
echo "Environment recorded under environment/:"
ls -1 "$OUT"
