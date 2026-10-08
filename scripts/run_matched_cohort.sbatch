#!/usr/bin/env bash
#
# MrOS matched-cohort secondary analysis.
#
# Refits every omics layer on the 332 matched men. Reads the re-fitted layer
# outputs under data/model_outputs/single_omics_rerun/ -- not the primary
# analysis outputs; see code/single_omics/00_setup.R on why the two are kept apart.
#
# Submit when those outputs are already in place:
#   sbatch scripts/run_matched_cohort.sbatch
#
# Set MROS_RENDER=true to render the release artifacts after the matched fit.

#SBATCH --job-name=mros-matched
#
# ---------------------------------------------------------------------------
# SITE SETTINGS
# ---------------------------------------------------------------------------
#SBATCH --partition=standard
#SBATCH --mail-user=hlin1239@umd.edu
# Runtime as D-HH:MM:SS. Keep comments on their own line: Slurm's #SBATCH
# parser is not documented to accept trailing comments, and a --time that
# fails to parse falls back to the partition default, which would kill a long
# run part-way through.
#SBATCH --time=2-00:00:00
#
# ---------------------------------------------------------------------------
# Allocation shape.
#
# This job runs a single R process that fans out to future workers, so one task
# with many CPUs is the correct shape here -- unlike the primary job, which
# runs three concurrent steps and therefore needs --ntasks=3.
#
# This job is NOT split by layer, unlike the primary one. The stacking models
# need the clinical and per-layer out-of-fold predictions together, in one
# process, on the same outer folds, and the stage writes a single
# provenance-validated artifact. Splitting it would mean inventing a
# partial-results contract and an assembly step -- restructuring model-fitting
# code for a scheduling gain that checkpointing already provides
# (MROS_MATCHED_RESUME=true resumes a validated interrupted run).
#
# --mem is deliberately not set: Zaratan grants 4 GiB per requested core, so 64
# cores already carries 256 GiB. If this job also sits in the queue, the lever
# is fewer cores -- MROS_MATCHED_WORKERS follows SLURM_CPUS_PER_TASK, so
# halving --cpus-per-task roughly doubles the run time but schedules far sooner.
# ---------------------------------------------------------------------------
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --requeue
#SBATCH --mail-type=BEGIN,TIME_LIMIT_90,END,FAIL
#SBATCH --output=mros_matched_%j.out
#SBATCH --error=mros_matched_%j.err

set -euo pipefail

# ---------------------------------------------------------------------------
# Runtime. A Slurm job starts a fresh shell, so conda must be activated here.
# MROS_R_MODULE now defaults to empty: the previous default (R/4.2.2) is too old
# for assert_modeling_api(), which requires dials::grid_space_filling() and
# tune::control_grid(parallel_over=).
# ---------------------------------------------------------------------------
CONDA_ROOT="${MROS_CONDA_ROOT:-$HOME/scratch/miniconda3}"
CONDA_ENV="${MROS_CONDA_ENV:-$HOME/scratch/conda/envs/r-4.4.3}"
if [[ -r "$CONDA_ROOT/etc/profile.d/conda.sh" ]]; then
  # shellcheck disable=SC1091
  source "$CONDA_ROOT/etc/profile.d/conda.sh"
  conda activate "$CONDA_ENV"
fi

R_MODULE="${MROS_R_MODULE-}"
if [[ -n "$R_MODULE" ]] && command -v module >/dev/null 2>&1; then
  module load "$R_MODULE"
fi
command -v Rscript >/dev/null || {
  echo "Rscript is not on PATH; check MROS_CONDA_ENV or set MROS_R_MODULE." >&2
  exit 127
}

REPO_ROOT="${MROS_REPO_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}}"
cd "$REPO_ROOT"
[[ -f code/single_omics/00_setup.R ]] || {
  echo "Not an MrOS repository root: $REPO_ROOT" >&2
  exit 2
}

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-64}"
MATCHED_WORKERS="${MROS_MATCHED_WORKERS:-$((TOTAL_CPUS - 4))}"
TUNE_PARALLEL_OVER="${MROS_TUNE_PARALLEL_OVER:-everything}"
SEED="${MROS_SEED:-2026}"
RESUME="${MROS_MATCHED_RESUME:-false}"
RENDER="${MROS_RENDER:-false}"

if [[ ! "$TOTAL_CPUS" =~ ^[0-9]+$ || ! "$MATCHED_WORKERS" =~ ^[0-9]+$ ]]; then
  echo "CPU settings must be positive integers." >&2
  exit 2
fi
if (( MATCHED_WORKERS < 1 || MATCHED_WORKERS >= TOTAL_CPUS )); then
  echo "MROS_MATCHED_WORKERS must be at least 1 and less than $TOTAL_CPUS." >&2
  exit 2
fi
if [[ "$TUNE_PARALLEL_OVER" != "resamples" && "$TUNE_PARALLEL_OVER" != "everything" ]]; then
  echo "MROS_TUNE_PARALLEL_OVER must be resamples or everything." >&2
  exit 2
fi
if [[ "$RESUME" != "true" && "$RESUME" != "false" ]]; then
  echo "MROS_MATCHED_RESUME must be true or false." >&2
  exit 2
fi
if [[ "$RENDER" != "true" && "$RENDER" != "false" ]]; then
  echo "MROS_RENDER must be true or false." >&2
  exit 2
fi

export MROS_REPO_ROOT="$REPO_ROOT"
export MROS_TUNE_PARALLEL_OVER="$TUNE_PARALLEL_OVER"
export MROS_MATCHED_WORKERS="$MATCHED_WORKERS"
export MROS_MATCHED_RESUME="$RESUME"
export MROS_MATCHED_SEED="$SEED"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_FUTURE_GLOBALS_MAXSIZE="${R_FUTURE_GLOBALS_MAXSIZE:-2000000000}"

# Keep large temporary files off the home filesystem.
# Temp files go to node-local storage, not the shared allocation. It is faster,
# it does not consume the project quota, and it keeps scratch I/O off a
# filesystem that has repeatedly produced intermittent errors on this cluster --
# "gzfile(): cannot open the connection" from inside tune being the latest.
export TMPDIR="${SLURM_TMPDIR:-/tmp}/mros_${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# Raise the open-file limit to the hard ceiling. Twenty forked workers each hold
# connections, and the default soft limit (commonly 1024) can be exhausted --
# which also surfaces as "cannot open the connection".
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
echo "TMPDIR: $TMPDIR | open-file limit: $(ulimit -n)"

echo "repository: $REPO_ROOT"
echo "git commit: $(git rev-parse HEAD 2>/dev/null || echo unavailable)"
echo "total CPUs: $TOTAL_CPUS"
echo "matched workers: $MATCHED_WORKERS"
echo "tune parallelization: $TUNE_PARALLEL_OVER"
echo "matched resume: $RESUME"
echo "render after fitting: $RENDER"
echo "seed: $SEED"
date

Rscript -e 'source("code/single_omics/00_setup.R"); print(utils::sessionInfo()); assert_packages(matched_packages, "MrOS matched workflow")'
Rscript code/single_omics/preflight.R --stage matched
Rscript -e 'source("code/single_omics/00_setup.R"); source("code/single_omics/06_fit_matched_cohort_models.R")'

if [[ "$RENDER" == "true" ]]; then
  Rscript code/single_omics/preflight.R --stage render
  Rscript code/single_omics/run_all.R
fi

echo "MrOS matched-cohort analysis completed successfully."
date
