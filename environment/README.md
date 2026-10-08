# Analysis environment

Everything needed to rebuild the R environment that produced the cluster
results: **R 4.4.3 with Bioconductor 3.20, installed through conda**. See
`../README.md` for how the analysis itself is run.

Three operations, easy to confuse:

| Script | Does | When |
|---|---|---|
| `restore_environment.sh` | Builds a conda environment from the lock | Fresh machine, or a deliberate reset |
| `restore_cran_packages.R` | Installs the required CRAN packages into the active environment | Repair: something is missing |
| `capture_environment.sh` | Verifies, then records | After any change; before submitting jobs |

Day to day you only need the last one.

## Build the environment

On a linux-64 machine with conda installed:

```bash
bash environment/restore_environment.sh "$HOME/scratch/conda/envs/r-4.4.3"
```

This recreates the conda environment from the hash-pinned lock, then installs
the CRAN-only packages at their recorded versions. 

On macOS or another platform the explicit lock does not apply. Solve
`conda-r-4.4.3.yml` instead, and expect different package versions and
therefore results that do not match to the last digit.

## Verify and record

```bash
conda activate "$HOME/scratch/conda/envs/r-4.4.3"
bash environment/capture_environment.sh
```

The capture verifies before it records, so a clean run is itself the check. It
confirms that R is the conda R rather than a cluster module, that the package
lists and modelling-API requirements in `../code/single_omics/00_setup.R` are
satisfied, and that every installed package is complete on disk. If anything
fails it stops with the specific problem and writes nothing.

Re-run it after any change to the environment, and commit the result alongside
whatever change prompted it.

## Repair an existing environment

`capture_environment.sh` verifies and records; it never installs. If it stops
with missing packages, put them back without rebuilding conda:

```bash
Rscript environment/restore_cran_packages.R
```

Then re-run the capture.

## Two things to keep in mind

Model fitting needs **forked** workers (`future::multicore`). The limma screen
is applied through a custom `recipes` step whose S3 methods live in the global
environment, and socket workers cannot dispatch them.
`analysis_parallel_plan()` in `../code/single_omics/00_setup.R` selects the
backend.

Keep BLAS and OpenMP single-threaded, and set them as **shell exports**:

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
```

## Files

| File | What it is |
|---|---|
| `conda-linux-64-r-4.4.3.lock` | Exact conda packages with URLs and hashes. Rebuilds the environment on linux-64. From `conda list --explicit`. |
| `conda-r-4.4.3.yml` | The explicitly requested packages, without pinned builds. The only option off linux-64. From `conda env export --from-history`. |
| `r-packages.csv` | Every installed R package with version, build, and whether it came from conda or CRAN. Base and recommended packages are omitted; those follow the R version. |
| `cran-packages.csv` | The packages the conda lock cannot rebuild *and* the analysis requires. Consumed by `restore_cran_packages.R`. |
| `r-session-info.txt` | `sessionInfo()`, including the BLAS and LAPACK in use, plus the thread-limit variables. |
| `MANIFEST.txt` | Capture time, host, kernel, distribution, R version, and the Git commit and working-tree state of the analysis code. |

## What this does not cover

The environment is half of reproducibility. The other half:

- **Analysis code** — the Git commit recorded in `MANIFEST.txt`.
- **Preprocessing** — runs on an approved local machine against
  controlled-access MrOS data, under a different R. Recorded separately in
  `../results/manifests/`.
- **Input data** — the three frozen modeling frames, with checksums in
  `../README.md` and provenance in `../results/manifests/`.
- **Model outputs** — each carries its own `*.manifest.csv`.

MrOS participant data is controlled-access and is not distributed with this
repository, nor deposited elsewhere. See "Data availability" in `../README.md`
for how investigators with study access reproduce the analysis.
