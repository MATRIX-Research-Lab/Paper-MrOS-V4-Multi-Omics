#!/usr/bin/env bash
#
# Build the LaTeX-generated manuscript figures.
#
#   bash code/latex_workflow/run_all.sh
#
# This is the LaTeX entry point, and it owns exactly the figures listed below.
# It is the counterpart of code/single_omics/run_all.R and of
# code/cross_omics/02_make_figures.R: each entry point regenerates only its
# own outputs and leaves the others alone. All three write into results/,
# which is the single source the manuscript is copied from.
#
# Requires a TeX distribution with TikZ and the standalone class.

set -euo pipefail

FIGURES=(
  Supplementary_Fig5_primary_modeling_workflow
  Supplementary_Fig6_matched_cohort_workflow
)

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
outdir="$repo/results/figures"

[[ -f "$repo/code/single_omics/00_setup.R" ]] || {
  echo "Not an MrOS repository root: $repo" >&2
  exit 2
}

engine=""
for candidate in latexmk pdflatex; do
  if command -v "$candidate" >/dev/null 2>&1; then engine="$candidate"; break; fi
done
if [[ -z "$engine" ]]; then
  echo "No LaTeX engine on PATH (looked for latexmk, pdflatex)." >&2
  echo "Install TeX Live, or leave the existing PDFs in results/figures in place:" >&2
  printf '  %s\n' "${FIGURES[@]/#/$outdir/}" >&2
  exit 127
fi

mkdir -p "$outdir"

# Build in a scratch directory so aux files never land in the repository, then
# move only the PDF into results/figures.
build="$(mktemp -d)"
trap 'rm -rf "$build"' EXIT

cd "$here"
for fig in "${FIGURES[@]}"; do
  echo "building $fig ..."
  if [[ "$engine" == "latexmk" ]]; then
    latexmk -pdf -interaction=nonstopmode -halt-on-error \
            -outdir="$build" "$fig.tex" >"$build/$fig.build.log" 2>&1 || {
      echo "LaTeX failed for $fig; last lines:" >&2
      tail -25 "$build/$fig.build.log" >&2
      exit 1
    }
  else
    pdflatex -interaction=nonstopmode -halt-on-error \
             -output-directory="$build" "$fig.tex" >"$build/$fig.build.log" 2>&1 || {
      echo "LaTeX failed for $fig; last lines:" >&2
      tail -25 "$build/$fig.build.log" >&2
      exit 1
    }
  fi
  [[ -s "$build/$fig.pdf" ]] || { echo "no PDF produced for $fig" >&2; exit 1; }
  mv "$build/$fig.pdf" "$outdir/$fig.pdf"
  printf '  -> results/figures/%s.pdf (%s bytes)\n' "$fig" "$(wc -c <"$outdir/$fig.pdf" | tr -d ' ')"
done

echo "Done. ${#FIGURES[@]} figure(s) written to results/figures/"
