#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MAIN_TEX="${1:-document.tex}"
JOB_NAME="${MAIN_TEX%.tex}"

cleanup_temp_files() {
  local base="$1"

  rm -f -- \
    "${base}.aux" \
    "${base}.bbl" \
    "${base}.bcf" \
    "${base}.blg" \
    "${base}.fdb_latexmk" \
    "${base}.fls" \
    "${base}.listing" \
    "${base}.log" \
    "${base}.out" \
    "${base}.run.xml" \
    "${base}.synctex.gz" \
    "${base}.toc"

  find pages -type f \( \
    -name '*.aux' -o \
    -name '*.listing' -o \
    -name '*.log' -o \
    -name '*.out' \
  \) -delete
}

if [[ ! -f "$MAIN_TEX" ]]; then
  echo "error: TeX file not found: $MAIN_TEX" >&2
  exit 1
fi

cleanup_temp_files "$JOB_NAME"

if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode -halt-on-error "$MAIN_TEX"
else
  if ! command -v pdflatex >/dev/null 2>&1; then
    echo "error: neither latexmk nor pdflatex is installed" >&2
    exit 1
  fi

  pdflatex -interaction=nonstopmode -halt-on-error "$MAIN_TEX"
  pdflatex -interaction=nonstopmode -halt-on-error "$MAIN_TEX"
  pdflatex -interaction=nonstopmode -halt-on-error "$MAIN_TEX"
fi

cleanup_temp_files "$JOB_NAME"
echo "built ${JOB_NAME}.pdf"
