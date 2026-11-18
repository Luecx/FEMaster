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

write_version_file() {
  local header="../src/core/version.h"
  local major minor patch

  major="$(sed -nE 's/^[[:space:]]*#[[:space:]]*define[[:space:]]+VERSION_MAJOR[[:space:]]+([0-9]+).*/\1/p' "$header")"
  minor="$(sed -nE 's/^[[:space:]]*#[[:space:]]*define[[:space:]]+VERSION_MINOR[[:space:]]+([0-9]+).*/\1/p' "$header")"
  patch="$(sed -nE 's/^[[:space:]]*#[[:space:]]*define[[:space:]]+VERSION_PATCH[[:space:]]+([0-9]+).*/\1/p' "$header")"

  if [[ -z "$major" || -z "$minor" || -z "$patch" ]]; then
    echo "error: could not read version from $header" >&2
    exit 1
  fi

  printf '\\newcommand{\\FEMasterVersion}{%s.%s.%s}\n' "$major" "$minor" "$patch" > version.tex
}

cleanup_temp_files "$JOB_NAME"
write_version_file

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
