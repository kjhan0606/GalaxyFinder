#!/usr/bin/env bash
set -euo pipefail

if [[ ! -f main.tex || ! -f RELEASE_SHA256SUMS ]]; then
  printf '%s\n' 'Run this script from the paper directory containing main.tex.' >&2
  exit 2
fi

sha256sum -c RELEASE_SHA256SUMS
supplement_check_output="$(cd supplement && sha256sum -c SHA256SUMS)" || {
  printf '%s\n' "$supplement_check_output" >&2
  exit 1
}
printf 'Verified %s packaged files with SHA-256.\n' "$(wc -l < supplement/SHA256SUMS)"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplement/extended_notes.tex
cp -- extended_notes.pdf supplement/extended_notes.pdf
