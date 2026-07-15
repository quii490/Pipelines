#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

declare -a entries=(
  "RNA-seq/rnaseq/run_auto_rnaseq.sh|docs/content/pipelines/rnaseq/parameters.md"
  "ATAC-seq/run_auto_atacseq.sh|docs/content/pipelines/atacseq/parameters.md"
  "CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh|docs/content/pipelines/cutrun/parameters.md"
)

for item in "${entries[@]}"; do
  entry="${item%%|*}"
  doc="${item#*|}"
  test -f "$entry" || { echo "ERROR: missing entry: $entry" >&2; exit 1; }
  test -f "$doc" || { echo "ERROR: missing parameter page: $doc" >&2; exit 1; }

  help_text="$(bash "$entry" --help 2>&1)"
  while IFS= read -r flag; do
    grep -F -q -- "$flag" "$doc" || {
      echo "ERROR: $entry documents $flag in --help, but $doc does not mention it" >&2
      exit 1
    }
  done < <(printf '%s\n' "$help_text" | grep -oE -- '--[A-Za-z0-9][A-Za-z0-9-]*' | sort -u)
done

echo "CLI documentation checks passed."
