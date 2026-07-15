#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

declare -a entries=(
  "RNA-seq/rnaseq/run_auto_rnaseq.sh"
  "ATAC-seq/run_auto_atacseq.sh"
  "CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh"
)

for entry in "${entries[@]}"; do
  test -f "$entry" || { echo "ERROR: missing entry: $entry" >&2; exit 1; }
  bash "$entry" --help >/dev/null
done

for flag in --fastq-dir --species --outdir --resume; do
  grep -R -F -q -- "$flag" docs/content/pipelines docs/content/reference/cli.md || {
    echo "ERROR: required documented flag missing: $flag" >&2
    exit 1
  }
done

echo "CLI documentation checks passed."
