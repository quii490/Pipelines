#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"

if [[ $# -eq 0 || "$1" == "-h" || "$1" == "--help" ]]; then
  cat <<'USAGE'
Usage:
  bash scripts/run_cross_contrast_heatmap.sh \
    --level-dir RESULTS/08_downstream/peak_level \
    --contrasts KO_vs_WT,Rescue_vs_KO [options]

Options:
  --top-n INT                 Number of regions to retain (default: 100)
  --outdir DIR                Output root (default: --level-dir)
  --output-prefix STR         Output filename prefix
  --annotation-mode MODE      gene_te | gene | none (default: gene_te)
USAGE
  exit 0
fi

CHIPSEQ_ENV_BIN="/path/to/.conda/envs/chipseq/bin"
if [[ -d "$CHIPSEQ_ENV_BIN" ]]; then
  export PATH="$CHIPSEQ_ENV_BIN:$PATH"
fi

exec Rscript "$ROOT_DIR/atacseq-downstream/run_cross_contrast_heatmap.R" "$@"
