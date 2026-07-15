#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<EOF
Usage: $0 --bam-glob <pattern> --outdir <dir> [options]

Generate nucleosome phasing graphs and compute nucleosome repeat length (NRL).

Required:
  --bam-glob PATTERN    Glob pattern for clean BAM files (e.g. "02_align/*.clean.bam")
  --outdir DIR          Output directory

Optional:
  --bam FILE            Single BAM file (alternative to --bam-glob)
  --bam-list FILE       File with one BAM path per line
  --labels LIST         Comma-separated sample labels (default: derived from filenames)
  --mapq INT            MAPQ filter (default: 30)
  --max-frag INT        Maximum fragment size in bp (default: 1000)
  --lspan NUM           First loess span for background fit (default: 0.35)
  --rspan NUM           Second loess span for residual fit (default: 0.1)
  --cores INT           Number of cores (default: 1)
  --conda-env NAME      Conda environment name (default: chipseq)
  --help                Show this help message
EOF
  exit 0
}

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BAM_GLOB=""
BAM_FILE=""
BAM_LIST=""
OUTDIR=""
MAPQ=30
MAX_FRAG=1000
LSPAN=0.35
RSPAN=0.1
CORES=1
CONDA_ENV="chipseq"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-glob)     BAM_GLOB="$2"; shift 2 ;;
    --bam)          BAM_FILE="$2"; shift 2 ;;
    --bam-list)     BAM_LIST="$2"; shift 2 ;;
    --outdir)       OUTDIR="$2"; shift 2 ;;
    --mapq)         MAPQ="$2"; shift 2 ;;
    --max-frag)     MAX_FRAG="$2"; shift 2 ;;
    --lspan)        LSPAN="$2"; shift 2 ;;
    --rspan)        RSPAN="$2"; shift 2 ;;
    --cores)        CORES="$2"; shift 2 ;;
    --conda-env)    CONDA_ENV="$2"; shift 2 ;;
    --help|-h)      usage ;;
    *) echo "Unknown option: $1"; usage ;;
  esac
done

if [[ -z "$OUTDIR" ]]; then
  echo "ERROR: --outdir is required" >&2
  usage
fi

mkdir -p "$OUTDIR"

RSCRIPT_BIN="/path/to/.conda/envs/${CONDA_ENV}/bin/Rscript"
NUC_PHASING_R="${PIPELINE_DIR}/atacseq-downstream/nuc_phasing.R"

RSCRIPT_ARGS=(
  --outdir "$OUTDIR"
  --mapq "$MAPQ"
  --max-frag "$MAX_FRAG"
  --lspan "$LSPAN"
  --rspan "$RSPAN"
  --cores "$CORES"
)

if [[ -n "$BAM_GLOB" ]]; then
  RSCRIPT_ARGS+=(--bam-glob "$BAM_GLOB")
fi
if [[ -n "$BAM_FILE" ]]; then
  RSCRIPT_ARGS+=(--bam "$BAM_FILE")
fi
if [[ -n "$BAM_LIST" ]]; then
  RSCRIPT_ARGS+=(--bam-list "$BAM_LIST")
fi

echo "[run_nuc_phasing] Starting nucleosome phasing analysis..."
echo "[run_nuc_phasing] Rscript: $RSCRIPT_BIN"
echo "[run_nuc_phasing] Script: $NUC_PHASING_R"
echo "[run_nuc_phasing] Outdir: $OUTDIR"

"$RSCRIPT_BIN" "$NUC_PHASING_R" "${RSCRIPT_ARGS[@]}"

echo "[run_nuc_phasing] Done. Outputs in $OUTDIR"
echo "  - *_nuc_phasing.pdf  : phasing graph per sample"
echo "  - NRL_summary.csv    : NRL table"
