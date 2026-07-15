#!/usr/bin/env bash
set -euo pipefail

BAM_DIR=""
OUT_PREFIX=""
THREADS=8
BIN_SIZE="${BAM_COR_BIN_SIZE:-5000}"

usage() { cat <<'USAGE'
Usage: cutrun_bam_cor.sh --bam-dir DIR --out-prefix PREFIX [--threads N] [--bin-size N]
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done
[[ -d "$BAM_DIR" ]] || { echo "ERROR: BAM directory missing: $BAM_DIR" >&2; exit 1; }
[[ -n "$OUT_PREFIX" ]] || { echo "ERROR: --out-prefix is required" >&2; exit 1; }
command -v multiBamSummary >/dev/null 2>&1 || { echo "ERROR: deepTools multiBamSummary not found" >&2; exit 1; }
command -v plotCorrelation >/dev/null 2>&1 || { echo "ERROR: deepTools plotCorrelation not found" >&2; exit 1; }
mkdir -p "$(dirname "$OUT_PREFIX")"
shopt -s nullglob
bams=("$BAM_DIR"/*.bam)
[[ ${#bams[@]} -gt 0 ]] || { echo "ERROR: no BAM files in $BAM_DIR" >&2; exit 1; }
matrix="${OUT_PREFIX}.npz"
raw="${OUT_PREFIX}.counts.tsv"
command=(multiBamSummary bins --bamfiles "${bams[@]}" --binSize "$BIN_SIZE" --skipZeros -p "$THREADS" --outFileName "$matrix" --outRawCounts "$raw")
printf '%q ' "${command[@]}" > "${OUT_PREFIX}.COMMANDS.sh"
printf '\n' >> "${OUT_PREFIX}.COMMANDS.sh"
"${command[@]}" > "${OUT_PREFIX}.multiBamSummary.log" 2>&1
plotCorrelation -in "$matrix" --corMethod pearson --whatToPlot heatmap --skipZeros \
  --colorMap RdYlBu_r --plotNumbers -o "${OUT_PREFIX}.pearson.heatmap.pdf" \
  --outFileCorMatrix "${OUT_PREFIX}.pearson.tsv" > "${OUT_PREFIX}.plotCorrelation.log" 2>&1
plotCorrelation -in "$matrix" --corMethod spearman --whatToPlot heatmap --skipZeros \
  --colorMap RdYlBu_r --plotNumbers -o "${OUT_PREFIX}.spearman.heatmap.pdf" \
  --outFileCorMatrix "${OUT_PREFIX}.spearman.tsv" >> "${OUT_PREFIX}.plotCorrelation.log" 2>&1
printf 'bam_count\t%s\nbin_size\t%s\n' "${#bams[@]}" "$BIN_SIZE" > "${OUT_PREFIX}.parameters.tsv"
