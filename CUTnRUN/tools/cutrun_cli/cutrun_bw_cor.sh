#!/usr/bin/env bash
set -euo pipefail

BW_DIR=""
OUT_PREFIX=""
THREADS=8
BIN_SIZE="${BIGWIG_COR_BIN_SIZE:-5000}"

usage() {
  cat <<'USAGE'
Usage: cutrun_bw_cor.sh --bw-dir DIR --out-prefix PREFIX [--threads N] [--bin-size N]
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-dir) BW_DIR="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done
[[ -d "$BW_DIR" ]] || { echo "ERROR: bigWig directory missing: $BW_DIR" >&2; exit 1; }
[[ -n "$OUT_PREFIX" ]] || { echo "ERROR: --out-prefix is required" >&2; exit 1; }
command -v multiBigwigSummary >/dev/null 2>&1 || { echo "ERROR: deepTools multiBigwigSummary not found" >&2; exit 1; }
command -v plotCorrelation >/dev/null 2>&1 || { echo "ERROR: deepTools plotCorrelation not found" >&2; exit 1; }
mkdir -p "$(dirname "$OUT_PREFIX")"
shopt -s nullglob
tracks=("$BW_DIR"/*.bw)
[[ ${#tracks[@]} -gt 0 ]] || { echo "ERROR: no .bw files in $BW_DIR" >&2; exit 1; }
matrix="${OUT_PREFIX}.npz"
raw="${OUT_PREFIX}.counts.tsv"
command=(multiBigwigSummary bins -b "${tracks[@]}" --binSize "$BIN_SIZE" --skipZeros -p "$THREADS" -out "$matrix" --outRawCounts "$raw")
printf '%q ' "${command[@]}" > "${OUT_PREFIX}.COMMANDS.sh"
printf '\n' >> "${OUT_PREFIX}.COMMANDS.sh"
"${command[@]}" > "${OUT_PREFIX}.multiBigwigSummary.log" 2>&1
plotCorrelation -in "$matrix" --corMethod pearson --whatToPlot heatmap --skipZeros \
  --colorMap RdYlBu_r --plotNumbers -o "${OUT_PREFIX}.pearson.heatmap.pdf" \
  --outFileCorMatrix "${OUT_PREFIX}.pearson.tsv" > "${OUT_PREFIX}.plotCorrelation.log" 2>&1
plotCorrelation -in "$matrix" --corMethod spearman --whatToPlot heatmap --skipZeros \
  --colorMap RdYlBu_r --plotNumbers -o "${OUT_PREFIX}.spearman.heatmap.pdf" \
  --outFileCorMatrix "${OUT_PREFIX}.spearman.tsv" >> "${OUT_PREFIX}.plotCorrelation.log" 2>&1
printf 'track_count\t%s\nbin_size\t%s\n' "${#tracks[@]}" "$BIN_SIZE" > "${OUT_PREFIX}.parameters.tsv"
