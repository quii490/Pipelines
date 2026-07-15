#!/usr/bin/env bash
set -euo pipefail

BW_DIR=""
BWS=()
OUT_PREFIX=""
LABELS=""
BINSIZE="10000"
THREADS="8"
METHODS="pearson,spearman"
SKIP_ZEROS="true"
EXTRA_SUMMARY_ARGS=""
EXTRA_PLOT_ARGS=""

usage() {
  cat <<'USAGE'
RNA-seq bigWig correlation wrapper using deepTools

示例：
  rnaseq-bw-cor \
    --bw-dir /path/to/bw \
    --out-prefix /path/qc/rnaseq_bw_cor \
    --binsize 10000 \
    --threads 8

输入：
  --bw-dir DIR             bigWig 目录，自动读取 *.bw/*.bigWig
  --bw FILE                指定 bigWig，可重复
  --labels "A B C"         可选，样本标签；数量需与 bigWig 一致

输出：
  <out-prefix>.multiBigwigSummary.npz
  <out-prefix>.multiBigwigSummary.tab
  <out-prefix>.<method>.heatmap.pdf/png
  <out-prefix>.<method>.scatter.pdf/png
  <out-prefix>.PCA.pdf/png

参数：
  --out-prefix PATH        输出前缀
  --binsize INT            bin size，默认 10000
  --threads INT            默认 8
  --methods LIST           pearson,spearman,kendall，默认 pearson,spearman
  --skip-zeros true|false  默认 true
  --extra-summary-args STR 额外传给 multiBigwigSummary bins
  --extra-plot-args STR    额外传给 plotCorrelation
  -h, --help               显示帮助
USAGE
}

die() {
  echo "[run_bw_correlation] ERROR: $*" >&2
  exit 1
}

log() {
  echo "[run_bw_correlation] $*" >&2
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-dir) BW_DIR="$2"; shift 2 ;;
    --bw|--bigwig) BWS+=("$2"); shift 2 ;;
    --labels) LABELS="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --binsize|--bin-size) BINSIZE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --methods) METHODS="$2"; shift 2 ;;
    --skip-zeros) SKIP_ZEROS="$2"; shift 2 ;;
    --extra-summary-args) EXTRA_SUMMARY_ARGS="$2"; shift 2 ;;
    --extra-plot-args) EXTRA_PLOT_ARGS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

command -v multiBigwigSummary >/dev/null 2>&1 || die "multiBigwigSummary not found in PATH"
command -v plotCorrelation >/dev/null 2>&1 || die "plotCorrelation not found in PATH"
command -v plotPCA >/dev/null 2>&1 || die "plotPCA not found in PATH"
[[ -n "${OUT_PREFIX}" ]] || die "--out-prefix is required"

if [[ -n "${BW_DIR}" ]]; then
  [[ -d "${BW_DIR}" ]] || die "--bw-dir does not exist: ${BW_DIR}"
  while IFS= read -r bw; do
    BWS+=("${bw}")
  done < <(find "${BW_DIR}" -maxdepth 2 -type f \( -name "*.bw" -o -name "*.bigWig" \) | sort)
fi
[[ "${#BWS[@]}" -ge 2 ]] || die "At least two bigWig files are required"

OUT_DIR="$(dirname "${OUT_PREFIX}")"
mkdir -p "${OUT_DIR}"
NPZ="${OUT_PREFIX}.multiBigwigSummary.npz"
TAB="${OUT_PREFIX}.multiBigwigSummary.tab"

summary_cmd=(multiBigwigSummary bins -b "${BWS[@]}" -o "${NPZ}" --outRawCounts "${TAB}" --binSize "${BINSIZE}" -p "${THREADS}")
if [[ -n "${LABELS}" ]]; then
  read -r -a labels_arr <<< "${LABELS}"
  [[ "${#labels_arr[@]}" -eq "${#BWS[@]}" ]] || die "--labels count must match bigWig count"
  summary_cmd+=(--labels "${labels_arr[@]}")
fi
if [[ -n "${EXTRA_SUMMARY_ARGS}" ]]; then
  read -r -a extra_summary <<< "${EXTRA_SUMMARY_ARGS}"
  summary_cmd+=("${extra_summary[@]}")
fi

log "running multiBigwigSummary on ${#BWS[@]} bigWig files"
"${summary_cmd[@]}"

IFS=',' read -r -a method_arr <<< "${METHODS}"
for method in "${method_arr[@]}"; do
  method="$(echo "${method}" | tr -d '[:space:]')"
  [[ -n "${method}" ]] || continue
  for fmt in pdf png; do
    heat="${OUT_PREFIX}.${method}.heatmap.${fmt}"
    scatter="${OUT_PREFIX}.${method}.scatter.${fmt}"
    cor_cmd=(plotCorrelation -in "${NPZ}" --corMethod "${method}" --whatToPlot heatmap --plotFile "${heat}" --outFileCorMatrix "${OUT_PREFIX}.${method}.correlation_matrix.tsv")
    scat_cmd=(plotCorrelation -in "${NPZ}" --corMethod "${method}" --whatToPlot scatterplot --plotFile "${scatter}")
    if [[ "${SKIP_ZEROS}" == "true" ]]; then
      cor_cmd+=(--skipZeros)
      scat_cmd+=(--skipZeros)
    fi
    if [[ -n "${EXTRA_PLOT_ARGS}" ]]; then
      read -r -a extra_plot <<< "${EXTRA_PLOT_ARGS}"
      cor_cmd+=("${extra_plot[@]}")
      scat_cmd+=("${extra_plot[@]}")
    fi
    "${cor_cmd[@]}"
    "${scat_cmd[@]}"
  done
done

plotPCA -in "${NPZ}" --plotFile "${OUT_PREFIX}.PCA.pdf"
plotPCA -in "${NPZ}" --plotFile "${OUT_PREFIX}.PCA.png"

log "DONE: ${OUT_PREFIX}"
