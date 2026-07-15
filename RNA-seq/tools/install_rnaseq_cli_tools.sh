#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PREFIX=""
MODE="symlink"
RSCRIPT_CMD=""

usage() {
  cat <<'USAGE'
安装 RNA-seq 独立 CLI 工具到 conda 环境 bin 目录

用法：
  bash install_rnaseq_cli_tools.sh --prefix /path/to/conda/env

可选：
  --prefix DIR      conda 环境前缀；默认读取当前 CONDA_PREFIX
  --mode STR        symlink | copy，默认 symlink
  --rscript-cmd CMD R 工具安装为 wrapper，并用 CMD 运行，例如：
                    'conda run -n downstream Rscript'
  -h, --help        显示帮助

安装后可直接调用：
  rnaseq-bam-to-counts
  rnaseq-diff-counts
  rnaseq-annotate-de
  rnaseq-pathway-de
  rnaseq-te-analysis
  rnaseq-counts-to-de
  rnaseq-de-visuals
  rnaseq-plot-counts
  rnaseq-gsea
  rnaseq-bw-cor
  rnaseq-two-sample-scatter
  rnaseq-qc-metrics-bam
  rnaseq-report
  rnaseq-clean-work
  rnaseq-publish
  rnaseq-rerun
  install_snakemake_922
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --prefix) PREFIX="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --rscript-cmd) RSCRIPT_CMD="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${PREFIX}" ]]; then
  PREFIX="${CONDA_PREFIX:-}"
fi
if [[ -z "${PREFIX}" ]]; then
  echo "[install_rnaseq_cli_tools] ERROR: 未提供 --prefix，且当前没有 CONDA_PREFIX" >&2
  exit 1
fi

BIN_DIR="${PREFIX}/bin"
mkdir -p "${BIN_DIR}"

echo "[install_rnaseq_cli_tools] prefix=${PREFIX}"
echo "[install_rnaseq_cli_tools] bin=${BIN_DIR}"
echo "[install_rnaseq_cli_tools] mode=${MODE}"
if [[ -n "${RSCRIPT_CMD}" ]]; then
  echo "[install_rnaseq_cli_tools] rscript_cmd=${RSCRIPT_CMD}"
fi

install_one() {
  local src="$1"
  local dest_name="$2"
  local dest="${BIN_DIR}/${dest_name}"
  if [[ ! -f "${src}" ]]; then
    echo "[install_rnaseq_cli_tools] ERROR: source not found: ${src}" >&2
    exit 1
  fi
  rm -f "${dest}"
  if [[ "${MODE}" == "copy" ]]; then
    cp "${src}" "${dest}"
  else
    ln -s "${src}" "${dest}"
  fi
  chmod +x "${dest}"
  echo "[install_rnaseq_cli_tools] installed ${dest_name} -> ${src}"
}

install_r_tool() {
  local src="$1"
  local dest_name="$2"
  local dest="${BIN_DIR}/${dest_name}"
  if [[ -z "${RSCRIPT_CMD}" ]]; then
    install_one "${src}" "${dest_name}"
    return
  fi
  if [[ ! -f "${src}" ]]; then
    echo "[install_rnaseq_cli_tools] ERROR: source not found: ${src}" >&2
    exit 1
  fi
  rm -f "${dest}"
  cat > "${dest}" <<EOF
#!/usr/bin/env bash
exec ${RSCRIPT_CMD} "${src}" "\$@"
EOF
  chmod +x "${dest}"
  echo "[install_rnaseq_cli_tools] installed ${dest_name} wrapper -> ${RSCRIPT_CMD} ${src}"
}

install_r_tool "${SCRIPT_DIR}/run_plot_from_counts.R" "rnaseq-plot-counts"
install_r_tool "${SCRIPT_DIR}/run_gsea_standalone.R" "rnaseq-gsea"
install_one "${SCRIPT_DIR}/run_bam_to_counts.sh" "rnaseq-bam-to-counts"
install_r_tool "${SCRIPT_DIR}/run_diff_from_counts.R" "rnaseq-diff-counts"
install_r_tool "${SCRIPT_DIR}/run_annotate_de.R" "rnaseq-annotate-de"
install_r_tool "${SCRIPT_DIR}/run_pathway_from_de.R" "rnaseq-pathway-de"
install_r_tool "${SCRIPT_DIR}/run_te_analysis.R" "rnaseq-te-analysis"
install_r_tool "${SCRIPT_DIR}/run_counts_to_de_matrix.R" "rnaseq-counts-to-de"
install_r_tool "${SCRIPT_DIR}/run_de_matrix_visuals.R" "rnaseq-de-visuals"
install_one "${SCRIPT_DIR}/run_bw_correlation.sh" "rnaseq-bw-cor"
install_one "${SCRIPT_DIR}/run_two_sample_scatter.py" "rnaseq-two-sample-scatter"
install_one "${SCRIPT_DIR}/run_qc_metrics_bam.sh" "rnaseq-qc-metrics-bam"
install_one "${SCRIPT_DIR}/run_report.sh" "rnaseq-report"
install_one "${SCRIPT_DIR}/rnaseq_clean_work.py" "rnaseq-clean-work"
install_one "${SCRIPT_DIR}/rnaseq_publish.py" "rnaseq-publish"
install_one "${SCRIPT_DIR}/rnaseq_rerun.sh" "rnaseq-rerun"
install_one "${SCRIPT_DIR}/install_snakemake_922.sh" "install_snakemake_922"
install_one "${SCRIPT_DIR}/../rnaseq-downstream/rnaseq-function.R" "rnaseq-function.R"

echo "[install_rnaseq_cli_tools] DONE"
