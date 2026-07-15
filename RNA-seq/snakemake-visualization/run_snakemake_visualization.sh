#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG=""
RESULTS_DIR=""
PLOT_DIR=""
CONFIG_OUT=""
SPECIES="hg38"
OUTDIR=""
RSCRIPT="conda run -n downstream Rscript"
CORES="8"
DRY_RUN="false"
EXTRA_ARGS=()
MAKE_CONFIG_ARGS=()

usage() {
  cat <<'USAGE'
Run RNA-seq Snakemake downstream visualization workflow

Usage:
  bash run_snakemake_visualization.sh --config my_project.yaml --cores 8
  bash run_snakemake_visualization.sh --results-dir /path/to/results --cores 8

Options:
  --config FILE             Existing Snakemake config YAML.
  --results-dir DIR         RNA-seq main pipeline results dir; auto-generate config.
  --plot-dir DIR            Plot/matrix dir inside results, e.g. plots or plot_with_replicates.
  --config-out FILE         Auto-generated config path.
  --outdir DIR              Snakemake output dir when auto-generating config.
  --species STR             hg38 | t2t | mm10 | mm39. Default: hg38.
  --rscript CMD             Rscript command in config. Default: conda run -n downstream Rscript.
  --include-matrix REGEX    Only include matrix names matching regex.
  --exclude-matrix REGEX    Exclude matrix names matching regex.
  --no-pathway              Disable gene pathway jobs in generated config.
  --no-go                   Disable GO in generated config.
  --no-gsea                 Disable GSEA in generated config.
  --no-te-analysis          Disable TE-specific jobs in generated config.
  --cores INT               CPU cores. Default: 8.
  --dry-run                 Run snakemake -n -p only.
  --extra "ARGS"            Extra raw args appended to snakemake. Can repeat.
  -h, --help                Show help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config) CONFIG="$2"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --plot-dir) PLOT_DIR="$2"; shift 2 ;;
    --config-out) CONFIG_OUT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --rscript) RSCRIPT="$2"; shift 2 ;;
    --include-matrix|--exclude-matrix|--sample-table|--contrast-file|--tx2gene-path|--te-annotation-tsv|--padj-cutoff|--lfc-cutoff|--base-mean-min|--label-top-n|--diff-threads)
      MAKE_CONFIG_ARGS+=("$1" "$2"); shift 2 ;;
    --no-pathway|--no-go|--no-gsea|--disable-gseaplot2|--no-te-analysis)
      MAKE_CONFIG_ARGS+=("$1"); shift 1 ;;
    --cores) CORES="$2"; shift 2 ;;
    --dry-run) DRY_RUN="true"; shift 1 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[run_snakemake_visualization] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${CONFIG}" && -z "${RESULTS_DIR}" ]]; then
  echo "[run_snakemake_visualization] ERROR: provide --config or --results-dir" >&2
  usage
  exit 1
fi
if [[ -n "${CONFIG}" && -n "${RESULTS_DIR}" ]]; then
  echo "[run_snakemake_visualization] ERROR: use only one of --config or --results-dir" >&2
  exit 1
fi
if [[ -z "${CONFIG}" ]]; then
  if [[ -z "${CONFIG_OUT}" ]]; then
    CONFIG_OUT="${RESULTS_DIR%/}/snakemake_visualization.config.yaml"
  fi
  make_cmd=(python3 "${SCRIPT_DIR}/make_config_from_results.py" --results-dir "${RESULTS_DIR}" --out "${CONFIG_OUT}" --species "${SPECIES}" --rscript "${RSCRIPT}")
  [[ -z "${PLOT_DIR}" ]] || make_cmd+=(--plot-dir "${PLOT_DIR}")
  [[ -z "${OUTDIR}" ]] || make_cmd+=(--outdir "${OUTDIR}")
  make_cmd+=("${MAKE_CONFIG_ARGS[@]}")
  echo "[run_snakemake_visualization] ${make_cmd[*]}"
  "${make_cmd[@]}"
  CONFIG="${CONFIG_OUT}"
fi
if [[ ! -f "${CONFIG}" ]]; then
  echo "[run_snakemake_visualization] ERROR: config not found: ${CONFIG}" >&2
  exit 1
fi
if ! command -v snakemake >/dev/null 2>&1; then
  echo "[run_snakemake_visualization] ERROR: snakemake not found in PATH" >&2
  echo "Activate the downstream/rnaseq conda environment first." >&2
  exit 1
fi

cmd=(snakemake -s "${SCRIPT_DIR}/visualization.smk" --configfile "${CONFIG}" --cores "${CORES}" -p --rerun-incomplete)
if [[ "${DRY_RUN}" == "true" ]]; then
  cmd+=(--dry-run)
fi
for x in "${EXTRA_ARGS[@]}"; do
  read -r -a parts <<< "${x}"
  cmd+=("${parts[@]}")
done

echo "[run_snakemake_visualization] ${cmd[*]}"
"${cmd[@]}"
