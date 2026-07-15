#!/usr/bin/env bash
set -euo pipefail
log(){ echo "[run_te_heatmap_batch] $*"; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
BW_GLOB=""
TE_BED=""
SPECIES="hg38"
GLOBAL_PEAK_BED=""
CONTRAST_BED_GLOB=""
OUTDIR="./te_heatmap_batch"
CORES=4
MAX_REGIONS=5000
MODE="center"
TE_LABEL_LEVEL="locus"
TE_CLASS_FILTER=""
TE_FAMILY_FILTER=""
TE_NAME_FILTER=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-glob) BW_GLOB="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --global-peak-bed) GLOBAL_PEAK_BED="$2"; shift 2 ;;
    --contrast-bed-glob) CONTRAST_BED_GLOB="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --max-regions) MAX_REGIONS="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --te-label-level) TE_LABEL_LEVEL="$2"; shift 2 ;;
    --te-class-filter) TE_CLASS_FILTER="$2"; shift 2 ;;
    --te-family-filter) TE_FAMILY_FILTER="$2"; shift 2 ;;
    --te-name-filter) TE_NAME_FILTER="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$BW_GLOB" ]] || { echo "--bw-glob is required" >&2; exit 1; }
if [[ -z "$TE_BED" && -f "$CONFIG_FILE" ]]; then TE_BED="$(get_species_param "$SPECIES" te_bed "$CONFIG_FILE" || true)"; fi
[[ -n "$TE_BED" ]] || { echo "--te-bed is required, or set params.te_bed for species $SPECIES in $CONFIG_FILE" >&2; exit 1; }
mkdir -p "$OUTDIR"
common_args=(--bw-glob "$BW_GLOB" --te-bed "$TE_BED" --species "$SPECIES" --cores "$CORES" --max-regions "$MAX_REGIONS" --mode "$MODE" --te-label-level "$TE_LABEL_LEVEL")
[[ -n "$TE_CLASS_FILTER" ]] && common_args+=(--te-class-filter "$TE_CLASS_FILTER")
[[ -n "$TE_FAMILY_FILTER" ]] && common_args+=(--te-family-filter "$TE_FAMILY_FILTER")
[[ -n "$TE_NAME_FILTER" ]] && common_args+=(--te-name-filter "$TE_NAME_FILTER")

if [[ -n "$GLOBAL_PEAK_BED" ]]; then
  log "run global accessible TE heatmap"
  bash "$SCRIPT_DIR/run_te_heatmap.sh" "${common_args[@]}" --peak-bed "$GLOBAL_PEAK_BED" --outdir "$OUTDIR/global" --purpose global
fi

if [[ -n "$CONTRAST_BED_GLOB" ]]; then
  shopt -s nullglob
  peak_files=( $CONTRAST_BED_GLOB )
  [[ ${#peak_files[@]} -gt 0 ]] || { echo "No contrast peak bed matched: $CONTRAST_BED_GLOB" >&2; exit 1; }
  for bed in "${peak_files[@]}"; do
    base=$(basename "$bed")
    parent=$(basename "$(dirname "$bed")")
    label="${parent}_${base%.bed}"
    subdir="$OUTDIR/contrast/$label"
    log "run contrast-specific TE heatmap: $label"
    bash "$SCRIPT_DIR/run_te_heatmap.sh" "${common_args[@]}" --peak-bed "$bed" --outdir "$subdir" --purpose contrast --label "$label"
  done
fi

log "done"
