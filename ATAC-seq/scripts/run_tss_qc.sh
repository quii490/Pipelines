#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[run_tss_qc] $*"; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
BW_GLOB=""
TSS_BED=""
SPECIES="hg38"
OUTDIR="./tss_qc"
CORES=4
UP=3000
DOWN=3000
usage() {
  cat <<'USAGE'
Usage:
  run_tss_qc.sh --bw-glob "04_bw/*.bw" [options]

Purpose:
  Draw a deepTools aggregate TSS signal profile and heatmap from existing
  bigWig tracks. This is not the insertion-based ENCODE TSS enrichment score.

Options:
  --bw-glob STR       Quoted bigWig glob (required).
  --tss-bed FILE      TSS BED; otherwise resolved from species config.
  --species STR       hg38 | mm10 | mm39. Default: hg38.
  --outdir DIR        Default: ./tss_qc.
  --cores INT         Default: 4.
  --upstream INT      Bases before TSS. Default: 3000.
  --downstream INT    Bases after TSS. Default: 3000.
  -h, --help          Show this help.
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-glob) BW_GLOB="$2"; shift 2 ;;
    --tss-bed) TSS_BED="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --upstream) UP="$2"; shift 2 ;;
    --downstream) DOWN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$BW_GLOB" ]] || { echo "--bw-glob is required" >&2; exit 1; }
if [[ -z "$TSS_BED" && -f "$CONFIG_FILE" ]]; then TSS_BED="$(get_species_param "$SPECIES" tss_bed "$CONFIG_FILE" || true)"; fi
[[ -n "$TSS_BED" ]] || { echo "--tss-bed is required, or set params.tss_bed for species $SPECIES in $CONFIG_FILE" >&2; exit 1; }
[[ -f "$TSS_BED" ]] || { echo "tss bed not found: $TSS_BED" >&2; exit 1; }
mkdir -p "$OUTDIR"
SANITIZED_BED="$OUTDIR/$(basename "${TSS_BED%.*}").valid.bed"
awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
  s=$2; e=$3;
  if (s < 0) s = 0;
  if (e <= s) e = s + 1;
  print $1, s, e, (NF>=4?$4:"."), (NF>=5?$5:"0"), (NF>=6?$6:".")
}' "$TSS_BED" | sort -k1,1 -k2,2n > "$SANITIZED_BED"
[[ -s "$SANITIZED_BED" ]] || { echo "No valid TSS intervals after sanitization: $TSS_BED" >&2; exit 1; }
shopt -s nullglob
bw_files=( $BW_GLOB )
[[ ${#bw_files[@]} -gt 0 ]] || { echo "No bigWig files matched: $BW_GLOB" >&2; exit 1; }
log "bigWig files: ${#bw_files[@]}"
log "sanitized tss bed: $SANITIZED_BED"
computeMatrix reference-point \
  --referencePoint TSS \
  -R "$SANITIZED_BED" \
  -S "${bw_files[@]}" \
  -b "$UP" -a "$DOWN" \
  --skipZeros \
  -p "$CORES" \
  -o "$OUTDIR/TSS_enrichment_matrix.gz"
plotProfile -m "$OUTDIR/TSS_enrichment_matrix.gz" -out "$OUTDIR/TSS_enrichment_profile.pdf" --perGroup
plotProfile -m "$OUTDIR/TSS_enrichment_matrix.gz" -out "$OUTDIR/TSS_enrichment_profile.png" --perGroup
plotHeatmap -m "$OUTDIR/TSS_enrichment_matrix.gz" -out "$OUTDIR/TSS_enrichment_heatmap.pdf" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
plotHeatmap -m "$OUTDIR/TSS_enrichment_matrix.gz" -out "$OUTDIR/TSS_enrichment_heatmap.png" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
log "done"
