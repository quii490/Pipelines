
#!/usr/bin/env bash
set -euo pipefail
log(){ echo "[run_gene_body_profile] $*"; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
BW_GLOB=""
GENE_BED=""
SPECIES="hg38"
OUTDIR="./gene_body_profile"
CORES=4
UP=3000
DOWN=3000
BODY_LEN=5000
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-glob) BW_GLOB="$2"; shift 2 ;;
    --gene-bed) GENE_BED="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --upstream) UP="$2"; shift 2 ;;
    --downstream) DOWN="$2"; shift 2 ;;
    --body-length) BODY_LEN="$2"; shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -n "$BW_GLOB" ]] || { echo "--bw-glob is required" >&2; exit 1; }
if [[ -z "$GENE_BED" && -f "$CONFIG_FILE" ]]; then GENE_BED="$(get_species_param "$SPECIES" gene_body_bed "$CONFIG_FILE" || true)"; fi
[[ -n "$GENE_BED" ]] || { echo "--gene-bed is required, or set params.gene_body_bed for species $SPECIES in $CONFIG_FILE" >&2; exit 1; }
[[ -f "$GENE_BED" ]] || { echo "gene body bed not found: $GENE_BED" >&2; exit 1; }
mkdir -p "$OUTDIR"
SANITIZED_BED="$OUTDIR/$(basename "${GENE_BED%.*}").valid.bed"
awk 'BEGIN{OFS="	"} !/^#/ && NF>=3 { s=$2; e=$3; if (s < 0) s = 0; if (e <= s) next; print $1, s, e, (NF>=4?$4:"."), (NF>=5?$5:"0"), (NF>=6?$6:".") }' "$GENE_BED" | sort -k1,1 -k2,2n > "$SANITIZED_BED"
[[ -s "$SANITIZED_BED" ]] || { echo "No valid gene body intervals after sanitization: $GENE_BED" >&2; exit 1; }
shopt -s nullglob
bw_files=( $BW_GLOB )
[[ ${#bw_files[@]} -gt 0 ]] || { echo "No bigWig files matched: $BW_GLOB" >&2; exit 1; }
log "bigWig files: ${#bw_files[@]}"
log "gene body bed: $SANITIZED_BED"
computeMatrix scale-regions   -R "$SANITIZED_BED"   -S "${bw_files[@]}"   -b "$UP" -a "$DOWN"   --regionBodyLength "$BODY_LEN"   --skipZeros   -p "$CORES"   -o "$OUTDIR/GeneBody_accessibility_matrix.gz"
plotProfile -m "$OUTDIR/GeneBody_accessibility_matrix.gz" -out "$OUTDIR/GeneBody_accessibility_profile.pdf" --perGroup
plotProfile -m "$OUTDIR/GeneBody_accessibility_matrix.gz" -out "$OUTDIR/GeneBody_accessibility_profile.png" --perGroup
plotHeatmap -m "$OUTDIR/GeneBody_accessibility_matrix.gz" -out "$OUTDIR/GeneBody_accessibility_heatmap.pdf" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
plotHeatmap -m "$OUTDIR/GeneBody_accessibility_matrix.gz" -out "$OUTDIR/GeneBody_accessibility_heatmap.png" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
log "done"
