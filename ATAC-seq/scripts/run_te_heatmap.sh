#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[run_te_heatmap] $*"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"

BW_GLOB=""
TE_BED=""
PEAK_BED=""
OUTDIR="./te_heatmap"
SPECIES="hg38"
PURPOSE="global"
CORES=4
UP=2000
DOWN=2000
BODY_LEN=3000
MAX_REGIONS=0
MODE="scale-regions"
LABEL=""
TE_LABEL_LEVEL="locus"
TE_CLASS_FILTER=""
TE_FAMILY_FILTER=""
TE_NAME_FILTER=""

usage() {
  cat <<'USAGE'
Usage:
  run_te_heatmap.sh --bw-glob "04_bw/*.bw" [options]

Required:
  --bw-glob PATTERN       bigWig 文件 glob

Reference options:
  --species STR           hg38 | mm10 | mm39; used to resolve TE annotation from conf/species_refs.config
  --te-bed FILE           Optional override. BED/GTF/GTF.GZ TE annotation
  --peak-bed FILE         Optional. If provided, also writes TE intervals overlapping peaks

Analysis options:
  --purpose STR           global | peak-overlap. global uses all TE intervals; peak-overlap uses TE intervals overlapping --peak-bed
  --outdir DIR
  --cores INT
  --upstream INT          default 2000
  --downstream INT        default 2000
  --body-length INT       default 3000
  --max-regions INT       default 0 (no sampling). If >0, sample at most N regions for computeMatrix
  --mode STR              accepted for compatibility; scale-regions behavior is used
  --label STR             accepted for compatibility with batch wrapper
  --te-label-level STR    accepted for compatibility with batch wrapper
  --te-class-filter STR   optional grep filter on TE annotation/name field
  --te-family-filter STR  optional grep filter on TE annotation/name field
  --te-name-filter STR    optional grep filter on TE annotation/name field
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-glob) BW_GLOB="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --peak-bed) PEAK_BED="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --purpose) PURPOSE="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --upstream) UP="$2"; shift 2 ;;
    --downstream) DOWN="$2"; shift 2 ;;
    --body-length) BODY_LEN="$2"; shift 2 ;;
    --max-regions) MAX_REGIONS="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --label) LABEL="$2"; shift 2 ;;
    --te-label-level) TE_LABEL_LEVEL="$2"; shift 2 ;;
    --te-class-filter) TE_CLASS_FILTER="$2"; shift 2 ;;
    --te-family-filter) TE_FAMILY_FILTER="$2"; shift 2 ;;
    --te-name-filter) TE_NAME_FILTER="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$BW_GLOB" ]] || { echo "--bw-glob is required" >&2; exit 1; }
if [[ -z "$TE_BED" && -f "$CONFIG_FILE" ]]; then
  TE_BED="$(get_species_param "$SPECIES" te_bed "$CONFIG_FILE" || true)"
fi
[[ -n "$TE_BED" ]] || { echo "--te-bed is required, or set params.te_bed for species $SPECIES in $CONFIG_FILE" >&2; exit 1; }
[[ -f "$TE_BED" ]] || { echo "TE annotation not found: $TE_BED" >&2; exit 1; }
if [[ -n "$PEAK_BED" && ! -f "$PEAK_BED" ]]; then
  echo "peak bed not found: $PEAK_BED" >&2
  exit 1
fi
[[ "$PURPOSE" == "global" || "$PURPOSE" == "peak-overlap" || "$PURPOSE" == "contrast" ]] || {
  echo "--purpose must be global, peak-overlap, or contrast" >&2
  exit 1
}

mkdir -p "$OUTDIR"
TE_LOCUS_BED="$OUTDIR/te.locus.bed"
TE_FILTERED_BED="$OUTDIR/te.filtered.bed"
TE_SAMPLED_BED="$OUTDIR/te.sampled.bed"

read_cmd=(cat "$TE_BED")
case "$TE_BED" in
  *.gz) read_cmd=(gzip -dc "$TE_BED") ;;
esac

if [[ "$TE_BED" == *.gtf || "$TE_BED" == *.gtf.gz || "$TE_BED" == *.gff || "$TE_BED" == *.gff3 || "$TE_BED" == *.gff.gz || "$TE_BED" == *.gff3.gz ]]; then
  log "Detected GTF/GFF format, converting to BED"
  "${read_cmd[@]}" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=5 {
    s=$4-1; e=$5; if (s<0) s=0; if (e<=s) next;
    name=$9; gsub(/[[:space:]]+/, "_", name);
    print $1,s,e,name,"0",(NF>=7?$7:".")
  }' | sort -k1,1 -k2,2n > "$TE_LOCUS_BED"
else
  log "Detected BED-like format"
  "${read_cmd[@]}" | awk 'BEGIN{FS=OFS="\t"} !/^#/ && NF>=3 {
    s=$2; e=$3; if (s<0) s=0; if (e<=s) next;
    print $1,s,e,(NF>=4?$4:"."),(NF>=5?$5:"0"),(NF>=6?$6:".")
  }' | sort -k1,1 -k2,2n > "$TE_LOCUS_BED"
fi
[[ -s "$TE_LOCUS_BED" ]] || { echo "No valid TE intervals after sanitization: $TE_BED" >&2; exit 1; }

for filter in "$TE_CLASS_FILTER" "$TE_FAMILY_FILTER" "$TE_NAME_FILTER"; do
  if [[ -n "$filter" ]]; then
    awk -v pat="$filter" 'BEGIN{IGNORECASE=1} $4 ~ pat' "$TE_LOCUS_BED" > "$TE_LOCUS_BED.tmp" || true
    if [[ -s "$TE_LOCUS_BED.tmp" ]]; then
      mv "$TE_LOCUS_BED.tmp" "$TE_LOCUS_BED"
    else
      rm -f "$TE_LOCUS_BED.tmp"
      echo "TE filter produced no intervals: $filter" >&2
      exit 1
    fi
  fi
done

cp "$TE_LOCUS_BED" "$TE_FILTERED_BED"
if [[ -n "$PEAK_BED" ]]; then
  bedtools intersect -u -a "$TE_LOCUS_BED" -b "$PEAK_BED" > "$TE_FILTERED_BED.tmp" || true
  if [[ -s "$TE_FILTERED_BED.tmp" ]]; then
    mv "$TE_FILTERED_BED.tmp" "$TE_FILTERED_BED"
  else
    rm -f "$TE_FILTERED_BED.tmp"
    log "No TE intervals overlap --peak-bed; te.filtered.bed falls back to all TE intervals"
    cp "$TE_LOCUS_BED" "$TE_FILTERED_BED"
  fi
fi

REGIONS="$TE_LOCUS_BED"
if [[ "$PURPOSE" == "peak-overlap" || "$PURPOSE" == "contrast" ]]; then
  REGIONS="$TE_FILTERED_BED"
fi

cp "$REGIONS" "$TE_SAMPLED_BED"
if [[ "$MAX_REGIONS" =~ ^[0-9]+$ && "$MAX_REGIONS" -gt 0 ]]; then
  n_regions="$(wc -l < "$REGIONS")"
  if [[ "$n_regions" -gt "$MAX_REGIONS" ]]; then
    awk -v seed=1 'BEGIN{srand(seed)} {print rand()"\t"$0}' "$REGIONS" | sort -k1,1n | cut -f2- | head -n "$MAX_REGIONS" | sort -k1,1 -k2,2n > "$TE_SAMPLED_BED"
    REGIONS="$TE_SAMPLED_BED"
    log "Sampled TE regions for computeMatrix: $MAX_REGIONS / $n_regions"
  fi
fi

shopt -s nullglob
bw_files=( $BW_GLOB )
[[ ${#bw_files[@]} -gt 0 ]] || { echo "No bigWig files matched: $BW_GLOB" >&2; exit 1; }

log "bigWig files: ${#bw_files[@]}"
log "TE annotation: $TE_BED"
log "TE regions for matrix: $REGIONS"
computeMatrix scale-regions \
  -R "$REGIONS" \
  -S "${bw_files[@]}" \
  -b "$UP" -a "$DOWN" \
  --regionBodyLength "$BODY_LEN" \
  --skipZeros \
  -p "$CORES" \
  -o "$OUTDIR/TE_accessibility_matrix.gz"

plotProfile -m "$OUTDIR/TE_accessibility_matrix.gz" -out "$OUTDIR/TE_accessibility_profile.pdf" --perGroup
plotProfile -m "$OUTDIR/TE_accessibility_matrix.gz" -out "$OUTDIR/TE_accessibility_profile.png" --perGroup
plotHeatmap -m "$OUTDIR/TE_accessibility_matrix.gz" -out "$OUTDIR/TE_accessibility_heatmap.pdf" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
plotHeatmap -m "$OUTDIR/TE_accessibility_matrix.gz" -out "$OUTDIR/TE_accessibility_heatmap.png" --sortRegions descend --whatToShow 'plot, heatmap and colorbar'
log "done"
