#!/usr/bin/env bash
set -euo pipefail
log(){ echo "[run_fixedbin_from_bam] $*"; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
source "$SCRIPT_DIR/species_config_lib.sh"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
BAM_GLOB=""
CHROM_SIZES=""
SAMPLE_META=""
CONTRAST_FILE=""
SPECIES="hg38"
OUTDIR="./fixedbin_from_bam"
BIN_SIZE=100000
CORES=4
TE_BED=""
GTF_FILE=""
TE_VIOLIN_CLASS=""
TE_VIOLIN_TOP_N="12"

usage() {
  cat <<'USAGE'
Usage:
  run_fixedbin_from_bam.sh --bam-glob PATTERN --sample-meta FILE --contrast-file FILE [options]

Required:
  --bam-glob PATTERN
  --sample-meta FILE       columns: sample,condition,replicate
  --contrast-file FILE     columns: case,control

Reference options:
  --species STR            hg38 | mm10 | mm39
  --chrom-sizes FILE       optional override; otherwise read from conf/species_refs.config
  --te-bed FILE            optional override
  --gtf FILE               optional override

Analysis options:
  --outdir DIR
  --bin-size INT           default 100000
  --cores INT              default 4
  --te-violin-class STR
  --te-violin-top-n INT
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-glob) BAM_GLOB="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --sample-meta) SAMPLE_META="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --gtf) GTF_FILE="$2"; shift 2 ;;
    --te-violin-class) TE_VIOLIN_CLASS="$2"; shift 2 ;;
    --te-violin-top-n) TE_VIOLIN_TOP_N="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done
[[ -n "$BAM_GLOB" ]] || { echo "--bam-glob is required" >&2; exit 1; }
[[ -n "$SAMPLE_META" ]] || { echo "--sample-meta is required" >&2; exit 1; }
[[ -n "$CONTRAST_FILE" ]] || { echo "--contrast-file is required" >&2; exit 1; }
if [[ -z "$CHROM_SIZES" && -f "$CONFIG_FILE" ]]; then CHROM_SIZES="$(get_species_param "$SPECIES" chrom_sizes "$CONFIG_FILE" || true)"; fi
if [[ -z "$TE_BED" && -f "$CONFIG_FILE" ]]; then TE_BED="$(get_species_param "$SPECIES" te_bed "$CONFIG_FILE" || true)"; fi
if [[ -z "$GTF_FILE" && -f "$CONFIG_FILE" ]]; then GTF_FILE="$(get_species_param "$SPECIES" gtf_genes "$CONFIG_FILE" || true)"; fi
[[ -n "$CHROM_SIZES" ]] || { echo "--chrom-sizes is required, or set params.chrom_sizes for species $SPECIES in $CONFIG_FILE" >&2; exit 1; }
[[ -f "$CHROM_SIZES" ]] || { echo "chrom sizes not found: $CHROM_SIZES" >&2; exit 1; }
mkdir -p "$OUTDIR" "$OUTDIR/counts"
shopt -s nullglob
bam_files=( $BAM_GLOB )
[[ ${#bam_files[@]} -gt 0 ]] || { echo "No BAM files matched: $BAM_GLOB" >&2; exit 1; }
log "bam files: ${#bam_files[@]}"

BINS_BED="$OUTDIR/counts/fixed_bins_${BIN_SIZE}.bed"
BINS_SAF="$OUTDIR/counts/fixed_bins_${BIN_SIZE}.saf"
COUNT_FILE="$OUTDIR/counts/fixedbin_counts.txt"

bedtools makewindows -g "$CHROM_SIZES" -w "$BIN_SIZE" \
  | awk 'BEGIN{OFS="\t"} $1 ~ /^chr([0-9]+|X|Y|M)$/ {print $0}' > "$BINS_BED"
[[ -s "$BINS_BED" ]] || { echo "No bins generated from chrom sizes: $CHROM_SIZES" >&2; exit 1; }
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print "bin_" NR, $1, $2+1, $3, "+"}' "$BINS_BED" > "$BINS_SAF"
log "bins generated: $(wc -l < "$BINS_BED")"

paired_flag=0
if samtools view -c -f 1 "${bam_files[0]}" | awk '{exit !($1>0)}'; then
  paired_flag=1
fi
log "paired-end inferred from first bam: $paired_flag"
if [[ "$paired_flag" -eq 1 ]]; then
  featureCounts -a "$BINS_SAF" -F SAF -o "$COUNT_FILE" -T "$CORES" -p --countReadPairs -B -C "${bam_files[@]}"
else
  featureCounts -a "$BINS_SAF" -F SAF -o "$COUNT_FILE" -T "$CORES" "${bam_files[@]}"
fi

Rscript "$ROOT_DIR/atacseq-downstream/run_downstream_atac.R" \
  --count-file "$COUNT_FILE" \
  --sample-meta "$SAMPLE_META" \
  --peak-bed "$BINS_BED" \
  --outdir "$OUTDIR" \
  --species "$SPECIES" \
  --contrast-file "$CONTRAST_FILE" \
  ${TE_BED:+--te-bed "$TE_BED"} \
  ${GTF_FILE:+--gtf "$GTF_FILE"} \
  ${TE_VIOLIN_CLASS:+--te-violin-class "$TE_VIOLIN_CLASS"} \
  ${TE_VIOLIN_TOP_N:+--te-violin-top-n "$TE_VIOLIN_TOP_N"}
log "done"
