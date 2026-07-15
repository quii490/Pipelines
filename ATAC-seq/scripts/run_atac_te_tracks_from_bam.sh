#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[run_atac_te_tracks_from_bam] $*"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"

BAM_GLOB=""
OUTDIR="./atac_te_tracks"
SPECIES="hg38"
BLACKLIST=""
MITO_CHR=""
EFFECTIVE_GENOME_SIZE=""
CORES=8
MAPQ=0
EXCLUDE_FLAGS=780
RUN_MARKDUP="false"
REMOVE_MITO="true"
REMOVE_BLACKLIST="true"
PROPER_PAIR_ONLY="true"
NORMALIZATION="CPM"
BINSIZE=10
IGNORE_FOR_NORMALIZATION=""

usage(){
  cat <<'USAGE'
Usage:
  run_atac_te_tracks_from_bam.sh --bam-glob "02_align/*.sorted.bam" --outdir results [options]

Purpose:
  Build TE/L1-relaxed ATAC BAM and bigWig tracks from pre-clean source BAMs.
  Use pre-clean sorted/raw BAMs when possible. Strict *.clean.bam cannot recover
  reads already removed by MAPQ, duplicate, or blacklist filters.

Required:
  --bam-glob STR                  Source BAM glob, quoted.

Options:
  --species STR                   hg38 | mm10 | mm39. Default: hg38
  --outdir DIR                    Default: ./atac_te_tracks
  --blacklist FILE                Optional override. Loaded from species config when possible
  --mito-chr STR                  Optional override. Loaded from species config when possible
  --effective-genome-size INT     Required only for --normalization RPGC
  --cores INT                     Default: 8
  --mapq INT                      Default: 0
  --exclude-flags INT             Default: 780; does not remove duplicate-marked reads
  --run-markdup BOOL              true | false. Default: false
  --remove-mito BOOL              true | false. Default: true
  --remove-blacklist BOOL         true | false. Default: true
  --proper-pair-only BOOL         true | false. Default: true for paired BAMs
  --normalization STR             CPM | RPGC | RPKM | BPM. Default: CPM
  --binsize INT                   Default: 10
  --ignore-for-normalization STR  Optional deepTools argument for RPGC
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-glob) BAM_GLOB="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --mito-chr) MITO_CHR="$2"; shift 2 ;;
    --effective-genome-size) EFFECTIVE_GENOME_SIZE="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --exclude-flags) EXCLUDE_FLAGS="$2"; shift 2 ;;
    --run-markdup) RUN_MARKDUP="$2"; shift 2 ;;
    --remove-mito) REMOVE_MITO="$2"; shift 2 ;;
    --remove-blacklist) REMOVE_BLACKLIST="$2"; shift 2 ;;
    --proper-pair-only) PROPER_PAIR_ONLY="$2"; shift 2 ;;
    --normalization) NORMALIZATION="$2"; shift 2 ;;
    --binsize) BINSIZE="$2"; shift 2 ;;
    --ignore-for-normalization) IGNORE_FOR_NORMALIZATION="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$BAM_GLOB" ]] || { echo "--bam-glob is required" >&2; usage >&2; exit 1; }
if [[ -f "$CONFIG_FILE" ]]; then
  [[ -n "$BLACKLIST" ]] || BLACKLIST="$(get_species_param "$SPECIES" blacklist "$CONFIG_FILE" || true)"
  [[ -n "$MITO_CHR" ]] || MITO_CHR="$(get_species_param "$SPECIES" mito_chr "$CONFIG_FILE" || true)"
  [[ -n "$EFFECTIVE_GENOME_SIZE" ]] || EFFECTIVE_GENOME_SIZE="$(get_species_param "$SPECIES" effective_genome_size "$CONFIG_FILE" || true)"
fi
[[ -n "$MITO_CHR" ]] || MITO_CHR="chrM"

case "${NORMALIZATION^^}" in
  CPM|RPGC|RPKM|BPM) NORMALIZATION="${NORMALIZATION^^}" ;;
  *) echo "--normalization must be CPM, RPGC, RPKM, or BPM" >&2; exit 1 ;;
esac
if [[ "$NORMALIZATION" == "RPGC" && -z "$EFFECTIVE_GENOME_SIZE" ]]; then
  echo "--effective-genome-size is required when --normalization RPGC" >&2
  exit 1
fi
if [[ "$REMOVE_BLACKLIST" == "true" ]]; then
  [[ -f "$BLACKLIST" ]] || { echo "blacklist not found: $BLACKLIST" >&2; exit 1; }
fi

mkdir -p "$OUTDIR/02_align_te" "$OUTDIR/04_bw_te" "$OUTDIR/logs"
shopt -s nullglob
bam_files=( $BAM_GLOB )
[[ ${#bam_files[@]} -gt 0 ]] || { echo "No BAM matched --bam-glob: $BAM_GLOB" >&2; exit 1; }
log "source BAM files: ${#bam_files[@]}"
log "TE/L1 settings: mapq=$MAPQ exclude_flags=$EXCLUDE_FLAGS markdup=$RUN_MARKDUP remove_mito=$REMOVE_MITO remove_blacklist=$REMOVE_BLACKLIST proper_pair_only=$PROPER_PAIR_ONLY normalization=$NORMALIZATION binsize=$BINSIZE"

for src_bam in "${bam_files[@]}"; do
  [[ -s "$src_bam" ]] || { echo "BAM not found or empty: $src_bam" >&2; exit 1; }
  base="$(basename "$src_bam")"
  sample="$base"
  sample="${sample%.bam}"
  sample="${sample%.sorted}"
  sample="${sample%.clean}"
  sample="${sample%.markdup}"

  work_prefix="$OUTDIR/02_align_te/${sample}.te"
  out_bam="$OUTDIR/02_align_te/${sample}.te.bam"
  out_bw="$OUTDIR/04_bw_te/${sample}.te.bw"
  counts_tsv="$OUTDIR/02_align_te/${sample}.te.clean_counts.tsv"
  flagstat_txt="$OUTDIR/02_align_te/${sample}.te.flagstat.txt"

  log "processing sample=$sample source=$src_bam"
  if [[ ! -s "${src_bam}.bai" && ! -s "${src_bam%.bam}.bai" ]]; then
    log "index source BAM for $sample"
    samtools index -@ "$CORES" "$src_bam"
  fi

  work_bam="$src_bam"
  if [[ "$RUN_MARKDUP" == "true" ]]; then
    log "remove duplicates for $sample"
    picard MarkDuplicates I="$src_bam" O="${work_prefix}.markdup.bam" M="${work_prefix}.markdup.metrics.txt" REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=false
    work_bam="${work_prefix}.markdup.bam"
  else
    log "keep duplicate-marked reads for $sample"
  fi

  paired_n="$(samtools view -c -f 1 "$work_bam")"
  pair_args=()
  if [[ "$PROPER_PAIR_ONLY" == "true" && "$paired_n" -gt 0 ]]; then
    pair_args=(-f 2)
  fi

  samtools view -@ "$CORES" -b -q "$MAPQ" "${pair_args[@]}" -F "$EXCLUDE_FLAGS" "$work_bam" > "${work_prefix}.mapq.bam"
  samtools index -@ "$CORES" "${work_prefix}.mapq.bam"

  if [[ "$REMOVE_MITO" == "true" ]]; then
    mapfile -t keep_chroms < <(samtools idxstats "${work_prefix}.mapq.bam" | awk -v mt="$MITO_CHR" '$1 != mt && $1 != "*" && $1 != "" {print $1}')
    printf '%s\n' "${keep_chroms[@]}" > "${work_prefix}.keep_chroms.txt"
    if [[ ${#keep_chroms[@]} -gt 0 ]]; then
      samtools view -@ "$CORES" -b "${work_prefix}.mapq.bam" "${keep_chroms[@]}" > "${work_prefix}.nomito.bam"
    else
      cp "${work_prefix}.mapq.bam" "${work_prefix}.nomito.bam"
    fi
  else
    cp "${work_prefix}.mapq.bam" "${work_prefix}.nomito.bam"
  fi

  if [[ "$REMOVE_BLACKLIST" == "true" ]]; then
    if [[ "$paired_n" -gt 0 ]]; then
      bedtools intersect -ubam -abam "${work_prefix}.nomito.bam" -b "$BLACKLIST" | \
        samtools view - | cut -f1 | sort -u > "${work_prefix}.blacklist.read_names.txt"
      if [[ -s "${work_prefix}.blacklist.read_names.txt" ]]; then
        samtools view -h "${work_prefix}.nomito.bam" | \
          awk -v bad="${work_prefix}.blacklist.read_names.txt" 'BEGIN{while((getline line < bad)>0) drop[line]=1} /^@/{print; next} !($1 in drop)' | \
          samtools view -@ "$CORES" -b - > "${work_prefix}.clean.unsorted.bam"
      else
        cp "${work_prefix}.nomito.bam" "${work_prefix}.clean.unsorted.bam"
      fi
    else
      bedtools intersect -v -abam "${work_prefix}.nomito.bam" -b "$BLACKLIST" > "${work_prefix}.clean.unsorted.bam"
    fi
  else
    cp "${work_prefix}.nomito.bam" "${work_prefix}.clean.unsorted.bam"
  fi

  samtools sort -@ "$CORES" -o "$out_bam" "${work_prefix}.clean.unsorted.bam"
  samtools index "$out_bam"
  samtools flagstat "$out_bam" > "$flagstat_txt"

  norm_args=(--normalizeUsing "$NORMALIZATION")
  if [[ "$NORMALIZATION" == "RPGC" ]]; then
    norm_args+=(--effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE")
    if [[ -n "$IGNORE_FOR_NORMALIZATION" ]]; then
      norm_args+=(--ignoreForNormalization "$IGNORE_FOR_NORMALIZATION")
    fi
  fi
  bamCoverage -b "$out_bam" -o "$out_bw" --binSize "$BINSIZE" -p "$CORES" "${norm_args[@]}"

  {
    echo -e "step\talignments"
    echo -e "input\t$(samtools view -c "$src_bam")"
    echo -e "post_markdup_step\t$(samtools view -c "$work_bam")"
    echo -e "post_mapq_flag_filter\t$(samtools view -c "${work_prefix}.mapq.bam")"
    echo -e "post_mito_step\t$(samtools view -c "${work_prefix}.nomito.bam")"
    echo -e "final\t$(samtools view -c "$out_bam")"
  } > "$counts_tsv"

  rm -f     "${work_prefix}.mapq.bam" "${work_prefix}.mapq.bam.bai"     "${work_prefix}.nomito.bam" "${work_prefix}.clean.unsorted.bam"     "${work_prefix}.blacklist.read_names.txt" "${work_prefix}.keep_chroms.txt"
  if [[ "$RUN_MARKDUP" == "true" ]]; then
    rm -f "${work_prefix}.markdup.bam" "${work_prefix}.markdup.bam.bai"
  fi

  log "done $sample: $out_bam ; $out_bw"
done

log "all TE/L1 relaxed tracks finished: $OUTDIR/04_bw_te"
