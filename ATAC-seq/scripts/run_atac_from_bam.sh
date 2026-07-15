#!/usr/bin/env bash
set -euo pipefail

resolve_script_dir() {
  local src="${BASH_SOURCE[0]}"
  while [ -L "$src" ]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    [[ "$src" != /* ]] && src="$dir/$src"
  done
  cd -P "$(dirname "$src")" && pwd
}
SCRIPT_DIR="$(resolve_script_dir)"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"

log(){ echo "[run_atac_from_bam] $*"; }

ORIG_ARGS=("$@")
BAM_DIR=""
SAMPLE_META=""
CONTRAST_FILE=""
SPECIES="hg38"
OUTDIR="./atac_from_bam"
TSS_BED=""
GENE_BODY_BED=""
TE_BED=""
MOTIF_GENOME=""
MOTIF_MEME=""
GENOME_FASTA=""
CHROM_SIZES=""
EFFECTIVE_GENOME_SIZE=""
MACS_GENOME_SIZE=""
CORES=8
CONS_HALF_WIDTH=250
RUN_BAMCOVERAGE="true"
RUN_TSS="false"
RUN_GENE_BODY="false"
RUN_TE_HEATMAP="false"
RUN_MOTIF="false"
RUN_TOBIAS="false"
RUN_NUC_PHASING="false"
RUN_FIXEDBIN="true"
RUN_BW_CORRELATION="true"
RUN_DIFF_PEAK_HEATMAP="false"
RUN_CONSENSUS_ANNOTATION="true"
FIXEDBIN_SIZE="100000"
TE_VIOLIN_CLASS=""
TE_VIOLIN_TOP_N="12"
TE_HEATMAP_PURPOSE="global"
BACKGROUND="false"
BACKGROUND_INTERNAL="false"

usage(){
  cat <<'USAGE'
Usage:
  run_atac_from_bam --bam-dir DIR --sample-meta FILE --contrast-file FILE [options]

Required:
  --bam-dir DIR              包含 *.clean.bam
  --sample-meta FILE         列: sample,condition,replicate
  --contrast-file FILE       列: case,control

Common options:
  --species STR              hg38 | mm10 | mm39
  --outdir DIR
  --cores INT
  --consensus-half-width INT
  --background

Optional refs (默认优先读 conf/species_refs.config):
  --genome FILE
  --chrom-sizes FILE
  --tss-bed FILE
  --gene-body-bed FILE
  --te-bed FILE
  --motif-genome STR
  --motif-meme FILE

Extra analysis switches:
  --run-bamcoverage BOOL
  --run-tss BOOL
  --run-gene-body BOOL
  --run-te-heatmap BOOL
  --run-motif BOOL
  --run-tobias BOOL
  --run-nuc-phasing BOOL
  --run-fixedbin BOOL
  --run-bw-correlation BOOL
  --run-diff-peak-heatmap BOOL
  --run-consensus-annotation BOOL
  --fixedbin-size INT
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --sample-meta) SAMPLE_META="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --genome) GENOME_FASTA="$2"; shift 2 ;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2 ;;
    --tss-bed) TSS_BED="$2"; shift 2 ;;
    --gene-body-bed) GENE_BODY_BED="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --motif-genome) MOTIF_GENOME="$2"; shift 2 ;;
    --motif-meme) MOTIF_MEME="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --consensus-half-width) CONS_HALF_WIDTH="$2"; shift 2 ;;
    --run-bamcoverage) RUN_BAMCOVERAGE="$2"; shift 2 ;;
    --run-tss) RUN_TSS="$2"; shift 2 ;;
    --run-gene-body) RUN_GENE_BODY="$2"; shift 2 ;;
    --run-te-heatmap) RUN_TE_HEATMAP="$2"; shift 2 ;;
    --run-motif) RUN_MOTIF="$2"; shift 2 ;;
    --run-tobias) RUN_TOBIAS="$2"; shift 2 ;;
    --run-nuc-phasing) RUN_NUC_PHASING="$2"; shift 2 ;;
    --run-fixedbin) RUN_FIXEDBIN="$2"; shift 2 ;;
    --run-bw-correlation) RUN_BW_CORRELATION="$2"; shift 2 ;;
    --run-diff-peak-heatmap) RUN_DIFF_PEAK_HEATMAP="$2"; shift 2 ;;
    --run-consensus-annotation) RUN_CONSENSUS_ANNOTATION="$2"; shift 2 ;;
    --fixedbin-size) FIXEDBIN_SIZE="$2"; shift 2 ;;
    --te-violin-class) TE_VIOLIN_CLASS="$2"; shift 2 ;;
    --te-violin-top-n) TE_VIOLIN_TOP_N="$2"; shift 2 ;;
    --te-heatmap-purpose) TE_HEATMAP_PURPOSE="$2"; shift 2 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --background-internal) BACKGROUND_INTERNAL="true"; BACKGROUND="false"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -d "$BAM_DIR" ]] || { echo "--bam-dir not found: $BAM_DIR" >&2; exit 1; }
[[ -f "$SAMPLE_META" ]] || { echo "--sample-meta not found: $SAMPLE_META" >&2; exit 1; }
[[ -f "$CONTRAST_FILE" ]] || { echo "--contrast-file not found: $CONTRAST_FILE" >&2; exit 1; }

if [[ "$BACKGROUND" == "true" && "$BACKGROUND_INTERNAL" != "true" && "${RUN_ATAC_FROM_BAM_BG_CHILD:-0}" != "1" ]]; then
  mkdir -p "$OUTDIR/_logs"
  TS="$(date +%Y%m%d_%H%M%S)"
  LOG_FILE="$OUTDIR/_logs/from_bam_${TS}.log"
  PID_FILE="$OUTDIR/_logs/from_bam_${TS}.pid"
  FILTERED_ARGS=()
  for x in "${ORIG_ARGS[@]}"; do
    [[ "$x" == "--background" ]] && continue
    [[ "$x" == "--background-internal" ]] && continue
    FILTERED_ARGS+=("$x")
  done
  log "后台运行，日志: $LOG_FILE"
  nohup env RUN_ATAC_FROM_BAM_BG_CHILD=1 bash "$0" "${FILTERED_ARGS[@]}" --background-internal > "$LOG_FILE" 2>&1 &
  echo $! > "$PID_FILE"
  log "PID: $(cat "$PID_FILE")"
  exit 0
fi

[[ -f "$CONFIG_FILE" ]] || { echo "config not found: $CONFIG_FILE" >&2; exit 1; }
[[ -n "$GENOME_FASTA" ]] || GENOME_FASTA="$(get_species_param "$SPECIES" genome_fasta "$CONFIG_FILE" || true)"
[[ -n "$CHROM_SIZES" ]] || CHROM_SIZES="$(get_species_param "$SPECIES" chrom_sizes "$CONFIG_FILE" || true)"
[[ -n "$TSS_BED" ]] || TSS_BED="$(get_species_param "$SPECIES" tss_bed "$CONFIG_FILE" || true)"
[[ -n "$GENE_BODY_BED" ]] || GENE_BODY_BED="$(get_species_param "$SPECIES" gene_body_bed "$CONFIG_FILE" || true)"
[[ -n "$TE_BED" ]] || TE_BED="$(get_species_param "$SPECIES" te_bed "$CONFIG_FILE" || true)"
[[ -n "$MOTIF_GENOME" ]] || MOTIF_GENOME="$(get_species_param "$SPECIES" motif_genome "$CONFIG_FILE" || true)"
[[ -n "$MOTIF_MEME" ]] || MOTIF_MEME="$(get_species_param "$SPECIES" motif_meme "$CONFIG_FILE" || true)"
EFFECTIVE_GENOME_SIZE="$(get_species_param "$SPECIES" effective_genome_size "$CONFIG_FILE" || true)"
MACS_GENOME_SIZE="$(get_species_param "$SPECIES" genome_size "$CONFIG_FILE" || true)"
GTF_FILE="$(get_species_param "$SPECIES" gtf_genes "$CONFIG_FILE" || true)"

[[ -f "$GENOME_FASTA" ]] || { echo "genome fasta not found: $GENOME_FASTA" >&2; exit 1; }
[[ -f "$CHROM_SIZES" ]] || { echo "chrom sizes not found: $CHROM_SIZES" >&2; exit 1; }
[[ -n "$MACS_GENOME_SIZE" ]] || { echo "genome_size missing for species $SPECIES in config" >&2; exit 1; }

mkdir -p "$OUTDIR" "$OUTDIR/03_qc" "$OUTDIR/04_bw" "$OUTDIR/05_peaks" "$OUTDIR/06_consensus_peaks" "$OUTDIR/07_counts" "$OUTDIR/08_downstream"
shopt -s nullglob
bam_files=( "$BAM_DIR"/*.clean.bam )
[[ ${#bam_files[@]} -gt 0 ]] || { echo "No *.clean.bam found in: $BAM_DIR" >&2; exit 1; }
log "clean bam files: ${#bam_files[@]}"

all_paired=true
qc_summary="$OUTDIR/03_qc/atac_qc_summary.tsv"
echo -e "sample\ttotal_raw\ttotal_clean\tmt_reads\ttotal_reads\tin_peak_reads\tfrip" > "$qc_summary"

for bam in "${bam_files[@]}"; do
  sample=$(basename "$bam" .clean.bam)
  if [[ ! -s "${bam}.bai" || "$bam" -nt "${bam}.bai" ]]; then
    log "index bam for $sample"
    samtools index -@ "$CORES" "$bam"
  fi
  paired_n=$(samtools view -c -f 1 "$bam")
  if [[ "$paired_n" -gt 0 ]]; then
    log "call peaks (PE) for $sample"
    macs3 callpeak -t "$bam" -f BAMPE -g "$MACS_GENOME_SIZE" -n "$sample" --nomodel --keep-dup all --call-summits -p 0.01 --outdir "$OUTDIR/05_peaks"
  else
    all_paired=false
    log "call peaks (SE) for $sample"
    macs3 callpeak -t "$bam" -f BAM -g "$MACS_GENOME_SIZE" -n "$sample" --nomodel --shift -100 --extsize 200 --keep-dup all --call-summits -p 0.01 --outdir "$OUTDIR/05_peaks"
  fi
  if [[ -s "$OUTDIR/05_peaks/${sample}_summits.bed" ]]; then
    cp "$OUTDIR/05_peaks/${sample}_summits.bed" "$OUTDIR/05_peaks/${sample}.summits.bed"
  else
    awk 'BEGIN{OFS="\t"} {c=int(($2+$3)/2); s=c; e=c+1; if(s<0)s=0; print $1,s,e}' "$OUTDIR/05_peaks/${sample}_peaks.narrowPeak" > "$OUTDIR/05_peaks/${sample}.summits.bed"
  fi

  total_clean=$(samtools view -c "$bam")
  mt_reads=$(samtools idxstats "$bam" | awk '/^chrM\t|^MT\t/{s+=$3} END{print s+0}')
  if [[ -s "$OUTDIR/05_peaks/${sample}_peaks.narrowPeak" ]]; then
    peak_file="$OUTDIR/05_peaks/${sample}_peaks.narrowPeak"
  elif [[ -s "$OUTDIR/05_peaks/${sample}_peaks.broadPeak" ]]; then
    peak_file="$OUTDIR/05_peaks/${sample}_peaks.broadPeak"
  else
    echo "Peak file missing for $sample" >&2
    exit 1
  fi
  in_peak=$(bedtools intersect -u -abam "$bam" -b "$peak_file" | samtools view -c)
  frip=$(awk -v tp="$total_clean" -v ip="$in_peak" 'BEGIN{print (tp==0 ? 0 : ip/tp)}')
  echo -e "$sample\t$total_clean\t$total_clean\t$mt_reads\t$total_clean\t$in_peak\t$frip" >> "$qc_summary"
done

Rscript "$ROOT_DIR/scripts/plot_atac_qc.R" --qc-summary "$qc_summary" --outdir "$OUTDIR/03_qc"

log "build consensus peaks"
cat "$OUTDIR"/05_peaks/*.summits.bed | cut -f1-3 | sort -k1,1 -k2,2n | \
  bedtools slop -i - -g "$CHROM_SIZES" -b "$CONS_HALF_WIDTH" | \
  bedtools sort -i - | bedtools merge -i - > "$OUTDIR/06_consensus_peaks/consensus_peaks.bed"
[[ -s "$OUTDIR/06_consensus_peaks/consensus_peaks.bed" ]] || { echo "Consensus peak file is empty" >&2; exit 1; }

awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print "peak_"NR,$1,$2+1,$3,"+"}' \
  "$OUTDIR/06_consensus_peaks/consensus_peaks.bed" > "$OUTDIR/07_counts/consensus_peaks.saf"

log "count consensus peaks"
if [[ "$all_paired" == "true" ]]; then
  featureCounts -a "$OUTDIR/07_counts/consensus_peaks.saf" -F SAF -o "$OUTDIR/07_counts/consensus_peak_counts.txt" -T "$CORES" -p --countReadPairs -B -C "$BAM_DIR"/*.clean.bam
else
  featureCounts -a "$OUTDIR/07_counts/consensus_peaks.saf" -F SAF -o "$OUTDIR/07_counts/consensus_peak_counts.txt" -T "$CORES" "$BAM_DIR"/*.clean.bam
fi
cp "$SAMPLE_META" "$OUTDIR/07_counts/sample_metadata.csv"

need_bw="false"
if [[ "$RUN_BAMCOVERAGE" == "true" || "$RUN_TSS" == "true" || "$RUN_GENE_BODY" == "true" || "$RUN_TE_HEATMAP" == "true" ]]; then
  need_bw="true"
fi

if [[ "$need_bw" == "true" ]]; then
  [[ -n "$EFFECTIVE_GENOME_SIZE" ]] || { echo "effective_genome_size missing for species $SPECIES" >&2; exit 1; }
  for bam in "${bam_files[@]}"; do
    sample=$(basename "$bam" .clean.bam)
    log "bamCoverage $sample"
    bamCoverage -b "$bam" -o "$OUTDIR/04_bw/${sample}.bw" --normalizeUsing RPGC --effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE" --binSize 10 -p "$CORES"
  done
fi

log "run downstream R"
CMD=(Rscript "$ROOT_DIR/atacseq-downstream/run_downstream_atac.R"
  --count-file "$OUTDIR/07_counts/consensus_peak_counts.txt"
  --sample-meta "$OUTDIR/07_counts/sample_metadata.csv"
  --peak-bed "$OUTDIR/06_consensus_peaks/consensus_peaks.bed"
  --outdir "$OUTDIR/08_downstream/peak_level"
  --species "$SPECIES"
  --contrast-file "$CONTRAST_FILE")
if [[ -n "$TE_BED" ]]; then CMD+=(--te-bed "$TE_BED"); fi
if [[ -n "$GTF_FILE" ]]; then CMD+=(--gtf "$GTF_FILE"); fi
if [[ -n "$TE_VIOLIN_CLASS" ]]; then CMD+=(--te-violin-class "$TE_VIOLIN_CLASS"); fi
if [[ -n "$TE_VIOLIN_TOP_N" ]]; then CMD+=(--te-violin-top-n "$TE_VIOLIN_TOP_N"); fi
"${CMD[@]}"

if [[ "$RUN_CONSENSUS_ANNOTATION" == "true" ]]; then
  log "annotate consensus peaks"
  bash "$SCRIPT_DIR/run_peak_annotation.sh" \
    --input "$OUTDIR/06_consensus_peaks/consensus_peaks.bed" \
    --species "$SPECIES" \
    --out-prefix "$OUTDIR/08_downstream/peak_level/annotation/consensus/consensus_peaks" \
    ${GTF_FILE:+--gtf "$GTF_FILE"} \
    ${TE_BED:+--te-bed "$TE_BED"}
fi

if [[ "$RUN_BW_CORRELATION" == "true" && "$need_bw" == "true" ]]; then
  log "run bigWig correlation/PCA"
  bash "$SCRIPT_DIR/run_bw_correlation.sh" \
    --bw-glob "$OUTDIR/04_bw/*.bw" \
    --out-prefix "$OUTDIR/03_qc/bw_correlation/all_samples" \
    --threads "$CORES"
fi

if [[ "$RUN_DIFF_PEAK_HEATMAP" == "true" && "$need_bw" == "true" ]]; then
  log "run differential peak heatmaps"
  bash "$SCRIPT_DIR/run_diff_peak_heatmap.sh" \
    --contrast-bed-dir "$OUTDIR/08_downstream/peak_level/contrast_beds" \
    --bw-glob "$OUTDIR/04_bw/*.bw" \
    --outdir "$OUTDIR/08_downstream/peak_level/region_heatmap/differential_peaks" \
    --threads "$CORES"
fi

if [[ "$RUN_FIXEDBIN" == "true" ]]; then
  log "run fixed-bin downstream (${FIXEDBIN_SIZE} bp)"
  bash "$SCRIPT_DIR/run_fixedbin_from_bam.sh" \
    --bam-glob "$BAM_DIR/*.clean.bam" \
    --chrom-sizes "$CHROM_SIZES" \
    --sample-meta "$OUTDIR/07_counts/sample_metadata.csv" \
    --contrast-file "$CONTRAST_FILE" \
    --species "$SPECIES" \
    --outdir "$OUTDIR/08_downstream/bin_level" \
    --bin-size "$FIXEDBIN_SIZE" \
    --cores "$CORES" \
    ${TE_BED:+--te-bed "$TE_BED"} \
    ${GTF_FILE:+--gtf "$GTF_FILE"} \
    ${TE_VIOLIN_CLASS:+--te-violin-class "$TE_VIOLIN_CLASS"} \
    ${TE_VIOLIN_TOP_N:+--te-violin-top-n "$TE_VIOLIN_TOP_N"}
fi

if [[ "$RUN_TSS" == "true" ]]; then
  [[ -f "$TSS_BED" ]] || { echo "tss bed not found: $TSS_BED" >&2; exit 1; }
  log "run TSS QC"
  bash "$SCRIPT_DIR/run_tss_qc.sh" --bw-glob "$OUTDIR/04_bw/*.bw" --species "$SPECIES" --tss-bed "$TSS_BED" --outdir "$OUTDIR/03_qc/tss_enrichment_manual" --cores "$CORES"
fi
if [[ "$RUN_GENE_BODY" == "true" ]]; then
  [[ -f "$GENE_BODY_BED" ]] || { echo "gene body bed not found: $GENE_BODY_BED" >&2; exit 1; }
  log "run gene body profile"
  bash "$SCRIPT_DIR/run_gene_body_profile.sh" --bw-glob "$OUTDIR/04_bw/*.bw" --species "$SPECIES" --gene-bed "$GENE_BODY_BED" --outdir "$OUTDIR/03_qc/gene_body_manual" --cores "$CORES"
fi
if [[ "$RUN_TE_HEATMAP" == "true" ]]; then
  [[ -f "$TE_BED" ]] || { echo "te annotation not found: $TE_BED" >&2; exit 1; }
  log "run TE heatmap"
  bash "$SCRIPT_DIR/run_te_heatmap.sh" --bw-glob "$OUTDIR/04_bw/*.bw" --species "$SPECIES" --te-bed "$TE_BED" --peak-bed "$OUTDIR/06_consensus_peaks/consensus_peaks.bed" --outdir "$OUTDIR/03_qc/te_heatmap_manual" --cores "$CORES" --purpose "$TE_HEATMAP_PURPOSE"
fi
if [[ "$RUN_MOTIF" == "true" ]]; then
  [[ -n "$MOTIF_GENOME" ]] || { echo "motif genome missing" >&2; exit 1; }
  log "run motif"
  bash "$SCRIPT_DIR/run_motif_homer.sh" --bed-dir "$OUTDIR/08_downstream/peak_level/contrast_beds" --genome "$MOTIF_GENOME" --outdir "$OUTDIR/09_motif_manual"
fi
if [[ "$RUN_TOBIAS" == "true" ]]; then
  [[ -f "$MOTIF_MEME" ]] || { echo "motif meme not found: $MOTIF_MEME" >&2; exit 1; }
  log "run TOBIAS"
  bash "$SCRIPT_DIR/run_tobias.sh" --bam-dir "$BAM_DIR" --peaks "$OUTDIR/06_consensus_peaks/consensus_peaks.bed" --sample-meta "$OUTDIR/07_counts/sample_metadata.csv" --contrast-file "$CONTRAST_FILE" --genome "$GENOME_FASTA" --motif-meme "$MOTIF_MEME" --outdir "$OUTDIR/10_footprinting_manual" --cores "$CORES"
fi

if [[ "$RUN_NUC_PHASING" == "true" ]]; then
  log "run nucleosome phasing"
  bash "$SCRIPT_DIR/run_nuc_phasing.sh" --bam-glob "$BAM_DIR/*.clean.bam" --outdir "$OUTDIR/11_nuc_phasing" --cores "$CORES"
fi

log "done"
