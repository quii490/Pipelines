#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[$(date '+%F %T')] [run_atac_te_tracks_from_fastq] $*"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"

SAMPLESHEET=""
SAMPLE=""
LAYOUT=""
FASTQ1=""
FASTQ2=""
OUTDIR="./atac_te_tracks_from_fastq"
SPECIES="hg38"
BOWTIE2_INDEX=""
BLACKLIST=""
MITO_CHR=""
EFFECTIVE_GENOME_SIZE=""
CORES=8
RUN_FASTP="true"
FASTP_TIMEOUT="120m"
MAX_INSERT=2000
MAPQ=0
EXCLUDE_FLAGS=780
RUN_MARKDUP="false"
REMOVE_MITO="true"
REMOVE_BLACKLIST="true"
PROPER_PAIR_ONLY="true"
NORMALIZATION="CPM"
BINSIZE=10
IGNORE_FOR_NORMALIZATION=""
KEEP_INTERMEDIATE="false"

usage(){
  cat <<'USAGE'
Usage:
  # Single sample, paired-end
  run_atac_te_tracks_from_fastq.sh \
    --sample SAMPLE --layout PE --fastq1 R1.fq.gz --fastq2 R2.fq.gz \
    --outdir results --species hg38 [options]

  # Single sample, single-end
  run_atac_te_tracks_from_fastq.sh \
    --sample SAMPLE --layout SE --fastq1 reads.fq.gz \
    --outdir results --species hg38 [options]

  # Batch mode: same CSV columns as main workflow: sample,layout,r1,r2,condition,replicate
  run_atac_te_tracks_from_fastq.sh \
    --samplesheet samplesheet.csv --outdir results --species hg38 [options]

Purpose:
  Build pre-clean sorted ATAC BAM + TE/L1-relaxed BAM + TE/L1 bigWig directly from FASTQ.
  This is a standalone diagnostic tool for repeat/TE/L1 tracks. It does not replace
  the standard ATAC peak-calling workflow.

Required, choose one mode:
  --samplesheet FILE              Batch CSV with sample,layout,r1,r2 columns
  OR
  --sample STR --layout PE|SE --fastq1 FILE [--fastq2 FILE]

Reference options:
  --species STR                   hg38 | mm10 | mm39. Default: hg38
  --bowtie2-index PREFIX          Optional override. Loaded from species config when possible
  --blacklist FILE                Optional override. Loaded from species config when possible
  --mito-chr STR                  Optional override. Loaded from species config when possible
  --effective-genome-size INT     Required only for --normalization RPGC

Processing options:
  --outdir DIR                    Default: ./atac_te_tracks_from_fastq
  --cores INT                     Default: 8
  --run-fastp BOOL                true | false. Default: true
  --fastp-timeout DURATION        Timeout for each fastp run. Default: 120m; use 0 to disable
  --max-insert INT                Bowtie2 -X for PE. Default: 2000
  --mapq INT                      TE/L1 relaxed MAPQ. Default: 0
  --exclude-flags INT             Default: 780; does not remove duplicate-marked reads
  --run-markdup BOOL              true | false. Default: false
  --remove-mito BOOL              true | false. Default: true
  --remove-blacklist BOOL         true | false. Default: true
  --proper-pair-only BOOL         true | false. Default: true for PE
  --normalization STR             CPM | RPGC | RPKM | BPM. Default: CPM
  --binsize INT                   bigWig bin size. Default: 10
  --ignore-for-normalization STR  Optional deepTools argument for RPGC
  --keep-intermediate BOOL        true | false. Default: false

Outputs:
  <outdir>/00_fastq_merged/*.merged.fastq(.gz)     merged FASTQ symlinks/files
  <outdir>/01_fastp/*                              fastp reports and clean FASTQs when enabled
  <outdir>/02_align_sorted/*.sorted.bam(.bai)      pre-clean sorted BAM, suitable as TE source BAM
  <outdir>/02_align_te/*.te.bam(.bai)              TE/L1-relaxed BAM
  <outdir>/02_align_te/*.te.clean_counts.tsv       stepwise alignment counts
  <outdir>/02_align_te/*.te.flagstat.txt           TE BAM flagstat
  <outdir>/04_bw_te/*.te.bw                        TE/L1-relaxed bigWig
  <outdir>/logs/*.bowtie2.log                      alignment logs
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --fastq1) FASTQ1="$2"; shift 2 ;;
    --fastq2) FASTQ2="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --bowtie2-index) BOWTIE2_INDEX="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --mito-chr) MITO_CHR="$2"; shift 2 ;;
    --effective-genome-size) EFFECTIVE_GENOME_SIZE="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --run-fastp) RUN_FASTP="$2"; shift 2 ;;
    --fastp-timeout) FASTP_TIMEOUT="$2"; shift 2 ;;
    --max-insert) MAX_INSERT="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --exclude-flags) EXCLUDE_FLAGS="$2"; shift 2 ;;
    --run-markdup) RUN_MARKDUP="$2"; shift 2 ;;
    --remove-mito) REMOVE_MITO="$2"; shift 2 ;;
    --remove-blacklist) REMOVE_BLACKLIST="$2"; shift 2 ;;
    --proper-pair-only) PROPER_PAIR_ONLY="$2"; shift 2 ;;
    --normalization) NORMALIZATION="$2"; shift 2 ;;
    --binsize) BINSIZE="$2"; shift 2 ;;
    --ignore-for-normalization) IGNORE_FOR_NORMALIZATION="$2"; shift 2 ;;
    --keep-intermediate) KEEP_INTERMEDIATE="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -f "$CONFIG_FILE" ]]; then
  [[ -n "$BOWTIE2_INDEX" ]] || BOWTIE2_INDEX="$(get_species_param "$SPECIES" bowtie2_index "$CONFIG_FILE" || true)"
  [[ -n "$BLACKLIST" ]] || BLACKLIST="$(get_species_param "$SPECIES" blacklist "$CONFIG_FILE" || true)"
  [[ -n "$MITO_CHR" ]] || MITO_CHR="$(get_species_param "$SPECIES" mito_chr "$CONFIG_FILE" || true)"
  [[ -n "$EFFECTIVE_GENOME_SIZE" ]] || EFFECTIVE_GENOME_SIZE="$(get_species_param "$SPECIES" effective_genome_size "$CONFIG_FILE" || true)"
fi
[[ -n "$BOWTIE2_INDEX" ]] || { echo "--bowtie2-index is required or must be available in species config" >&2; exit 1; }
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

require_cmd() {
  local cmd="$1"
  command -v "$cmd" >/dev/null 2>&1 || {
    echo "Required command not found in PATH: $cmd" >&2
    echo "Activate the correct conda environment first, for example: conda activate chipseq" >&2
    exit 1
  }
}
[[ "$RUN_FASTP" == "true" ]] && require_cmd fastp
require_cmd bowtie2
require_cmd samtools
require_cmd bamCoverage
[[ "$REMOVE_BLACKLIST" == "true" ]] && require_cmd bedtools
[[ "$RUN_MARKDUP" == "true" ]] && require_cmd picard

if [[ -n "$SAMPLESHEET" ]]; then
  [[ -f "$SAMPLESHEET" ]] || { echo "samplesheet not found: $SAMPLESHEET" >&2; exit 1; }
else
  [[ -n "$SAMPLE" && -n "$LAYOUT" && -n "$FASTQ1" ]] || { echo "Use --samplesheet or provide --sample --layout --fastq1 [--fastq2]" >&2; usage >&2; exit 1; }
  LAYOUT="${LAYOUT^^}"
  [[ "$LAYOUT" == "PE" || "$LAYOUT" == "SE" ]] || { echo "--layout must be PE or SE" >&2; exit 1; }
  if [[ "$LAYOUT" == "PE" && -z "$FASTQ2" ]]; then
    echo "--fastq2 is required for PE" >&2
    exit 1
  fi
fi

mkdir -p "$OUTDIR/00_fastq_merged" "$OUTDIR/01_fastp" "$OUTDIR/02_align_sorted" "$OUTDIR/02_align_te" "$OUTDIR/04_bw_te" "$OUTDIR/logs"
SAMPLE_TABLE="$OUTDIR/logs/samples.resolved.tsv"

if [[ -n "$SAMPLESHEET" ]]; then
  python - "$SAMPLESHEET" > "$SAMPLE_TABLE" <<'PY'
import csv, sys
path = sys.argv[1]
with open(path, newline='') as fh:
    reader = csv.DictReader(fh)
    required = {'sample', 'layout', 'r1'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"samplesheet missing columns: {','.join(sorted(missing))}")
    print("sample\tlayout\tr1\tr2")
    for row in reader:
        sample = (row.get('sample') or '').strip()
        layout = (row.get('layout') or '').strip().upper()
        r1 = (row.get('r1') or '').strip()
        r2 = (row.get('r2') or '').strip()
        if not sample or not layout or not r1:
            raise SystemExit(f"invalid row: {row}")
        if layout == 'PE' and not r2:
            raise SystemExit(f"PE sample {sample} missing r2")
        if layout not in {'PE','SE'}:
            raise SystemExit(f"invalid layout for {sample}: {layout}")
        print(f"{sample}\t{layout}\t{r1}\t{r2}")
PY
else
  {
    echo -e "sample\tlayout\tr1\tr2"
    echo -e "${SAMPLE}\t${LAYOUT}\t${FASTQ1}\t${FASTQ2}"
  } > "$SAMPLE_TABLE"
fi

split_csv_to_array() {
  local input="$1"
  local -n out_arr="$2"
  IFS=',' read -r -a out_arr <<< "$input"
}

make_merged_fastq() {
  local sample="$1"
  local mate="$2"
  local files_csv="$3"
  local -n out_var="$4"
  local files=()
  split_csv_to_array "$files_csv" files
  [[ ${#files[@]} -gt 0 && -n "${files[0]}" ]] || { echo "No FASTQ files for $sample $mate" >&2; exit 1; }
  for fq in "${files[@]}"; do
    [[ -s "$fq" ]] || { echo "FASTQ not found or empty: $fq" >&2; exit 1; }
  done
  local ext="fastq"
  [[ "${files[0]}" == *.gz ]] && ext="fastq.gz"
  local out="$OUTDIR/00_fastq_merged/${sample}.${mate}.merged.${ext}"
  if [[ ${#files[@]} -eq 1 ]]; then
    ln -sf "$(readlink -f "${files[0]}")" "$out"
  else
    cat "${files[@]}" > "$out"
  fi
  out_var="$out"
}

build_norm_args() {
  NORM_ARGS=(--normalizeUsing "$NORMALIZATION")
  if [[ "$NORMALIZATION" == "RPGC" ]]; then
    NORM_ARGS+=(--effectiveGenomeSize "$EFFECTIVE_GENOME_SIZE")
    if [[ -n "$IGNORE_FOR_NORMALIZATION" ]]; then
      NORM_ARGS+=(--ignoreForNormalization "$IGNORE_FOR_NORMALIZATION")
    fi
  fi
}

run_with_optional_timeout() {
  local timeout_value="$1"
  local log_file="$2"
  shift 2
  if [[ "$timeout_value" == "0" || "$timeout_value" == "none" || "$timeout_value" == "NONE" ]]; then
    "$@" > "$log_file" 2>&1
  elif command -v timeout >/dev/null 2>&1; then
    timeout "$timeout_value" "$@" > "$log_file" 2>&1
  else
    log "WARN: timeout command not found; running without timeout"
    "$@" > "$log_file" 2>&1
  fi
}

process_one_sample() {
  local sample="$1" layout="$2" r1_csv="$3" r2_csv="$4"
  local merged_r1="" merged_r2="" align_r1="" align_r2=""
  log "processing sample=$sample layout=$layout"

  make_merged_fastq "$sample" "R1" "$r1_csv" merged_r1
  log "FASTQ R1 ready for $sample: $merged_r1"
  if [[ "$layout" == "PE" ]]; then
    make_merged_fastq "$sample" "R2" "$r2_csv" merged_r2
    log "FASTQ R2 ready for $sample: $merged_r2"
  fi

  if [[ "$RUN_FASTP" == "true" ]]; then
    local fastp_log="$OUTDIR/logs/${sample}.fastp.log"
    if [[ "$layout" == "PE" ]]; then
      align_r1="$OUTDIR/01_fastp/${sample}.R1.clean.fastq.gz"
      align_r2="$OUTDIR/01_fastp/${sample}.R2.clean.fastq.gz"
      if [[ -s "$align_r1" && -s "$align_r2" && -s "$OUTDIR/01_fastp/${sample}.fastp.json" ]]; then
        log "fastp outputs exist for $sample; skip fastp"
      else
        rm -f "$align_r1" "$align_r2" "$OUTDIR/01_fastp/${sample}.fastp.html" "$OUTDIR/01_fastp/${sample}.fastp.json" "$fastp_log"
        log "running fastp for $sample; log=$fastp_log; timeout=$FASTP_TIMEOUT"
        set +e
        run_with_optional_timeout "$FASTP_TIMEOUT" "$fastp_log" \
          fastp --thread "$CORES" -i "$merged_r1" -I "$merged_r2" -o "$align_r1" -O "$align_r2" \
            -h "$OUTDIR/01_fastp/${sample}.fastp.html" -j "$OUTDIR/01_fastp/${sample}.fastp.json"
        status=$?
        set -e
        if [[ "$status" -ne 0 ]]; then
          log "ERROR: fastp failed or timed out for $sample with status=$status. Last log lines:"
          tail -n 60 "$fastp_log" >&2 || true
          exit "$status"
        fi
        [[ -s "$align_r1" && -s "$align_r2" ]] || { echo "fastp did not create clean FASTQ for $sample; see $fastp_log" >&2; exit 1; }
      fi
    else
      align_r1="$OUTDIR/01_fastp/${sample}.SE.clean.fastq.gz"
      if [[ -s "$align_r1" && -s "$OUTDIR/01_fastp/${sample}.fastp.json" ]]; then
        log "fastp outputs exist for $sample; skip fastp"
      else
        rm -f "$align_r1" "$OUTDIR/01_fastp/${sample}.fastp.html" "$OUTDIR/01_fastp/${sample}.fastp.json" "$fastp_log"
        log "running fastp for $sample; log=$fastp_log; timeout=$FASTP_TIMEOUT"
        set +e
        run_with_optional_timeout "$FASTP_TIMEOUT" "$fastp_log" \
          fastp --thread "$CORES" -i "$merged_r1" -o "$align_r1" \
            -h "$OUTDIR/01_fastp/${sample}.fastp.html" -j "$OUTDIR/01_fastp/${sample}.fastp.json"
        status=$?
        set -e
        if [[ "$status" -ne 0 ]]; then
          log "ERROR: fastp failed or timed out for $sample with status=$status. Last log lines:"
          tail -n 60 "$fastp_log" >&2 || true
          exit "$status"
        fi
        [[ -s "$align_r1" ]] || { echo "fastp did not create clean FASTQ for $sample; see $fastp_log" >&2; exit 1; }
      fi
    fi
    log "fastp finished for $sample"
  else
    log "fastp skipped for $sample"
    align_r1="$merged_r1"
    align_r2="$merged_r2"
  fi

  local raw_bam="$OUTDIR/02_align_sorted/${sample}.raw.bam"
  local sorted_bam="$OUTDIR/02_align_sorted/${sample}.sorted.bam"
  local bowtie_log="$OUTDIR/logs/${sample}.bowtie2.log"

  log "running bowtie2 for $sample; log=$bowtie_log; raw BAM will grow at $raw_bam"
  if [[ "$layout" == "PE" ]]; then
    bowtie2 -p "$CORES" --very-sensitive --no-mixed --no-discordant -X "$MAX_INSERT" \
      -x "$BOWTIE2_INDEX" \
      --rg-id "$sample" --rg "SM:${sample}" --rg "LB:lib1" --rg "PL:ILLUMINA" \
      -1 "$align_r1" -2 "$align_r2" 2> "$bowtie_log" | \
      samtools view -@ "$CORES" -bS - > "$raw_bam"
  else
    bowtie2 -p "$CORES" --very-sensitive \
      -x "$BOWTIE2_INDEX" \
      --rg-id "$sample" --rg "SM:${sample}" --rg "LB:lib1" --rg "PL:ILLUMINA" \
      -U "$align_r1" 2> "$bowtie_log" | \
      samtools view -@ "$CORES" -bS - > "$raw_bam"
  fi
  log "bowtie2/samtools view finished for $sample; raw_bam=$raw_bam"
  log "sorting BAM for $sample"
  samtools sort -@ "$CORES" -o "$sorted_bam" "$raw_bam"
  samtools index -@ "$CORES" "$sorted_bam"
  log "wrote sorted BAM: $sorted_bam"

  local work_prefix="$OUTDIR/02_align_te/${sample}.te"
  local work_bam="$sorted_bam"
  local out_bam="$OUTDIR/02_align_te/${sample}.te.bam"
  local out_bw="$OUTDIR/04_bw_te/${sample}.te.bw"
  local counts_tsv="$OUTDIR/02_align_te/${sample}.te.clean_counts.tsv"
  local flagstat_txt="$OUTDIR/02_align_te/${sample}.te.flagstat.txt"

  if [[ "$RUN_MARKDUP" == "true" ]]; then
    log "remove duplicates for $sample"
    picard MarkDuplicates I="$sorted_bam" O="${work_prefix}.markdup.bam" M="${work_prefix}.markdup.metrics.txt" REMOVE_DUPLICATES=true ASSUME_SORTED=true CREATE_INDEX=false
    work_bam="${work_prefix}.markdup.bam"
  else
    log "keep duplicate-marked reads for $sample"
  fi

  pair_args=()
  if [[ "$layout" == "PE" && "$PROPER_PAIR_ONLY" == "true" ]]; then
    pair_args=(-f 2)
  fi
  log "filtering TE/L1 BAM for $sample: MAPQ=$MAPQ exclude_flags=$EXCLUDE_FLAGS proper_pair_only=$PROPER_PAIR_ONLY"
  samtools view -@ "$CORES" -b -q "$MAPQ" "${pair_args[@]}" -F "$EXCLUDE_FLAGS" "$work_bam" > "${work_prefix}.mapq.bam"
  samtools index -@ "$CORES" "${work_prefix}.mapq.bam"

  if [[ "$REMOVE_MITO" == "true" ]]; then
    log "removing mitochondrial contig for $sample: $MITO_CHR"
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
    log "removing blacklist regions for $sample: $BLACKLIST"
    bedtools intersect -v -abam "${work_prefix}.nomito.bam" -b "$BLACKLIST" > "${work_prefix}.clean.unsorted.bam"
  else
    cp "${work_prefix}.nomito.bam" "${work_prefix}.clean.unsorted.bam"
  fi

  log "sorting/indexing final TE/L1 BAM for $sample"
  samtools sort -@ "$CORES" -o "$out_bam" "${work_prefix}.clean.unsorted.bam"
  samtools index -@ "$CORES" "$out_bam"
  samtools flagstat "$out_bam" > "$flagstat_txt"

  build_norm_args
  log "running bamCoverage for $sample: normalization=$NORMALIZATION binsize=$BINSIZE out=$out_bw"
  bamCoverage -b "$out_bam" -o "$out_bw" --binSize "$BINSIZE" -p "$CORES" "${NORM_ARGS[@]}"
  log "bamCoverage finished for $sample"

  log "writing stepwise alignment counts for $sample"
  {
    echo -e "step\talignments"
    echo -e "raw_bam\t$(samtools view -c "$raw_bam")"
    echo -e "sorted_bam\t$(samtools view -c "$sorted_bam")"
    echo -e "post_markdup_step\t$(samtools view -c "$work_bam")"
    echo -e "post_mapq_flag_filter\t$(samtools view -c "${work_prefix}.mapq.bam")"
    echo -e "post_mito_step\t$(samtools view -c "${work_prefix}.nomito.bam")"
    echo -e "final_te_bam\t$(samtools view -c "$out_bam")"
  } > "$counts_tsv"

  if [[ "$KEEP_INTERMEDIATE" != "true" ]]; then
    rm -f "$raw_bam" "${raw_bam}.bai" \
      "${work_prefix}.markdup.bam" "${work_prefix}.markdup.bam.bai" \
      "${work_prefix}.mapq.bam" "${work_prefix}.mapq.bam.bai" \
      "${work_prefix}.nomito.bam" "${work_prefix}.nomito.bam.bai" \
      "${work_prefix}.clean.unsorted.bam" "${work_prefix}.keep_chroms.txt"
  fi
  log "done $sample: sorted=$sorted_bam ; te_bam=$out_bam ; te_bw=$out_bw"
}

tail -n +2 "$SAMPLE_TABLE" | while IFS=$'\t' read -r sample layout r1 r2; do
  [[ -n "$sample" ]] || continue
  process_one_sample "$sample" "${layout^^}" "$r1" "$r2"
done

log "all FASTQ -> TE/L1 BAM/bigWig jobs finished: $OUTDIR"
