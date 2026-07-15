#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAM=""
BAM_DIR=""
GTF=""
REF_FLAT=""
OUTDIR="qc_metrics"
STRAND="unstranded"
RIBOSOMAL_INTERVALS=""
THREADS=4
RUN_MARKDUP=true
RUN_RNASEQ_METRICS=true
PICARD_CMD=""
SAMTOOLS_CMD=""

usage() {
  cat <<'USAGE'
从已有 BAM 补做 RNA-seq QC metrics，不修改 BAM，也不改变 counts。

用法：
  rnaseq-qc-metrics-bam --bam FILE --gtf genes.gtf --outdir qc_metrics
  rnaseq-qc-metrics-bam --bam-dir DIR --gtf genes.gtf --strand reverse --outdir qc_metrics

参数：
  --bam FILE                   单个坐标排序 BAM
  --bam-dir DIR                批量扫描目录内的 *.bam
  --gtf FILE                   gene GTF；未给 --ref-flat 时用于自动生成
  --ref-flat FILE              已有 Picard refFlat，可代替 --gtf
  --strand STR                 unstranded|forward|reverse，默认 unstranded
  --ribosomal-intervals FILE   可选 Picard interval_list
  --threads INT                每个样本 Picard 使用的线程提示，默认 4
  --run-markdup BOOL           只统计重复率，默认 true
  --run-rnaseq-metrics BOOL    统计 exonic/intronic/intergenic 和 5'/3' bias，默认 true
  --picard PATH                Picard 可执行文件；默认自动查找
  --outdir DIR                 输出目录，默认 qc_metrics
  -h, --help                   显示帮助
USAGE
}

is_true() {
  case "${1,,}" in true|t|1|yes|y) return 0 ;; false|f|0|no|n) return 1 ;; *) return 2 ;; esac
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --ref-flat|--refFlat) REF_FLAT="$2"; shift 2 ;;
    --strand) STRAND="$2"; shift 2 ;;
    --ribosomal-intervals) RIBOSOMAL_INTERVALS="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --run-markdup) RUN_MARKDUP="$2"; shift 2 ;;
    --run-rnaseq-metrics) RUN_RNASEQ_METRICS="$2"; shift 2 ;;
    --picard) PICARD_CMD="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[rnaseq-qc-metrics-bam] 未知参数: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$BAM" || -n "$BAM_DIR" ]] || { echo "必须提供 --bam 或 --bam-dir" >&2; exit 2; }
[[ -z "$BAM" || -z "$BAM_DIR" ]] || { echo "--bam 与 --bam-dir 只能选择一个" >&2; exit 2; }
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "--threads 必须是正整数" >&2; exit 2; }
is_true "$RUN_MARKDUP" || [[ $? -eq 1 ]] || { echo "--run-markdup 必须是 true/false" >&2; exit 2; }
is_true "$RUN_RNASEQ_METRICS" || [[ $? -eq 1 ]] || { echo "--run-rnaseq-metrics 必须是 true/false" >&2; exit 2; }

case "$STRAND" in
  unstranded|none|NONE) METRIC_STRAND="NONE" ;;
  forward) METRIC_STRAND="FIRST_READ_TRANSCRIPTION_STRAND" ;;
  reverse) METRIC_STRAND="SECOND_READ_TRANSCRIPTION_STRAND" ;;
  *) echo "--strand 仅支持 unstranded、forward、reverse" >&2; exit 2 ;;
esac

if [[ -z "$PICARD_CMD" ]]; then
  PICARD_CMD="$(command -v picard || true)"
  [[ -n "$PICARD_CMD" ]] || PICARD_CMD="/path/to/.conda/envs/rnaseq/bin/picard"
fi
[[ -x "$PICARD_CMD" ]] || { echo "找不到 Picard: $PICARD_CMD" >&2; exit 127; }
SAMTOOLS_CMD="$(command -v samtools || true)"
[[ -n "$SAMTOOLS_CMD" ]] || SAMTOOLS_CMD="/path/to/.conda/envs/rnaseq/bin/samtools"
[[ -x "$SAMTOOLS_CMD" ]] || { echo "找不到 samtools: $SAMTOOLS_CMD" >&2; exit 127; }

mkdir -p "$OUTDIR/markdup" "$OUTDIR/rnaseq_metrics" "$OUTDIR/logs" "$OUTDIR/reference" "$OUTDIR/tmp"

if is_true "$RUN_RNASEQ_METRICS"; then
  if [[ -z "$REF_FLAT" ]]; then
    [[ -s "$GTF" ]] || { echo "RNA-seq metrics 需要 --gtf 或 --ref-flat" >&2; exit 2; }
    REF_FLAT="$OUTDIR/reference/generated.refFlat.txt"
    if [[ ! -s "$REF_FLAT" || "$GTF" -nt "$REF_FLAT" ]]; then
      python3 "$SCRIPT_DIR/../rnaseq/scripts/gtf_to_refflat.py" --input "$GTF" --output "$REF_FLAT"
    fi
  fi
  [[ -s "$REF_FLAT" ]] || { echo "refFlat 不存在或为空: $REF_FLAT" >&2; exit 2; }
fi

bams=()
if [[ -n "$BAM" ]]; then
  [[ -s "$BAM" ]] || { echo "BAM 不存在或为空: $BAM" >&2; exit 2; }
  bams+=("$BAM")
else
  while IFS= read -r -d '' path; do bams+=("$path"); done < <(find "$BAM_DIR" -maxdepth 1 -type f -name '*.bam' -print0 | sort -z)
  (( ${#bams[@]} > 0 )) || { echo "目录中未找到 BAM: $BAM_DIR" >&2; exit 2; }
fi

rrna_args=()
[[ -n "$RIBOSOMAL_INTERVALS" ]] && rrna_args+=("RIBOSOMAL_INTERVALS=$RIBOSOMAL_INTERVALS")

for bam_path in "${bams[@]}"; do
  sample="$(basename "$bam_path")"
  sample="${sample%.bam}"
  log="$OUTDIR/logs/${sample}.qc_metrics.log"
  echo "[qc_metrics] sample=$sample bam=$bam_path" | tee "$log"

  if is_true "$RUN_MARKDUP"; then
    picard_input="$OUTDIR/tmp/${sample}.picard_input.bam"
    bash "$SCRIPT_DIR/../rnaseq/scripts/prepare_picard_bam.sh" \
      --bam "$bam_path" \
      --output "$picard_input" \
      --sample "$sample" \
      --library RNAseq_gene \
      --picard "$PICARD_CMD" \
      --samtools "$SAMTOOLS_CMD" \
      2>&1 | tee -a "$log"
    "$PICARD_CMD" MarkDuplicates \
      "I=$picard_input" \
      "O=$OUTDIR/tmp/${sample}.markdup_qc.tmp.bam" \
      "M=$OUTDIR/markdup/${sample}.markdup.metrics.txt" \
      REMOVE_DUPLICATES=false ASSUME_SORTED=true CREATE_INDEX=false \
      VALIDATION_STRINGENCY=LENIENT "TMP_DIR=$OUTDIR/tmp" \
      2>&1 | tee -a "$log"
    rm -f "$OUTDIR/tmp/${sample}.markdup_qc.tmp.bam" "$picard_input"
    [[ -s "$OUTDIR/markdup/${sample}.markdup.metrics.txt" ]]
  fi

  if is_true "$RUN_RNASEQ_METRICS"; then
    "$PICARD_CMD" CollectRnaSeqMetrics \
      "I=$bam_path" \
      "O=$OUTDIR/rnaseq_metrics/${sample}.rnaseq.metrics.txt" \
      "REF_FLAT=$REF_FLAT" \
      "STRAND_SPECIFICITY=$METRIC_STRAND" \
      "${rrna_args[@]}" \
      VALIDATION_STRINGENCY=LENIENT "TMP_DIR=$OUTDIR/tmp" \
      2>&1 | tee -a "$log"
    [[ -s "$OUTDIR/rnaseq_metrics/${sample}.rnaseq.metrics.txt" ]]
  fi
done

echo "[qc_metrics] 完成，样本数=${#bams[@]}，输出=$OUTDIR"
