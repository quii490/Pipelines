#!/usr/bin/env bash
set -euo pipefail

SCRIPT_NAME="$(basename "$0")"

BAM_DIR=""
BAMS=()
OUTDIR=""
PREFIX="rnaseq"
GENE_GTF=""
TE_GTF=""
GENE_SAF=""
TE_SAF=""
FEATURE_TYPE="exon"
GENE_ID_ATTR="gene_id"
TE_ID_ATTR="gene_id"
STRANDEDNESS="0"
LAYOUT="auto"
THREADS="8"
MIN_MQ="10"
EXTRA_FEATURECOUNTS_ARGS=""
RUN_GENE="true"
RUN_TE="true"

usage() {
  cat <<'USAGE'
RNA-seq BAM -> raw count matrix wrapper

用途：
  对已经比对好的 BAM 重新做 gene / TE featureCounts，并导出可直接给
  rnaseq-counts-to-de 或 rnaseq-plot-counts 使用的 count matrix。

常用示例：
  rnaseq-bam-to-counts \
    --bam-dir /path/to/02_align \
    --gene-gtf /path/gencode.gtf \
    --te-gtf /path/rmsk_TE.gtf \
    --outdir /path/counts \
    --prefix project \
    --layout PE \
    --strandedness 2 \
    --threads 16

输入：
  --bam-dir DIR                  BAM 所在目录，自动读取 *.bam
  --bam FILE                     指定单个 BAM，可重复；与 --bam-dir 可同时使用
  --gene-gtf FILE                gene GTF/GFF
  --te-gtf FILE                  TE GTF/GFF
  --gene-saf FILE                gene SAF；提供时优先于 --gene-gtf
  --te-saf FILE                  TE SAF；提供时优先于 --te-gtf

输出与开关：
  --outdir DIR                   输出目录
  --prefix STR                   输出前缀，默认 rnaseq
  --gene-only                    只跑 gene counts
  --te-only                      只跑 TE counts

featureCounts 参数：
  --layout auto|PE|SE            默认 auto；PE 会加 -p --countReadPairs
  --strandedness 0|1|2           0 unstranded, 1 forward, 2 reverse；默认 0
  --feature-type STR             GTF feature type，默认 exon
  --gene-id-attr STR             gene GTF ID 属性，默认 gene_id
  --te-id-attr STR               TE GTF ID 属性，默认 gene_id
  --min-mapq INT                 -Q，默认 10
  --threads INT                  默认 8
  --extra-featurecounts-args STR 额外原样传给 featureCounts，例如 "--primary -M --fraction"

输出文件：
  <outdir>/<prefix>.gene.featureCounts.txt
  <outdir>/<prefix>.gene.count_matrix.csv
  <outdir>/<prefix>.TE.featureCounts.txt
  <outdir>/<prefix>.TE.count_matrix.csv
  <outdir>/<prefix>.bam_to_counts.summary.txt
USAGE
}

log() {
  echo "[${SCRIPT_NAME}] $*" >&2
}

die() {
  echo "[${SCRIPT_NAME}] ERROR: $*" >&2
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --bam) BAMS+=("$2"); shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --gene-gtf) GENE_GTF="$2"; shift 2 ;;
    --te-gtf) TE_GTF="$2"; shift 2 ;;
    --gene-saf) GENE_SAF="$2"; shift 2 ;;
    --te-saf) TE_SAF="$2"; shift 2 ;;
    --feature-type) FEATURE_TYPE="$2"; shift 2 ;;
    --gene-id-attr) GENE_ID_ATTR="$2"; shift 2 ;;
    --te-id-attr) TE_ID_ATTR="$2"; shift 2 ;;
    --strandedness) STRANDEDNESS="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --min-mapq) MIN_MQ="$2"; shift 2 ;;
    --extra-featurecounts-args) EXTRA_FEATURECOUNTS_ARGS="$2"; shift 2 ;;
    --gene-only) RUN_TE="false"; shift 1 ;;
    --te-only) RUN_GENE="false"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

command -v featureCounts >/dev/null 2>&1 || die "featureCounts not found in PATH"
[[ -n "${OUTDIR}" ]] || die "--outdir is required"
mkdir -p "${OUTDIR}"

if [[ -n "${BAM_DIR}" ]]; then
  [[ -d "${BAM_DIR}" ]] || die "--bam-dir does not exist: ${BAM_DIR}"
  while IFS= read -r bam; do
    BAMS+=("${bam}")
  done < <(find "${BAM_DIR}" -maxdepth 2 -type f -name "*.bam" | sort)
fi

[[ "${#BAMS[@]}" -gt 0 ]] || die "No BAM files found. Use --bam-dir or --bam"

case "${LAYOUT}" in
  auto|PE|SE) ;;
  *) die "--layout must be auto, PE or SE" ;;
esac
case "${STRANDEDNESS}" in
  0|1|2|unstranded|forward|reverse)
    if [[ "${STRANDEDNESS}" == "unstranded" ]]; then STRANDEDNESS="0"; fi
    if [[ "${STRANDEDNESS}" == "forward" ]]; then STRANDEDNESS="1"; fi
    if [[ "${STRANDEDNESS}" == "reverse" ]]; then STRANDEDNESS="2"; fi
    ;;
  *) die "--strandedness must be 0/1/2 or unstranded/forward/reverse" ;;
esac

if [[ "${RUN_GENE}" == "true" ]]; then
  [[ -n "${GENE_GTF}${GENE_SAF}" ]] || die "gene counting requires --gene-gtf or --gene-saf"
fi
if [[ "${RUN_TE}" == "true" ]]; then
  [[ -n "${TE_GTF}${TE_SAF}" ]] || die "TE counting requires --te-gtf or --te-saf"
fi

detect_layout() {
  local first_bam="$1"
  if [[ "${LAYOUT}" != "auto" ]]; then
    echo "${LAYOUT}"
    return
  fi
  if command -v samtools >/dev/null 2>&1; then
    local paired
    paired="$(samtools view -c -f 1 "${first_bam}" 2>/dev/null || echo 0)"
    if [[ "${paired}" =~ ^[0-9]+$ && "${paired}" -gt 0 ]]; then
      echo "PE"
    else
      echo "SE"
    fi
  else
    echo "SE"
  fi
}

write_matrix_from_featurecounts() {
  local fc_file="$1"
  local matrix_file="$2"
  python3 - "$fc_file" "$matrix_file" <<'PY'
import csv
import os
import sys

fc_file, matrix_file = sys.argv[1], sys.argv[2]
with open(fc_file, newline="") as handle:
    rows = [line for line in handle if not line.startswith("#")]
reader = csv.reader(rows, delimiter="\t")
header = next(reader)
if len(header) < 7:
    raise SystemExit("featureCounts output has fewer than 7 columns")
sample_cols = header[6:]
sample_names = []
for x in sample_cols:
    name = os.path.basename(x)
    if name.endswith(".bam"):
        name = name[:-4]
    sample_names.append(name)
with open(matrix_file, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["feature_id"] + sample_names)
    for row in reader:
        if len(row) < 7:
            continue
        writer.writerow([row[0]] + row[6:])
PY
}

run_featurecounts() {
  local tag="$1"
  local annot="$2"
  local annot_format="$3"
  local attr="$4"
  local out_txt="${OUTDIR}/${PREFIX}.${tag}.featureCounts.txt"
  local out_matrix="${OUTDIR}/${PREFIX}.${tag}.count_matrix.csv"
  local fc_layout
  fc_layout="$(detect_layout "${BAMS[0]}")"

  log "run ${tag} featureCounts; layout=${fc_layout}; strandedness=${STRANDEDNESS}"
  local cmd=(featureCounts -T "${THREADS}" -s "${STRANDEDNESS}" -Q "${MIN_MQ}" -o "${out_txt}")
  if [[ "${fc_layout}" == "PE" ]]; then
    cmd+=(-p --countReadPairs)
  fi
  if [[ "${annot_format}" == "SAF" ]]; then
    cmd+=(-F SAF -a "${annot}")
  else
    cmd+=(-t "${FEATURE_TYPE}" -g "${attr}" -a "${annot}")
  fi
  if [[ -n "${EXTRA_FEATURECOUNTS_ARGS}" ]]; then
    read -r -a extra_args <<< "${EXTRA_FEATURECOUNTS_ARGS}"
    cmd+=("${extra_args[@]}")
  fi
  cmd+=("${BAMS[@]}")
  "${cmd[@]}"
  write_matrix_from_featurecounts "${out_txt}" "${out_matrix}"
  log "matrix written: ${out_matrix}"
}

if [[ "${RUN_GENE}" == "true" ]]; then
  if [[ -n "${GENE_SAF}" ]]; then
    run_featurecounts "gene" "${GENE_SAF}" "SAF" "${GENE_ID_ATTR}"
  else
    run_featurecounts "gene" "${GENE_GTF}" "GTF" "${GENE_ID_ATTR}"
  fi
fi

if [[ "${RUN_TE}" == "true" ]]; then
  if [[ -n "${TE_SAF}" ]]; then
    run_featurecounts "TE" "${TE_SAF}" "SAF" "${TE_ID_ATTR}"
  else
    run_featurecounts "TE" "${TE_GTF}" "GTF" "${TE_ID_ATTR}"
  fi
fi

{
  echo "bam_count=${#BAMS[@]}"
  echo "layout=$(detect_layout "${BAMS[0]}")"
  echo "strandedness=${STRANDEDNESS}"
  echo "run_gene=${RUN_GENE}"
  echo "run_te=${RUN_TE}"
  echo "prefix=${PREFIX}"
} > "${OUTDIR}/${PREFIX}.bam_to_counts.summary.txt"

log "DONE"
