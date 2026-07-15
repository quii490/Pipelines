#!/usr/bin/env bash
set -euo pipefail

BAM=""
OUTPUT=""
SAMPLE=""
LIBRARY="RNAseq"
PICARD_CMD=""
SAMTOOLS_CMD=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam) BAM="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --sample) SAMPLE="$2"; shift 2 ;;
    --library) LIBRARY="$2"; shift 2 ;;
    --picard) PICARD_CMD="$2"; shift 2 ;;
    --samtools) SAMTOOLS_CMD="$2"; shift 2 ;;
    *) echo "[prepare_picard_bam] unknown argument: $1" >&2; exit 2 ;;
  esac
done

[[ -s "$BAM" ]] || { echo "[prepare_picard_bam] BAM missing or empty: $BAM" >&2; exit 2; }
[[ -n "$OUTPUT" && -n "$SAMPLE" ]] || { echo "[prepare_picard_bam] --output and --sample are required" >&2; exit 2; }

[[ -n "$PICARD_CMD" ]] || PICARD_CMD="$(command -v picard || true)"
[[ -n "$SAMTOOLS_CMD" ]] || SAMTOOLS_CMD="$(command -v samtools || true)"
[[ -x "$PICARD_CMD" ]] || { echo "[prepare_picard_bam] picard not found" >&2; exit 127; }
[[ -x "$SAMTOOLS_CMD" ]] || { echo "[prepare_picard_bam] samtools not found" >&2; exit 127; }

if "$SAMTOOLS_CMD" view -H "$BAM" | grep '^@RG' >/dev/null; then
  ln -sf "$(realpath "$BAM")" "$OUTPUT"
  echo "[prepare_picard_bam] existing read group retained: $SAMPLE"
else
  "$PICARD_CMD" AddOrReplaceReadGroups \
    "I=$BAM" \
    "O=$OUTPUT" \
    "RGID=$SAMPLE" \
    "RGLB=$LIBRARY" \
    RGPL=ILLUMINA \
    "RGPU=$SAMPLE" \
    "RGSM=$SAMPLE" \
    SORT_ORDER=coordinate \
    CREATE_INDEX=false \
    VALIDATION_STRINGENCY=LENIENT
  echo "[prepare_picard_bam] temporary read group added: $SAMPLE"
fi

[[ -s "$OUTPUT" ]]
