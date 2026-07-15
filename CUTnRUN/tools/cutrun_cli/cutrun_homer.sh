#!/usr/bin/env bash
set -euo pipefail

INPUT=""
REF="hg38"
OUT_DIR=""
THREADS=8

usage() {
  cat <<'USAGE'
Usage: cutrun_homer.sh --input PEAK --ref hg38|mm39 --out-dir DIR [--threads N]

Run HOMER annotation and de-novo motif discovery. The wrapper records the
exact commands and writes a SKIP marker when HOMER is not installed.
USAGE
}
while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --ref) REF="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done
[[ -s "$INPUT" ]] || { echo "ERROR: peak file not found: $INPUT" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: --out-dir is required" >&2; exit 1; }
mkdir -p "$OUT_DIR"
case "$REF" in hg38) HOMER_REF=hg38 ;; mm39) HOMER_REF=mm39 ;; *) echo "ERROR: unsupported ref: $REF" >&2; exit 1 ;; esac

FIND_MOTIFS="${FIND_MOTIFS_GENOME:-$(command -v findMotifsGenome.pl || true)}"
ANNOTATE="${ANNOTATE_PEAKS:-$(command -v annotatePeaks.pl || true)}"
if [[ -z "$FIND_MOTIFS" || -z "$ANNOTATE" ]]; then
  printf 'status\tSKIP\nreason\tHOMER findMotifsGenome.pl/annotatePeaks.pl not found\n' > "$OUT_DIR/STATUS.tsv"
  exit 0
fi

BED="$OUT_DIR/input.bed"
awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {print $1,$2,$3,"peak_"NR}' "$INPUT" > "$BED"
printf '%q ' "$FIND_MOTIFS" "$BED" "$HOMER_REF" "$OUT_DIR" -p "$THREADS" > "$OUT_DIR/COMMANDS.sh"
printf '\n%q ' "$ANNOTATE" "$BED" "$HOMER_REF" >> "$OUT_DIR/COMMANDS.sh"
printf '\n' >> "$OUT_DIR/COMMANDS.sh"

"$ANNOTATE" "$BED" "$HOMER_REF" > "$OUT_DIR/annotatePeaks.txt" 2> "$OUT_DIR/annotatePeaks.log"
"$FIND_MOTIFS" "$BED" "$HOMER_REF" "$OUT_DIR" -p "$THREADS" > "$OUT_DIR/findMotifs.stdout.log" 2> "$OUT_DIR/findMotifs.stderr.log"
[[ -s "$OUT_DIR/annotatePeaks.txt" ]] || { echo "HOMER annotation output is empty" >&2; exit 1; }
printf 'status\tPASS\ninput\t%s\nreference\t%s\n' "$INPUT" "$HOMER_REF" > "$OUT_DIR/STATUS.tsv"
