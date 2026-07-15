#!/usr/bin/env bash
set -euo pipefail
log(){ echo "[run_tobias] $*"; }
BACKGROUND="false"
BAM_DIR=""
PEAKS=""
SAMPLE_META=""
CONTRAST_FILE=""
GENOME=""
MOTIF_MEME=""
OUTDIR="./tobias"
CORES=4

usage() {
  cat <<'USAGE'
Usage:
  run_tobias.sh --bam-dir DIR --peaks BED --sample-meta CSV \
    --genome FASTA --motif-meme FILE [options]

Required:
  --bam-dir DIR          Directory containing *.clean.bam files.
  --peaks BED            Consensus/accessibility regions for TOBIAS.
  --sample-meta CSV      Columns: sample,condition,replicate.
  --genome FASTA         Reference FASTA matching BAM and peaks.
  --motif-meme FILE      Motifs in MEME format.

Optional:
  --contrast-file CSV    Columns: case,control; enables BINDetect.
  --outdir DIR           Default: ./tobias.
  --cores INT            Default: 4.
  -h, --help             Show this help.

The wrapper runs ATACorrect and ScoreBigwig per sample. When a contrast file is
provided, BINDetect currently chooses the first sample listed for each condition;
review bindetect_sample_choice.tsv before interpreting differential binding.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bam-dir) BAM_DIR="$2"; shift 2 ;;
    --peaks) PEAKS="$2"; shift 2 ;;
    --sample-meta) SAMPLE_META="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --motif-meme) MOTIF_MEME="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -d "$BAM_DIR" ]] || { echo "--bam-dir not found: $BAM_DIR" >&2; exit 1; }
[[ -f "$PEAKS" ]] || { echo "--peaks not found: $PEAKS" >&2; exit 1; }
[[ -f "$SAMPLE_META" ]] || { echo "--sample-meta not found: $SAMPLE_META" >&2; exit 1; }
[[ -f "$GENOME" ]] || { echo "--genome not found: $GENOME" >&2; exit 1; }
[[ -f "$MOTIF_MEME" ]] || { echo "--motif-meme not found: $MOTIF_MEME" >&2; exit 1; }
command -v TOBIAS >/dev/null 2>&1 || { echo 'TOBIAS not found in PATH' >&2; exit 1; }
mkdir -p "$OUTDIR/corrected" "$OUTDIR/footprints" "$OUTDIR/bindetect"
echo -e "case\tcontrol\tcase_sample\tcontrol_sample" > "$OUTDIR/bindetect/bindetect_sample_choice.tsv"
shopt -s nullglob
bam_files=( "$BAM_DIR"/*.clean.bam )
[[ ${#bam_files[@]} -gt 0 ]] || { echo "No *.clean.bam found in: $BAM_DIR" >&2; exit 1; }

for bam in "${bam_files[@]}"; do
  sample=$(basename "$bam" .clean.bam)
  log "ATACorrect $sample"
  TOBIAS ATACorrect --bam "$bam" --genome "$GENOME" --peaks "$PEAKS" --outdir "$OUTDIR/corrected" --prefix "$sample" --cores "$CORES"
  log "ScoreBigwig $sample"
  TOBIAS ScoreBigwig --signal "$OUTDIR/corrected/${sample}_corrected.bw" --regions "$PEAKS" --output "$OUTDIR/footprints/${sample}_footprints.bw" --cores "$CORES"
done

if [[ -f "$CONTRAST_FILE" ]]; then
  tail -n +2 "$CONTRAST_FILE" | while IFS=',' read -r case control; do
    [[ -n "$case" && -n "$control" ]] || continue
    case_sample=$(awk -F',' -v c="$case" 'NR>1 && $2==c {print $1; exit}' "$SAMPLE_META" || true)
    control_sample=$(awk -F',' -v c="$control" 'NR>1 && $2==c {print $1; exit}' "$SAMPLE_META" || true)
    if [[ -z "$case_sample" || -z "$control_sample" ]]; then
      log "skip BINDetect for $case vs $control: cannot identify single representative sample per condition"
      continue
    fi
    echo -e "$case\t$control\t$case_sample\t$control_sample" >> "$OUTDIR/bindetect/bindetect_sample_choice.tsv"
    sig1="$OUTDIR/footprints/${control_sample}_footprints.bw"
    sig2="$OUTDIR/footprints/${case_sample}_footprints.bw"
    [[ -f "$sig1" && -f "$sig2" ]] || { log "skip BINDetect for $case vs $control: footprint bw missing"; continue; }
    cmp_out="$OUTDIR/bindetect/${case}_vs_${control}"
    mkdir -p "$cmp_out"
    log "BINDetect $case vs $control"
    TOBIAS BINDetect --motifs "$MOTIF_MEME" --signals "$sig1" "$sig2" --genome "$GENOME" --peaks "$PEAKS" --outdir "$cmp_out" --cond_names "$control_sample" "$case_sample" --cores "$CORES" || true
  done
fi
log "done"
