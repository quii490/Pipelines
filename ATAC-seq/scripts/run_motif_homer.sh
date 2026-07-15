#!/usr/bin/env bash
set -euo pipefail
log(){ echo "[run_motif_homer] $*"; }
BACKGROUND="false"
BED_DIR=""
GENOME=""
OUTDIR="./motif"
SIZE="200"
RUN_ANNOTATION="true"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed-dir) BED_DIR="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --size) SIZE="$2"; shift 2 ;;
    --run-annotation) RUN_ANNOTATION="$2"; shift 2 ;;
    -h|--help)
      cat <<'USAGE'
Usage:
  run_motif_homer.sh --bed-dir DIR --genome hg38 [options]

Input BED layouts:
  Organized: DIR/<contrast>.up.bed and DIR/<contrast>.down.bed
  Legacy:    DIR/<contrast>/up.bed and DIR/<contrast>/down.bed

Options:
  --outdir DIR
  --size STR              HOMER motif window. Default: 200 bp around peak center
  --run-annotation BOOL   Also run annotatePeaks.pl. Default: true
USAGE
      exit 0
      ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done
[[ -d "$BED_DIR" ]] || { echo "--bed-dir not found: $BED_DIR" >&2; exit 1; }
[[ -n "$GENOME" ]] || { echo "--genome is required" >&2; exit 1; }
command -v findMotifsGenome.pl >/dev/null 2>&1 || { echo 'HOMER findMotifsGenome.pl not found' >&2; exit 1; }
if [[ "$RUN_ANNOTATION" == "true" ]]; then
  command -v annotatePeaks.pl >/dev/null 2>&1 || { echo 'HOMER annotatePeaks.pl not found' >&2; exit 1; }
fi
mkdir -p "$OUTDIR"
shopt -s nullglob
beds=( "$BED_DIR"/*.bed "$BED_DIR"/*/*.bed )
[[ ${#beds[@]} -gt 0 ]] || { echo "No bed files found under: $BED_DIR" >&2; exit 1; }
summary="$OUTDIR/motif_summary.tsv"
echo -e "contrast\tdirection\tmotif_name\tlog_pvalue\ttarget_pct\tbg_pct" > "$summary"
failures=0
for bed in "${beds[@]}"; do
  [[ -s "$bed" ]] || { log "skip empty BED: $bed"; continue; }
  parent="$(dirname "$bed")"
  stem="$(basename "$bed" .bed)"
  if [[ "$parent" == "$BED_DIR" ]]; then
    direction="${stem##*.}"
    contrast="${stem%.$direction}"
  else
    contrast="$(basename "$parent")"
    direction="$stem"
  fi
  subdir="$OUTDIR/$contrast/$direction"
  mkdir -p "$subdir"
  log "run $contrast / $direction"
  if [[ "$RUN_ANNOTATION" == "true" ]]; then
    if ! annotatePeaks.pl "$bed" "$GENOME" > "$subdir/annotated_peaks.tsv" 2> "$subdir/annotatePeaks.log"; then
      log "WARN: annotatePeaks failed for $contrast / $direction"
      failures=$((failures + 1))
    fi
  fi
  if ! findMotifsGenome.pl "$bed" "$GENOME" "$subdir" -size "$SIZE" -mask > "$subdir/homer.log" 2>&1; then
    log "WARN: motif enrichment failed for $contrast / $direction"
    failures=$((failures + 1))
    continue
  fi
  if [[ -f "$subdir/knownResults.txt" ]]; then
    awk -F'\t' -v c="$contrast" -v d="$direction" 'NR>1 && NR<=21 {print c"\t"d"\t"$1"\t"$3"\t"$8"\t"$9}' "$subdir/knownResults.txt" >> "$summary"
  fi
done
if [[ "$failures" -gt 0 ]]; then
  log "completed with $failures failed HOMER jobs; inspect per-contrast logs"
  exit 1
fi
log "done"
