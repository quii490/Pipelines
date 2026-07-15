#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR=""
MANIFEST=""
SPECIES="hg38"

usage() {
  cat <<'USAGE'
Usage: cutrun_call_peaks_published.sh --results-dir DIR --manifest manifest.csv [--species hg38|mm39]

Recovery/helper mode: call MACS3 narrow and broad peaks from already published
04_clean_bam/*.bam files.  This is useful after an interrupted Nextflow run;
normal new runs should call CALL_PEAKS inside main.nf.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "$RESULTS_DIR/04_clean_bam" ]] || { echo "ERROR: 04_clean_bam not found" >&2; exit 1; }
[[ -s "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
command -v macs3 >/dev/null || { echo "ERROR: macs3 not found in PATH" >&2; exit 1; }

case "$SPECIES" in
  hg38) MACS_GENOME="hs" ;;
  mm39) MACS_GENOME="mm" ;;
  *) echo "ERROR: unsupported species: $SPECIES" >&2; exit 1 ;;
esac
MACS_CUTOFF_ANALYSIS="${MACS_CUTOFF_ANALYSIS:-false}"
cutoff_arg=()
[[ "$MACS_CUTOFF_ANALYSIS" == true ]] && cutoff_arg+=(--cutoff-analysis)

PEAK_ROOT="$RESULTS_DIR/06_peaks"
mkdir -p "$PEAK_ROOT"
status=0

while IFS=$'\t' read -r sample control layout; do
  [[ -n "$sample" ]] || continue
  target_bam="$RESULTS_DIR/04_clean_bam/${sample}_clean.bam"
  [[ -s "$target_bam" ]] || { echo "ERROR: target BAM missing: $target_bam" >&2; status=1; continue; }
  control_arg=()
  if [[ -n "$control" && -s "$RESULTS_DIR/04_clean_bam/${control}_clean.bam" ]]; then
    control_arg=(-c "$RESULTS_DIR/04_clean_bam/${control}_clean.bam")
  fi
  format="BAM"
  [[ "${layout^^}" == "PE" ]] && format="BAMPE"
  sample_dir="$PEAK_ROOT/$sample"
  mkdir -p "$sample_dir/narrow" "$sample_dir/broad"
  echo "[cutrun_call_peaks] sample=$sample control=${control:-NONE} format=$format"

  if ! macs3 callpeak -t "$target_bam" "${control_arg[@]}" -f "$format" \
      -n "$sample" -g "$MACS_GENOME" --keep-dup 1 "${cutoff_arg[@]}" \
      --outdir "$sample_dir/narrow" >"$sample_dir/narrow/${sample}.macs3.log" 2>&1; then
    echo "[cutrun_call_peaks] narrow failed: $sample" >&2
    status=1
  fi
  if ! macs3 callpeak -t "$target_bam" "${control_arg[@]}" -f "$format" \
      -n "$sample" -g "$MACS_GENOME" --keep-dup 1 "${cutoff_arg[@]}" \
      --broad --outdir "$sample_dir/broad" >"$sample_dir/broad/${sample}.macs3.log" 2>&1; then
    echo "[cutrun_call_peaks] broad failed: $sample" >&2
    status=1
  fi
done < <(python - "$MANIFEST" <<'PY'
import csv, sys
with open(sys.argv[1], newline="") as handle:
    for row in csv.DictReader(handle):
        truthy = str(row.get("is_igg", "")).strip().lower() in {"1", "true", "yes", "y"}
        if not truthy:
            print("\t".join((row.get("sample", "").strip(), row.get("igg", "").strip(), row.get("layout", "").strip())))
PY
)

echo "[cutrun_call_peaks] outputs under $PEAK_ROOT"
exit "$status"
