#!/usr/bin/env bash
set -euo pipefail

REGIONS=""
ANCHOR=""
SIGNALS=""
LABELS=""
OUT_PREFIX=""
MODE="both"
BEFORE=3000
AFTER=3000
BODY_LENGTH=6000
THREADS=8
MAX_REGIONS="${TE_LOCUS_MAX_REGIONS:-100000}"
BIN_SIZE="${TE_LOCUS_BIN_SIZE:-25}"
WRITE_VALUES="${TE_HEATMAP_WRITE_VALUES:-false}"

usage() {
  cat <<'USAGE'
Usage:
  cutrun_te_locus_heatmap.sh \
    --regions L1.bed \
    --anchor anchor_te_locus_best_5bp_cpm.bw \
    --signals "anchor.bw comparison.bw" \
    --labels "Anchor Comparison" \
    --out-prefix results/L1

The input BED must contain six columns with a valid strand. Rows are sorted once
by mean anchor signal around the strand-aware TE 5-prime end, then that exact
order is reused for all signal tracks and both matrix views.

Options:
  --mode reference-point|scale-regions|both  default both
  --before INT                              default 3000
  --after INT                               default 3000
  --body-length INT                         default 6000
  --threads INT                             default 8
  --max-regions INT                         deterministic display cap (default 100000; 0=all)
  --bin-size INT                            matrix bin size in bp (default 25)
  --write-values                            additionally export large raw TSV matrices
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --regions) REGIONS="$2"; shift 2 ;;
    --anchor) ANCHOR="$2"; shift 2 ;;
    --signals) SIGNALS="$2"; shift 2 ;;
    --labels) LABELS="$2"; shift 2 ;;
    --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --before) BEFORE="$2"; shift 2 ;;
    --after) AFTER="$2"; shift 2 ;;
    --body-length) BODY_LENGTH="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --max-regions) MAX_REGIONS="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    --write-values) WRITE_VALUES=true; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 1 ;;
  esac
done

[[ -s "$REGIONS" ]] || { echo "ERROR: --regions BED6 is required" >&2; exit 1; }
[[ -s "$ANCHOR" ]] || { echo "ERROR: --anchor bigWig is required" >&2; exit 1; }
[[ -n "$SIGNALS" ]] || { echo "ERROR: --signals is required" >&2; exit 1; }
[[ -n "$LABELS" ]] || { echo "ERROR: --labels is required" >&2; exit 1; }
[[ -n "$OUT_PREFIX" ]] || { echo "ERROR: --out-prefix is required" >&2; exit 1; }
[[ "$MODE" =~ ^(reference-point|scale-regions|both)$ ]] || { echo "ERROR: invalid --mode $MODE" >&2; exit 1; }
[[ "$MAX_REGIONS" =~ ^[0-9]+$ ]] || { echo "ERROR: --max-regions must be a non-negative integer" >&2; exit 1; }
[[ "$BIN_SIZE" =~ ^[0-9]+$ && "$BIN_SIZE" -gt 0 ]] || { echo "ERROR: --bin-size must be a positive integer" >&2; exit 1; }

read -r -a signal_array <<< "$SIGNALS"
read -r -a label_array <<< "$LABELS"
[[ ${#signal_array[@]} -eq ${#label_array[@]} ]] || {
  echo "ERROR: signal and label counts differ" >&2
  exit 1
}
for signal in "${signal_array[@]}"; do
  [[ -s "$signal" ]] || { echo "ERROR: signal not found: $signal" >&2; exit 1; }
done

mkdir -p "$(dirname "$OUT_PREFIX")"
sorted_regions="${OUT_PREFIX}.anchor_sorted.bed"
matrix_regions="$REGIONS"
if (( MAX_REGIONS > 0 )); then
  region_count="$(wc -l < "$REGIONS")"
  if (( region_count > MAX_REGIONS )); then
    stride=$(( (region_count + MAX_REGIONS - 1) / MAX_REGIONS ))
    matrix_regions="${OUT_PREFIX}.display_regions.bed"
    awk -v stride="$stride" 'NR == 1 || (NR - 1) % stride == 0' "$REGIONS" > "$matrix_regions"
    echo "[cutrun_te_locus_heatmap] capping regions: ${region_count} -> $(wc -l < "$matrix_regions") (stride=${stride})"
  fi
fi

computeMatrix reference-point \
  --referencePoint TSS \
  --regionsFileName "$matrix_regions" \
  --scoreFileName "$ANCHOR" \
  --beforeRegionStartLength "$BEFORE" \
  --afterRegionStartLength "$AFTER" \
  --missingDataAsZero \
  --binSize "$BIN_SIZE" \
  --sortRegions descend \
  --sortUsing mean \
  --outFileSortedRegions "$sorted_regions" \
  --numberOfProcessors "$THREADS" \
  --outFileName "${OUT_PREFIX}.anchor_sort.matrix.gz"

if [[ "$MODE" == "reference-point" || "$MODE" == "both" ]]; then
  values_args=()
  [[ "$WRITE_VALUES" == true ]] && values_args+=(--outFileNameMatrix "${OUT_PREFIX}.5prime.values.tsv")
  computeMatrix reference-point \
    --referencePoint TSS \
    --regionsFileName "$sorted_regions" \
    --scoreFileName "${signal_array[@]}" \
    --samplesLabel "${label_array[@]}" \
    --beforeRegionStartLength "$BEFORE" \
    --afterRegionStartLength "$AFTER" \
    --missingDataAsZero \
    --binSize "$BIN_SIZE" \
    --sortRegions keep \
    --numberOfProcessors "$THREADS" \
    "${values_args[@]}" \
    --outFileName "${OUT_PREFIX}.5prime.matrix.gz"

  plotHeatmap \
    --matrixFile "${OUT_PREFIX}.5prime.matrix.gz" \
    --outFileName "${OUT_PREFIX}.5prime.heatmap.pdf" \
    --sortRegions no \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar'
  plotProfile \
    --matrixFile "${OUT_PREFIX}.5prime.matrix.gz" \
    --outFileName "${OUT_PREFIX}.5prime.profile.pdf" \
    --plotType lines \
    --perGroup
fi

if [[ "$MODE" == "scale-regions" || "$MODE" == "both" ]]; then
  values_args=()
  [[ "$WRITE_VALUES" == true ]] && values_args+=(--outFileNameMatrix "${OUT_PREFIX}.scaled_body.values.tsv")
  computeMatrix scale-regions \
    --regionsFileName "$sorted_regions" \
    --scoreFileName "${signal_array[@]}" \
    --samplesLabel "${label_array[@]}" \
    --beforeRegionStartLength "$BEFORE" \
    --afterRegionStartLength "$AFTER" \
    --regionBodyLength "$BODY_LENGTH" \
    --missingDataAsZero \
    --binSize "$BIN_SIZE" \
    --sortRegions keep \
    --numberOfProcessors "$THREADS" \
    "${values_args[@]}" \
    --outFileName "${OUT_PREFIX}.scaled_body.matrix.gz"

  plotHeatmap \
    --matrixFile "${OUT_PREFIX}.scaled_body.matrix.gz" \
    --outFileName "${OUT_PREFIX}.scaled_body.heatmap.pdf" \
    --sortRegions no \
    --colorMap RdBu_r \
    --whatToShow 'heatmap and colorbar'
  plotProfile \
    --matrixFile "${OUT_PREFIX}.scaled_body.matrix.gz" \
    --outFileName "${OUT_PREFIX}.scaled_body.profile.pdf" \
    --plotType lines \
    --perGroup
fi

printf 'regions\t%s\nregions_used\t%s\nmax_regions\t%s\nbin_size\t%s\nwrite_values\t%s\nanchor\t%s\nmode\t%s\nbefore\t%s\nafter\t%s\nbody_length\t%s\n' \
  "$REGIONS" "$matrix_regions" "$MAX_REGIONS" "$BIN_SIZE" "$WRITE_VALUES" "$ANCHOR" "$MODE" "$BEFORE" "$AFTER" "$BODY_LENGTH" \
  > "${OUT_PREFIX}.parameters.tsv"

echo "[cutrun_te_locus_heatmap] done: ${OUT_PREFIX}"
