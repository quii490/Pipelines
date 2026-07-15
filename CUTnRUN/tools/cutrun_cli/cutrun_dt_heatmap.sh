#!/usr/bin/env bash
set -euo pipefail

REGIONS=""
SIGNALS=""
LABELS=""
OUT_PREFIX=""
MODE="reference-point"
REFERENCE_POINT="center"
BEFORE=3000
AFTER=3000
BODY_LENGTH=6000
BIN_SIZE=50
MAX_REGIONS="${MAX_HEATMAP_REGIONS:-0}"
MIN_REGION_LENGTH=0
THREADS=8
REGION_LABEL="Regions"
COLOR_MAP="viridis"
X_AXIS_LABEL=""
Y_AXIS_LABEL="Signal"
SORT_REGIONS="descend"
BLACKLIST=""
SKIP_ZEROS=false
MISSING_DATA_AS_ZERO=false
WRITE_MATRIX_TSV=false

usage() {
  cat <<'USAGE'
Usage:
  cutrun_dt_heatmap.sh --regions BED --signals "A.bw B.bw" \
    --labels "A B" --out-prefix PREFIX [options]

Matrix options:
  --mode reference-point|scale-regions  default: reference-point
  --reference-point TSS|TES|center      default: center
  --before INT --after INT              default: 3000, 3000
  --body-length INT                     default: 6000
  --bin-size INT                        default: 50
  --min-region-length INT               discard shorter regions; default: 0
  --max-regions INT                     0 means all regions; default: 0
  --blacklist BED                       optional deepTools blacklist
  --skip-zeros                          opt in to outcome-dependent zero filtering
  --missing-data-as-zero                opt in to replacing missing values with 0
  --sort-regions keep|ascend|descend    default: descend

Plot/output options:
  --region-label TEXT                   default: Regions
  --color-map NAME                      default: viridis
  --x-axis-label TEXT                   heatmap x-axis label
  --y-axis-label TEXT                   default: Signal
  --write-matrix-tsv                    also write the uncompressed numeric matrix
  --threads INT                         default: 8
USAGE
}

die() { echo "ERROR: $*" >&2; exit 1; }
is_uint() { [[ "$1" =~ ^[0-9]+$ ]]; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --regions|-R) REGIONS="$2"; shift 2 ;;
    --signals|-S) SIGNALS="$2"; shift 2 ;;
    --labels) LABELS="$2"; shift 2 ;;
    --out-prefix|-o) OUT_PREFIX="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --reference-point) REFERENCE_POINT="$2"; shift 2 ;;
    --before) BEFORE="$2"; shift 2 ;;
    --after) AFTER="$2"; shift 2 ;;
    --body-length) BODY_LENGTH="$2"; shift 2 ;;
    --bin-size) BIN_SIZE="$2"; shift 2 ;;
    --min-region-length) MIN_REGION_LENGTH="$2"; shift 2 ;;
    --max-regions) MAX_REGIONS="$2"; shift 2 ;;
    --blacklist) BLACKLIST="$2"; shift 2 ;;
    --skip-zeros) SKIP_ZEROS=true; shift ;;
    --missing-data-as-zero) MISSING_DATA_AS_ZERO=true; shift ;;
    --sort-regions) SORT_REGIONS="$2"; shift 2 ;;
    --region-label) REGION_LABEL="$2"; shift 2 ;;
    --color-map) COLOR_MAP="$2"; shift 2 ;;
    --x-axis-label) X_AXIS_LABEL="$2"; shift 2 ;;
    --y-axis-label) Y_AXIS_LABEL="$2"; shift 2 ;;
    --write-matrix-tsv) WRITE_MATRIX_TSV=true; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) die "unknown option $1" ;;
  esac
done

[[ -n "$REGIONS" && -s "$REGIONS" ]] || die "regions missing or empty: $REGIONS"
[[ -n "$SIGNALS" && -n "$LABELS" && -n "$OUT_PREFIX" ]] || die "signals, labels, and out-prefix are required"
[[ "$MODE" =~ ^(reference-point|scale-regions)$ ]] || die "invalid --mode: $MODE"
[[ "$REFERENCE_POINT" =~ ^(TSS|TES|center)$ ]] || die "invalid --reference-point: $REFERENCE_POINT"
[[ "$SORT_REGIONS" =~ ^(keep|ascend|descend)$ ]] || die "invalid --sort-regions: $SORT_REGIONS"
for value in "$BEFORE" "$AFTER" "$BODY_LENGTH" "$BIN_SIZE" "$MIN_REGION_LENGTH" "$MAX_REGIONS" "$THREADS"; do
  is_uint "$value" || die "numeric options must be non-negative integers: $value"
done
(( BIN_SIZE > 0 && THREADS > 0 )) || die "--bin-size and --threads must be greater than zero"
[[ -z "$BLACKLIST" || -s "$BLACKLIST" ]] || die "blacklist missing or empty: $BLACKLIST"
command -v computeMatrix >/dev/null 2>&1 || die "deepTools computeMatrix not found"
command -v plotHeatmap >/dev/null 2>&1 || die "deepTools plotHeatmap not found"
command -v plotProfile >/dev/null 2>&1 || die "deepTools plotProfile not found"

read -r -a signal_array <<< "$SIGNALS"
read -r -a label_array <<< "$LABELS"
(( ${#signal_array[@]} > 0 )) || die "no signal files"
(( ${#signal_array[@]} == ${#label_array[@]} )) || \
  die "signal count (${#signal_array[@]}) differs from label count (${#label_array[@]})"
for signal in "${signal_array[@]}"; do
  [[ -s "$signal" ]] || die "signal missing or empty: $signal"
done

mkdir -p "$(dirname "$OUT_PREFIX")"
valid_regions="${OUT_PREFIX}.valid_regions.bed"
display_regions="${OUT_PREFIX}.display_regions.bed"

# Normalize to BED6, reject malformed coordinates, and optionally remove regions
# that are too short for an interpretable scale-regions metaprofile.
awk -v min_len="$MIN_REGION_LENGTH" 'BEGIN{FS=OFS="\t"}
  !/^#/ && NF>=3 && $2~/^[0-9]+$/ && $3~/^[0-9]+$/ && $2>=0 && $3>$2 && ($3-$2)>=min_len {
    name=(NF>=4 ? $4 : "."); score=(NF>=5 ? $5 : 0); strand=(NF>=6 ? $6 : ".");
    print $1,$2,$3,name,score,strand
  }' "$REGIONS" > "$valid_regions"
[[ -s "$valid_regions" ]] || die "no valid regions remain after BED validation/filtering"

input_regions="$(awk '!/^#/ && NF>=3 {n++} END{print n+0}' "$REGIONS")"
valid_count="$(wc -l < "$valid_regions")"
bad_strands="$(awk '$6!="+" && $6!="-" {n++} END{print n+0}' "$valid_regions")"
if [[ "$MODE" == "scale-regions" && "$bad_strands" -gt 0 ]]; then
  die "scale-regions requires BED6 strand (+/-); invalid or missing strand rows: $bad_strands"
fi

if (( MAX_REGIONS > 0 && valid_count > MAX_REGIONS )); then
  # Deterministic sampling prevents chromosome/input-order bias.
  awk 'BEGIN{srand(1)} {print rand(),$0}' "$valid_regions" \
    | sort -k1,1n \
    | awk -v n="$MAX_REGIONS" 'NR<=n {sub(/^[^ ]+ /, ""); print}' \
    | sort -k1,1 -k2,2n > "$display_regions"
else
  cp -f "$valid_regions" "$display_regions"
fi
display_count="$(wc -l < "$display_regions")"

matrix="${OUT_PREFIX}.matrix.gz"
sorted_regions="${OUT_PREFIX}.matrix_regions.bed"
command=(computeMatrix)
if [[ "$MODE" == "scale-regions" ]]; then
  command+=(scale-regions -S "${signal_array[@]}" -R "$display_regions" \
    -b "$BEFORE" -a "$AFTER" -m "$BODY_LENGTH")
else
  command+=(reference-point -S "${signal_array[@]}" -R "$display_regions" \
    --referencePoint "$REFERENCE_POINT" -b "$BEFORE" -a "$AFTER")
fi
command+=(--binSize "$BIN_SIZE" --samplesLabel "${label_array[@]}" \
  --sortRegions "$SORT_REGIONS" --sortUsing mean \
  --outFileSortedRegions "$sorted_regions" -p "$THREADS" -o "$matrix")
[[ -n "$BLACKLIST" ]] && command+=(--blackListFileName "$BLACKLIST")
[[ "$SKIP_ZEROS" == true ]] && command+=(--skipZeros)
[[ "$MISSING_DATA_AS_ZERO" == true ]] && command+=(--missingDataAsZero)
[[ "$WRITE_MATRIX_TSV" == true ]] && command+=(--outFileNameMatrix "${OUT_PREFIX}.matrix.tsv")

printf '%q ' "${command[@]}" > "${OUT_PREFIX}.COMMANDS.sh"
printf '\n' >> "${OUT_PREFIX}.COMMANDS.sh"
"${command[@]}" > "${OUT_PREFIX}.computeMatrix.log" 2>&1

plot_labels=()
if [[ "$MODE" == "scale-regions" ]]; then
  plot_labels=(--startLabel TSS --endLabel TES)
  [[ -n "$X_AXIS_LABEL" ]] || X_AXIS_LABEL="Relative position (scaled region body)"
else
  plot_labels=(--refPointLabel "$REFERENCE_POINT")
  [[ -n "$X_AXIS_LABEL" ]] || X_AXIS_LABEL="Distance from $REFERENCE_POINT (bp)"
fi

plotHeatmap -m "$matrix" -out "${OUT_PREFIX}.heatmap.pdf" \
  --colorMap "$COLOR_MAP" --whatToShow 'plot, heatmap and colorbar' \
  --samplesLabel "${label_array[@]}" --regionsLabel "$REGION_LABEL" \
  --sortRegions keep --heatmapHeight 12 --heatmapWidth 10 \
  --xAxisLabel "$X_AXIS_LABEL" --yAxisLabel "$Y_AXIS_LABEL" \
  "${plot_labels[@]}" > "${OUT_PREFIX}.plotHeatmap.log" 2>&1

plotProfile -m "$matrix" -out "${OUT_PREFIX}.profile.pdf" \
  --samplesLabel "${label_array[@]}" --regionsLabel "$REGION_LABEL" \
  --plotHeight 5 --plotWidth 8 --yAxisLabel "$Y_AXIS_LABEL" "${plot_labels[@]}" \
  > "${OUT_PREFIX}.plotProfile.log" 2>&1

matrix_count="$(awk '!/^#/ && NF>=3 {n++} END{print n+0}' "$sorted_regions")"
{
  printf 'parameter\tvalue\n'
  printf 'regions_input\t%s\n' "$REGIONS"
  printf 'regions_input_count\t%s\n' "$input_regions"
  printf 'regions_valid_count\t%s\n' "$valid_count"
  printf 'regions_displayed_count\t%s\n' "$display_count"
  printf 'regions_in_matrix_count\t%s\n' "$matrix_count"
  printf 'mode\t%s\n' "$MODE"
  printf 'reference_point\t%s\n' "$REFERENCE_POINT"
  printf 'before\t%s\n' "$BEFORE"
  printf 'after\t%s\n' "$AFTER"
  printf 'body_length\t%s\n' "$BODY_LENGTH"
  printf 'bin_size\t%s\n' "$BIN_SIZE"
  printf 'min_region_length\t%s\n' "$MIN_REGION_LENGTH"
  printf 'max_regions\t%s\n' "$MAX_REGIONS"
  printf 'skip_zeros\t%s\n' "$SKIP_ZEROS"
  printf 'missing_data_as_zero\t%s\n' "$MISSING_DATA_AS_ZERO"
  printf 'sort_regions\t%s\n' "$SORT_REGIONS"
  printf 'blacklist\t%s\n' "${BLACKLIST:-none}"
  printf 'region_label\t%s\n' "$REGION_LABEL"
  printf 'color_map\t%s\n' "$COLOR_MAP"
  printf 'x_axis_label\t%s\n' "$X_AXIS_LABEL"
  printf 'y_axis_label\t%s\n' "$Y_AXIS_LABEL"
} > "${OUT_PREFIX}.parameters.tsv"

echo "DONE: matrix=$matrix regions=$matrix_count heatmap=${OUT_PREFIX}.heatmap.pdf profile=${OUT_PREFIX}.profile.pdf"
