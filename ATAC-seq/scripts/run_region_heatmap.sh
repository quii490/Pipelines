#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_region_heatmap.sh --regions regions.bed --signals "A.bw B.bw" --out-prefix out/name [options]
  run_region_heatmap.sh --regions regions.bed --bw-glob "04_bw/*.bw" --outdir out/region_heatmap --name targets [options]

Purpose:
  Draw ATAC accessibility heatmap/profile over any BED regions using deepTools.
  Suitable for target peaks, promoters, enhancers, TE subsets, motif peaks, or differential peaks.

Required:
  --regions FILE                 BED/narrowPeak/broadPeak regions.
  --signals STR                  Space-separated bigWig files, quoted.
  --bw-glob STR                  Alternative to --signals.
  --out-prefix PREFIX            Output prefix.
  --outdir DIR --name STR        Alternative output layout.

Options:
  --labels STR                   Space-separated sample labels.
  --mode STR                     reference-point | scale-regions. Default: reference-point.
  --reference-point STR          center | TSS | TES. Default: center.
  --before INT                   Upstream bp. Default: 3000.
  --after INT                    Downstream bp. Default: 3000.
  --body-length INT              scale-regions body length. Default: 5000.
  --sort-regions STR             descend | ascend | no | keep. Default: descend.
  --sort-using STR               mean | median | max | sum | region_length. Default: mean.
  --kmeans INT                   Optional k-means clusters for plotHeatmap.
  --z-min FLOAT                  Optional heatmap zMin.
  --z-max FLOAT                  Optional heatmap zMax.
  --threads INT                  Default: SLURM_CPUS_PER_TASK or 8.
USAGE
}

regions=""
signals=""
bw_glob=""
labels=""
outprefix=""
outdir=""
name=""
mode="reference-point"
reference_point="center"
before=3000
after=3000
body_length=5000
sort_regions="descend"
sort_using="mean"
kmeans=""
z_min=""
z_max=""
threads="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    -R|--regions) regions="$2"; shift 2 ;;
    -S|--signals) signals="$2"; shift 2 ;;
    --bw-glob) bw_glob="$2"; shift 2 ;;
    --labels|--samples-label) labels="$2"; shift 2 ;;
    -o|--out-prefix|--outprefix) outprefix="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --name) name="$2"; shift 2 ;;
    --mode) mode="$2"; shift 2 ;;
    --reference-point) reference_point="$2"; shift 2 ;;
    --before|-b) before="$2"; shift 2 ;;
    --after|-a) after="$2"; shift 2 ;;
    --body-length) body_length="$2"; shift 2 ;;
    --sort-regions) sort_regions="$2"; shift 2 ;;
    --sort-using) sort_using="$2"; shift 2 ;;
    --kmeans) kmeans="$2"; shift 2 ;;
    --z-min) z_min="$2"; shift 2 ;;
    --z-max) z_max="$2"; shift 2 ;;
    --threads|--cores) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$regions" ]] || { usage >&2; exit 1; }
[[ -s "$regions" ]] || { echo "ERROR: regions file not found: $regions" >&2; exit 1; }
if [[ -z "$outprefix" ]]; then
  [[ -n "$outdir" && -n "$name" ]] || { echo "ERROR: set --out-prefix or both --outdir/--name" >&2; exit 1; }
  outprefix="$outdir/$name"
fi
if [[ -z "$signals" && -n "$bw_glob" ]]; then
  shopt -s nullglob
  bw_files=( $bw_glob )
  [[ ${#bw_files[@]} -gt 0 ]] || { echo "ERROR: no bigWig matched --bw-glob: $bw_glob" >&2; exit 1; }
else
  read -r -a bw_files <<< "$signals"
fi
[[ ${#bw_files[@]} -gt 0 ]] || { echo "ERROR: set --signals or --bw-glob" >&2; exit 1; }

case "$mode" in
  reference-point|scale-regions) ;;
  *) echo "ERROR: --mode must be reference-point or scale-regions" >&2; exit 1 ;;
esac

command -v computeMatrix >/dev/null || { echo "ERROR: computeMatrix not found in PATH" >&2; exit 1; }
command -v plotProfile >/dev/null || { echo "ERROR: plotProfile not found in PATH" >&2; exit 1; }
command -v plotHeatmap >/dev/null || { echo "ERROR: plotHeatmap not found in PATH" >&2; exit 1; }

for bw in "${bw_files[@]}"; do
  [[ -s "$bw" ]] || { echo "ERROR: bigWig not found: $bw" >&2; exit 1; }
done

base_dir="$(dirname "$outprefix")"
matrix_dir="$base_dir/matrix"
plot_dir="$base_dir/plots"
log_dir="$base_dir/logs"
mkdir -p "$matrix_dir" "$plot_dir" "$log_dir"
prefix_name="$(basename "$outprefix")"

matrix="$matrix_dir/${prefix_name}.matrix.gz"
matrix_tab="$matrix_dir/${prefix_name}.matrix.tab"
regions_sorted="$matrix_dir/${prefix_name}.regions.sorted.bed"

label_args=()
if [[ -n "$labels" ]]; then
  read -r -a label_values <<< "$labels"
  label_args=(--samplesLabel "${label_values[@]}")
fi

compute_args=(-R "$regions" -S "${bw_files[@]}" -b "$before" -a "$after"
  --numberOfProcessors "$threads" "${label_args[@]}"
  --sortRegions "$sort_regions" --sortUsing "$sort_using"
  -o "$matrix" --outFileNameMatrix "$matrix_tab" --outFileSortedRegions "$regions_sorted")

if [[ "$mode" == "scale-regions" ]]; then
  computeMatrix scale-regions --regionBodyLength "$body_length" "${compute_args[@]}" 2>&1 | tee "$log_dir/${prefix_name}.computeMatrix.log"
else
  computeMatrix reference-point --referencePoint "$reference_point" "${compute_args[@]}" 2>&1 | tee "$log_dir/${prefix_name}.computeMatrix.log"
fi

heatmap_args=(-m "$matrix" --plotTitle "${prefix_name} heatmap")
profile_args=(-m "$matrix" --plotTitle "${prefix_name} profile")
[[ -n "$kmeans" ]] && heatmap_args+=(--kmeans "$kmeans")
[[ -n "$z_min" ]] && heatmap_args+=(--zMin "$z_min")
[[ -n "$z_max" ]] && heatmap_args+=(--zMax "$z_max")

plotHeatmap "${heatmap_args[@]}" -out "$plot_dir/${prefix_name}.heatmap.pdf"
plotHeatmap "${heatmap_args[@]}" -out "$plot_dir/${prefix_name}.heatmap.png"
plotProfile "${profile_args[@]}" -out "$plot_dir/${prefix_name}.profile.pdf"
plotProfile "${profile_args[@]}" -out "$plot_dir/${prefix_name}.profile.png"

printf "Finished region heatmap: %s\n" "$base_dir"
