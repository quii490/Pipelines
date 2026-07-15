#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_bw_correlation.sh --bw-dir DIR --out-prefix out/name [options]
  run_bw_correlation.sh --bw-glob "04_bw/*.bw" --out-prefix out/name [options]

Purpose:
  Run deepTools multiBigwigSummary, correlation heatmap/scatterplot and PCA.

Options:
  --bw-dir DIR                   Directory containing .bw/.bigWig files.
  --bw-glob STR                  Quoted glob for bigWig files.
  --out-prefix PREFIX            Output prefix.
  --bin-size INT                 Default: 10000.
  --cor-method STR               pearson | spearman. Default: pearson.
  --threads INT                  Default: SLURM_CPUS_PER_TASK or 8.
USAGE
}

bw_dir=""
bw_glob=""
outprefix=""
bin_size=10000
cor_method="pearson"
threads="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bw-dir) bw_dir="$2"; shift 2 ;;
    --bw-glob) bw_glob="$2"; shift 2 ;;
    --out-prefix|--outprefix) outprefix="$2"; shift 2 ;;
    --bin-size|--binsize) bin_size="$2"; shift 2 ;;
    --cor-method) cor_method="$2"; shift 2 ;;
    --threads|--cores) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$outprefix" ]] || { usage >&2; exit 1; }
command -v multiBigwigSummary >/dev/null || { echo "ERROR: multiBigwigSummary not found in PATH" >&2; exit 1; }
command -v plotCorrelation >/dev/null || { echo "ERROR: plotCorrelation not found in PATH" >&2; exit 1; }
command -v plotPCA >/dev/null || { echo "ERROR: plotPCA not found in PATH" >&2; exit 1; }

shopt -s nullglob
if [[ -n "$bw_glob" ]]; then
  bw_files=( $bw_glob )
else
  [[ -d "$bw_dir" ]] || { echo "ERROR: bw dir not found: $bw_dir" >&2; exit 1; }
  bw_files=( "$bw_dir"/*.bw "$bw_dir"/*.bigWig "$bw_dir"/*.bigwig )
fi
[[ ${#bw_files[@]} -gt 1 ]] || { echo "ERROR: at least two bigWig files are required" >&2; exit 1; }

base_dir="$(dirname "$outprefix")"
matrix_dir="$base_dir/matrix"
plot_dir="$base_dir/plots"
mkdir -p "$matrix_dir" "$plot_dir"
name="$(basename "$outprefix")"

npz="$matrix_dir/${name}.bins.npz"
raw="$matrix_dir/${name}.raw_counts.tsv"
cor_tsv="$matrix_dir/${name}.${cor_method}_correlation.tsv"
pca_tsv="$matrix_dir/${name}.pca.tsv"

multiBigwigSummary bins \
  --bwfiles "${bw_files[@]}" \
  --binSize "$bin_size" \
  --numberOfProcessors "$threads" \
  --outFileName "$npz" \
  --outRawCounts "$raw"

for ext in pdf png; do
  plotCorrelation \
    --corData "$npz" \
    --corMethod "$cor_method" \
    --whatToPlot heatmap \
    --skipZeros \
    --plotNumbers \
    --removeOutliers \
    --outFileCorMatrix "$cor_tsv" \
    --plotFile "$plot_dir/${name}.${cor_method}_heatmap.${ext}"

  plotCorrelation \
    --corData "$npz" \
    --corMethod "$cor_method" \
    --whatToPlot scatterplot \
    --skipZeros \
    --removeOutliers \
    --plotFile "$plot_dir/${name}.${cor_method}_scatter.${ext}"

  plotPCA \
    --corData "$npz" \
    --plotFile "$plot_dir/${name}.pca.${ext}" \
    --outFileNameData "$pca_tsv"
done

printf "Finished bigWig correlation: %s\n" "$base_dir"
