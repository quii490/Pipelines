#!/usr/bin/env bash
set -euo pipefail

resolve_script_dir() {
  local src="${BASH_SOURCE[0]}"
  while [ -L "$src" ]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    [[ "$src" != /* ]] && src="$dir/$src"
  done
  cd -P "$(dirname "$src")" && pwd
}

SCRIPT_DIR="$(resolve_script_dir)"

usage() {
  cat <<'USAGE'
Usage:
  run_diff_peak_heatmap.sh --contrast-bed-dir DIR --bw-glob "04_bw/*.bw" --outdir DIR [options]

Purpose:
  Draw ATAC bigWig heatmaps/profiles over differential peak BED files.
  It scans each contrast subdirectory for up/down/sig BED files.

Options:
  --contrast-bed-dir DIR         Usually 08_downstream/peak_level/contrast_beds.
  --bw-glob STR                  BigWig glob, quoted.
  --signals STR                  Alternative space-separated bigWig list.
  --outdir DIR                   Output directory. Default: region_heatmap/differential_peaks.
  --labels STR                   Optional sample labels.
  --before INT                   Default: 3000.
  --after INT                    Default: 3000.
  --threads INT                  Default: SLURM_CPUS_PER_TASK or 8.
USAGE
}

contrast_bed_dir=""
bw_glob=""
signals=""
outdir="region_heatmap/differential_peaks"
labels=""
before=3000
after=3000
threads="${SLURM_CPUS_PER_TASK:-8}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --contrast-bed-dir|--bed-dir) contrast_bed_dir="$2"; shift 2 ;;
    --bw-glob) bw_glob="$2"; shift 2 ;;
    --signals) signals="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --labels) labels="$2"; shift 2 ;;
    --before|-b) before="$2"; shift 2 ;;
    --after|-a) after="$2"; shift 2 ;;
    --threads|--cores) threads="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "$contrast_bed_dir" ]] || { echo "ERROR: contrast bed dir not found: $contrast_bed_dir" >&2; usage >&2; exit 1; }
[[ -n "$bw_glob" || -n "$signals" ]] || { echo "ERROR: set --bw-glob or --signals" >&2; exit 1; }

shopt -s nullglob
bed_files=( "$contrast_bed_dir"/*/*.bed "$contrast_bed_dir"/*.bed )
[[ ${#bed_files[@]} -gt 0 ]] || { echo "ERROR: no BED files found under $contrast_bed_dir" >&2; exit 1; }

for bed in "${bed_files[@]}"; do
  [[ -s "$bed" ]] || continue
  contrast="$(basename "$(dirname "$bed")")"
  [[ "$contrast" == "$(basename "$contrast_bed_dir")" ]] && contrast="all_contrasts"
  stem="$(basename "$bed")"
  stem="${stem%.bed}"
  out_prefix="$outdir/$contrast/$stem"
  cmd=(bash "$SCRIPT_DIR/run_region_heatmap.sh" --regions "$bed" --out-prefix "$out_prefix"
    --before "$before" --after "$after" --threads "$threads")
  [[ -n "$signals" ]] && cmd+=(--signals "$signals")
  [[ -n "$bw_glob" ]] && cmd+=(--bw-glob "$bw_glob")
  [[ -n "$labels" ]] && cmd+=(--labels "$labels")
  "${cmd[@]}"
done

printf "Finished differential peak heatmaps: %s\n" "$outdir"
