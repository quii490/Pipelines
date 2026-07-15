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
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"

usage() {
  cat <<'USAGE'
Usage:
  run_peak_annotation.sh --input peaks.narrowPeak --species hg38 --out-prefix out/sample [options]

Purpose:
  Annotate open chromatin peaks using ChIPseeker gene annotation and optional TE overlap.

Options:
  --input FILE                   BED/narrowPeak/broadPeak.
  --species STR                  hg38 | mm10 | mm39. Default: hg38.
  --out-prefix PREFIX            Output prefix.
  --gtf FILE                     Override GTF. Default from conf/species_refs.config.
  --te-bed FILE                  Override TE BED/GTF. Default from config.
  --promoter-up INT              Default: 3000.
  --promoter-down INT            Default: 3000.
  --txdb-cache FILE              Optional TxDb sqlite cache.
  --no-plots                     Only write tables.
USAGE
}

input=""
species="hg38"
out_prefix=""
gtf=""
te_bed=""
promoter_up=3000
promoter_down=3000
txdb_cache=""
no_plots="false"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input|-i) input="$2"; shift 2 ;;
    --species) species="$2"; shift 2 ;;
    --out-prefix|--outprefix|-o) out_prefix="$2"; shift 2 ;;
    --gtf) gtf="$2"; shift 2 ;;
    --te-bed|--te-anno) te_bed="$2"; shift 2 ;;
    --promoter-up) promoter_up="$2"; shift 2 ;;
    --promoter-down) promoter_down="$2"; shift 2 ;;
    --txdb-cache) txdb_cache="$2"; shift 2 ;;
    --no-plots) no_plots="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -s "$input" && -n "$out_prefix" ]] || { usage >&2; exit 1; }
[[ -f "$CONFIG_FILE" ]] || { echo "ERROR: config not found: $CONFIG_FILE" >&2; exit 1; }
[[ -n "$gtf" ]] || gtf="$(get_species_param "$species" gtf_genes "$CONFIG_FILE" || true)"
[[ -n "$te_bed" ]] || te_bed="$(get_species_param "$species" te_bed "$CONFIG_FILE" || true)"
[[ -s "$gtf" ]] || { echo "ERROR: GTF not found: $gtf" >&2; exit 1; }

cmd=(Rscript "$SCRIPT_DIR/run_peak_annotation_chipseeker.R"
  --input "$input"
  --species "$species"
  --out-prefix "$out_prefix"
  --gtf "$gtf"
  --promoter-up "$promoter_up"
  --promoter-down "$promoter_down")
[[ -n "$te_bed" && -s "$te_bed" ]] && cmd+=(--te-bed "$te_bed")
[[ -n "$txdb_cache" ]] && cmd+=(--txdb-cache "$txdb_cache")
[[ "$no_plots" == "true" ]] && cmd+=(--no-plots)

"${cmd[@]}"
