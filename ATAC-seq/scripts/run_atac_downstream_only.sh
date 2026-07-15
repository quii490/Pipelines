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
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG_FILE="$ROOT_DIR/conf/species_refs.config"
ORGANIZER="$SCRIPT_DIR/organize_atac_downstream_outputs.py"
PROFILE_HEATMAP_RUNNER="$SCRIPT_DIR/run_atac_profile_heatmaps.sh"
# shellcheck source=/dev/null
source "$SCRIPT_DIR/species_config_lib.sh"

# Prefer the project conda environment when present, so downstream can find Rscript
# even in non-interactive/background runs.
CHIPSEQ_ENV_BIN="/path/to/.conda/envs/chipseq/bin"
if [[ -d "$CHIPSEQ_ENV_BIN" ]]; then
  export PATH="$CHIPSEQ_ENV_BIN:$PATH"
fi

RESULT_DIR=""
SPECIES="hg38"
CONTRAST_FILE=""
OUTDIR=""
LEVELS="both"
SAMPLESHEET=""
TE_BED=""
GTF_FILE=""
ANNOTATION_MODE="gene_te"
TE_VIOLIN_CLASS=""
TE_VIOLIN_TOP_N="12"
TE_FAMILY_FILTER=""
TE_NAME_FILTER=""
LABEL_TOP_N="40"
BACKGROUND="false"
OUTPUT_LAYOUT="standard"
RUN_PROFILE_HEATMAPS="false"
PROFILE_CORES="8"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_atac_downstream_only.sh --result-dir RESULTS [options]

Purpose:
  Re-run downstream analysis from an existing ATAC result directory.
  Normal users only need to edit:
    RESULTS/_automation/inputs/samplesheet.csv
    RESULTS/_automation/inputs/contrasts.csv
  This script will automatically create the metadata used by downstream R.

Required:
  --result-dir DIR                  Existing ATAC result directory.

Common options:
  --species STR                     hg38 | mm10 | mm39. Default: hg38
  --levels STR                      peak | bin | both. Default: both
  --outdir DIR                      Default: RESULTS/08_downstream
  --samplesheet FILE                Default: RESULTS/_automation/inputs/samplesheet.csv
  --contrast-file FILE              Default: RESULTS/_automation/inputs/contrasts.csv
  --background                      Run in background and write a log under RESULTS/_automation/logs
  --output-layout STR               standard | legacy. Default: standard
  --run-profile-heatmaps true|false Default: false. Profile is a standalone, compute-heavy tool
  --skip-profile-heatmaps           Shortcut for --run-profile-heatmaps false
  --profile-cores INT               Cores for computeMatrix. Default: 8

Annotation options:
  --te-bed FILE                     Override TE annotation. Default from conf/species_refs.config
  --gtf FILE                        Override gene GTF. Default from conf/species_refs.config
  --annotation-mode STR             gene_te | gene | none. Default: gene_te
  --te-violin-class STR
  --te-violin-top-n INT             Default: 12
  --te-family-filter STR
  --te-name-filter STR
  --label-top-n INT                Labels on MA/volcano/scatter plots. Default: 40

Outputs:
  OUTDIR/peak_level                 Peak-level downstream results
  OUTDIR/bin_level                  Fixed-bin downstream results
  OUTDIR/profile_heatmaps           deepTools TSS/gene body/TE/L1 profile heatmaps
  OUTDIR/_inputs/sample_metadata.csv Synced metadata generated from samplesheet
  OUTDIR/<level>/legacy/raw_r_output Raw R output used to build the organized view
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --result-dir) RESULT_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --level|--levels) LEVELS="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --gtf) GTF_FILE="$2"; shift 2 ;;
    --annotation-mode) ANNOTATION_MODE="$2"; shift 2 ;;
    --te-violin-class) TE_VIOLIN_CLASS="$2"; shift 2 ;;
    --te-violin-top-n) TE_VIOLIN_TOP_N="$2"; shift 2 ;;
    --te-family-filter) TE_FAMILY_FILTER="$2"; shift 2 ;;
    --te-name-filter) TE_NAME_FILTER="$2"; shift 2 ;;
    --label-top-n) LABEL_TOP_N="$2"; shift 2 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --output-layout) OUTPUT_LAYOUT="$2"; shift 2 ;;
    --run-profile-heatmaps) RUN_PROFILE_HEATMAPS="$2"; shift 2 ;;
    --skip-profile-heatmaps) RUN_PROFILE_HEATMAPS="false"; shift 1 ;;
    --profile-cores) PROFILE_CORES="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -d "$RESULT_DIR" ]] || { echo "--result-dir not found" >&2; usage >&2; exit 1; }
[[ -n "$OUTDIR" ]] || OUTDIR="$RESULT_DIR/08_downstream"
[[ -n "$SAMPLESHEET" ]] || SAMPLESHEET="$RESULT_DIR/_automation/inputs/samplesheet.csv"
[[ -n "$CONTRAST_FILE" ]] || CONTRAST_FILE="$RESULT_DIR/_automation/inputs/contrasts.csv"

case "$LEVELS" in
  peak|bin|both) ;;
  *) echo "--levels must be peak, bin, or both" >&2; exit 1 ;;
esac
case "$OUTPUT_LAYOUT" in
  standard|legacy) ;;
  *) echo "--output-layout must be standard or legacy" >&2; exit 1 ;;
esac
case "$RUN_PROFILE_HEATMAPS" in
  true|false) ;;
  *) echo "--run-profile-heatmaps must be true or false" >&2; exit 1 ;;
esac
case "$ANNOTATION_MODE" in
  gene_te|gene|none) ;;
  *) echo "--annotation-mode must be gene_te, gene, or none" >&2; exit 1 ;;
esac

if [[ "$BACKGROUND" == "true" && "${RUN_ATAC_DOWNSTREAM_BG_CHILD:-0}" != "1" ]]; then
  log_dir="$RESULT_DIR/_automation/logs"
  mkdir -p "$log_dir"
  ts="$(date +%Y%m%d_%H%M%S)"
  log_file="$log_dir/downstream_${LEVELS}_${ts}.log"
  echo "[run_atac_downstream_only] background log: $log_file"
  nohup env RUN_ATAC_DOWNSTREAM_BG_CHILD=1 bash "$0" \
    --result-dir "$RESULT_DIR" \
    --species "$SPECIES" \
    --levels "$LEVELS" \
    --outdir "$OUTDIR" \
    --samplesheet "$SAMPLESHEET" \
    --contrast-file "$CONTRAST_FILE" \
    --output-layout "$OUTPUT_LAYOUT" \
    --run-profile-heatmaps "$RUN_PROFILE_HEATMAPS" \
    --profile-cores "$PROFILE_CORES" \
    --annotation-mode "$ANNOTATION_MODE" \
    ${TE_BED:+--te-bed "$TE_BED"} \
    ${GTF_FILE:+--gtf "$GTF_FILE"} \
    ${TE_VIOLIN_CLASS:+--te-violin-class "$TE_VIOLIN_CLASS"} \
    --te-violin-top-n "$TE_VIOLIN_TOP_N" \
    ${TE_FAMILY_FILTER:+--te-family-filter "$TE_FAMILY_FILTER"} \
    ${TE_NAME_FILTER:+--te-name-filter "$TE_NAME_FILTER"} \
    --label-top-n "$LABEL_TOP_N" \
    > "$log_file" 2>&1 &
  echo "[run_atac_downstream_only] PID: $!"
  exit 0
fi

mkdir -p "$OUTDIR/_inputs"
SYNC_META="$OUTDIR/_inputs/sample_metadata.csv"

if [[ -f "$SAMPLESHEET" ]]; then
  python3 - "$SAMPLESHEET" "$SYNC_META" <<'PY'
import csv, sys
inp, outp = sys.argv[1], sys.argv[2]
with open(inp, newline='') as f:
    reader = csv.DictReader(f)
    required = {'sample', 'condition', 'replicate'}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f'samplesheet missing columns: {sorted(missing)}')
    rows = []
    for row in reader:
        rows.append({
            'sample': str(row.get('sample', '')).strip(),
            'condition': str(row.get('condition', 'NA')).strip(),
            'replicate': str(row.get('replicate', 'NA')).strip(),
        })
if not rows:
    raise SystemExit('samplesheet has no sample rows')
with open(outp, 'w', newline='') as wf:
    writer = csv.DictWriter(wf, fieldnames=['sample', 'condition', 'replicate'])
    writer.writeheader()
    writer.writerows(rows)
PY
elif [[ -f "$RESULT_DIR/07_counts/sample_metadata.csv" ]]; then
  cp "$RESULT_DIR/07_counts/sample_metadata.csv" "$SYNC_META"
else
  echo "sample metadata not found. Please edit $SAMPLESHEET first." >&2
  exit 1
fi

if [[ -f "$CONTRAST_FILE" ]]; then
  python3 - "$SYNC_META" "$CONTRAST_FILE" <<'PY'
import csv, sys
meta_file, contrast_file = sys.argv[1], sys.argv[2]
conds = {r['condition'] for r in csv.DictReader(open(meta_file, newline='')) if r.get('condition') not in ('', 'NA', None)}
bad = []
rows = list(csv.DictReader(open(contrast_file, newline='')))
if not rows:
    raise SystemExit(f'contrast file has no contrasts: {contrast_file}')
pairs = []
for i, row in enumerate(rows, start=2):
    case = row.get('case', '').strip()
    control = row.get('control', '').strip()
    pairs.append((case, control, i))
    for key in ('case', 'control'):
        val = row.get(key, '').strip()
        if not val or val not in conds:
            bad.append(f'line {i}: {key}={val!r}')
if bad:
    raise SystemExit('contrast groups not found in sample metadata: ' + '; '.join(bad))
seen = {}
for case, control, line in pairs:
    pair = (case, control)
    if pair in seen:
        raise SystemExit(
            f'duplicate contrast: lines {seen[pair]} and {line} both contain {case} vs {control}'
        )
    seen[pair] = line
warned = set()
for case, control, line in pairs:
    reverse = (control, case)
    if reverse in seen and reverse not in warned and (case, control) not in warned:
        print(
            f'WARNING: reciprocal contrasts detected at lines {line} and {seen[reverse]}: '
            f'{case} vs {control} and {control} vs {case}. Both will run, but their effects are mirror directions.',
            file=sys.stderr,
        )
        warned.add((case, control))
        warned.add(reverse)
PY
else
  echo "contrast file not found: $CONTRAST_FILE" >&2
  exit 1
fi

if [[ -z "$TE_BED" && -f "$CONFIG_FILE" ]]; then TE_BED="$(get_species_param "$SPECIES" te_bed "$CONFIG_FILE" || true)"; fi
if [[ -z "$GTF_FILE" && -f "$CONFIG_FILE" ]]; then GTF_FILE="$(get_species_param "$SPECIES" gtf_genes "$CONFIG_FILE" || true)"; fi

run_level() {
  local level="$1"
  local count_file="$2"
  local region_bed="$3"
  local level_out="$4"
  local r_out="$level_out"
  [[ -f "$count_file" ]] || { echo "count file not found: $count_file" >&2; exit 1; }
  [[ -f "$region_bed" ]] || { echo "region BED not found: $region_bed" >&2; exit 1; }
  mkdir -p "$level_out"
  if [[ "$OUTPUT_LAYOUT" == "standard" ]]; then
    r_out="$level_out/legacy/raw_r_output"
    mkdir -p "$r_out"
  fi
  CMD=(Rscript "$ROOT_DIR/atacseq-downstream/run_downstream_atac.R"
    --count-file "$count_file"
    --sample-meta "$SYNC_META"
    --peak-bed "$region_bed"
    --outdir "$r_out"
    --analysis-level "$level"
    --species "$SPECIES"
    --annotation-mode "$ANNOTATION_MODE"
    --label-top-n "$LABEL_TOP_N"
    --contrast-file "$CONTRAST_FILE")
  if [[ -n "$TE_BED" ]]; then CMD+=(--te-bed "$TE_BED"); fi
  if [[ -n "$GTF_FILE" ]]; then CMD+=(--gtf "$GTF_FILE"); fi
  if [[ -n "$TE_VIOLIN_CLASS" ]]; then CMD+=(--te-violin-class "$TE_VIOLIN_CLASS"); fi
  if [[ -n "$TE_VIOLIN_TOP_N" ]]; then CMD+=(--te-violin-top-n "$TE_VIOLIN_TOP_N"); fi
  if [[ -n "$TE_FAMILY_FILTER" ]]; then CMD+=(--te-family-filter "$TE_FAMILY_FILTER"); fi
  if [[ -n "$TE_NAME_FILTER" ]]; then CMD+=(--te-name-filter "$TE_NAME_FILTER"); fi
  echo "[run_atac_downstream_only] level=$level"
  echo "[run_atac_downstream_only] ${CMD[*]}"
  "${CMD[@]}"
  if [[ "$OUTPUT_LAYOUT" == "standard" ]]; then
    if [[ ! -f "$ORGANIZER" ]]; then
      echo "organizer not found: $ORGANIZER" >&2
      exit 1
    fi
    python3 "$ORGANIZER" --raw-dir "$r_out" --outdir "$level_out" --clean
  fi
}

if [[ "$LEVELS" == "peak" || "$LEVELS" == "both" ]]; then
  run_level peak \
    "$RESULT_DIR/07_counts/consensus_peak_counts.txt" \
    "$RESULT_DIR/06_consensus_peaks/consensus_peaks.bed" \
    "$OUTDIR/peak_level"
fi

if [[ "$LEVELS" == "bin" || "$LEVELS" == "both" ]]; then
  bin_count="$RESULT_DIR/07_counts/bin_level/fixedbin_counts.txt"
  bin_bed="$(ls "$RESULT_DIR"/07_counts/bin_level/fixed_bins_*.bed 2>/dev/null | head -1 || true)"
  if [[ -z "$bin_bed" ]]; then
    if [[ "$LEVELS" == "bin" ]]; then
      echo "fixed-bin BED not found under $RESULT_DIR/07_counts/bin_level" >&2
      exit 1
    fi
    echo "[run_atac_downstream_only] skip bin level: fixed-bin BED not found" >&2
  else
    run_level bin "$bin_count" "$bin_bed" "$OUTDIR/bin_level"
  fi
fi

if [[ "$RUN_PROFILE_HEATMAPS" == "true" ]]; then
  if [[ ! -x "$PROFILE_HEATMAP_RUNNER" ]]; then
    echo "profile heatmap runner not found or not executable: $PROFILE_HEATMAP_RUNNER" >&2
    exit 1
  fi
  echo "[run_atac_downstream_only] run profile heatmaps"
  bash "$PROFILE_HEATMAP_RUNNER" \
    --result-dir "$RESULT_DIR" \
    --species "$SPECIES" \
    --outdir "$OUTDIR/profile_heatmaps" \
    --samplesheet "$SAMPLESHEET" \
    --cores "$PROFILE_CORES"
fi
