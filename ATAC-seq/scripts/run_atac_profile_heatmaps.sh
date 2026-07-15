#!/usr/bin/env bash
set -euo pipefail

log(){ echo "[run_atac_profile_heatmaps] $*"; }

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RESULT_DIR=""
OUTDIR=""
SPECIES="hg38"
BW_GLOB=""
SAMPLESHEET=""
CORES=8
PRESET="standard"
FEATURES="tss,gene_body,te,l1"
TSS_BED=""
GENE_BODY_BED=""
TE_BED=""
L1_BED=""
UP=3000
DOWN=3000
BODY_LEN=5000
TE_FLANK=15000
MAX_TE_REGIONS=""
MAX_L1_REGIONS=""
MAX_TSS_REGIONS=""
MAX_GENE_REGIONS=""
WRITE_MATRIX_TAB="false"

usage() {
  cat <<'USAGE'
Usage:
  bash scripts/run_atac_profile_heatmaps.sh --result-dir RESULTS [options]

Purpose:
  Standalone deepTools profile tool. The main FASTQ workflow calls quick/standard
  automatically; downstream-only does not rerun it.
  Gene body is scaled from TSS to TES. TE/L1 are anchored at the strand-aware
  5' end and displayed with 15 kb flanks by default.

Defaults:
  bigWig:      RESULTS/04_bw/*.bw
  output:      RESULTS/08_downstream/profile_heatmaps
  TSS BED:     /path/to/reference/TSS/<species>_tss.bed
  gene body:   /path/to/reference/TSS/<species>_gene_body.bed
  TE BED:      /path/to/reference/TE_bed/<species>_TE.bed
  L1 BED:      /path/to/reference/TE_bed/L1.bed

Options:
  --species STR
  --outdir DIR
  --bw-glob STR
  --samplesheet FILE
  --cores INT
  --preset STR                   quick | standard | full. Default: standard
  --features STR                 Comma-separated: tss,gene_body,te,l1
  --tss-bed FILE
  --gene-body-bed FILE
  --te-bed FILE
  --l1-bed FILE
  --upstream INT
  --downstream INT
  --body-length INT
  --te-flank INT                 TE/L1 5' flank on each side. Default: 15000
  --max-te-regions INT
  --max-l1-regions INT
  --max-tss-regions INT
  --max-gene-regions INT
  --write-matrix-tab BOOL        Write large plain-text matrices. Default: false

Preset region caps:
  quick:     1,000 per feature
  standard:  3,000 per feature
  full:     10,000 per feature
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --result-dir) RESULT_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --bw-glob) BW_GLOB="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --cores) CORES="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    --features) FEATURES="$2"; shift 2 ;;
    --tss-bed) TSS_BED="$2"; shift 2 ;;
    --gene-body-bed) GENE_BODY_BED="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --l1-bed) L1_BED="$2"; shift 2 ;;
    --upstream) UP="$2"; shift 2 ;;
    --downstream) DOWN="$2"; shift 2 ;;
    --body-length) BODY_LEN="$2"; shift 2 ;;
    --te-flank) TE_FLANK="$2"; shift 2 ;;
    --max-te-regions) MAX_TE_REGIONS="$2"; shift 2 ;;
    --max-l1-regions) MAX_L1_REGIONS="$2"; shift 2 ;;
    --max-tss-regions) MAX_TSS_REGIONS="$2"; shift 2 ;;
    --max-gene-regions) MAX_GENE_REGIONS="$2"; shift 2 ;;
    --write-matrix-tab) WRITE_MATRIX_TAB="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$RESULT_DIR" ]] || { usage >&2; exit 1; }
[[ -d "$RESULT_DIR" ]] || { echo "result dir not found: $RESULT_DIR" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || OUTDIR="$RESULT_DIR/08_downstream/profile_heatmaps"
[[ -n "$BW_GLOB" ]] || BW_GLOB="$RESULT_DIR/04_bw/*.bw"
[[ -n "$SAMPLESHEET" ]] || SAMPLESHEET="$RESULT_DIR/_automation/inputs/samplesheet.csv"
[[ -n "$TSS_BED" ]] || TSS_BED="/path/to/reference/TSS/${SPECIES}_tss.bed"
[[ -n "$GENE_BODY_BED" ]] || GENE_BODY_BED="/path/to/reference/TSS/${SPECIES}_gene_body.bed"
[[ -n "$TE_BED" ]] || TE_BED="/path/to/reference/TE_bed/${SPECIES}_TE.bed"
[[ -n "$L1_BED" ]] || L1_BED="/path/to/reference/TE_bed/L1.bed"

case "$PRESET" in
  quick)
    : "${MAX_TSS_REGIONS:=1000}"
    : "${MAX_GENE_REGIONS:=1000}"
    : "${MAX_TE_REGIONS:=1000}"
    : "${MAX_L1_REGIONS:=1000}"
    ;;
  standard)
    : "${MAX_TSS_REGIONS:=3000}"
    : "${MAX_GENE_REGIONS:=3000}"
    : "${MAX_TE_REGIONS:=3000}"
    : "${MAX_L1_REGIONS:=3000}"
    ;;
  full)
    : "${MAX_TSS_REGIONS:=10000}"
    : "${MAX_GENE_REGIONS:=10000}"
    : "${MAX_TE_REGIONS:=10000}"
    : "${MAX_L1_REGIONS:=10000}"
    ;;
  *) echo "--preset must be quick, standard, or full" >&2; exit 1 ;;
esac
case "$WRITE_MATRIX_TAB" in
  true|false) ;;
  *) echo "--write-matrix-tab must be true or false" >&2; exit 1 ;;
esac

export PATH="/path/to/.conda/envs/chipseq/bin:$PATH"
command -v computeMatrix >/dev/null || { echo "computeMatrix not found in PATH" >&2; exit 1; }
command -v plotHeatmap >/dev/null || { echo "plotHeatmap not found in PATH" >&2; exit 1; }
command -v plotProfile >/dev/null || { echo "plotProfile not found in PATH" >&2; exit 1; }

shopt -s nullglob
bw_files=( $BW_GLOB )
[[ ${#bw_files[@]} -gt 0 ]] || { echo "No bigWig matched: $BW_GLOB" >&2; exit 1; }

mkdir -p "$OUTDIR"/{figures,matrix,logs,regions}

labels_file="$OUTDIR/sample_labels.txt"
python3 - "$SAMPLESHEET" "$labels_file" "${bw_files[@]}" <<'PY'
import csv, os, sys
samplesheet, outp, *bws = sys.argv[1:]
mapping = {}
if os.path.exists(samplesheet):
    with open(samplesheet, newline='') as fh:
        for row in csv.DictReader(fh):
            sample = str(row.get('sample', '')).strip()
            cond = str(row.get('condition', '')).strip()
            rep = str(row.get('replicate', '')).strip()
            label = cond if cond and cond != 'NA' else sample
            if rep and rep != 'NA':
                label = f'{label}_rep{rep}' if label else f'rep{rep}'
            if sample:
                mapping[sample] = label
labels = []
seen = {}
for bw in bws:
    stem = os.path.basename(bw)
    for suffix in ('.bw', '.bigWig', '.clean', '.normalized'):
        stem = stem.replace(suffix, '')
    label = mapping.get(stem, stem)
    n = seen.get(label, 0) + 1
    seen[label] = n
    if n > 1:
        label = f'{label}_{n}'
    labels.append(label)
with open(outp, 'w') as fh:
    fh.write('\n'.join(labels) + '\n')
PY
mapfile -t labels < "$labels_file"

feature_enabled() {
  [[ ",${FEATURES}," == *",$1,"* ]]
}

sanitize_bed() {
  local in_bed="$1"
  local out_bed="$2"
  [[ -f "$in_bed" ]] || return 1
  awk 'BEGIN{OFS="\t"} !/^#/ && NF>=3 {
    s=$2; e=$3; if (s < 0) s = 0; if (e <= s) next;
    print $1, s, e, (NF>=4?$4:"."), (NF>=5?$5:"0"), (NF>=6?$6:".")
  }' "$in_bed" | sort -k1,1 -k2,2n > "$out_bed"
  [[ -s "$out_bed" ]]
}

sample_bed() {
  local in_bed="$1"
  local out_bed="$2"
  local max_n="$3"
  cp "$in_bed" "$out_bed"
  if [[ "$max_n" =~ ^[0-9]+$ && "$max_n" -gt 0 ]]; then
    local n
    n="$(wc -l < "$in_bed")"
    if [[ "$n" -gt "$max_n" ]]; then
      set +o pipefail
      awk 'BEGIN{srand(1)} {print rand()"\t"$0}' "$in_bed" | sort -k1,1n | cut -f2- | head -n "$max_n" | sort -k1,1 -k2,2n > "$out_bed"
      set -o pipefail
      log "sampled $(basename "$in_bed"): $max_n / $n"
    fi
  fi
}

run_matrix() {
  local name="$1"
  local mode="$2"
  local bed="$3"
  local body_len="$4"
  local before="$5"
  local after="$6"
  local matrix="$OUTDIR/matrix/${name}.matrix.gz"
  local matrix_tab="$OUTDIR/matrix/${name}.matrix.tab"
  local sorted="$OUTDIR/regions/${name}.sorted_regions.bed"
  local log_file="$OUTDIR/logs/${name}.computeMatrix.log"
  local common=(-R "$bed" -S "${bw_files[@]}" -b "$before" -a "$after" --samplesLabel "${labels[@]}" --numberOfProcessors "$CORES" --skipZeros --sortRegions descend --sortUsing mean -o "$matrix" --outFileSortedRegions "$sorted")
  if [[ "$WRITE_MATRIX_TAB" == "true" ]]; then
    common+=(--outFileNameMatrix "$matrix_tab")
  fi
  log "computeMatrix $name"
  if [[ "$mode" == "reference-point" ]]; then
    computeMatrix reference-point --referencePoint TSS "${common[@]}" > "$log_file" 2>&1
  else
    computeMatrix scale-regions --regionBodyLength "$body_len" "${common[@]}" > "$log_file" 2>&1
  fi
  mkdir -p "$OUTDIR/figures/$name"
  label_args=()
  display_name="$name"
  case "$name" in
    tss) display_name="TSS" ;;
    gene_body) display_name="Gene body"; label_args=(--startLabel "TSS" --endLabel "TES") ;;
    te) display_name="TE"; label_args=(--refPointLabel "5' end") ;;
    l1) display_name="L1"; label_args=(--refPointLabel "5' end") ;;
  esac
  if [[ "$name" == "te" || "$name" == "l1" ]]; then
    x_label="Position relative to $display_name 5' end (bp)"
  else
    x_label="Position relative to $display_name (bp)"
  fi
  plotHeatmap -m "$matrix" -out "$OUTDIR/figures/$name/${name}.profile_heatmap.pdf" --sortRegions descend --whatToShow 'plot, heatmap and colorbar' --heatmapWidth 7 --regionsLabel "$display_name" --xAxisLabel "$x_label" --plotTitle "$display_name ATAC accessibility" "${label_args[@]}"
  plotHeatmap -m "$matrix" -out "$OUTDIR/figures/$name/${name}.profile_heatmap.png" --sortRegions descend --whatToShow 'plot, heatmap and colorbar' --heatmapWidth 7 --regionsLabel "$display_name" --xAxisLabel "$x_label" --plotTitle "$display_name ATAC accessibility" "${label_args[@]}"
  plotProfile -m "$matrix" -out "$OUTDIR/figures/$name/${name}.profile.pdf" --perGroup --plotWidth 16 --regionsLabel "$display_name" --plotTitle "$display_name ATAC profile" "${label_args[@]}"
  plotProfile -m "$matrix" -out "$OUTDIR/figures/$name/${name}.profile.png" --perGroup --plotWidth 16 --regionsLabel "$display_name" --plotTitle "$display_name ATAC profile" "${label_args[@]}"
}

run_feature() {
  local name="$1" display="$2" bed="$3" mode="$4" cap="$5" before="$6" after="$7"
  if ! feature_enabled "$name"; then
    log "skip $display: disabled by --features"
    return
  fi
  if ! sanitize_bed "$bed" "$OUTDIR/regions/${name}.valid.bed"; then
    log "skip $display: BED not found or empty: $bed"
    return
  fi
  sample_bed "$OUTDIR/regions/${name}.valid.bed" "$OUTDIR/regions/${name}.sampled.bed" "$cap"
  run_matrix "$name" "$mode" "$OUTDIR/regions/${name}.sampled.bed" "$BODY_LEN" "$before" "$after"
}

run_feature tss "TSS" "$TSS_BED" reference-point "$MAX_TSS_REGIONS" "$UP" "$DOWN"
run_feature gene_body "gene body" "$GENE_BODY_BED" scale-regions "$MAX_GENE_REGIONS" "$UP" "$DOWN"
run_feature te "TE" "$TE_BED" reference-point "$MAX_TE_REGIONS" "$TE_FLANK" "$TE_FLANK"
run_feature l1 "L1" "$L1_BED" reference-point "$MAX_L1_REGIONS" "$TE_FLANK" "$TE_FLANK"

cat > "$OUTDIR/README_profile_heatmaps.md" <<EOF
# ATAC profile heatmaps

This directory contains deepTools signal profile heatmaps generated from bigWig files.

- \`figures/tss\`: reference-point heatmap/profile around TSS using \`$TSS_BED\`.
- \`figures/gene_body\`: TSS-to-TES scale-regions profile using \`$GENE_BODY_BED\`.
- \`figures/te\`: strand-aware TE 5' end profile with ±${TE_FLANK} bp using \`$TE_BED\`.
- \`figures/l1\`: strand-aware L1 5' end profile with ±${TE_FLANK} bp using \`$L1_BED\`.
- preset: \`$PRESET\`; selected features: \`$FEATURES\`.
- \`matrix\`: computeMatrix output files.
- \`regions\`: sanitized and sampled BED files actually used.
- \`logs\`: computeMatrix logs.

These are global signal-profile plots, not differential accessibility tests. Differential tables remain under \`peak_level/results/differential\` and \`bin_level/results/differential\`.
EOF

log "done: $OUTDIR"
