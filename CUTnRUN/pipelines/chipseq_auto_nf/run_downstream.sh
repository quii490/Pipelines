#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR=""
SPECIES="hg38"
THREADS=8
PEAK_JOBS="${PEAK_JOBS:-1}"
RESOURCE_TIER="${CUTRUN_RESOURCE_TIER:-standard}"
# A 100k deterministic cap keeps the default L1 matrices practical while
# retaining enough loci for stable aggregate profiles.  Set 0 for all loci.
TE_LOCUS_MAX_REGIONS="${TE_LOCUS_MAX_REGIONS:-}"
TE_LOCUS_BIN_SIZE="${TE_LOCUS_BIN_SIZE:-}"
PEAK_HEATMAP_MAX_REGIONS="${PEAK_HEATMAP_MAX_REGIONS:-}"
SKIP_BW_COR=false
SKIP_COUNT_DRAW=false
TE_METHODS=""
TE_METHODS_EXECUTE=false
RESUME=true
TE_REPEAT_BED=""
TE_SAF=""
TE_ANNO=""
T3E_DIR=""
T3E_PYTHON=""
T3E_MAX_BED_READS="${T3E_MAX_BED_READS:-1000000}"
T3E_ITERATIONS="${T3E_ITERATIONS:-100}"
ALLO_COMMAND=""
REPENTOOLS_COMMAND=""
REPENTOOLS_INDEX_DIR="${REPENTOOLS_INDEX_DIR:-/path/to/CUTnRUN/resources/RepEnTools/chm13v2/indexes}"
REPENTOOLS_GTF="${REPENTOOLS_GTF:-/path/to/CUTnRUN/resources/RepEnTools/chm13v2/annotation/rmsk.gtf}"
REPENTOOLS_RET="${REPENTOOLS_RET:-/path/to/.local/share/CUTnRUN-tools/RepEnTools/ret}"
REPENTOOLS_TARGET_GROUPS="${REPENTOOLS_TARGET_GROUPS:-}"
REPENTOOLS_INPUT_SAMPLES="${REPENTOOLS_INPUT_SAMPLES:-}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --resource-tier) RESOURCE_TIER="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --peak-jobs) PEAK_JOBS="$2"; shift 2 ;;
    --skip-bw-cor) SKIP_BW_COR=true; shift ;;
    --skip-count-draw) SKIP_COUNT_DRAW=true; shift ;;
    --te-methods) TE_METHODS="$2"; shift 2 ;;
    --te-methods-execute) TE_METHODS_EXECUTE=true; shift ;;
    --resume) RESUME=true; shift ;;
    --no-resume|--force) RESUME=false; shift ;;
    --te-repeat-bed) TE_REPEAT_BED="$2"; shift 2 ;;
    --te-saf) TE_SAF="$2"; shift 2 ;;
    --te-anno) TE_ANNO="$2"; shift 2 ;;
    --t3e-dir) T3E_DIR="$2"; shift 2 ;;
    --t3e-python) T3E_PYTHON="$2"; shift 2 ;;
    --t3e-max-bed-reads) T3E_MAX_BED_READS="$2"; shift 2 ;;
    --t3e-iterations) T3E_ITERATIONS="$2"; shift 2 ;;
    --allo-command) ALLO_COMMAND="$2"; shift 2 ;;
    --repentools-command) REPENTOOLS_COMMAND="$2"; shift 2 ;;
    --repentools-index-dir) REPENTOOLS_INDEX_DIR="$2"; shift 2 ;;
    --repentools-gtf) REPENTOOLS_GTF="$2"; shift 2 ;;
    --repentools-ret) REPENTOOLS_RET="$2"; shift 2 ;;
    --repentools-target-groups) REPENTOOLS_TARGET_GROUPS="$2"; shift 2 ;;
    --repentools-input-samples) REPENTOOLS_INPUT_SAMPLES="$2"; shift 2 ;;
    -h|--help)
      echo "Usage: run_downstream.sh --results-dir DIR [--species hg38|mm39] [--threads N] [--peak-jobs N] [--resource-tier small|standard|full] [--resume|--no-resume] [--te-methods LIST] [--repentools-ret FILE]"
      echo "Peak annotation: ChIPseeker by default; set PEAK_ANNOTATOR=chippeakanno for legacy ChIPpeakAnno"
      exit 0
      ;;
    *) echo "ERROR: unknown option $1" >&2; exit 1 ;;
  esac
done

case "$RESOURCE_TIER" in
  small) DEFAULT_TE_LOCUS_MAX=25000; DEFAULT_TE_BIN=50; DEFAULT_PEAK_HEATMAP_MAX=25000 ;;
  standard) DEFAULT_TE_LOCUS_MAX=100000; DEFAULT_TE_BIN=25; DEFAULT_PEAK_HEATMAP_MAX=100000 ;;
  full) DEFAULT_TE_LOCUS_MAX=0; DEFAULT_TE_BIN=5; DEFAULT_PEAK_HEATMAP_MAX=0 ;;
  *) echo "ERROR: --resource-tier must be small, standard or full" >&2; exit 1 ;;
esac
[[ -n "$TE_LOCUS_MAX_REGIONS" ]] || TE_LOCUS_MAX_REGIONS="$DEFAULT_TE_LOCUS_MAX"
[[ -n "$TE_LOCUS_BIN_SIZE" ]] || TE_LOCUS_BIN_SIZE="$DEFAULT_TE_BIN"
[[ -n "$PEAK_HEATMAP_MAX_REGIONS" ]] || PEAK_HEATMAP_MAX_REGIONS="$DEFAULT_PEAK_HEATMAP_MAX"

[[ -d "$RESULTS_DIR" ]] || { echo "ERROR: results directory not found: $RESULTS_DIR" >&2; exit 1; }
[[ "$PEAK_JOBS" =~ ^[0-9]+$ && "$PEAK_JOBS" -ge 1 ]] || { echo "ERROR: --peak-jobs must be a positive integer" >&2; exit 1; }
export CUTRUN_RESOURCE_TIER="$RESOURCE_TIER"
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
downstream_dir="${RESULTS_DIR}/09_downstream"
chipseq_bin="${CHIPSEQ_BIN:-${HOME}/.conda/envs/chipseq/bin}"
if [[ ! -x "${chipseq_bin}/python" && -x "/path/to/.conda/envs/chipseq/bin/python" ]]; then
  chipseq_bin="/path/to/.conda/envs/chipseq/bin"
fi
downstream_bin="${DOWNSTREAM_BIN:-${HOME}/.conda/envs/downstream/bin}"
if [[ ! -x "${downstream_bin}/Rscript" && -x "/path/to/.conda/envs/downstream/bin/Rscript" ]]; then
  downstream_bin="/path/to/.conda/envs/downstream/bin"
fi
if [[ ! -x "${downstream_bin}/Rscript" && -x "${chipseq_bin}/Rscript" ]]; then
  downstream_bin="$chipseq_bin"
fi
R_BIN="${R_BIN:-${downstream_bin}/Rscript}"
if [[ ! -x "$R_BIN" ]]; then
  R_BIN="$(command -v Rscript || true)"
fi
export PATH="${chipseq_bin}:${HOME}/.local/bin:${PATH}"

# ChIPseeker is the default annotation implementation. ChIPpeakAnno remains
# available only when explicitly requested with PEAK_ANNOTATOR=chippeakanno.
ANNOTATION_SCRIPT="${repo_root}/tools/cutrun_cli/cutrun_annotate_peaks_chipseeker.R"
if [[ "${PEAK_ANNOTATOR:-chipseeker}" == "chippeakanno" ]]; then
  if [[ -x "$R_BIN" ]] && "$R_BIN" -e 'quit(status = if (requireNamespace("ChIPpeakAnno", quietly = TRUE)) 0 else 1)' >/dev/null 2>&1; then
    ANNOTATION_SCRIPT="${repo_root}/tools/cutrun_cli/cutrun_annotate_peaks_chippeakanno.R"
  else
    echo "[run_downstream] PEAK_ANNOTATOR=chippeakanno requested but ChIPpeakAnno is unavailable; using ChIPseeker" >&2
  fi
fi
HOMER_AVAILABLE=false
if command -v findMotifsGenome.pl >/dev/null 2>&1 && command -v annotatePeaks.pl >/dev/null 2>&1; then
  HOMER_AVAILABLE=true
fi

mkdir -p "$downstream_dir"
RESOLVED_MANIFEST="${RESULTS_DIR}/manifest/resolved_manifest.csv"
if [[ ! -s "$RESOLVED_MANIFEST" ]]; then
  RESOLVED_MANIFEST="${RESULTS_DIR}/_automation/inputs/manifest.csv"
fi
STATUS_FILE="${downstream_dir}/module_status.tsv"
STATUS_PY="${repo_root}/tools/cutrun_cli/cutrun_status.py"
RUN_MANIFEST_PY="${repo_root}/tools/cutrun_cli/cutrun_run_manifest.py"
TE_QC_PY="${repo_root}/tools/cutrun_cli/cutrun_te_qc.py"
TE_MULTIMAP_PY="${repo_root}/tools/cutrun_cli/cutrun_te_multimap.py"
CONSENSUS_PY="${repo_root}/tools/cutrun_cli/cutrun_consensus_peaks.py"
PATHWAY_SCRIPT="${repo_root}/tools/cutrun_cli/cutrun_pathway_enrichment.R"
REPORT_PY="${repo_root}/tools/cutrun_cli/cutrun_report.py"
TE_VIS_PY="${repo_root}/tools/te_methods/render_te_method_visuals.py"
PYTHON_BIN="${chipseq_bin}/python"
if [[ ! -x "$PYTHON_BIN" ]]; then
  PYTHON_BIN="$(command -v python3 || command -v python)"
fi
RUN_ID="${CUTRUN_RUN_ID:-run_$(date -u +%Y%m%dT%H%M%SZ)}"
export CUTRUN_RUN_ID="$RUN_ID"
if [[ -x "$STATUS_PY" || -f "$STATUS_PY" ]]; then
  "$PYTHON_BIN" "$STATUS_PY" init \
    --status-file "$STATUS_FILE" \
    --run-id "$RUN_ID" \
    --meta "${downstream_dir}/status_ledger.json"
fi

record_status() {
  local label="$1" status="$2" started="$3" finished="$4" reason="${5:-}" outputs_ok="${6:-false}" output_paths="${7:-}"
  [[ -f "$STATUS_PY" ]] || return 0
  local args=("$PYTHON_BIN" "$STATUS_PY" record
    --status-file "$STATUS_FILE" --run-id "$RUN_ID" --module "$label" --status "$status"
    --started "$started" --finished "$finished" --reason "$reason" --output-paths "$output_paths")
  [[ "$outputs_ok" == true ]] && args+=(--outputs-ok)
  if command -v flock >/dev/null 2>&1; then
    (
      flock 9
      "${args[@]}"
    ) 9>"${STATUS_FILE}.lock"
  else
    "${args[@]}"
  fi
}

run_optional() {
  local label="$1"
  shift
  local started finished
  started="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "[run_downstream] START ${label}"
  if "$@"; then
    finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    record_status "optional:${label}" PASS "$started" "$finished" "command completed" true
    echo "[run_downstream] DONE ${label}"
  else
    local exit_code=$?
    finished="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    if [[ "$exit_code" -eq 42 ]]; then
      record_status "optional:${label}" SKIP "$started" "$finished" "module reported scientifically unsupported or unavailable" false
      echo "[run_downstream] SKIP ${label}"
    else
      record_status "optional:${label}" FAIL "$started" "$finished" "command failed (rc=${exit_code}); inspect module log" false
      echo "[run_downstream] WARNING ${label} failed; continuing" >&2
    fi
  fi
}

run_cached() {
  local label="$1" sentinel="$2"
  shift 2
  if [[ "$RESUME" == true && -s "$sentinel" ]]; then
    local stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    local cached_status=PASS cached_reason="cached output reused"
    if grep -q $'^status\tSKIP' "$sentinel" 2>/dev/null; then
      cached_status=SKIP; cached_reason="cached module status is SKIP"
    elif grep -q $'^status\tFAIL' "$sentinel" 2>/dev/null; then
      cached_status=FAIL; cached_reason="cached module status is FAIL"
    fi
    if [[ "$cached_status" == FAIL ]]; then
      echo "[run_downstream] CACHE INVALID ${label}: ${sentinel}; rerunning"
    else
      record_status "optional:${label}" "$cached_status" "$stamp" "$stamp" "$cached_reason" true "$sentinel"
      echo "[run_downstream] CACHE ${label}: ${sentinel} (${cached_status})"
      return 0
    fi
  fi
  run_optional "$label" "$@"
}

if [[ -f "$RUN_MANIFEST_PY" ]]; then
  run_optional "run manifest" "$PYTHON_BIN" "$RUN_MANIFEST_PY" \
    --results-dir "$RESULTS_DIR" \
    --manifest "$RESOLVED_MANIFEST" \
    --config "${repo_root}/pipelines/chipseq_auto_nf/nextflow.config" \
    --run-id "$RUN_ID"
fi

stage_tracks() {
  local source_dir="$1"
  local pattern="$2"
  local stage_dir="$3"
  local track
  mkdir -p "$stage_dir"
  find "$stage_dir" -maxdepth 1 -type l -delete
  shopt -s nullglob
  for track in "$source_dir"/$pattern; do
    ln -s "$track" "$stage_dir/$(basename "$track")"
  done
}

if [[ "$SKIP_BW_COR" != true ]]; then
  stage_tracks "${RESULTS_DIR}/05_tracks" "*_100bp_rpkm.bw" "${downstream_dir}/_inputs/standard_rpkm"
  stage_tracks "${RESULTS_DIR}/05_tracks" "*_10bp_rpgc.bw" "${downstream_dir}/_inputs/standard_rpgc"
  stage_tracks "${RESULTS_DIR}/05_tracks_te" "*_te_*bp_rpgc.bw" "${downstream_dir}/_inputs/te_rpgc"
  stage_tracks "${RESULTS_DIR}/05_tracks_te_locus_best" "*_te_locus_best_*bp_cpm.bw" "${downstream_dir}/_inputs/te_locus_best_cpm"

  for mode in standard_rpkm standard_rpgc te_rpgc te_locus_best_cpm; do
    input_dir="${downstream_dir}/_inputs/${mode}"
    if compgen -G "${input_dir}/*.bw" >/dev/null; then
      run_cached "bigWig correlation ${mode}" "${downstream_dir}/correlation/bigwig/${mode}/${mode}.pearson.heatmap.pdf" \
        bash "${repo_root}/tools/cutrun_cli/cutrun_bw_cor.sh" \
          --bw-dir "$input_dir" \
          --out-prefix "${downstream_dir}/correlation/bigwig/${mode}/${mode}" \
          --threads "$THREADS"
    fi
  done
fi

if compgen -G "${RESULTS_DIR}/04_clean_bam/*.bam" >/dev/null; then
  run_cached "standard BAM correlation" "${downstream_dir}/correlation/bam/standard/standard.pearson.heatmap.pdf" \
    bash "${repo_root}/tools/cutrun_cli/cutrun_bam_cor.sh" \
      --bam-dir "${RESULTS_DIR}/04_clean_bam" \
      --out-prefix "${downstream_dir}/correlation/bam/standard/standard" \
      --threads "$THREADS"
fi

if compgen -G "${RESULTS_DIR}/04_te_bam/*.bam" >/dev/null; then
  run_cached "TE BAM correlation" "${downstream_dir}/correlation/bam/te/te.pearson.heatmap.pdf" \
    bash "${repo_root}/tools/cutrun_cli/cutrun_bam_cor.sh" \
      --bam-dir "${RESULTS_DIR}/04_te_bam" \
      --out-prefix "${downstream_dir}/correlation/bam/te/te" \
      --threads "$THREADS"
fi

if [[ "$SKIP_COUNT_DRAW" != true ]]; then
  gene_counts="${RESULTS_DIR}/07_featurecounts/featurecounts_gene.txt"
  te_counts="${RESULTS_DIR}/07_featurecounts/featurecounts_te.txt"
  if [[ -s "$gene_counts" && -s "$te_counts" ]]; then
    run_cached "unified gene/TE count plots" "${RESULTS_DIR}/08_analysis/results/intermediate_data/02_plot_data.rds" \
      bash "${repo_root}/pipelines/count_draw/count_draw.sh" \
        --reference "$SPECIES" \
        --manifest "${RESULTS_DIR}/manifest/resolved_manifest.csv" \
        --gene-counts "$gene_counts" \
        --te-counts "$te_counts" \
        --gene-bam-dir "${RESULTS_DIR}/04_clean_bam" \
        --te-bam-dir "${RESULTS_DIR}/04_te_bam" \
        --output-dir "${RESULTS_DIR}/08_analysis" \
        --threads "$THREADS" \
        --skip-featurecounts
  fi
fi

if [[ -n "$TE_METHODS" ]]; then
  te_manifest="${RESULTS_DIR}/manifest/resolved_manifest.csv"
  te_adapter="${repo_root}/tools/te_methods/run_te_methods.sh"
  # RepEnTools is a separate CHM13 FASTQ workflow.  Do not pass it to the
  # BAM-oriented adapters (which would only report a misleading hg38 skip).
  te_methods_bam="$(printf '%s' "$TE_METHODS" | tr ',' '\n' | awk 'tolower($0)!="repentools" && tolower($0)!="repen-tools" && NF' | paste -sd, -)"
  if [[ -n "$te_methods_bam" && -s "$te_manifest" && -x "$te_adapter" ]]; then
    te_cmd=(bash "$te_adapter"
      --manifest "$te_manifest"
      --bam-dir "${RESULTS_DIR}/04_te_bam"
      --out-dir "${downstream_dir}/te_methods"
      --species "$SPECIES"
      --methods "$te_methods_bam"
      --threads "$THREADS")
    [[ -n "$TE_REPEAT_BED" ]] && te_cmd+=(--repeat-bed "$TE_REPEAT_BED")
    [[ -n "$TE_SAF" ]] && te_cmd+=(--te-saf "$TE_SAF")
    [[ -n "$TE_ANNO" ]] && te_cmd+=(--te-anno "$TE_ANNO")
    [[ -n "$T3E_DIR" ]] && te_cmd+=(--t3e-dir "$T3E_DIR")
    [[ -n "$T3E_PYTHON" ]] && te_cmd+=(--t3e-python "$T3E_PYTHON")
    te_cmd+=(--t3e-max-bed-reads "$T3E_MAX_BED_READS" --t3e-iterations "$T3E_ITERATIONS")
    [[ -n "$ALLO_COMMAND" ]] && te_cmd+=(--allo-command "$ALLO_COMMAND")
    [[ -n "$REPENTOOLS_COMMAND" ]] && te_cmd+=(--repentools-command "$REPENTOOLS_COMMAND")
    [[ "$TE_METHODS_EXECUTE" == true ]] && te_cmd+=(--execute)
    run_optional "TE method adapters (${TE_METHODS})" "${te_cmd[@]}"
  elif [[ -n "$te_methods_bam" ]]; then
    now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    record_status "optional:TE method adapters (${TE_METHODS})" SKIP "$now" "$now" "resolved manifest or adapter missing" false
    echo "[run_downstream] TE method adapters skipped: resolved manifest or adapter missing"
  fi

  if printf '%s' "$TE_METHODS" | tr ',' '\n' | awk 'tolower($0)=="repentools" || tolower($0)=="repen-tools"' | grep -q .; then
    repentools_adapter="${repo_root}/tools/te_methods/run_repentools_fastq.sh"
    repentools_dir="${downstream_dir}/te_methods/repentools_fastq"
    if [[ -s "$te_manifest" && -x "$repentools_adapter" ]]; then
      repentools_cmd=(bash "$repentools_adapter"
        --manifest "$te_manifest"
        --out-dir "$repentools_dir"
        --index-dir "$REPENTOOLS_INDEX_DIR"
        --gtf "$REPENTOOLS_GTF"
        --ret "$REPENTOOLS_RET"
        --target-groups "$REPENTOOLS_TARGET_GROUPS"
        --input-samples "$REPENTOOLS_INPUT_SAMPLES"
        --threads "$THREADS")
      [[ "$RESUME" == false ]] && repentools_cmd+=(--force)
      [[ "$TE_METHODS_EXECUTE" == true ]] && repentools_cmd+=(--execute)
      run_optional "RepEnTools FASTQ adapter" "${repentools_cmd[@]}"
    else
      now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
      record_status "optional:RepEnTools FASTQ adapter" SKIP "$now" "$now" "manifest or adapter missing" false
      echo "[run_downstream] RepEnTools FASTQ adapter skipped: manifest or adapter missing"
    fi
  fi

  # The three upstream tools have incompatible result schemas and none of
  # them provides a common figure contract.  Render validated, method-specific
  # summaries after both adapters have run.  This step is intentionally cheap
  # and is rerun on every invocation so a newly completed adapter cannot be
  # hidden by a stale visualization sentinel.
  if [[ -f "$TE_VIS_PY" ]]; then
    run_optional "TE method visualizations" "$PYTHON_BIN" "$TE_VIS_PY" \
      --methods-dir "${downstream_dir}/te_methods" \
      --output-dir "${downstream_dir}/te_methods/visuals"
  fi
fi

shopt -s nullglob
# MACS3 always emits both narrow and broad calls in the core pipeline.  Keep
# both peak classes in downstream annotation, motif, overlap and signal plots;
# include the class in every output name so broad/narrow results never
# overwrite one another.
peak_files=(
  "${RESULTS_DIR}"/06_peaks/*/broad/*_peaks.broadPeak
  "${RESULTS_DIR}"/06_peaks/*/narrow/*_peaks.narrowPeak
)
if [[ ${#peak_files[@]} -eq 0 ]]; then
  echo "[run_downstream] No narrowPeak or broadPeak files found; peak-derived analyses skipped"
  now="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  record_status "optional:peak-derived analyses" SKIP "$now" "$now" "no narrowPeak or broadPeak files found" false
  # Keep the run inventory useful even when peak calling failed or was
  # intentionally skipped.  This makes a partial run auditable instead of
  # silently omitting the manifest/QC summary.
  summary_script="${repo_root}/tools/cutrun_cli/cutrun_results_summary.py"
  summary_manifest="${RESULTS_DIR}/manifest/resolved_manifest.csv"
  if [[ ! -s "$summary_manifest" ]]; then
    summary_manifest="${RESULTS_DIR}/_automation/inputs/manifest.csv"
  fi
  if [[ -f "$summary_script" && -s "$summary_manifest" ]]; then
    run_optional "results inventory and QC summary" \
      "$PYTHON_BIN" "$summary_script" \
        --results-dir "$RESULTS_DIR" \
        --manifest "$summary_manifest"
  fi
  if [[ -f "$TE_QC_PY" && -s "$summary_manifest" ]]; then
    run_optional "TE and classical QC metrics" "$PYTHON_BIN" "$TE_QC_PY" \
      --results-dir "$RESULTS_DIR" --manifest "$summary_manifest"
  fi
  if [[ -f "$TE_MULTIMAP_PY" && -d "${RESULTS_DIR}/04_te_bam" ]]; then
    run_optional "TE multimapping policy summary" "$PYTHON_BIN" "$TE_MULTIMAP_PY" \
      --bam-dir "${RESULTS_DIR}/04_te_bam" \
      --output "${downstream_dir}/te_multimap_sensitivity.tsv"
  fi
  if [[ -f "$RUN_MANIFEST_PY" && -s "$summary_manifest" ]]; then
    run_optional "final run manifest" "$PYTHON_BIN" "$RUN_MANIFEST_PY" \
      --results-dir "$RESULTS_DIR" --manifest "$summary_manifest" \
      --config "${repo_root}/pipelines/chipseq_auto_nf/nextflow.config" --run-id "$RUN_ID"
  fi
  "$PYTHON_BIN" "$STATUS_PY" finalize --status-file "$STATUS_FILE" --run-id "$RUN_ID" \
    --output "${downstream_dir}/module_status_summary.json" || true
  if [[ -f "$REPORT_PY" ]]; then
    run_optional "unified HTML/PDF report" "$PYTHON_BIN" "$REPORT_PY" --results-dir "$RESULTS_DIR"
  fi
  echo "[run_downstream] Module status: ${STATUS_FILE}"
  exit 0
fi

summary_manifest="${RESULTS_DIR}/manifest/resolved_manifest.csv"
if [[ ! -s "$summary_manifest" ]]; then
  summary_manifest="${RESULTS_DIR}/_automation/inputs/manifest.csv"
fi
if [[ -f "$CONSENSUS_PY" && -s "$summary_manifest" ]]; then
  run_cached "replicate consensus peaks" "${downstream_dir}/consensus_peaks/consensus_summary.tsv" \
    "$PYTHON_BIN" "$CONSENSUS_PY" \
      --results-dir "$RESULTS_DIR" \
      --manifest "$summary_manifest" \
      --out-dir "${downstream_dir}/consensus_peaks"
fi

if [[ ${#peak_files[@]} -gt 1 ]]; then
  peak_args=()
  for peak in "${peak_files[@]}"; do
    sample="$(basename "$(dirname "$(dirname "$peak")")")"
    peak_type="$(basename "$(dirname "$peak")")"
    peak_args+=(--peak "${sample}_${peak_type}=${peak}")
  done
  run_cached "peak overlap" "${downstream_dir}/peak_overlap/pairwise_overlap.tsv" \
    "$PYTHON_BIN" "${repo_root}/tools/cutrun_cli/cutrun_peak_overlap.py" \
      "${peak_args[@]}" \
      --out-dir "${downstream_dir}/peak_overlap"
fi

standard_signals=()
te_signals=()
locus_signals=()
for signal in "${RESULTS_DIR}"/05_tracks/*_10bp_rpgc.bw; do standard_signals+=("$signal"); done
for signal in "${RESULTS_DIR}"/05_tracks_te/*_te_*bp_rpgc.bw; do te_signals+=("$signal"); done
for signal in "${RESULTS_DIR}"/05_tracks_te_locus_best/*_te_locus_best_*bp_cpm.bw; do locus_signals+=("$signal"); done
standard_labels=()
te_labels=()
locus_labels=()
if [[ -n "${standard_signals[*]-}" ]]; then
  for signal in "${standard_signals[@]}"; do
    standard_labels+=("$(basename "$signal" _10bp_rpgc.bw)")
  done
fi
if [[ -n "${te_signals[*]-}" ]]; then
  for signal in "${te_signals[@]}"; do
    label="$(basename "$signal")"
    te_labels+=("${label%%_te_*}")
  done
fi
if [[ -n "${locus_signals[*]-}" ]]; then
  for signal in "${locus_signals[@]}"; do
    label="$(basename "$signal")"
    locus_labels+=("${label%%_te_locus_best_*}")
  done
fi

l1_regions="${repo_root}/resources/CUTRUN_analysis/anno/L1_${SPECIES}.bed"
if [[ -s "$l1_regions" && -n "${locus_signals[*]-}" ]]; then
  run_cached "TE/L1 locus heatmap" "${downstream_dir}/heatmaps/te_locus_best/L1.profile.pdf" \
    bash "${repo_root}/tools/cutrun_cli/cutrun_te_locus_heatmap.sh" \
      --regions "$l1_regions" \
      --anchor "${locus_signals[0]}" \
      --signals "${locus_signals[*]}" \
      --labels "${locus_labels[*]}" \
      --out-prefix "${downstream_dir}/heatmaps/te_locus_best/L1" \
      --max-regions "$TE_LOCUS_MAX_REGIONS" \
      --bin-size "$TE_LOCUS_BIN_SIZE" \
      --threads "$THREADS"
fi

process_peak() {
  local peak="$1"
  sample="$(basename "$(dirname "$(dirname "$peak")")")"
  peak_type="$(basename "$(dirname "$peak")")"
  peak_label="${sample}.${peak_type}"

  if [[ -x "$R_BIN" ]]; then
    run_cached "peak annotation ${peak_label}" "${downstream_dir}/annotation/${sample}/${peak_label}.gene_structure.tsv" \
      "$R_BIN" "$ANNOTATION_SCRIPT" \
        --input "$peak" \
        --ref "$SPECIES" \
        --out-prefix "${downstream_dir}/annotation/${sample}/${peak_label}"
  else
    stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    record_status "optional:peak annotation ${peak_label}" SKIP "$stamp" "$stamp" "Rscript not found" false
  fi

  if [[ -x "$R_BIN" && -f "$PATHWAY_SCRIPT" ]]; then
    run_cached "pathway enrichment ${peak_label}" "${downstream_dir}/annotation/${sample}/${peak_label}.pathway.STATUS.tsv" \
      "$R_BIN" "$PATHWAY_SCRIPT" \
        --annotation "${downstream_dir}/annotation/${sample}/${peak_label}.gene_structure.tsv" \
        --out-prefix "${downstream_dir}/annotation/${sample}/${peak_label}.pathway" \
        --organism "$SPECIES"
  else
    stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    record_status "optional:pathway enrichment ${peak_label}" SKIP "$stamp" "$stamp" "Rscript or pathway script unavailable" false
  fi

  if [[ "$HOMER_AVAILABLE" == true ]]; then
    run_cached "HOMER ${peak_label}" "${downstream_dir}/homer/${sample}/${peak_type}/STATUS.tsv" \
      bash "${repo_root}/tools/cutrun_cli/cutrun_homer.sh" \
        --input "$peak" \
        --ref "$SPECIES" \
        --out-dir "${downstream_dir}/homer/${sample}/${peak_type}" \
        --threads "$THREADS"
  else
    stamp="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    record_status "optional:HOMER ${peak_label}" SKIP "$stamp" "$stamp" "findMotifsGenome.pl/annotatePeaks.pl not found" false
  fi

  if [[ -n "${standard_signals[*]-}" ]]; then
    run_cached "standard signal heatmap ${peak_label}" "${downstream_dir}/heatmaps/standard/${peak_label}.heatmap.pdf" \
      bash "${repo_root}/tools/cutrun_cli/cutrun_dt_heatmap.sh" \
        --regions "$peak" \
        --signals "${standard_signals[*]}" \
        --labels "${standard_labels[*]}" \
        --out-prefix "${downstream_dir}/heatmaps/standard/${peak_label}" \
        --max-regions "$PEAK_HEATMAP_MAX_REGIONS" \
        --threads "$THREADS"
  fi

  if [[ -n "${te_signals[*]-}" ]]; then
    run_cached "TE signal heatmap ${peak_label}" "${downstream_dir}/heatmaps/te/${peak_label}.heatmap.pdf" \
      bash "${repo_root}/tools/cutrun_cli/cutrun_dt_heatmap.sh" \
        --regions "$peak" \
        --signals "${te_signals[*]}" \
        --labels "${te_labels[*]}" \
        --out-prefix "${downstream_dir}/heatmaps/te/${peak_label}" \
        --max-regions "$PEAK_HEATMAP_MAX_REGIONS" \
        --threads "$THREADS"
  fi
}

peak_pids=()
for peak in "${peak_files[@]}"; do
  process_peak "$peak" &
  peak_pids+=("$!")
  if (( ${#peak_pids[@]} >= PEAK_JOBS )); then
    wait "${peak_pids[0]}" || true
    peak_pids=("${peak_pids[@]:1}")
  fi
done
for peak_pid in "${peak_pids[@]}"; do
  wait "$peak_pid" || true
done

summary_script="${repo_root}/tools/cutrun_cli/cutrun_results_summary.py"
summary_manifest="${RESULTS_DIR}/manifest/resolved_manifest.csv"
if [[ ! -s "$summary_manifest" ]]; then
  summary_manifest="${RESULTS_DIR}/_automation/inputs/manifest.csv"
fi
if [[ -f "$summary_script" && -s "$summary_manifest" ]]; then
  run_optional "results inventory and QC summary" \
    "$PYTHON_BIN" "$summary_script" \
      --results-dir "$RESULTS_DIR" \
        --manifest "$summary_manifest"
fi

if [[ -f "$TE_QC_PY" && -s "$summary_manifest" ]]; then
  run_optional "TE and classical QC metrics" \
    "$PYTHON_BIN" "$TE_QC_PY" \
      --results-dir "$RESULTS_DIR" \
      --manifest "$summary_manifest"
fi

if [[ -f "$TE_MULTIMAP_PY" && -d "${RESULTS_DIR}/04_te_bam" ]]; then
  run_optional "TE multimapping policy summary" "$PYTHON_BIN" "$TE_MULTIMAP_PY" \
    --bam-dir "${RESULTS_DIR}/04_te_bam" \
    --output "${downstream_dir}/te_multimap_sensitivity.tsv"
fi

if [[ -f "$RUN_MANIFEST_PY" && -s "$summary_manifest" ]]; then
  run_optional "final run manifest" "$PYTHON_BIN" "$RUN_MANIFEST_PY" \
    --results-dir "$RESULTS_DIR" \
    --manifest "$summary_manifest" \
    --config "${repo_root}/pipelines/chipseq_auto_nf/nextflow.config" \
    --run-id "$RUN_ID"
fi

if [[ -f "$STATUS_PY" ]]; then
  "$PYTHON_BIN" "$STATUS_PY" finalize \
    --status-file "$STATUS_FILE" \
    --run-id "$RUN_ID" \
    --output "${downstream_dir}/module_status_summary.json" || true
fi

if [[ -f "$REPORT_PY" ]]; then
  run_optional "unified HTML/PDF report" "$PYTHON_BIN" "$REPORT_PY" \
    --results-dir "$RESULTS_DIR"
fi

echo "[run_downstream] All available downstream analyses attempted: ${downstream_dir}"
failed_modules="$(awk -F '\t' 'NR > 1 && $4 == "FAIL" {n++} END {print n+0}' "$STATUS_FILE")"
echo "[run_downstream] Module status: ${STATUS_FILE} (failed=${failed_modules})"
