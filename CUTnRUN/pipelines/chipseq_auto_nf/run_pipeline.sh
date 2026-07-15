#!/usr/bin/env bash
set -euo pipefail

ORIG_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"
MAIN_NF="${PROJECT_ROOT}/main.nf"

MANIFEST=""
OUTDIR="${PWD}/chipseq_auto_results"
WORK_DIR=""
PROFILE="conda"
BACKGROUND="false"
BACKGROUND_INTERNAL="false"
RESUME="true"
DRY_RUN="false"
PREVIEW="false"
SKIP_PREFLIGHT="false"
TRIM="true"
RUN_ANALYSIS="true"
MAKE_RPGC_TRACK="true"
MAKE_TE_TRACK="true"
MAKE_TE_LOCUS_BEST_TRACK="true"
TE_K="25"
TE_DUPLICATE_POLICY="mark"
TE_REMOVE_BLACKLIST="false"
RUN_DOWNSTREAM="false"
TE_METHODS=""
TE_METHODS_EXECUTE="false"
TE_REPEAT_BED=""
TE_SAF=""
TE_ANNO=""
T3E_DIR=""
T3E_PYTHON=""
ALLO_COMMAND=""
REPENTOOLS_COMMAND=""
REPENTOOLS_INDEX_DIR=""
REPENTOOLS_GTF=""
REPENTOOLS_RET=""
REPENTOOLS_TARGET_GROUPS=""
REPENTOOLS_INPUT_SAMPLES=""
SPECIES="hg38"
RESOURCE_TIER="standard"
MIN_MAPQ="30"
MIN_FRAG="30"
MAX_FRAG="1200"
EXTEND_READS_SE="250"
IGNORE_FOR_NORMALIZATION="chrX"
DOWNLOAD_THREADS="8"
MACS_QVALUE=""
MACS_PVALUE=""
MACS_BROAD_CUTOFF=""
EXTRA_ARGS=()

usage() {
  cat <<'USAGE'
Usage:
  run_pipeline --manifest manifest.csv --outdir results [options]

Required:
  --manifest FILE

Core:
  --outdir DIR                  default ./chipseq_auto_results
  --work-dir DIR                default <outdir>/_automation/work
  --profile STR                 default conda
  --background
  --run-downstream              run downstream tools after Nextflow succeeds
  --te-methods LIST             optional adapters: t3e,allo,repentools
  --te-methods-execute          execute available external TE adapters
  --te-repeat-bed FILE          RepeatMasker BED for T3E/Allo
  --te-saf FILE                 TE SAF used to build RepeatMasker BED
  --te-anno FILE                TE annotation TSV used to build RepeatMasker BED
  --t3e-dir DIR                 T3E repository directory
  --t3e-python FILE             isolated T3E Python executable
  --t3e-max-bed-reads INT       deterministic T3E BED cap (default 1000000; 0=all)
  --t3e-iterations INT          T3E iterations (default 100)
  --allo-command CMD             reviewed Allo command template
  --repentools-command CMD       reviewed RepEnTools command template
  --repentools-index-dir DIR     CHM13 index directory for FASTQ adapter
  --repentools-gtf FILE          RepeatMasker GTF for FASTQ adapter
  --repentools-ret FILE          RepEnTools ret executable
  --repentools-target-groups SPEC  group=chip1,chip2;... mapping
  --repentools-input-samples SPEC  input1,input2 mapping
  --species STR                 hg38 | mm39, default hg38
  --resource-tier STR           downstream resources: small | standard | full
  --resume / --no-resume
  --dry-run                     print command only
  --preview                     nextflow -preview, checks graph without running jobs
  --skip-preflight              skip manifest/reference/FASTQ/disk checks (not recommended)

Pipeline params:
  --trim BOOL                   default true
  --run-analysis BOOL           default true
  --make-rpgc-track BOOL        default true
  --make-te-tracks BOOL         default true; make TE/L1 relaxed bigWig tracks
  --make-te-locus-best-track BOOL
                                default true; make reproducible one-best-location CPM track
  --te-k INT                    TE multi-mapping alignments, default 25
  --te-duplicate-policy STR     mark | keep | remove, default mark
  --te-remove-blacklist BOOL    default false
  --min-mapq INT                default 30
  --min-frag INT                default 30
  --max-frag INT                default 1200
  --extend-reads-se INT         default 250
  --ignore-for-normalization STR default chrX
  --download-threads INT        default 8
  --macs-qvalue FLOAT           MACS -q cutoff; lower=stricter, higher=more relaxed
  --macs-pvalue FLOAT           MACS -p cutoff; if set, MACS uses p-value instead of q-value
  --macs-broad-cutoff FLOAT     MACS --broad-cutoff; lower=stricter, higher=more relaxed
  --extra "ARGS"                append raw Nextflow args/params
USAGE
}

log_msg() { echo "[run_pipeline] $*"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --run-downstream) RUN_DOWNSTREAM="true"; shift 1 ;;
    --te-methods) TE_METHODS="$2"; shift 2 ;;
    --te-methods-execute) TE_METHODS_EXECUTE=true; shift ;;
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
    --species) SPECIES="$2"; shift 2 ;;
    --resource-tier) RESOURCE_TIER="$2"; shift 2 ;;
    --background-internal) BACKGROUND_INTERNAL="true"; BACKGROUND="false"; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --dry-run) DRY_RUN="true"; shift 1 ;;
    --preview) PREVIEW="true"; shift 1 ;;
    --skip-preflight) SKIP_PREFLIGHT="true"; shift 1 ;;
    --trim) TRIM="$2"; shift 2 ;;
    --run-analysis) RUN_ANALYSIS="$2"; shift 2 ;;
    --make-rpgc-track) MAKE_RPGC_TRACK="$2"; shift 2 ;;
    --make-te-tracks|--make-te-track) MAKE_TE_TRACK="$2"; shift 2 ;;
    --make-te-locus-best-track|--make-te-locus-best-tracks) MAKE_TE_LOCUS_BEST_TRACK="$2"; shift 2 ;;
    --te-k) TE_K="$2"; shift 2 ;;
    --te-duplicate-policy) TE_DUPLICATE_POLICY="$2"; shift 2 ;;
    --te-remove-blacklist) TE_REMOVE_BLACKLIST="$2"; shift 2 ;;
    --min-mapq) MIN_MAPQ="$2"; shift 2 ;;
    --min-frag) MIN_FRAG="$2"; shift 2 ;;
    --max-frag) MAX_FRAG="$2"; shift 2 ;;
    --extend-reads-se) EXTEND_READS_SE="$2"; shift 2 ;;
    --ignore-for-normalization) IGNORE_FOR_NORMALIZATION="$2"; shift 2 ;;
    --download-threads) DOWNLOAD_THREADS="$2"; shift 2 ;;
    --macs-qvalue) MACS_QVALUE="$2"; shift 2 ;;
    --macs-pvalue) MACS_PVALUE="$2"; shift 2 ;;
    --macs-broad-cutoff) MACS_BROAD_CUTOFF="$2"; shift 2 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[run_pipeline] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$MANIFEST" ]] || { echo "[run_pipeline] ERROR: --manifest is required" >&2; exit 1; }
[[ -f "$MANIFEST" ]] || { echo "[run_pipeline] ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
[[ -f "$MAIN_NF" ]] || { echo "[run_pipeline] ERROR: main.nf not found: $MAIN_NF" >&2; exit 1; }

# The pod exports a legacy Java 8 through JAVA_HOME/JAVA_CMD, while the
# chipseq environment contains the Java runtime required by Nextflow/Picard.
if [[ -z "${CHIPSEQ_ENV:-}" ]]; then
  for candidate in "${HOME}/.conda/envs/chipseq" "/path/to/.conda/envs/chipseq"; do
    if [[ -x "${candidate}/bin/java" ]]; then
      CHIPSEQ_ENV="$candidate"
      break
    fi
  done
fi
CHIPSEQ_ENV="${CHIPSEQ_ENV:-${HOME}/.conda/envs/chipseq}"
if [[ -x "${CHIPSEQ_ENV}/bin/java" ]]; then
  unset JAVA_CMD
  export JAVA_HOME="${CHIPSEQ_ENV}/lib/jvm"
  export PATH="${CHIPSEQ_ENV}/bin:${PATH}"
fi
command -v nextflow >/dev/null 2>&1 || { echo "[run_pipeline] ERROR: nextflow not found in PATH" >&2; exit 1; }
command -v python >/dev/null 2>&1 || { echo "[run_pipeline] ERROR: python not found in PATH" >&2; exit 1; }
JAVA_MAJOR="$(java -version 2>&1 | awk -F '[\".]' '/version/ {print $2; exit}')"
[[ "$JAVA_MAJOR" =~ ^[0-9]+$ && "$JAVA_MAJOR" -ge 17 ]] || {
  echo "[run_pipeline] ERROR: Java 17+ is required; detected: $(java -version 2>&1 | head -1)" >&2
  exit 1
}
[[ -n "$WORK_DIR" ]] || WORK_DIR="${OUTDIR}/_automation/work"

AUTOMATION_DIR="${OUTDIR}/_automation"
LOG_DIR="${AUTOMATION_DIR}/logs"
TS="$(date +%Y%m%d_%H%M%S)"
RUN_ID="${CUTRUN_RUN_ID:-run_${TS}Z}"
export CUTRUN_RUN_ID="$RUN_ID"
# Nextflow obtains a history.lock during startup.  The project directory is a
# network filesystem in the Pod and its advisory locks can block indefinitely;
# keep only the small launch/history state on local /tmp.  The task work
# directory and all published outputs remain under OUTDIR, so --resume still
# reuses the durable task cache.  Override NEXTFLOW_RUNTIME_ROOT when local
# scratch is unavailable.
NEXTFLOW_RUNTIME_ROOT="${NEXTFLOW_RUNTIME_ROOT:-/tmp/cutrun-nextflow}"
NXF_HOME_DIR="${NEXTFLOW_RUNTIME_ROOT}/home"
NF_LAUNCH_DIR="${NEXTFLOW_RUNTIME_ROOT}/launch_${TS}"
CONDA_CACHE_DIR="${AUTOMATION_DIR}/conda_cache"
TMP_DIR="${AUTOMATION_DIR}/tmp"
LOG_FILE="${LOG_DIR}/pipeline_${TS}.log"
NXF_LOG_FILE="${LOG_DIR}/nextflow_${TS}.log"
PID_FILE="${LOG_DIR}/pipeline_${TS}.pid"
mkdir -p "$OUTDIR" "$WORK_DIR" "$LOG_DIR" "$NXF_HOME_DIR" "$NF_LAUNCH_DIR" "$CONDA_CACHE_DIR" "$TMP_DIR"

if [[ "$BACKGROUND" == "true" && "$BACKGROUND_INTERNAL" != "true" && "${RUN_PIPELINE_BG_CHILD:-0}" != "1" ]]; then
  FILTERED_ARGS=()
  for x in "${ORIG_ARGS[@]}"; do
    [[ "$x" == "--background" ]] && continue
    [[ "$x" == "--background-internal" ]] && continue
    FILTERED_ARGS+=("$x")
  done
  log_msg "后台运行，日志: ${LOG_FILE}"
  nohup env RUN_PIPELINE_BG_CHILD=1 bash "$0" "${FILTERED_ARGS[@]}" --background-internal > "$LOG_FILE" 2>&1 &
  echo $! > "$PID_FILE"
  log_msg "PID: $(cat "$PID_FILE")"
  exit 0
fi

export NXF_HOME="$NXF_HOME_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
export TMPDIR="$TMP_DIR"
export TMP="$TMP_DIR"
export TEMP="$TMP_DIR"

CMD=(nextflow -log "$NXF_LOG_FILE" run "$MAIN_NF" -profile "$PROFILE" -work-dir "$WORK_DIR"
  --manifest "$MANIFEST"
  --outdir "$OUTDIR"
  --download_threads "$DOWNLOAD_THREADS"
  --trim "$TRIM"
  --run_analysis "$RUN_ANALYSIS"
  --make_rpgc_track "$MAKE_RPGC_TRACK"
  --make_te_tracks "$MAKE_TE_TRACK"
  --make_te_locus_best_track "$MAKE_TE_LOCUS_BEST_TRACK"
  --te_k "$TE_K"
  --te_duplicate_policy "$TE_DUPLICATE_POLICY"
  --te_remove_blacklist "$TE_REMOVE_BLACKLIST"
  --min_mapq "$MIN_MAPQ"
  --min_frag "$MIN_FRAG"
  --max_frag "$MAX_FRAG"
  --extend_reads_se "$EXTEND_READS_SE"
  --ignore_for_normalization "$IGNORE_FOR_NORMALIZATION")

if [[ -n "$MACS_QVALUE" ]]; then CMD+=(--macs_qvalue "$MACS_QVALUE"); fi
if [[ -n "$MACS_PVALUE" ]]; then CMD+=(--macs_pvalue "$MACS_PVALUE"); fi
if [[ -n "$MACS_BROAD_CUTOFF" ]]; then CMD+=(--macs_broad_cutoff "$MACS_BROAD_CUTOFF"); fi

if [[ "$RESUME" == "true" ]]; then CMD+=(-resume); fi
if [[ "$PREVIEW" == "true" ]]; then CMD+=(-preview); fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for x in "${EXTRA_ARGS[@]}"; do
    # shellcheck disable=SC2206
    EXTRA_SPLIT=( $x )
    CMD+=("${EXTRA_SPLIT[@]}")
  done
fi

log_msg "OUTDIR=$OUTDIR"
log_msg "WORK_DIR=$WORK_DIR"
log_msg "MANIFEST=$MANIFEST"
log_msg "RUN_ID=$RUN_ID"
log_msg "JAVA=$(java -version 2>&1 | head -1)"
log_msg "NEXTFLOW=$(nextflow -version 2>&1 | awk '/version/ {print $2; exit}')"
log_msg "CMD=${CMD[*]}"

PREFLIGHT_SCRIPT="${PROJECT_ROOT}/../../tools/cutrun_cli/cutrun_preflight.py"
if [[ "$SKIP_PREFLIGHT" != "true" ]]; then
  [[ -f "$PREFLIGHT_SCRIPT" ]] || {
    echo "[run_pipeline] ERROR: preflight script not found: $PREFLIGHT_SCRIPT" >&2
    exit 1
  }
  PREFLIGHT_REPORT="${AUTOMATION_DIR}/preflight_$(date +%Y%m%d_%H%M%S).json"
  PREFLIGHT_ARGS=(
    --manifest "$MANIFEST" \
    --species "$SPECIES" \
    --config "${PROJECT_ROOT}/nextflow.config" \
    --outdir "$OUTDIR" \
    --json-out "$PREFLIGHT_REPORT"
  )
  # A resume after a successful run should not spend several minutes
  # re-reading multi-gigabyte FASTQ files.  Structural/reference checks still
  # run; set CUTRUN_PREFLIGHT_RECHECK_GZIP=1 to force full integrity checks.
  if compgen -G "${AUTOMATION_DIR}/preflight_*.json" >/dev/null 2>&1 \
      && [[ "${CUTRUN_PREFLIGHT_RECHECK_GZIP:-0}" != "1" ]]; then
    PREFLIGHT_ARGS+=(--skip-gzip)
    log_msg "Existing preflight report found; skipping repeated gzip scan"
  fi
  if printf '%s' "$TE_METHODS" | tr ',' '\n' | awk 'tolower($0)=="repentools" || tolower($0)=="repen-tools"' | grep -q .; then
    [[ -n "$REPENTOOLS_INDEX_DIR" ]] && PREFLIGHT_ARGS+=(--repentools-index-dir "$REPENTOOLS_INDEX_DIR")
    [[ -n "$REPENTOOLS_GTF" ]] && PREFLIGHT_ARGS+=(--repentools-gtf "$REPENTOOLS_GTF")
  fi
  python "$PREFLIGHT_SCRIPT" "${PREFLIGHT_ARGS[@]}"
else
  log_msg "WARNING: preflight skipped by --skip-preflight"
fi

if [[ "$DRY_RUN" == "true" ]]; then
  exit 0
fi

# Preserve the complete wrapper invocation (including downstream/TE adapter
# switches) alongside the generated Nextflow command.
export CUTRUN_WRAPPER_ARGS="$(printf '%q ' "${ORIG_ARGS[@]}")"
python "${PROJECT_ROOT}/bin/write_run_metadata.py" \
  --output "${AUTOMATION_DIR}/run_metadata_${TS}.json" \
  --run-id "$RUN_ID" \
  --manifest "$MANIFEST" \
  --main-nf "$MAIN_NF" \
  --nextflow-config "${PROJECT_ROOT}/nextflow.config" \
  -- "${CMD[@]}"

cd "$NF_LAUNCH_DIR"
"${CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if [[ "$RUN_DOWNSTREAM" == "true" && "$PREVIEW" != "true" ]]; then
  log_msg "Starting downstream analyses"
  downstream_cmd=(bash "${SCRIPT_DIR}/run_downstream.sh" \
    --results-dir "$OUTDIR" \
    --species "$SPECIES" \
    --resource-tier "$RESOURCE_TIER" \
    --threads 8)
  if [[ "$RESUME" == true ]]; then downstream_cmd+=(--resume); else downstream_cmd+=(--no-resume); fi
  export CUTRUN_RUN_ID="$RUN_ID"
  if [[ -n "$TE_METHODS" ]]; then
    downstream_cmd+=(--te-methods "$TE_METHODS")
  fi
  if [[ "$TE_METHODS_EXECUTE" == true ]]; then
    downstream_cmd+=(--te-methods-execute)
  fi
  [[ -n "$TE_REPEAT_BED" ]] && downstream_cmd+=(--te-repeat-bed "$TE_REPEAT_BED")
  [[ -n "$TE_SAF" ]] && downstream_cmd+=(--te-saf "$TE_SAF")
  [[ -n "$TE_ANNO" ]] && downstream_cmd+=(--te-anno "$TE_ANNO")
  [[ -n "$T3E_DIR" ]] && downstream_cmd+=(--t3e-dir "$T3E_DIR")
  [[ -n "$T3E_PYTHON" ]] && downstream_cmd+=(--t3e-python "$T3E_PYTHON")
  downstream_cmd+=(--t3e-max-bed-reads "${T3E_MAX_BED_READS:-1000000}" --t3e-iterations "${T3E_ITERATIONS:-100}")
  [[ -n "$ALLO_COMMAND" ]] && downstream_cmd+=(--allo-command "$ALLO_COMMAND")
  [[ -n "$REPENTOOLS_COMMAND" ]] && downstream_cmd+=(--repentools-command "$REPENTOOLS_COMMAND")
  [[ -n "$REPENTOOLS_INDEX_DIR" ]] && downstream_cmd+=(--repentools-index-dir "$REPENTOOLS_INDEX_DIR")
  [[ -n "$REPENTOOLS_GTF" ]] && downstream_cmd+=(--repentools-gtf "$REPENTOOLS_GTF")
  [[ -n "$REPENTOOLS_RET" ]] && downstream_cmd+=(--repentools-ret "$REPENTOOLS_RET")
  [[ -n "$REPENTOOLS_TARGET_GROUPS" ]] && downstream_cmd+=(--repentools-target-groups "$REPENTOOLS_TARGET_GROUPS")
  [[ -n "$REPENTOOLS_INPUT_SAMPLES" ]] && downstream_cmd+=(--repentools-input-samples "$REPENTOOLS_INPUT_SAMPLES")
  "${downstream_cmd[@]}"
elif [[ "$RUN_DOWNSTREAM" == "true" ]]; then
  log_msg "Preview mode: downstream analyses skipped"
fi
