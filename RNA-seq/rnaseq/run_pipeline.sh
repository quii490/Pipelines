#!/usr/bin/env bash
set -euo pipefail

# ==============================
# RNA-seq pipeline launcher
# - configurable output/work/log dirs
# - safe background mode
# - resume / no-resume
# - disables Nextflow ANSI TTY output by default
# ==============================
export NXF_ANSI_LOG=false
export NXF_OPTS="-Djdk.lang.Process.launchMechanism=FORK"

# Use the user-installed Temurin JDK for Nextflow and Picard.  This avoids
# conda-bundled Java/NFS lock issues and keeps the JVM version consistent.
TEMURIN_JDK=/path/to/softwares/jdk-17.0.8.1+1
if [ -x "$TEMURIN_JDK/bin/java" ]; then
  export JAVA_HOME="$TEMURIN_JDK"
  export PATH="$TEMURIN_JDK/bin:$PATH"
  unset JAVA_LD_LIBRARY_PATH
  export LD_LIBRARY_PATH="$TEMURIN_JDK/lib/server:$TEMURIN_JDK/lib:${LD_LIBRARY_PATH:-}"
fi
# ---------- defaults ----------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_DIR="$SCRIPT_DIR"
MAIN_NF="$PIPELINE_DIR/main.nf"
MODULE_DEFAULTS="${PIPELINE_DIR}/nextflow.config"
DEFAULT_READER="${PIPELINE_DIR}/scripts/read_module_default.py"
module_default() { python3 "${DEFAULT_READER}" "${MODULE_DEFAULTS}" "$1" "$2"; }
OUTDIR="${PWD}/results"
WORK_DIR="${PWD}/work"
LOG_DIR="${OUTDIR}/logs"
RUN_TS="$(date +%Y%m%d_%H%M%S)"
LOG_FILE="${LOG_DIR}/pipeline_${RUN_TS}.log"
AUTOMATION_DIR="${OUTDIR}/_automation"
NXF_HOME_DIR="${NXF_HOME:-/tmp/${USER:-$(id -un)}/nextflow/home}"
NF_LAUNCH_DIR="${NF_LAUNCH_DIR:-/tmp/${USER:-$(id -un)}/nextflow/launch/${RUN_NAME:-default}}"
CONDA_CACHE_DIR="${AUTOMATION_DIR}/conda_cache"
TMP_DIR="${AUTOMATION_DIR}/tmp"
NXF_LOG_FILE="${LOG_DIR}/nextflow_${RUN_TS}.log"
NEXTFLOW_PID_FILE="${LOG_DIR}/nextflow_${RUN_TS}.pid"
NXF_TRACE_FILE="${LOG_DIR}/nextflow_${RUN_TS}.trace.tsv"
NXF_REPORT_FILE="${LOG_DIR}/nextflow_${RUN_TS}.report.html"
NXF_TIMELINE_FILE="${LOG_DIR}/nextflow_${RUN_TS}.timeline.html"
PROGRESS_SCRIPT="${PIPELINE_DIR}/scripts/rnaseq_progress.sh"
PROGRESS_FILE="${LOG_DIR}/progress_${RUN_TS}.txt"
PROGRESS_PID_FILE="${LOG_DIR}/progress_${RUN_TS}.pid"
PROGRESS_INTERVAL=30
USER_SET_WORK_DIR="false"
USER_SET_NXF_HOME="false"
USER_SET_NF_LAUNCH="false"
USER_SET_CONDA_CACHE="false"
USER_SET_TMP_DIR="false"

SPECIES="hg38"
HUMAN_REF="hg38"
ALIGNER="star"
STRANDEDNESS="reverse"
SAMPLESHEET=""

RESUME="true"
RESUME_SESSION=""
BACKGROUND="false"

RUN_FASTP="$(module_default run_fastp true)"
RUN_STAR_FC="$(module_default run_star_fc true)"
RUN_SALMON="$(module_default run_salmon true)"
RUN_STRINGTIE="$(module_default run_stringtie false)"
RUN_TECOUNT="$(module_default run_tecount false)"
RUN_TELOCAL="$(module_default run_telocal false)"
RUN_TETRANSCRIPTS="$(module_default run_tetranscripts true)"
RUN_DEDUP="$(module_default run_dedup false)"
RUN_MARKDUP_QC="$(module_default run_markdup_qc true)"
RUN_SALMONTE="$(module_default run_salmonte false)"
RUN_MULTIQC="$(module_default run_multiqc true)"
MULTIQC_CMD="$(module_default multiqc_cmd /path/to/.conda/envs/emseq/bin/multiqc)"
# Reserved compatibility switches. They are forwarded by run_auto_rnaseq.sh;
# explicit defaults prevent strict-shell failures while their workflow branches
# are being maintained separately.
RUN_TELESCOPE="$(module_default run_telescope true)"
RUN_FASTQC="$(module_default run_fastqc true)"
RUN_RNASEQ_METRICS="$(module_default run_rnaseq_metrics true)"
RUN_REDISCOVERTE_ROLLUP="$(module_default run_rediscoverte_rollup true)"
REF_FLAT=""
RIBOSOMAL_INTERVALS=""
REDISCOVERTE_ROLLUP_CONDA_PREFIX="${REDISCOVERTE_ROLLUP_CONDA_PREFIX:-/path/to/.conda/envs/downstream}"
RNASEQ_METRICS_MEM="${RNASEQ_METRICS_MEM:-12 GB}"
RNASEQ_METRICS_JAVA_HEAP="${RNASEQ_METRICS_JAVA_HEAP:-8g}"
PICARD_MARKDUP_JAVA_HEAP="${PICARD_MARKDUP_JAVA_HEAP:-16g}"
RUN_SESSION_ISOLATED="${RUN_SESSION_ISOLATED:-false}"

# REdiscoverTE: default auto behavior
# hg38 -> true ; mm10/mm39 -> false
RUN_REDISCOVERTE="auto"
DRY_RUN="false"
FAILURE_POLICY="core"

ts_now() {
  date +"%Y-%m-%d %H:%M:%S"
}

log_msg() {
  echo "[$(ts_now)] [run_pipeline] $*"
}

# Nextflow may leave a Conda child behind after a failed environment create.
# Keep cleanup scoped to this run's process tree and cache directory.
terminate_descendants() {
  local root="$1"
  local child
  for child in $(ps --no-headers --ppid "$root" -o pid= 2>/dev/null || true); do
    terminate_descendants "$child"
    kill -TERM "$child" 2>/dev/null || true
  done
}

terminate_process_group() {
  local root="$1"
  local pgid
  pgid="$(ps -o pgid= -p "$root" 2>/dev/null | tr -d ' ' || true)"
  if [[ "$RUN_SESSION_ISOLATED" == "true" && "$pgid" =~ ^[0-9]+$ && "$pgid" != "0" ]]; then
    trap - TERM INT HUP
    kill -TERM -- "-${pgid}" 2>/dev/null || true
    for _ in {1..20}; do
      kill -0 -- "-${pgid}" 2>/dev/null || return 0
      sleep 0.25
    done
    kill -KILL -- "-${pgid}" 2>/dev/null || true
  else
    terminate_descendants "$root"
  fi
}

terminate_conda_for_cache() {
  local pid
  local user_name="${USER:-$(id -un)}"
  for pid in $(pgrep -u "$user_name" -f "conda env create --prefix ${CONDA_CACHE_DIR}" 2>/dev/null || true); do
    [[ "$pid" == "$$" ]] && continue
    kill -TERM "$pid" 2>/dev/null || true
  done
  sleep 1
  for pid in $(pgrep -u "$user_name" -f "conda env create --prefix ${CONDA_CACHE_DIR}" 2>/dev/null || true); do
    [[ "$pid" == "$$" ]] && continue
    kill -KILL "$pid" 2>/dev/null || true
  done
}

handle_pipeline_signal() {
  local signal_code=143
  log_msg "收到终止信号，清理 Nextflow 子进程和本运行的 Conda 创建进程。"
  terminate_process_group "$$"
  terminate_conda_for_cache
  exit "$signal_code"
}

append_engine_log() {
  local rc="$1"
  if [[ "${LOG_FILE}" == "${NXF_LOG_FILE}" ]]; then
    return 0
  fi
  if [[ ! -f "${NXF_LOG_FILE}" ]]; then
    return 0
  fi
  if grep -q "BEGIN NEXTFLOW ENGINE LOG" "${LOG_FILE}" 2>/dev/null; then
    return 0
  fi
  {
    echo
    echo "========== BEGIN NEXTFLOW ENGINE LOG: ${NXF_LOG_FILE} =========="
    cat "${NXF_LOG_FILE}"
    echo "========== END NEXTFLOW ENGINE LOG: exit_code=${rc} =========="
  } >> "${LOG_FILE}"
}

start_progress_monitor() {
  local watched_pid="$1"
  [[ -f "${PROGRESS_SCRIPT}" ]] || return 0

  if [[ "${BACKGROUND}" == "true" ]]; then
    nohup bash "${PROGRESS_SCRIPT}" \
      --results-dir "${OUTDIR}" \
      --trace "${NXF_TRACE_FILE}" \
      --engine-log "${NXF_LOG_FILE}" \
      --pid "${watched_pid}" \
      --output "${PROGRESS_FILE}" \
      --watch \
      --interval "${PROGRESS_INTERVAL}" \
      > /dev/null 2>&1 < /dev/null &
  else
    bash "${PROGRESS_SCRIPT}" \
      --results-dir "${OUTDIR}" \
      --trace "${NXF_TRACE_FILE}" \
      --engine-log "${NXF_LOG_FILE}" \
      --pid "${watched_pid}" \
      --output "${PROGRESS_FILE}" \
      --watch \
      --interval "${PROGRESS_INTERVAL}" \
      > /dev/null 2>&1 < /dev/null &
  fi
  echo "$!" > "${PROGRESS_PID_FILE}"
}

# optional runtime controls
PROFILE=""
QUEUE_SIZE=""
MAX_CPUS=""
MAX_MEMORY=""
EXTRA_ARGS=()

usage() {
  cat <<'USAGE'
Usage:
  bash run_pipeline.sh --samplesheet samplesheet.csv [options]

Required:
  --samplesheet FILE          CSV with columns: sample,layout,condition,replicate,r1,r2

Core options:
  --pipeline-dir DIR          Pipeline directory (default: script directory)
  --outdir DIR                Results directory (default: ./results)
  --work-dir DIR              Nextflow work directory (default: ./work)
  --log-dir DIR               Log directory (default: <outdir>/logs)
  --log-file FILE             Pipeline log file (default: <log-dir>/pipeline_YYYYmmdd_HHMMSS.log)
  --species STR               hg38 | mm10 | mm39
  --human-ref STR             hg38 | t2t，默认 hg38，仅 species=hg38 生效
  --aligner STR               star | hisat2，默认 star
  --strand STR                unstranded | forward | reverse
  --strandedness STR          deprecated alias of --strand

Run control:
  --background                Run with nohup in background
  --resume                    Enable -resume (default)
  --no-resume                 Disable -resume
  --resume-session ID         Resume a specific Nextflow session ID
  --profile NAME              Nextflow profile name
  --queue-size INT            Override params.queue_size
  --max-cpus INT              Override params.max_cpus
  --max-memory STR            Override params.max_memory (example: "300 GB")
  --progress-interval INT     Seconds between progress updates (default: 30)
  --dry-run                  Preview the resolved Nextflow workflow without running tasks
  --failure-policy STR       core|strict; default core, optional modules may fail without stopping core outputs

Pipeline switches:
  --run-fastp BOOL            true|false
  --run-gene-count-branch BOOL          true|false
  --run-salmon BOOL           true|false
  --run-stringtie BOOL        true|false
  --run-tecount BOOL          true|false
  --run-telocal BOOL          true|false
  --run-tetranscripts BOOL    true|false
  --run-dedup BOOL            true|false  Advanced: use deduplicated BAMs for counts; not recommended for ordinary RNA-seq
  --run-markdup-qc BOOL       true|false  Picard duplicate metrics only; does not change counts
  --run-rediscoverte BOOL     true|false|auto  (default: auto)
  --run-rediscoverte-rollup BOOL  true|false; default true (required for layered REdiscoverTE plots)
  --run-salmonte BOOL         true|false
  --run-multiqc BOOL          true|false
  --multiqc-cmd PATH          MultiQC executable (default: configured multiqc_cmd)
  --run-telescope BOOL       true|false
  --run-fastqc BOOL          true|false
  --run-rnaseq-metrics BOOL  true|false; default true
  --ref-flat FILE            Optional Picard refFlat; generated from gene GTF if omitted
  --refFlat FILE             Alias of --ref-flat
  --ribosomal-intervals FILE Optional Picard ribosomal interval list
  --rediscoverte-rollup-conda-prefix DIR
                              Existing environment for REdiscoverTE rollup
  --rnaseq-metrics-mem STR    Picard task memory (default: 12 GB)
  --rnaseq-metrics-java-heap STR
                              Picard JVM heap (default: 8g)
  --picard-markdup-java-heap STR
                              MarkDuplicates JVM heap (default: 16g; task memory: 24 GB)

Other:
  --extra "ARGS"              Extra raw args appended to nextflow run
  -h, --help                  Show help

Examples:
  bash run_pipeline.sh \
    --samplesheet samplesheet.csv \
    --species hg38 \
    --outdir results

  bash run_pipeline.sh \
    --samplesheet samplesheet.csv \
    --species hg38 \
    --outdir results \
    --background

  bash run_pipeline.sh \
    --samplesheet samplesheet.csv \
    --species mm10 \
    --outdir results \
    --no-resume \
    --max-cpus 30 \
    --max-memory "300 GB"
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --pipeline-dir) PIPELINE_DIR="$2"; MAIN_NF="$2/main.nf"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; USER_SET_WORK_DIR="true"; shift 2 ;;
    --log-dir) LOG_DIR="$2"; shift 2 ;;
    --log-file) LOG_FILE="$2"; shift 2 ;;
    --nxf-home) NXF_HOME_DIR="$2"; USER_SET_NXF_HOME="true"; shift 2 ;;
    --nf-launch-dir) NF_LAUNCH_DIR="$2"; USER_SET_NF_LAUNCH="true"; shift 2 ;;
    --conda-cache-dir) CONDA_CACHE_DIR="$2"; USER_SET_CONDA_CACHE="true"; shift 2 ;;
    --tmp-dir) TMP_DIR="$2"; USER_SET_TMP_DIR="true"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --human-ref) HUMAN_REF="$2"; shift 2 ;;
    --aligner) ALIGNER="$2"; shift 2 ;;
    --strand|--strandedness) STRANDEDNESS="$2"; shift 2 ;;

    --background) BACKGROUND="true"; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --resume-session) RESUME_SESSION="$2"; RESUME="true"; shift 2 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --queue-size) QUEUE_SIZE="$2"; shift 2 ;;
    --max-cpus) MAX_CPUS="$2"; shift 2 ;;
    --max-memory) MAX_MEMORY="$2"; shift 2 ;;
    --progress-interval) PROGRESS_INTERVAL="$2"; shift 2 ;;
    --dry-run) DRY_RUN="true"; shift ;;
    --failure-policy) FAILURE_POLICY="$2"; shift 2 ;;

    --run-fastp) RUN_FASTP="$2"; shift 2 ;;
    --run-gene-count-branch) RUN_STAR_FC="$2"; shift 2 ;;
    --run-salmon) RUN_SALMON="$2"; shift 2 ;;
    --run-stringtie) RUN_STRINGTIE="$2"; shift 2 ;;
    --run-tecount) RUN_TECOUNT="$2"; shift 2 ;;
    --run-telocal) RUN_TELOCAL="$2"; shift 2 ;;
    --run-tetranscripts) RUN_TETRANSCRIPTS="$2"; shift 2 ;;
    --run-dedup) RUN_DEDUP="$2"; shift 2 ;;
    --run-markdup-qc) RUN_MARKDUP_QC="$2"; shift 2 ;;
    --run-rediscoverte) RUN_REDISCOVERTE="$2"; shift 2 ;;
    --run-rediscoverte-rollup) RUN_REDISCOVERTE_ROLLUP="$2"; shift 2 ;;
    --run-salmonte) RUN_SALMONTE="$2"; shift 2 ;;
    --run-multiqc) RUN_MULTIQC="$2"; shift 2 ;;
    --multiqc-cmd) MULTIQC_CMD="$2"; shift 2 ;;
    --run-telescope) RUN_TELESCOPE="$2"; shift 2 ;;
    --run-fastqc) RUN_FASTQC="$2"; shift 2 ;;
    --run-rnaseq-metrics) RUN_RNASEQ_METRICS="$2"; shift 2 ;;
    --ref-flat|--refFlat) REF_FLAT="$2"; shift 2 ;;
    --ribosomal-intervals|--ribosomalIntervals) RIBOSOMAL_INTERVALS="$2"; shift 2 ;;
    --rediscoverte-rollup-conda-prefix) REDISCOVERTE_ROLLUP_CONDA_PREFIX="$2"; shift 2 ;;
    --rnaseq-metrics-mem) RNASEQ_METRICS_MEM="$2"; shift 2 ;;
    --rnaseq-metrics-java-heap) RNASEQ_METRICS_JAVA_HEAP="$2"; shift 2 ;;
    --picard-markdup-java-heap) PICARD_MARKDUP_JAVA_HEAP="$2"; shift 2 ;;

    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ "$DRY_RUN" == "true" ]]; then
  BACKGROUND="false"
  RESUME="false"
fi
case "$FAILURE_POLICY" in core|strict) ;; *) echo "Error: --failure-policy must be core or strict" >&2; exit 1 ;; esac

# ---------- normalize paths ----------
[[ -n "$SAMPLESHEET" ]] || { echo "Error: --samplesheet is required" >&2; usage; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "Error: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }
[[ -f "$MAIN_NF" ]] || { echo "Error: main.nf not found: $MAIN_NF" >&2; exit 1; }
command -v nextflow >/dev/null 2>&1 || { echo "Error: nextflow not found in PATH" >&2; exit 1; }

# Make log defaults follow outdir unless user explicitly changed them after --outdir
if [[ "$LOG_DIR" == "${OUTDIR}/logs" || "$LOG_DIR" == "${PWD}/results/logs" ]]; then
  LOG_DIR="${OUTDIR}/logs"
fi
if [[ "$LOG_FILE" == "${LOG_DIR}/pipeline.log" || "$LOG_FILE" == "${PWD}/results/logs/pipeline.log" || "$LOG_FILE" == "${LOG_DIR}/pipeline_${RUN_TS}.log" || "$LOG_FILE" == "${PWD}/results/logs/pipeline_${RUN_TS}.log" ]]; then
  LOG_FILE="${LOG_DIR}/pipeline_${RUN_TS}.log"
fi
AUTOMATION_DIR="${OUTDIR}/_automation"
[[ "${USER_SET_WORK_DIR}" == "true" ]] || WORK_DIR="${AUTOMATION_DIR}/work"
[[ "${USER_SET_NXF_HOME}" == "true" ]] || NXF_HOME_DIR="${NXF_HOME:-/tmp/${USER:-$(id -un)}/nextflow/home}"
if [[ "${USER_SET_NF_LAUNCH}" != "true" ]]; then
  # Keep Nextflow's launcher/history lock on local scratch.  The result and
  # work directories may be on NFS, but FileChannel.lock() in the launch
  # directory can block indefinitely on some NFS mounts.
  NF_LAUNCH_DIR="/tmp/${USER:-$(id -un)}/nextflow/launch/$(basename "${OUTDIR}")"
fi
[[ "${USER_SET_CONDA_CACHE}" == "true" ]] || CONDA_CACHE_DIR="${AUTOMATION_DIR}/conda_cache"
[[ "${USER_SET_TMP_DIR}" == "true" ]] || TMP_DIR="${AUTOMATION_DIR}/tmp"
NXF_LOG_FILE="${LOG_DIR}/nextflow_${RUN_TS}.log"
NEXTFLOW_PID_FILE="${LOG_DIR}/nextflow_${RUN_TS}.pid"
NXF_TRACE_FILE="${LOG_DIR}/nextflow_${RUN_TS}.trace.tsv"
NXF_REPORT_FILE="${LOG_DIR}/nextflow_${RUN_TS}.report.html"
NXF_TIMELINE_FILE="${LOG_DIR}/nextflow_${RUN_TS}.timeline.html"
PROGRESS_SCRIPT="${PIPELINE_DIR}/scripts/rnaseq_progress.sh"
PROGRESS_FILE="${LOG_DIR}/progress_${RUN_TS}.txt"
PROGRESS_PID_FILE="${LOG_DIR}/progress_${RUN_TS}.pid"

[[ "${PROGRESS_INTERVAL}" =~ ^[0-9]+$ ]] && (( PROGRESS_INTERVAL > 0 )) || {
  echo "Error: --progress-interval must be a positive integer" >&2
  exit 1
}

mkdir -p "$OUTDIR" "$WORK_DIR" "$LOG_DIR" "$NXF_HOME_DIR" "$NF_LAUNCH_DIR" "$CONDA_CACHE_DIR" "$TMP_DIR"

log_msg "RUN_TS=$RUN_TS"
log_msg "NXF_HOME=$NXF_HOME_DIR"
log_msg "NF_LAUNCH_DIR=$NF_LAUNCH_DIR"
log_msg "CONDA_CACHE_DIR=$CONDA_CACHE_DIR"
log_msg "TMP_DIR=$TMP_DIR"
log_msg "NXF_LOG_FILE=$NXF_LOG_FILE"
log_msg "NEXTFLOW_PID_FILE=$NEXTFLOW_PID_FILE"
log_msg "NXF_TRACE_FILE=$NXF_TRACE_FILE"
log_msg "NXF_REPORT_FILE=$NXF_REPORT_FILE"
log_msg "NXF_TIMELINE_FILE=$NXF_TIMELINE_FILE"
log_msg "PROGRESS_FILE=$PROGRESS_FILE"

# ---------- validate / auto-resolve ----------
SPECIES="$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]')"
case "$SPECIES" in
  hg38|mm10|mm39) ;;
  *) echo "Error: --species must be hg38 | mm10 | mm39" >&2; exit 1 ;;
esac
HUMAN_REF="$(echo "$HUMAN_REF" | tr '[:upper:]' '[:lower:]')"
case "$HUMAN_REF" in
  hg38|t2t) ;;
  *) echo "Error: --human-ref must be hg38 | t2t" >&2; exit 1 ;;
esac
if [[ "$SPECIES" != "hg38" && "$HUMAN_REF" != "hg38" ]]; then
  echo "Error: --human-ref t2t only supports --species hg38" >&2
  exit 1
fi
ALIGNER="$(echo "$ALIGNER" | tr '[:upper:]' '[:lower:]')"
case "$ALIGNER" in
  star|hisat2) ;;
  *) echo "Error: --aligner must be star | hisat2" >&2; exit 1 ;;
esac

STRANDEDNESS="$(echo "$STRANDEDNESS" | tr '[:upper:]' '[:lower:]')"
case "$STRANDEDNESS" in
  unstranded|forward|reverse) ;;
  *) echo "Error: --strand must be unstranded | forward | reverse" >&2; exit 1 ;;
esac

RUN_NAME="${RUN_NAME:-rnaseq_${SPECIES}_${RUN_TS}}"
RUN_REDISCOVERTE="$(echo "$RUN_REDISCOVERTE" | tr '[:upper:]' '[:lower:]')"
if [[ "$RUN_REDISCOVERTE" == "auto" ]]; then
  if [[ "$SPECIES" == "hg38" && "$HUMAN_REF" == "hg38" ]]; then
    RUN_REDISCOVERTE="true"
  else
    RUN_REDISCOVERTE="false"
  fi
fi

if [[ "${RUN_REDISCOVERTE_ROLLUP,,}" == "true" && ! -d "${REDISCOVERTE_ROLLUP_CONDA_PREFIX}" ]]; then
  echo "Error: REdiscoverTE rollup environment not found: ${REDISCOVERTE_ROLLUP_CONDA_PREFIX}" >&2
  echo "       Use --rediscoverte-rollup-conda-prefix DIR or disable --run-rediscoverte-rollup false." >&2
  exit 1
fi

trap 'handle_pipeline_signal' INT TERM HUP
export NXF_HOME="$NXF_HOME_DIR"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
export TMPDIR="$TMP_DIR"
export TMP="$TMP_DIR"
export TEMP="$TMP_DIR"
# `-resume` needs Nextflow history.  Older wrappers disabled it by default,
# which makes Nextflow 26 abort before workflow parsing.
unset NXF_IGNORE_RESUME_HISTORY
export NXF_DISABLE_CHECK_LATEST="${NXF_DISABLE_CHECK_LATEST:-true}"
# ---------- nextflow command ----------

CMD=(nextflow -log "$NXF_LOG_FILE" run "$MAIN_NF" -name "$RUN_NAME"
  -with-trace "$NXF_TRACE_FILE"
  -with-report "$NXF_REPORT_FILE"
  -with-timeline "$NXF_TIMELINE_FILE"
  -work-dir "$WORK_DIR"
  --samplesheet "$SAMPLESHEET"
  --outdir "$OUTDIR"
  --conda_cache_dir "$CONDA_CACHE_DIR"
  --species "$SPECIES"
  --human_ref "$HUMAN_REF"
  --aligner "$ALIGNER"
  --strandedness "$STRANDEDNESS"
  --failure_policy "$FAILURE_POLICY"
  --run_fastp "$RUN_FASTP"
  --run_star_fc "$RUN_STAR_FC"
  --run_salmon "$RUN_SALMON"
  --run_stringtie "$RUN_STRINGTIE"
  --run_tecount "$RUN_TECOUNT"
  --run_telocal "$RUN_TELOCAL"
  --run_tetranscripts "$RUN_TETRANSCRIPTS"
  --run_dedup "$RUN_DEDUP"
  --run_markdup_qc "$RUN_MARKDUP_QC"
  --run_rediscoverte "$RUN_REDISCOVERTE"
  --run_rediscoverte_rollup "$RUN_REDISCOVERTE_ROLLUP"
  --rediscoverte_rollup_conda_prefix "$REDISCOVERTE_ROLLUP_CONDA_PREFIX"
  --rnaseq_metrics_mem "$RNASEQ_METRICS_MEM"
  --rnaseq_metrics_java_heap "$RNASEQ_METRICS_JAVA_HEAP"
  --picard_markdup_java_heap "$PICARD_MARKDUP_JAVA_HEAP"
  --run_salmonte "$RUN_SALMONTE"
  --run_telescope "$RUN_TELESCOPE"
  --run_fastqc "$RUN_FASTQC"
  --run_rnaseq_metrics "$RUN_RNASEQ_METRICS"
  --run_multiqc "$RUN_MULTIQC"
  --multiqc_cmd "$MULTIQC_CMD")

if [[ "$DRY_RUN" == "true" ]]; then
  CMD+=(-preview)
fi

[[ -n "$REF_FLAT" ]] && CMD+=(--ref_flat "$REF_FLAT")
[[ -n "$RIBOSOMAL_INTERVALS" ]] && CMD+=(--ribosomal_intervals "$RIBOSOMAL_INTERVALS")

if [[ "$RESUME" == "true" ]]; then
  if [[ -n "$RESUME_SESSION" ]]; then
    CMD+=(-resume "$RESUME_SESSION")
  else
    CMD+=(-resume)
  fi
fi

if [[ -n "$PROFILE" ]]; then
  CMD+=(-profile "$PROFILE")
fi

if [[ -n "$QUEUE_SIZE" ]]; then
  CMD+=(--queue_size "$QUEUE_SIZE")
fi

if [[ -n "$MAX_CPUS" ]]; then
  CMD+=(--max_cpus "$MAX_CPUS")
fi

if [[ -n "$MAX_MEMORY" ]]; then
  CMD+=(--max_memory "$MAX_MEMORY")
fi

if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for x in "${EXTRA_ARGS[@]}"; do
    # shellcheck disable=SC2206
    EXTRA_SPLIT=( $x )
    CMD+=("${EXTRA_SPLIT[@]}")
  done
fi

# ---------- print summary ----------
echo "========================================"
echo "Pipeline dir    : $PIPELINE_DIR"
echo "Main script     : $MAIN_NF"
echo "Samplesheet     : $SAMPLESHEET"
echo "Species         : $SPECIES"
echo "Human ref       : $HUMAN_REF"
echo "Gene aligner    : $ALIGNER"
echo "Strandedness    : $STRANDEDNESS"
echo "Outdir          : $OUTDIR"
echo "Work dir        : $WORK_DIR"
echo "Log file        : $LOG_FILE"
echo "Background      : $BACKGROUND"
echo "Resume          : $RESUME"
echo "Dry run         : $DRY_RUN"
echo "Failure policy  : $FAILURE_POLICY"
echo "TEtranscripts   : $RUN_TETRANSCRIPTS"
echo "Dedup           : $RUN_DEDUP"
echo "REdiscoverTE    : $RUN_REDISCOVERTE"
echo "REdiscoverTE rollup env: $REDISCOVERTE_ROLLUP_CONDA_PREFIX"
echo "RNA-seq metrics : mem=$RNASEQ_METRICS_MEM, java_heap=$RNASEQ_METRICS_JAVA_HEAP"
[[ -n "$PROFILE" ]] && echo "Profile         : $PROFILE"
[[ -n "$QUEUE_SIZE" ]] && echo "Queue size      : $QUEUE_SIZE"
[[ -n "$MAX_CPUS" ]] && echo "Max CPUs        : $MAX_CPUS"
[[ -n "$MAX_MEMORY" ]] && echo "Max memory      : $MAX_MEMORY"
echo "----------------------------------------"
printf 'Command: '
printf '%q ' "${CMD[@]}"
printf '\n'
echo "========================================"

# ---------- run ----------
# Disable interactive ANSI terminal behavior to avoid "Stopped (tty output)"
export NXF_ANSI_LOG=false
CMD_STR=$(printf '%q ' "${CMD[@]}")
python3 "${PIPELINE_DIR}/scripts/write_run_manifest.py" \
  --out "${AUTOMATION_DIR}/run_manifest_${RUN_TS}.json" \
  --repo "${PIPELINE_DIR}/.." \
  --samplesheet "${SAMPLESHEET}" \
  --command "${CMD_STR}" \
  --param "species=${SPECIES}" \
  --param "human_ref=${HUMAN_REF}" \
  --param "aligner=${ALIGNER}" \
  --param "strand=${STRANDEDNESS}" \
  --param "failure_policy=${FAILURE_POLICY}" \
  --param "run_telocal=${RUN_TELOCAL}" \
  --param "run_tetranscripts=${RUN_TETRANSCRIPTS}" \
  --param "run_telescope=${RUN_TELESCOPE}" \
  --param "run_rediscoverte=${RUN_REDISCOVERTE}" \
  --param "run_rediscoverte_rollup=${RUN_REDISCOVERTE_ROLLUP}" \
  --param "rediscoverte_rollup_conda_prefix=${REDISCOVERTE_ROLLUP_CONDA_PREFIX}" \
  --param "rnaseq_metrics_mem=${RNASEQ_METRICS_MEM}" \
  --param "rnaseq_metrics_java_heap=${RNASEQ_METRICS_JAVA_HEAP}" \
  --param "picard_markdup_java_heap=${PICARD_MARKDUP_JAVA_HEAP}" \
  --param "run_salmonte=${RUN_SALMONTE}"
cp -f "${AUTOMATION_DIR}/run_manifest_${RUN_TS}.json" "${AUTOMATION_DIR}/run_manifest_latest.json"
if [[ "$BACKGROUND" == "true" ]]; then
  pushd "$NF_LAUNCH_DIR" >/dev/null
  BACKGROUND_SCRIPT=$(cat <<EOF
    set +e
    RUN_SESSION_ISOLATED=true
    terminate_descendants() {
      local root="\$1"
      local child
      for child in \$(ps --no-headers --ppid "\$root" -o pid= 2>/dev/null || true); do
        terminate_descendants "\$child"
        kill -TERM "\$child" 2>/dev/null || true
      done
    }
    terminate_process_group() {
      local root="\$1"
      local pgid
      pgid="\$(ps -o pgid= -p "\$root" 2>/dev/null | tr -d ' ' || true)"
      if [[ "\$pgid" =~ ^[0-9]+$ && "\$pgid" != "0" ]]; then
        trap - TERM INT HUP
        kill -TERM -- "-\$pgid" 2>/dev/null || true
        for _ in {1..20}; do
          kill -0 -- "-\$pgid" 2>/dev/null || return 0
          sleep 0.25
        done
        kill -KILL -- "-\$pgid" 2>/dev/null || true
      else
        terminate_descendants "\$root"
      fi
    }
    terminate_conda_for_cache() {
      local pid
      for pid in \$(pgrep -u "\$(id -un)" -f "conda env create --prefix $CONDA_CACHE_DIR" 2>/dev/null || true); do
        [[ "\$pid" == "\$\$" ]] && continue
        kill -TERM "\$pid" 2>/dev/null || true
      done
    }
    handle_signal() {
      terminate_process_group "\$\$"
      terminate_conda_for_cache
      exit 143
    }
    trap handle_signal TERM INT HUP
    ts_now() { date +'%Y-%m-%d %H:%M:%S'; }
    echo "[\$(ts_now)] [run_pipeline] Nextflow started (background)."
    echo "[\$(ts_now)] [run_pipeline] Console log: $LOG_FILE"
    echo "[\$(ts_now)] [run_pipeline] Engine log: $NXF_LOG_FILE"
    echo "[\$(ts_now)] [run_pipeline] Trace: $NXF_TRACE_FILE"
    echo "[\$(ts_now)] [run_pipeline] Report: $NXF_REPORT_FILE"
    echo "[\$(ts_now)] [run_pipeline] Timeline: $NXF_TIMELINE_FILE"
    $CMD_STR
    rc=\$?
    if [[ \$rc -ne 0 ]]; then
      terminate_process_group "\$\$"
      terminate_conda_for_cache
    fi
    echo "[\$(ts_now)] [run_pipeline] Nextflow finished with exit code \$rc"
    if [[ "$LOG_FILE" != "$NXF_LOG_FILE" && -f "$NXF_LOG_FILE" ]] && ! grep -q "BEGIN NEXTFLOW ENGINE LOG" "$LOG_FILE" 2>/dev/null; then
      {
        echo
        echo "========== BEGIN NEXTFLOW ENGINE LOG: $NXF_LOG_FILE =========="
        cat "$NXF_LOG_FILE"
        echo "========== END NEXTFLOW ENGINE LOG: exit_code=\$rc =========="
      } >> "$LOG_FILE"
    fi
    exit \$rc
EOF
)
  nohup setsid bash -lc "$BACKGROUND_SCRIPT" > "$LOG_FILE" 2>&1 < /dev/null &
  pid=$!
  disown || true
  popd >/dev/null

  log_msg "Started in background."
  log_msg "PID: $pid"
  log_msg "Log: $LOG_FILE"
  log_msg "Engine log: $NXF_LOG_FILE"
  log_msg "Trace: $NXF_TRACE_FILE"
  log_msg "Report: $NXF_REPORT_FILE"
  log_msg "Timeline: $NXF_TIMELINE_FILE"
  log_msg "Progress: $PROGRESS_FILE"
  log_msg "PID file: $NEXTFLOW_PID_FILE"
  echo "$pid" > "$NEXTFLOW_PID_FILE"
  start_progress_monitor "$pid"
  log_msg "Progress PID file: $PROGRESS_PID_FILE"
else
  pushd "$NF_LAUNCH_DIR" >/dev/null
  start_progress_monitor "$$"
  set +e
  {
    log_msg "Nextflow started (foreground)."
    log_msg "Console log: $LOG_FILE"
    log_msg "Engine log: $NXF_LOG_FILE"
    log_msg "Trace: $NXF_TRACE_FILE"
    log_msg "Report: $NXF_REPORT_FILE"
    log_msg "Timeline: $NXF_TIMELINE_FILE"
    log_msg "Progress: $PROGRESS_FILE"
    "${CMD[@]}"
    cmd_rc=$?
    log_msg "Nextflow finished with exit code $cmd_rc"
    append_engine_log "$cmd_rc"
    exit "$cmd_rc"
  } 2>&1 | tee -a "$LOG_FILE"
  rc=${PIPESTATUS[0]}
  set -e
  if [[ "$rc" -ne 0 ]]; then
    terminate_descendants "$$"
    terminate_conda_for_cache
  fi
  popd >/dev/null
  exit "$rc"
fi
