#!/usr/bin/env bash
set -euo pipefail

ORIG_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}" && pwd)"
PIPELINE_RUNNER="${PROJECT_ROOT}/run_pipeline.sh"
SAMPLESHEET_GEN="${PROJECT_ROOT}/generate_samplesheet.sh"
DOWNSTREAM_RUNNER="${PROJECT_ROOT}/scripts/run_atac_downstream_only.sh"

MODE="auto"
FASTQ_DIR=""
SPECIES="hg38"
LAYOUT="auto"
BACKGROUND="false"
BACKGROUND_INTERNAL="false"
RESUME="true"
RUN_PIPELINE="true"
INIT_ONLY="false"
OUTDIR=""
WORK_DIR=""
SAMPLESHEET=""
CONTRAST_FILE=""
QUEUE_SIZE=""
MAX_CPUS=""
MAX_MEMORY=""
FASTP_MAX_FORKS=""
FASTP_TIMEOUT=""
METADATA_CSV=""
OVERWRITE_INPUTS="false"
PRESET="standard"
PROFILE_PRESET=""
PROFILE_CORES="8"
LEVELS="both"
RUN_TE_RELAXED_TRACKS="false"
TE_MAPQ="0"
TE_TRACK_NORMALIZATION="CPM"
TE_BW_BINSIZE="10"
EXTRA_ARGS=()
CONTRAST_SPECS=()

usage() {
  cat <<'USAGE'
ATAC-seq 自动化入口：自动扫描 FASTQ -> 生成 samplesheet/contrast -> 编辑后运行 pipeline

常用：
  1) 初始化模板：
     bash run_auto_atacseq.sh --fastq-dir FASTQ_DIR --init-only

  2) 模板编辑后直接运行：
     bash run_auto_atacseq.sh --fastq-dir FASTQ_DIR --species hg38

  3) 上游完成后只重跑下游：
     bash run_auto_atacseq.sh --mode downstream --outdir RESULTS --species hg38 --levels both

参数：
  --mode STR               auto | init | upstream | downstream，默认 auto
  --fastq-dir DIR
  --species STR             hg38 | mm10 | mm39
  --layout STR              auto | PE | SE，默认 auto
  --outdir DIR              默认 <fastq-dir同级>/atacseq_results
  --work-dir DIR            默认 <outdir>/_automation/work
  --samplesheet FILE        默认 <outdir>/_automation/inputs/samplesheet.csv
  --contrast-file FILE      默认 <outdir>/_automation/inputs/contrasts.csv
  --init-only               只生成模板，不启动流程
  --background              后台运行
  --resume / --no-resume
  --queue-size INT
  --max-cpus INT
  --max-memory STR
  --fastp-max-forks INT
  --fastp-timeout STR
  --metadata-csv FILE       可选，列包含 sample,condition,replicate；自动回填 samplesheet
  --contrast CASE,CONTROL   可重复；直接写入 contrasts.csv，例如 --contrast KO,WT
  --overwrite-inputs        重新生成 samplesheet/contrasts，覆盖旧模板
  --preset STR              quick | standard | full，默认 standard
  --profile-preset STR      off | quick | standard | full；默认跟随主 preset
  --profile-cores INT       profile heatmap 线程数，默认 8
  --levels STR              downstream 模式：peak | bin | both，默认 both
  --downstream-only         等价于 --mode downstream
  --run-te-relaxed-tracks BOOL  是否额外生成 TE/L1 relaxed BAM+bigWig，默认 false
  --te-mapq INT                 TE relaxed track 的 MAPQ，默认 0
  --te-track-normalization STR  CPM | RPGC | RPKM | BPM，默认 CPM
  --te-bw-binsize INT           TE relaxed bigWig bin size，默认 10
  --extra "ARGS"
USAGE
}

log_msg() {
  echo "[run_auto_atacseq] $*"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --mode) MODE="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --init-only) INIT_ONLY="true"; RUN_PIPELINE="false"; shift 1 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --background-internal) BACKGROUND_INTERNAL="true"; BACKGROUND="false"; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --queue-size) QUEUE_SIZE="$2"; shift 2 ;;
    --max-cpus) MAX_CPUS="$2"; shift 2 ;;
    --max-memory) MAX_MEMORY="$2"; shift 2 ;;
    --fastp-max-forks) FASTP_MAX_FORKS="$2"; shift 2 ;;
    --fastp-timeout) FASTP_TIMEOUT="$2"; shift 2 ;;
    --metadata-csv) METADATA_CSV="$2"; shift 2 ;;
    --contrast) CONTRAST_SPECS+=("$2"); shift 2 ;;
    --overwrite-inputs) OVERWRITE_INPUTS="true"; shift 1 ;;
    --preset) PRESET="$2"; shift 2 ;;
    --profile-preset) PROFILE_PRESET="$2"; shift 2 ;;
    --profile-cores) PROFILE_CORES="$2"; shift 2 ;;
    --levels|--level) LEVELS="$2"; shift 2 ;;
    --downstream-only) MODE="downstream"; shift 1 ;;
    --run-te-relaxed-tracks) RUN_TE_RELAXED_TRACKS="$2"; shift 2 ;;
    --te-mapq) TE_MAPQ="$2"; shift 2 ;;
    --te-track-normalization) TE_TRACK_NORMALIZATION="$2"; shift 2 ;;
    --te-bw-binsize) TE_BW_BINSIZE="$2"; shift 2 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[run_auto_atacseq] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

case "$MODE" in
  auto|init|upstream|downstream) ;;
  *) echo "[run_auto_atacseq] ERROR: --mode must be auto, init, upstream, or downstream" >&2; exit 1 ;;
esac

if [[ "$MODE" == "init" ]]; then
  INIT_ONLY="true"
  RUN_PIPELINE="false"
fi

if [[ "$MODE" == "downstream" ]]; then
  [[ -n "$OUTDIR" ]] || { echo "[run_auto_atacseq] ERROR: downstream mode requires --outdir RESULTS" >&2; exit 1; }
  [[ -d "$OUTDIR" ]] || { echo "[run_auto_atacseq] ERROR: result dir not found: $OUTDIR" >&2; exit 1; }
  INPUTS_DIR="${OUTDIR}/_automation/inputs"
  LOG_DIR="${OUTDIR}/_automation/logs"
  mkdir -p "$INPUTS_DIR" "$LOG_DIR"
  [[ -n "$SAMPLESHEET" ]] || SAMPLESHEET="${INPUTS_DIR}/samplesheet.csv"
  [[ -n "$CONTRAST_FILE" ]] || CONTRAST_FILE="${INPUTS_DIR}/contrasts.csv"
  CMD=(bash "$DOWNSTREAM_RUNNER"
    --result-dir "$OUTDIR"
    --species "$SPECIES"
    --levels "$LEVELS"
    --samplesheet "$SAMPLESHEET"
    --contrast-file "$CONTRAST_FILE")
  if [[ "$BACKGROUND" == "true" ]]; then CMD+=(--background); fi
  log_msg "启动下游分析"
  log_msg "CMD=${CMD[*]}"
  "${CMD[@]}"
  exit 0
fi

[[ -n "$FASTQ_DIR" ]] || { echo "[run_auto_atacseq] ERROR: --fastq-dir is required unless --mode downstream" >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "[run_auto_atacseq] ERROR: fastq dir not found: $FASTQ_DIR" >&2; exit 1; }

BASE_DIR="$(cd "$(dirname "$FASTQ_DIR")" && pwd)"
[[ -n "$OUTDIR" ]] || OUTDIR="${BASE_DIR}/atacseq_results"
[[ -n "$WORK_DIR" ]] || WORK_DIR="${OUTDIR}/_automation/work"
INPUTS_DIR="${OUTDIR}/_automation/inputs"
LOG_DIR="${OUTDIR}/_automation/logs"
mkdir -p "$INPUTS_DIR" "$LOG_DIR"
[[ -n "$SAMPLESHEET" ]] || SAMPLESHEET="${INPUTS_DIR}/samplesheet.csv"
[[ -n "$CONTRAST_FILE" ]] || CONTRAST_FILE="${INPUTS_DIR}/contrasts.csv"
AUTOMATION_LOG="${LOG_DIR}/automation_$(date +%Y%m%d_%H%M%S).log"
PID_FILE="${LOG_DIR}/automation.pid"

if [[ "$BACKGROUND" == "true" && "$BACKGROUND_INTERNAL" != "true" && "${RUN_AUTO_ATACSEQ_BG_CHILD:-0}" != "1" ]]; then
  FILTERED_ARGS=()
  for x in "${ORIG_ARGS[@]}"; do
    [[ "$x" == "--background" ]] && continue
    [[ "$x" == "--background-internal" ]] && continue
    FILTERED_ARGS+=("$x")
  done
  log_msg "后台运行，日志: ${AUTOMATION_LOG}"
  nohup env RUN_AUTO_ATACSEQ_BG_CHILD=1 bash "$0" "${FILTERED_ARGS[@]}" --background-internal > "$AUTOMATION_LOG" 2>&1 &
  echo $! > "$PID_FILE"
  log_msg "PID: $(cat "$PID_FILE")"
  exit 0
fi

log_msg "扫描 FASTQ 并生成模板"
if [[ ! -f "$SAMPLESHEET" || ! -f "$CONTRAST_FILE" ]]; then
  GEN_CMD=(bash "$SAMPLESHEET_GEN"
    --input "$FASTQ_DIR"
    --samplesheet "$SAMPLESHEET"
    --contrasts "$CONTRAST_FILE"
    --layout "$LAYOUT")
  if [[ -n "$METADATA_CSV" ]]; then GEN_CMD+=(--metadata-csv "$METADATA_CSV"); fi
  "${GEN_CMD[@]}"
  log_msg "请直接编辑以下文件："
  log_msg "  $SAMPLESHEET"
  log_msg "  $CONTRAST_FILE"
elif [[ "$OVERWRITE_INPUTS" == "true" ]]; then
  log_msg "根据 --overwrite-inputs 重新生成设计文件"
  GEN_CMD=(bash "$SAMPLESHEET_GEN"
    --input "$FASTQ_DIR"
    --samplesheet "$SAMPLESHEET"
    --contrasts "$CONTRAST_FILE"
    --layout "$LAYOUT")
  if [[ -n "$METADATA_CSV" ]]; then GEN_CMD+=(--metadata-csv "$METADATA_CSV"); fi
  "${GEN_CMD[@]}"
else
  log_msg "检测到已存在的设计文件，保留现有内容，不重新生成："
  log_msg "  $SAMPLESHEET"
  log_msg "  $CONTRAST_FILE"
fi

if [[ ${#CONTRAST_SPECS[@]} -gt 0 ]]; then
  log_msg "根据 --contrast 写入 contrast 文件: $CONTRAST_FILE"
  {
    echo "case,control"
    for spec in "${CONTRAST_SPECS[@]}"; do
      spec="${spec/:/,}"
      IFS=',' read -r case_name control_name <<< "$spec"
      case_name="${case_name//[[:space:]]/}"
      control_name="${control_name//[[:space:]]/}"
      [[ -n "$case_name" && -n "$control_name" ]] || { echo "[run_auto_atacseq] ERROR: invalid --contrast '$spec', use CASE,CONTROL" >&2; exit 1; }
      echo "${case_name},${control_name}"
    done
  } > "$CONTRAST_FILE"
fi

if [[ "$INIT_ONLY" == "true" ]]; then
  log_msg "初始化完成。请编辑以下文件后再启动："
  log_msg "  $SAMPLESHEET"
  log_msg "  $CONTRAST_FILE"
  exit 0
fi

HAS_NON_NA_CONDITION="$(awk -F',' 'NR>1 && $3!="" && $3!="NA" {print "yes"; exit}' "$SAMPLESHEET" || true)"
HAS_CONTRAST_ROW="$(awk -F',' 'NR>1 {gsub(/\r/,"",$1); gsub(/\r/,"",$2); if($1!="" && $2!="") {print "yes"; exit}}' "$CONTRAST_FILE" || true)"

if [[ "$HAS_NON_NA_CONDITION" != "yes" ]]; then
  log_msg "samplesheet 里的 condition 仍然是 NA；已停止。请先手动填写后再运行。"
  log_msg "文件: $SAMPLESHEET"
  exit 0
fi
if [[ "$HAS_CONTRAST_ROW" != "yes" ]]; then
  log_msg "contrast 文件还没有有效比较；已停止。请先手动填写 case/control 后再运行。"
  log_msg "文件: $CONTRAST_FILE"
  exit 0
fi

CMD=(bash "$PIPELINE_RUNNER"
  --samplesheet "$SAMPLESHEET"
  --contrast-file "$CONTRAST_FILE"
  --species "$SPECIES"
  --outdir "$OUTDIR"
  --work-dir "$WORK_DIR"
  --preset "$PRESET"
  --profile-cores "$PROFILE_CORES"
  --run-te-relaxed-tracks "$RUN_TE_RELAXED_TRACKS"
  --te-mapq "$TE_MAPQ"
  --te-track-normalization "$TE_TRACK_NORMALIZATION"
  --te-bw-binsize "$TE_BW_BINSIZE")
if [[ -n "$PROFILE_PRESET" ]]; then CMD+=(--profile-preset "$PROFILE_PRESET"); fi

if [[ "$RESUME" == "true" ]]; then CMD+=(--resume); else CMD+=(--no-resume); fi
if [[ -n "$QUEUE_SIZE" ]]; then CMD+=(--queue-size "$QUEUE_SIZE"); fi
if [[ -n "$MAX_CPUS" ]]; then CMD+=(--max-cpus "$MAX_CPUS"); fi
if [[ -n "$MAX_MEMORY" ]]; then CMD+=(--max-memory "$MAX_MEMORY"); fi
if [[ -n "$FASTP_MAX_FORKS" ]]; then CMD+=(--fastp-max-forks "$FASTP_MAX_FORKS"); fi
if [[ -n "$FASTP_TIMEOUT" ]]; then CMD+=(--fastp-timeout "$FASTP_TIMEOUT"); fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for x in "${EXTRA_ARGS[@]}"; do
    CMD+=(--extra "$x")
  done
fi

log_msg "启动流程"
log_msg "CMD=${CMD[*]}"
"${CMD[@]}"
