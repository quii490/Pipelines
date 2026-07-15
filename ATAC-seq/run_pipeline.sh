
#!/usr/bin/env bash
set -euo pipefail

ORIG_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"
MAIN_NF="${PROJECT_ROOT}/main.nf"

SAMPLESHEET=""
CONTRAST_FILE=""
OUTDIR="${PWD}/atacseq_results"
WORK_DIR="${PWD}/atacseq_results/_automation/work"
SPECIES="hg38"
BACKGROUND="false"
RESUME="true"
QUEUE_SIZE=""
MAX_CPUS=""
MAX_MEMORY=""
RUN_FASTP="true"
RUN_MARKDUP="true"
REMOVE_MITO="true"
REMOVE_BLACKLIST="true"
RUN_BAMCOVERAGE="true"
RUN_PEAK_CALLING="true"
RUN_DOWNSTREAM="true"
CALL_BROAD="false"
TRACK_NORMALIZATION="RPGC"
MAPQ="30"
CONSENSUS_HALF_WIDTH="250"
EXTRA_ARGS=()
RUN_TSS_ENRICH="false"
TSS_BED=""
RUN_GENE_BODY_PROFILE="false"
GENE_BODY_BED=""
RUN_TE_HEATMAP="false"
TE_BED=""
RUN_TE_RELAXED_TRACKS="false"
TE_MAPQ="0"
TE_EXCLUDE_FLAGS="780"
TE_RUN_MARKDUP="false"
TE_REMOVE_MITO="true"
TE_REMOVE_BLACKLIST="true"
TE_PROPER_PAIR_ONLY="true"
TE_TRACK_NORMALIZATION="CPM"
TE_BW_BINSIZE="10"
RUN_FIXEDBIN="true"
FIXEDBIN_SIZE="100000"
RUN_MOTIF="false"
MOTIF_GENOME=""
RUN_FOOTPRINTING="false"
MOTIF_MEME=""
RUN_NUC_PHASING="true"
FASTP_MAX_FORKS=""
FASTP_TIMEOUT=""
PRESET="standard"
PROFILE_PRESET=""
PROFILE_CORES="8"

usage() {
  cat <<'USAGE'
Usage:
  bash run_pipeline.sh --samplesheet samplesheet.csv [options]

Required:
  --samplesheet FILE

Core options:
  --contrast-file FILE
  --outdir DIR
  --work-dir DIR
  --species STR              hg38 | mm10 | mm39

Run control:
  --background
  --resume
  --no-resume
  --queue-size INT
  --max-cpus INT
  --max-memory STR

Pipeline switches:
  --run-fastp BOOL
  --run-markdup BOOL
  --remove-mito BOOL
  --remove-blacklist BOOL
  --run-bamcoverage BOOL
  --run-peak-calling BOOL
  --run-downstream BOOL
  --call-broad BOOL
  --track-normalization STR  RPGC | CPM
  --mapq INT                 默认 30
  --consensus-half-width INT 默认 250 (summit ±250bp)
  --run-tss-enrich BOOL
  --tss-bed FILE
  --run-gene-body-profile BOOL
  --gene-body-bed FILE
  --run-te-heatmap BOOL
  --te-bed FILE
  --run-te-relaxed-tracks BOOL      生成 TE/L1 relaxed BAM+bigWig，供 TE heatmap 使用
  --te-mapq INT                     默认 0
  --te-exclude-flags INT            默认 780，不去掉 duplicate-marked reads
  --te-run-markdup BOOL             默认 false
  --te-remove-mito BOOL             默认 true
  --te-remove-blacklist BOOL        默认 true
  --te-proper-pair-only BOOL        默认 true
  --te-track-normalization STR      CPM | RPGC | RPKM | BPM，默认 CPM
  --te-bw-binsize INT               默认 10
  --run-fixedbin BOOL
  --fixedbin-size INT        默认 100000
  --run-motif BOOL
  --motif-genome STR
  --run-footprinting BOOL
  --motif-meme FILE
  --run-nuc-phasing BOOL
  --fastp-max-forks INT      默认读取 nextflow.config (当前 3)
  --fastp-timeout STR        默认读取 nextflow.config (当前 30m)
  --preset STR               quick | standard | full，默认 standard
  --profile-preset STR       off | quick | standard | full；默认跟随主 preset
  --profile-cores INT        独立 profile 工具线程数，默认 8

Other:
  --extra "ARGS"
  -h, --help
USAGE
}

log_msg() {
  echo "[run_pipeline] $*"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --samplesheet) SAMPLESHEET="$2"; shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --background-internal) BACKGROUND="false"; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --queue-size) QUEUE_SIZE="$2"; shift 2 ;;
    --max-cpus) MAX_CPUS="$2"; shift 2 ;;
    --max-memory) MAX_MEMORY="$2"; shift 2 ;;
    --run-fastp) RUN_FASTP="$2"; shift 2 ;;
    --run-markdup) RUN_MARKDUP="$2"; shift 2 ;;
    --remove-mito) REMOVE_MITO="$2"; shift 2 ;;
    --remove-blacklist) REMOVE_BLACKLIST="$2"; shift 2 ;;
    --run-bamcoverage) RUN_BAMCOVERAGE="$2"; shift 2 ;;
    --run-peak-calling) RUN_PEAK_CALLING="$2"; shift 2 ;;
    --run-downstream) RUN_DOWNSTREAM="$2"; shift 2 ;;
    --call-broad) CALL_BROAD="$2"; shift 2 ;;
    --track-normalization) TRACK_NORMALIZATION="$2"; shift 2 ;;
    --mapq) MAPQ="$2"; shift 2 ;;
    --consensus-half-width) CONSENSUS_HALF_WIDTH="$2"; shift 2 ;;
    --run-tss-enrich) RUN_TSS_ENRICH="$2"; shift 2 ;;
    --tss-bed) TSS_BED="$2"; shift 2 ;;
    --run-gene-body-profile) RUN_GENE_BODY_PROFILE="$2"; shift 2 ;;
    --gene-body-bed) GENE_BODY_BED="$2"; shift 2 ;;
    --run-te-heatmap) RUN_TE_HEATMAP="$2"; shift 2 ;;
    --te-bed) TE_BED="$2"; shift 2 ;;
    --run-te-relaxed-tracks) RUN_TE_RELAXED_TRACKS="$2"; shift 2 ;;
    --te-mapq) TE_MAPQ="$2"; shift 2 ;;
    --te-exclude-flags) TE_EXCLUDE_FLAGS="$2"; shift 2 ;;
    --te-run-markdup) TE_RUN_MARKDUP="$2"; shift 2 ;;
    --te-remove-mito) TE_REMOVE_MITO="$2"; shift 2 ;;
    --te-remove-blacklist) TE_REMOVE_BLACKLIST="$2"; shift 2 ;;
    --te-proper-pair-only) TE_PROPER_PAIR_ONLY="$2"; shift 2 ;;
    --te-track-normalization) TE_TRACK_NORMALIZATION="$2"; shift 2 ;;
    --te-bw-binsize) TE_BW_BINSIZE="$2"; shift 2 ;;
    --run-fixedbin) RUN_FIXEDBIN="$2"; shift 2 ;;
    --fixedbin-size) FIXEDBIN_SIZE="$2"; shift 2 ;;
    --run-motif) RUN_MOTIF="$2"; shift 2 ;;
    --motif-genome) MOTIF_GENOME="$2"; shift 2 ;;
    --run-footprinting) RUN_FOOTPRINTING="$2"; shift 2 ;;
    --motif-meme) MOTIF_MEME="$2"; shift 2 ;;
    --run-nuc-phasing) RUN_NUC_PHASING="$2"; shift 2 ;;
    --fastp-max-forks) FASTP_MAX_FORKS="$2"; shift 2 ;;
    --fastp-timeout) FASTP_TIMEOUT="$2"; shift 2 ;;
    --preset) PRESET="$2"; shift 2 ;;
    --profile-preset) PROFILE_PRESET="$2"; shift 2 ;;
    --profile-cores) PROFILE_CORES="$2"; shift 2 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$SAMPLESHEET" ]] || { echo "[run_pipeline] ERROR: --samplesheet is required" >&2; exit 1; }
[[ -f "$SAMPLESHEET" ]] || { echo "[run_pipeline] ERROR: samplesheet not found: $SAMPLESHEET" >&2; exit 1; }
if [[ -n "$CONTRAST_FILE" && ! -f "$CONTRAST_FILE" ]]; then
  echo "[run_pipeline] ERROR: contrast file not found: $CONTRAST_FILE" >&2
  exit 1
fi
[[ -f "$MAIN_NF" ]] || { echo "[run_pipeline] ERROR: main.nf not found: $MAIN_NF" >&2; exit 1; }
command -v nextflow >/dev/null 2>&1 || { echo "[run_pipeline] ERROR: nextflow not found in PATH" >&2; exit 1; }

case "$PRESET" in
  quick)
    : "${PROFILE_PRESET:=quick}"
    RUN_TSS_ENRICH="false"
    RUN_GENE_BODY_PROFILE="false"
    RUN_TE_HEATMAP="false"
    RUN_FIXEDBIN="false"
    RUN_MOTIF="false"
    RUN_FOOTPRINTING="false"
    RUN_NUC_PHASING="false"
    ;;
  standard)
    : "${PROFILE_PRESET:=standard}"
    ;;
  full)
    : "${PROFILE_PRESET:=standard}"
    RUN_FIXEDBIN="true"
    RUN_MOTIF="true"
    RUN_FOOTPRINTING="true"
    RUN_NUC_PHASING="true"
    ;;
  *)
    echo "[run_pipeline] ERROR: --preset must be quick, standard, or full" >&2
    exit 1
    ;;
esac
case "$PROFILE_PRESET" in
  off|quick|standard|full) ;;
  *) echo "[run_pipeline] ERROR: --profile-preset must be off, quick, standard, or full" >&2; exit 1 ;;
esac

AUTOMATION_DIR="${OUTDIR}/_automation"
LOG_DIR="${AUTOMATION_DIR}/logs"
CONDA_CACHE_DIR="${AUTOMATION_DIR}/conda_cache"
TMP_DIR="${AUTOMATION_DIR}/tmp"
TS="$(date +%Y%m%d_%H%M%S)"
RUN_NAME="${NXF_RUN_NAME:-atac_${SPECIES}_${TS}}"
NXF_BASE="${NXF_BASE:-/tmp/${USER:-$(id -un)}/nextflow}"
NXF_HOME_DIR="${NXF_HOME:-${NXF_BASE}/home}"
NF_LAUNCH_DIR="${NF_LAUNCH_DIR:-${NXF_BASE}/launch/${RUN_NAME}}"
LOG_FILE="${LOG_DIR}/pipeline_${TS}.log"
NXF_LOG_FILE="${LOG_DIR}/nextflow_${TS}.log"
PID_FILE="${LOG_DIR}/pipeline_${TS}.pid"
mkdir -p "$OUTDIR" "$WORK_DIR" "$LOG_DIR" "$NXF_HOME_DIR" "$NF_LAUNCH_DIR" "$CONDA_CACHE_DIR" "$TMP_DIR"

if [[ "$BACKGROUND" == "true" && "${RUN_PIPELINE_BG_CHILD:-0}" != "1" ]]; then
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

log_msg "OUTDIR=$OUTDIR"
log_msg "WORK_DIR=$WORK_DIR"
log_msg "NXF_BASE=$NXF_BASE"
log_msg "NF_LAUNCH_DIR=$NF_LAUNCH_DIR"
log_msg "NXF_HOME=$NXF_HOME_DIR"
log_msg "RUN_NAME=$RUN_NAME"
log_msg "SPECIES=$SPECIES"
log_msg "SAMPLESHEET=$SAMPLESHEET"
[[ -n "$CONTRAST_FILE" ]] && log_msg "CONTRAST_FILE=$CONTRAST_FILE"

export NXF_HOME="$NXF_HOME_DIR"
# Avoid Nextflow history/cache file locks on /home NFS (mounted with local_lock=none).
# Defaults keep Nextflow launch/cache under /tmp; result files and work-dir still stay in OUTDIR.
export NXF_IGNORE_RESUME_HISTORY="${NXF_IGNORE_RESUME_HISTORY:-true}"
export NXF_DISABLE_CHECK_LATEST="${NXF_DISABLE_CHECK_LATEST:-true}"
export NXF_CONDA_CACHEDIR="$CONDA_CACHE_DIR"
export TMPDIR="$TMP_DIR"
export TMP="$TMP_DIR"
export TEMP="$TMP_DIR"

CMD=(nextflow -log "$NXF_LOG_FILE" run "$MAIN_NF"
  -name "$RUN_NAME"
  -work-dir "$WORK_DIR"
  --samplesheet "$SAMPLESHEET"
  --outdir "$OUTDIR"
  --conda_cache_dir "$CONDA_CACHE_DIR"
  --species "$SPECIES"
  --run_fastp "$RUN_FASTP"
  --run_markdup "$RUN_MARKDUP"
  --remove_mito "$REMOVE_MITO"
  --remove_blacklist "$REMOVE_BLACKLIST"
  --run_bamcoverage "$RUN_BAMCOVERAGE"
  --run_peak_calling "$RUN_PEAK_CALLING"
  --run_downstream "$RUN_DOWNSTREAM"
  --call_broad "$CALL_BROAD"
  --track_normalization "$TRACK_NORMALIZATION"
  --mapq "$MAPQ"
  --consensus_half_width "$CONSENSUS_HALF_WIDTH"
  --run_tss_enrich "$RUN_TSS_ENRICH"
  --run_gene_body_profile "$RUN_GENE_BODY_PROFILE"
  --run_te_heatmap "$RUN_TE_HEATMAP"
  --run_te_relaxed_tracks "$RUN_TE_RELAXED_TRACKS"
  --te_mapq "$TE_MAPQ"
  --te_exclude_flags "$TE_EXCLUDE_FLAGS"
  --te_run_markdup "$TE_RUN_MARKDUP"
  --te_remove_mito "$TE_REMOVE_MITO"
  --te_remove_blacklist "$TE_REMOVE_BLACKLIST"
  --te_proper_pair_only "$TE_PROPER_PAIR_ONLY"
  --te_track_normalization "$TE_TRACK_NORMALIZATION"
  --te_bw_binsize "$TE_BW_BINSIZE"
  --run_fixedbin "$RUN_FIXEDBIN"
  --fixedbin_size "$FIXEDBIN_SIZE"
  --run_motif "$RUN_MOTIF"
  --run_footprinting "$RUN_FOOTPRINTING"
  --run_nuc_phasing "$RUN_NUC_PHASING")
if [[ -n "$TSS_BED" ]]; then CMD+=(--tss_bed "$TSS_BED"); fi
if [[ -n "$GENE_BODY_BED" ]]; then CMD+=(--gene_body_bed "$GENE_BODY_BED"); fi
if [[ -n "$TE_BED" ]]; then CMD+=(--te_bed "$TE_BED"); fi
if [[ -n "$MOTIF_GENOME" ]]; then CMD+=(--motif_genome "$MOTIF_GENOME"); fi
if [[ -n "$MOTIF_MEME" ]]; then CMD+=(--motif_meme "$MOTIF_MEME"); fi
if [[ -n "$CONTRAST_FILE" ]]; then CMD+=(--contrast_file "$CONTRAST_FILE"); fi
if [[ "$RESUME" == "true" ]]; then
  if [[ "${NXF_IGNORE_RESUME_HISTORY:-}" == "true" && -z "${NXF_RESUME_ID:-}" ]]; then
    log_msg "NXF_IGNORE_RESUME_HISTORY=true；未设置 NXF_RESUME_ID，自动跳过 -resume 以避免 Nextflow history lock 报错。"
  elif [[ -n "${NXF_RESUME_ID:-}" ]]; then
    CMD+=(-resume "$NXF_RESUME_ID")
  else
    CMD+=(-resume)
  fi
fi
if [[ -n "$QUEUE_SIZE" ]]; then CMD+=(--queue_size "$QUEUE_SIZE"); fi
if [[ -n "$MAX_CPUS" ]]; then CMD+=(--max_cpus "$MAX_CPUS"); fi
if [[ -n "$MAX_MEMORY" ]]; then CMD+=(--max_memory "$MAX_MEMORY"); fi
if [[ -n "$FASTP_MAX_FORKS" ]]; then CMD+=(--fastp_max_forks "$FASTP_MAX_FORKS"); fi
if [[ -n "$FASTP_TIMEOUT" ]]; then CMD+=(--fastp_timeout "$FASTP_TIMEOUT"); fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for x in "${EXTRA_ARGS[@]}"; do
    EXTRA_SPLIT=( $x )
    CMD+=("${EXTRA_SPLIT[@]}")
  done
fi
log_msg "CMD=${CMD[*]}"
cd "$NF_LAUNCH_DIR"
"${CMD[@]}" 2>&1 | tee -a "$LOG_FILE"

if command -v python3 >/dev/null 2>&1 && [[ -f "${PROJECT_ROOT}/scripts/generate_qc_report.py" ]]; then
  log_msg "生成中文 QC 报告"
  python3 "${PROJECT_ROOT}/scripts/generate_qc_report.py" --outdir "$OUTDIR" --species "$SPECIES" | tee -a "$LOG_FILE"
fi

if [[ "$PROFILE_PRESET" != "off" && "$RUN_BAMCOVERAGE" == "true" ]]; then
  log_msg "生成 ${PROFILE_PRESET} profile heatmaps（独立工具，可用 --profile-preset off 关闭）"
  bash "${PROJECT_ROOT}/scripts/run_atac_profile_heatmaps.sh" \
    --result-dir "$OUTDIR" \
    --species "$SPECIES" \
    --samplesheet "$SAMPLESHEET" \
    --preset "$PROFILE_PRESET" \
    --cores "$PROFILE_CORES" | tee -a "$LOG_FILE"
fi
