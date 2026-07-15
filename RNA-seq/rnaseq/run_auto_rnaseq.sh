#!/usr/bin/env bash
set -euo pipefail

ORIG_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOWNSTREAM_DIR="${PROJECT_ROOT}/rnaseq-downstream"
PIPELINE_RUNNER="${SCRIPT_DIR}/run_pipeline.sh"
SAMPLESHEET_GEN="${SCRIPT_DIR}/generate_samplesheet.sh"
MANIFEST_HELPER="${SCRIPT_DIR}/scripts/prepare_manifest_fastq.py"
DOWNSTREAM_RUNNER="${DOWNSTREAM_DIR}/run_downstream.R"
MODULE_DEFAULTS="${SCRIPT_DIR}/nextflow.config"
DEFAULT_READER="${SCRIPT_DIR}/scripts/read_module_default.py"

module_default() {
  python3 "${DEFAULT_READER}" "${MODULE_DEFAULTS}" "$1" "$2"
}

MANIFEST=""
FASTQ_DIR=""
SPECIES="hg38"
HUMAN_REF="hg38"
ALIGNER="star"
STRANDEDNESS="reverse"
STRAND_SOURCE="default"
LAYOUT="auto"
R1_PATTERN="_R1"
R2_PATTERN="_R2"
TECOUNT_STRAND=""
BACKGROUND="false"
RUN_UPSTREAM="true"
RUN_DOWNSTREAM="false"
RESUME="true"
RESUME_SESSION=""
PROFILE=""
QUEUE_SIZE=""
MAX_CPUS=""
MAX_MEMORY=""
EXTRA_ARGS=()

WORK_ROOT="${PWD}/rnaseq_auto_run"
RESULTS_DIR=""
PLOT_OUTDIR=""
DOWNLOAD_DIR=""
INPUTS_DIR=""
USER_SET_WORK_ROOT="false"
USER_SET_RESULTS_DIR="false"
USER_SET_PLOT_OUTDIR="false"
USER_SET_DOWNLOAD_DIR="false"
AUTO_LAYOUT_BASE=""
RUN_TS="$(date +%Y%m%d_%H%M%S)"

SAMPLE_METADATA=""
CONTRAST_FILE=""
SAMPLE_META_ARGS=()
CONDITION_MAP_ARGS=()
CONTRAST_ARGS=()
DEFAULT_CONDITION="NA"
DEFAULT_REPLICATE="NA"

RUN_FASTP="$(module_default run_fastp true)"
RUN_STAR_FC="$(module_default run_star_fc true)"
RUN_SALMON="$(module_default run_salmon true)"
RUN_STRINGTIE="$(module_default run_stringtie false)"
RUN_TECOUNT="$(module_default run_tecount false)"
RUN_TELOCAL="$(module_default run_telocal false)"
RUN_TETRANSCRIPTS="$(module_default run_tetranscripts true)"
RUN_DEDUP="$(module_default run_dedup false)"
RUN_MARKDUP_QC="$(module_default run_markdup_qc true)"
RUN_REDISCOVERTE="auto"
RUN_REDISCOVERTE_ROLLUP="$(module_default run_rediscoverte_rollup true)"
RUN_SALMONTE="$(module_default run_salmonte false)"
RUN_MULTIQC="$(module_default run_multiqc true)"
MULTIQC_CMD="$(module_default multiqc_cmd /path/to/.conda/envs/emseq/bin/multiqc)"
RUN_TELESCOPE="$(module_default run_telescope true)"
RUN_FASTQC="$(module_default run_fastqc true)"
RUN_RNASEQ_METRICS="$(module_default run_rnaseq_metrics true)"
REF_FLAT=""
RIBOSOMAL_INTERVALS=""

PADJ_CUTOFF="0.05"
LFC_CUTOFF="0.58"
BASEMEAN_MIN="5"
LABEL_TOP_N="40"
HEATMAP_TOP_N="40"
VOLCANO_ORIENTATION="classic"
GRAY_NONSIG="true"
EXPLORATORY_METHOD="logCPM_diff"
EXPLORATORY_FIXED_BCV="0.4"
PLOT_THREADS=""
ONLY_TOOLS=""
SKIP_TOOLS=""
PARTIAL_INPUT_POLICY="skip"
DRY_RUN="false"
FAILURE_POLICY="core"
REPLACE_DESIGN="false"
RNASEQ_ENV="${RNASEQ_ENV:-/path/to/.conda/envs/rnaseq}"
DOWNSTREAM_ENV="${DOWNSTREAM_ENV:-/path/to/.conda/envs/downstream}"
ROLLUP_CONDA_PREFIX=""
PICARD_MARKDUP_JAVA_HEAP="${PICARD_MARKDUP_JAVA_HEAP:-16g}"
SRA_THREADS="8"
PROGRESS_INTERVAL="30"
RUN_SESSION_ISOLATED="false"
UPSTREAM_STATUS="NOT_RUN"

usage() {
  cat <<'USAGE'
统一 RNA-seq 自动化入口（manifest 下载 / 直接 FASTQ / 仅重跑下游绘图）

两步常用模式：
  1) fastq-dir / manifest -> 自动识别与上游分析 -> 生成 results-dir 下的 condition.csv / contrast.csv 模板
  2) 手动修改 condition.csv / contrast.csv 后，再用 --downstream 只跑下游差异分析与出图

输入模式（至少满足一种）：
  --manifest FILE           manifest 表；支持本地 fastq、fastq URL、SRR/ERR/DRR、部分 GSM
  --fastq-dir DIR           已有 FASTQ 根目录
  --downstream              仅重跑下游，默认读取 results-dir 下的 condition.csv / contrast.csv
  --downstream-only         同 --downstream

设计输入：
  --sample-metadata FILE    可选，CSV 至少包含 sample,condition,replicate
  --sample-meta STR         命令行逐样本指定；格式 sample:condition:replicate；可重复
  --condition-map STR       命令行逐组指定；格式 condition:sample1,sample2,...；replicate 自动按顺序 1..N
  --contrast-file FILE      CSV，列名: case,control（兼容旧格式 group_col,case,control）
  --contrast STR            命令行指定 contrast；推荐格式 CASE:CONTROL；兼容旧格式 group_col:CASE:CONTROL；可重复

路径：
  --work-root DIR           工作根目录；默认 <fastq-dir同级>/rnaseq_auto_run 或 <manifest同级>/rnaseq_auto_run
  --results-dir DIR         上游结果目录；默认 <fastq-dir同级>/rnaseq_results 或 <manifest同级>/rnaseq_results
  --plot-outdir DIR         下游绘图输出目录；默认 <results-dir>/plots
  --plot-dir DIR            同 --plot-outdir
  --download-dir DIR        manifest 下载出的 FASTQ 目录；默认 <work-root>/fastq

核心参数：
  --species STR             hg38 | mm10 | mm39
  --human-ref STR           hg38 | t2t，默认 hg38；仅 species=hg38 时生效
  --aligner STR             star | hisat2，默认 star（基因比对）
  --strand STR              unstranded | forward | reverse；同时用于 gene 和 TE
  --strandedness STR        兼容旧参数，同 --strand
  --layout STR              auto | PE | SE，默认 auto
  --r1-pattern STR          默认 _R1
  --r2-pattern STR          默认 _R2
  --tecount-strand STR      兼容旧参数，同 --strand；不建议再单独使用
  --sra-threads INT         manifest 中跑 fasterq-dump 的线程数，默认 8

运行控制：
  --background              后台运行整个自动化流程
  --resume                  上游启用 -resume（默认）
  --no-resume               上游禁用 -resume
  --resume-session ID       恢复指定 Nextflow session（高级选项）
  --upstream-only           仅跑上游（默认行为）
  --downstream              仅跑下游
  --downstream-only         同 --downstream
  --profile NAME            Nextflow profile
  --queue-size INT          Nextflow queue_size
  --max-cpus INT            最大 CPU
  --max-memory STR          最大内存，例如 "300 GB"
  --progress-interval INT   进度文件刷新间隔（秒），默认 30
  --dry-run                只预检输入、参数和计划模块，不提交分析任务
  --failure-policy STR     core|strict，默认 core；可选模块失败不阻断核心结果

上游开关：
  --run-fastp BOOL
  --run-gene-count-branch BOOL
  --run-salmon BOOL
  --run-stringtie BOOL      默认 false
  --run-tecount BOOL        默认 false
  --run-telocal BOOL        默认 false；需要时显式开启
  --run-tetranscripts BOOL
  --run-dedup BOOL         默认 false；高级选项：用 dedup BAM 做 counts，普通 RNA-seq 不推荐
  --run-markdup-qc BOOL    默认 true；只输出 Picard duplicate metrics，不改变 counts
  --run-rediscoverte BOOL    true|false|auto，默认 auto；hg38 跑，鼠跳过
  --run-rediscoverte-rollup BOOL
                            是否调用 REdiscoverTE rollup，默认 true；用于默认生成分层 TE 可视化
  --run-salmonte BOOL       默认 false；实验性兼容模块
  --run-multiqc BOOL       默认 true；汇总 FastQC、比对、计数和 Picard QC
  --multiqc-cmd PATH       MultiQC 可执行文件；默认使用 config 中的 multiqc_cmd
  --run-telescope BOOL    默认 true；基于贝叶斯EM模型重分配multi-mapped TE reads
  --run-fastqc BOOL       默认 true；对原始 FASTQ 运行 FastQC（fastp 同时提供 clean FASTQ QC）
  --run-rnaseq-metrics BOOL  默认 true；Picard RNA-seq 区域分布和偏倚指标
  --ref-flat FILE         可选；不提供时自动由基因 GTF 生成
  --refFlat FILE          --ref-flat 的兼容别名
  --ribosomal-intervals FILE
                          可选 Picard rRNA interval list

下游绘图参数：
  --padj-cutoff NUM         默认 0.05
  --lfc-cutoff NUM          默认 0.58
  --baseMean-min NUM        默认 5
  --label-top-n NUM         默认 40
  --heatmap-top-n NUM       默认 40
  --volcano-orientation STR classic|horizontal，默认 classic
  --gray-nonsig BOOL        非显著点是否统一灰色，默认 true
  --exploratory-method STR  无重复差异算法：logCPM_diff|edgeR_fixedBCV，默认 logCPM_diff
  --exploratory-fixed-bcv NUM
                            edgeR_fixedBCV 的 BCV；human/TE 可用 0.4，默认 0.4
  --plot-threads NUM        下游 contrast 并行数；默认 min(--max-cpus, 8)，未设置 max-cpus 时为 1
  --only-tools STR          仅运行指定下游模块，逗号分隔；例如 TE_TEtranscripts
  --skip-tools STR          跳过指定下游模块，逗号分隔；例如 TE_TElocal,TE_TEcount
  --partial-input-policy STR  skip|error|allow，默认 skip；缺样本时默认整模块跳过
  --allow-partial-inputs BOOL  旧兼容参数；true 等价于 policy=allow

其它：
  --default-condition STR   没有命中 metadata 时的默认 condition，默认 NA
  --default-replicate STR   没有命中 metadata 时的默认 replicate，默认 NA
  --replace-design          允许覆盖 condition.csv/contrast.csv；覆盖前自动备份
  --rnaseq-env DIR          上游 conda 环境路径
  --downstream-env DIR      下游 conda 环境路径
  --rediscoverte-rollup-conda-prefix DIR
                            REdiscoverTE rollup 使用的现有环境；默认跟随 --downstream-env
  --picard-markdup-java-heap STR
                            MarkDuplicates JVM heap；默认 16g（任务内存默认 24 GB）
  --extra "ARGS"           追加到 nextflow run 的原始参数
  -h, --help                显示帮助
USAGE
}

log_msg() {
  echo "[run_auto_rnaseq] $*"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --work-root) WORK_ROOT="$2"; USER_SET_WORK_ROOT="true"; shift 2 ;;
    --results-dir) RESULTS_DIR="$2"; USER_SET_RESULTS_DIR="true"; shift 2 ;;
    --plot-outdir|--plot-dir) PLOT_OUTDIR="$2"; USER_SET_PLOT_OUTDIR="true"; shift 2 ;;
    --download-dir) DOWNLOAD_DIR="$2"; USER_SET_DOWNLOAD_DIR="true"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --human-ref) HUMAN_REF="$2"; shift 2 ;;
    --aligner) ALIGNER="$2"; shift 2 ;;
    --strand|--strandedness)
      if [[ "${STRAND_SOURCE}" != "default" && "${STRANDEDNESS}" != "$2" ]]; then
        echo "[run_auto_rnaseq] 错误：链特异性参数冲突：已有 ${STRAND_SOURCE}=${STRANDEDNESS}，又收到 $1=$2" >&2
        exit 1
      fi
      STRANDEDNESS="$2"
      STRAND_SOURCE="$1"
      shift 2
      ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --r1-pattern) R1_PATTERN="$2"; shift 2 ;;
    --r2-pattern) R2_PATTERN="$2"; shift 2 ;;
    --sample-metadata) SAMPLE_METADATA="$2"; shift 2 ;;
    --sample-meta) SAMPLE_META_ARGS+=("$2"); shift 2 ;;
    --condition-map) CONDITION_MAP_ARGS+=("$2"); shift 2 ;;
    --contrast-file) CONTRAST_FILE="$2"; shift 2 ;;
    --contrast) CONTRAST_ARGS+=("$2"); shift 2 ;;
    --tecount-strand)
      if [[ "${STRAND_SOURCE}" != "default" && "${STRANDEDNESS}" != "$2" ]]; then
        echo "[run_auto_rnaseq] 错误：--tecount-strand 现在与 --strand 合并，不能和 ${STRAND_SOURCE}=${STRANDEDNESS} 不一致" >&2
        exit 1
      fi
      STRANDEDNESS="$2"
      STRAND_SOURCE="$1"
      shift 2
      ;;
    --sra-threads) SRA_THREADS="$2"; shift 2 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --background-internal) RUN_SESSION_ISOLATED="true"; export RUN_SESSION_ISOLATED; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --resume-session) RESUME_SESSION="$2"; RESUME="true"; shift 2 ;;
    --upstream-only) RUN_DOWNSTREAM="false"; shift 1 ;;
    --downstream|--downstream-only) RUN_UPSTREAM="false"; RUN_DOWNSTREAM="true"; shift 1 ;;
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

    --padj-cutoff) PADJ_CUTOFF="$2"; shift 2 ;;
    --lfc-cutoff) LFC_CUTOFF="$2"; shift 2 ;;
    --baseMean-min) BASEMEAN_MIN="$2"; shift 2 ;;
    --label-top-n) LABEL_TOP_N="$2"; shift 2 ;;
    --heatmap-top-n) HEATMAP_TOP_N="$2"; shift 2 ;;
    --volcano-orientation) VOLCANO_ORIENTATION="$2"; shift 2 ;;
    --gray-nonsig) GRAY_NONSIG="$2"; shift 2 ;;
    --exploratory-method) EXPLORATORY_METHOD="$2"; shift 2 ;;
    --exploratory-fixed-bcv) EXPLORATORY_FIXED_BCV="$2"; shift 2 ;;
    --plot-threads) PLOT_THREADS="$2"; shift 2 ;;
    --only-tools) ONLY_TOOLS="$2"; shift 2 ;;
    --skip-tools) SKIP_TOOLS="$2"; shift 2 ;;
    --partial-input-policy) PARTIAL_INPUT_POLICY="$2"; shift 2 ;;
    --allow-partial-inputs)
      if [[ "${2,,}" == "true" ]]; then PARTIAL_INPUT_POLICY="allow"; else PARTIAL_INPUT_POLICY="skip"; fi
      shift 2 ;;

    --default-condition) DEFAULT_CONDITION="$2"; shift 2 ;;
    --default-replicate) DEFAULT_REPLICATE="$2"; shift 2 ;;
    --replace-design) REPLACE_DESIGN="true"; shift ;;
    --rnaseq-env) RNASEQ_ENV="$2"; shift 2 ;;
    --downstream-env) DOWNSTREAM_ENV="$2"; shift 2 ;;
    --rediscoverte-rollup-conda-prefix) ROLLUP_CONDA_PREFIX="$2"; shift 2 ;;
    --picard-markdup-java-heap) PICARD_MARKDUP_JAVA_HEAP="$2"; shift 2 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[run_auto_rnaseq] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "${ROLLUP_CONDA_PREFIX}" ]] || ROLLUP_CONDA_PREFIX="${DOWNSTREAM_ENV}"

resolve_default_layout_base() {
  if [[ -n "${FASTQ_DIR}" ]]; then
    local abs_fastq
    abs_fastq="$(cd "$(dirname "${FASTQ_DIR}")" && pwd)/$(basename "${FASTQ_DIR}")"
    AUTO_LAYOUT_BASE="$(dirname "${abs_fastq}")"
  elif [[ -n "${MANIFEST}" ]]; then
    local abs_manifest
    abs_manifest="$(cd "$(dirname "${MANIFEST}")" && pwd)/$(basename "${MANIFEST}")"
    AUTO_LAYOUT_BASE="$(dirname "${abs_manifest}")"
  else
    AUTO_LAYOUT_BASE="${PWD}"
  fi
}

resolve_default_layout_base

if [[ "${USER_SET_RESULTS_DIR}" == "true" && "${USER_SET_WORK_ROOT}" == "false" ]]; then
  RESULTS_PARENT="$(dirname "${RESULTS_DIR}")"
  if ! mkdir -p "${RESULTS_PARENT}" 2>/dev/null; then
    echo "[run_auto_rnaseq] 错误：无法创建 results-dir 的父目录: ${RESULTS_PARENT}" >&2
    echo "[run_auto_rnaseq] 提示：请把 --results-dir 放到当前用户可写目录。" >&2
    exit 1
  fi
  WORK_ROOT="$(cd "${RESULTS_PARENT}" && pwd)"
  RESULTS_DIR="${WORK_ROOT}/$(basename "${RESULTS_DIR}")"
elif [[ "${USER_SET_WORK_ROOT}" == "false" ]]; then
  WORK_ROOT="${AUTO_LAYOUT_BASE}/rnaseq_auto_run"
fi
[[ -n "${RESULTS_DIR}" ]] || RESULTS_DIR="${AUTO_LAYOUT_BASE}/rnaseq_results"
[[ -n "${PLOT_OUTDIR}" ]] || PLOT_OUTDIR="${RESULTS_DIR}/plots"

AUTOMATION_DIR="${RESULTS_DIR}/_automation"
WORK_DIR="${AUTOMATION_DIR}/work"
INPUTS_DIR="${AUTOMATION_DIR}/inputs"
LOG_DIR="${AUTOMATION_DIR}/logs"
NXF_HOME_DIR="/tmp/${USER:-$(id -un)}/nextflow/home"
# Keep Nextflow's launcher/history lock on local scratch.  The result and work
# directories may be on NFS, but FileChannel.lock() in the launch directory can
# block indefinitely on some NFS mounts.
NF_LAUNCH_DIR="/tmp/${USER:-$(id -un)}/nextflow/launch/$(basename "${RESULTS_DIR}")"
CONDA_CACHE_DIR="${AUTOMATION_DIR}/conda_cache"
TMP_DIR="${AUTOMATION_DIR}/tmp"

[[ -n "${DOWNLOAD_DIR}" ]] || DOWNLOAD_DIR="${AUTOMATION_DIR}/fastq"

AUTOMATION_LOG="${LOG_DIR}/automation_${RUN_TS}.log"
AUTOMATION_PID_FILE="${LOG_DIR}/automation_${RUN_TS}.pid"
PIPELINE_LOG_FILE="${LOG_DIR}/pipeline_${RUN_TS}.log"
SAMPLESHEET="${INPUTS_DIR}/samplesheet.csv"
CONTRAST_OUT="${RESULTS_DIR}/contrast.csv"
LEGACY_CONTRAST_OUT="${INPUTS_DIR}/contrasts.csv"
SAMPLE_METADATA_AUTO="${RESULTS_DIR}/condition.csv"
LEGACY_SAMPLE_METADATA_AUTO="${INPUTS_DIR}/sample_metadata.auto.csv"
RESOLVED_MANIFEST_OUT="${INPUTS_DIR}/manifest.resolved.csv"

ensure_writable_dir() {
  local dir="$1"
  local label="$2"
  if ! mkdir -p "${dir}" 2>/dev/null; then
    echo "[run_auto_rnaseq] 错误：无法创建 ${label}: ${dir}" >&2
    echo "[run_auto_rnaseq] 提示：如果 FASTQ 在别人目录下，请显式设置 --results-dir 到当前用户可写目录。" >&2
    exit 1
  fi
  if [[ ! -w "${dir}" ]]; then
    echo "[run_auto_rnaseq] 错误：当前用户不能写 ${label}: ${dir}" >&2
    echo "[run_auto_rnaseq] 提示：请换一个 --results-dir/--work-root，或让目录所有者开放写权限。" >&2
    exit 1
  fi
}

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

handle_auto_signal() {
  log_msg "收到终止信号，清理自动化流程的子进程。"
  terminate_process_group "$$"
  terminate_conda_for_cache
  exit 143
}

trap 'handle_auto_signal' INT TERM HUP

summarize_nextflow_failures() {
  local trace_file summary_file status_col process_col name_col exit_col
  trace_file="$(ls -1t "${LOG_DIR}"/nextflow_*.trace.tsv 2>/dev/null | head -n 1 || true)"
  summary_file="${AUTOMATION_DIR}/failed_tasks_${RUN_TS}.tsv"
  if [[ -z "${trace_file}" || ! -s "${trace_file}" ]]; then
    return 1
  fi
  awk -F '\t' -v OFS='\t' -v out="${summary_file}" '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "status") s=i
        if ($i == "process") p=i
        if ($i == "name") n=i
        if ($i == "exit") e=i
      }
      print "process", "status", "exit", "name" > out
      next
    }
    s && $s != "" && $s !~ /^(COMPLETED|CACHED)$/ {
      proc = (p ? $p : "")
      name = (n ? $n : "")
      code = (e ? $e : "")
      print proc, $s, code, name >> out
      count++
    }
    END { print count + 0 }
  ' "${trace_file}" > "${summary_file}.count"
  local failed_count
  failed_count="$(tail -n 1 "${summary_file}.count")"
  rm -f "${summary_file}.count"
  if [[ "${failed_count}" =~ ^[1-9][0-9]*$ ]]; then
    log_msg "Nextflow 存在 ${failed_count} 个非 COMPLETED/CACHED 任务，状态记为 PARTIAL_SUCCESS。失败清单: ${summary_file}" | tee -a "${AUTOMATION_LOG}"
    return 0
  fi
  rm -f "${summary_file}"
  return 1
}

if [[ "${BACKGROUND}" == "true" ]]; then
  if [[ "${DRY_RUN}" == "true" ]]; then
    BACKGROUND="false"
    log_msg "dry-run 使用前台输出，忽略 --background"
  else
  ensure_writable_dir "${LOG_DIR}" "automation log dir"
  BG_LOG="${AUTOMATION_LOG}"
  if ! touch "${BG_LOG}" "${AUTOMATION_PID_FILE}" 2>/dev/null; then
    echo "[run_auto_rnaseq] 错误：无法写后台日志或 PID 文件：" >&2
    echo "[run_auto_rnaseq]   ${BG_LOG}" >&2
    echo "[run_auto_rnaseq]   ${AUTOMATION_PID_FILE}" >&2
    echo "[run_auto_rnaseq] 提示：请添加 --results-dir 到当前用户可写目录，例如 /path/to/lab-data/RNAseq/<project>/results。" >&2
    exit 1
  fi
  FILTERED_ARGS=()
  for x in "${ORIG_ARGS[@]}"; do
    [[ "$x" == "--background" ]] && continue
    FILTERED_ARGS+=("$x")
  done
  log_msg "后台运行，日志: ${BG_LOG}"
  nohup setsid bash "$0" "${FILTERED_ARGS[@]}" --background-internal > "${BG_LOG}" 2>&1 &
  echo $! > "${AUTOMATION_PID_FILE}"
  log_msg "PID 文件: ${AUTOMATION_PID_FILE}"
  log_msg "PID: $(cat "${AUTOMATION_PID_FILE}")"
  log_msg "PGID: $(ps -o pgid= -p "$(cat "${AUTOMATION_PID_FILE}")" 2>/dev/null | tr -d ' ' || true)"
  exit 0
  fi
fi

SPECIES="$(echo "${SPECIES}" | tr '[:upper:]' '[:lower:]')"
HUMAN_REF="$(echo "${HUMAN_REF}" | tr '[:upper:]' '[:lower:]')"
ALIGNER="$(echo "${ALIGNER}" | tr '[:upper:]' '[:lower:]')"
STRANDEDNESS="$(echo "${STRANDEDNESS}" | tr '[:upper:]' '[:lower:]')"
TECOUNT_STRAND="${STRANDEDNESS}"
case "${SPECIES}" in
  hg38|mm10|mm39) ;;
  *) echo "[run_auto_rnaseq] 错误：--species 只能是 hg38 | mm10 | mm39" >&2; exit 1 ;;
esac
case "${HUMAN_REF}" in
  hg38|t2t) ;;
  *) echo "[run_auto_rnaseq] 错误：--human-ref 只能是 hg38 | t2t" >&2; exit 1 ;;
esac
case "${ALIGNER}" in
  star|hisat2) ;;
  *) echo "[run_auto_rnaseq] 错误：--aligner 只能是 star | hisat2" >&2; exit 1 ;;
esac
case "${STRANDEDNESS}" in
  unstranded|forward|reverse) ;;
  *) echo "[run_auto_rnaseq] 错误：--strand 只能是 unstranded | forward | reverse" >&2; exit 1 ;;
esac
if [[ "${SPECIES}" != "hg38" && "${HUMAN_REF}" != "hg38" ]]; then
  echo "[run_auto_rnaseq] 错误：--human-ref 仅在 --species hg38 时可设为 t2t" >&2
  exit 1
fi

ensure_writable_dir "${RESULTS_DIR}" "results-dir"
ensure_writable_dir "${PLOT_OUTDIR}" "plot-outdir"
ensure_writable_dir "${INPUTS_DIR}" "automation inputs dir"
ensure_writable_dir "${LOG_DIR}" "automation log dir"
ensure_writable_dir "${WORK_DIR}" "Nextflow work dir"
ensure_writable_dir "${NXF_HOME_DIR}" "Nextflow home dir"
ensure_writable_dir "${NF_LAUNCH_DIR}" "Nextflow launch dir"
ensure_writable_dir "${CONDA_CACHE_DIR}" "conda cache dir"
ensure_writable_dir "${TMP_DIR}" "tmp dir"
ensure_writable_dir "${DOWNLOAD_DIR}" "download dir"

[[ "${PROGRESS_INTERVAL}" =~ ^[0-9]+$ ]] && (( PROGRESS_INTERVAL > 0 )) || {
  echo "[run_auto_rnaseq] 错误：--progress-interval 必须是正整数" >&2
  exit 1
}
PARTIAL_INPUT_POLICY="${PARTIAL_INPUT_POLICY,,}"
case "${PARTIAL_INPUT_POLICY}" in
  skip|error|allow) ;;
  *) echo "[run_auto_rnaseq] 错误：--partial-input-policy 只能是 skip | error | allow" >&2; exit 1 ;;
esac
case "${FAILURE_POLICY}" in core|strict) ;; *) echo "[run_auto_rnaseq] 错误：--failure-policy 只能是 core | strict" >&2; exit 1 ;; esac
if [[ "${RUN_REDISCOVERTE_ROLLUP,,}" == "true" && ! -d "${ROLLUP_CONDA_PREFIX}" ]]; then
  echo "[run_auto_rnaseq] 错误：REdiscoverTE rollup 环境不存在: ${ROLLUP_CONDA_PREFIX}" >&2
  echo "[run_auto_rnaseq] 提示：用 --rediscoverte-rollup-conda-prefix DIR 指定环境，或关闭 --run-rediscoverte-rollup false。" >&2
  exit 1
fi

log_msg "run timestamp: ${RUN_TS}"
log_msg "automation log: ${AUTOMATION_LOG}"
log_msg "pipeline log: ${PIPELINE_LOG_FILE}"
log_msg "automation pid file: ${AUTOMATION_PID_FILE}"
log_msg "automation dir: ${AUTOMATION_DIR}"
log_msg "nextflow work dir: ${WORK_DIR}"
log_msg "REdiscoverTE rollup env: ${ROLLUP_CONDA_PREFIX}"
log_msg "nextflow launch dir: ${NF_LAUNCH_DIR}"
log_msg "gene aligner: ${ALIGNER}; human ref: ${HUMAN_REF}"
log_msg "进度文件: ${LOG_DIR}/progress_latest.txt（每 ${PROGRESS_INTERVAL} 秒刷新）"
[[ "${USER_SET_PLOT_OUTDIR}" == "false" ]] || mkdir -p "${PLOT_OUTDIR}"

resolve_refs() {
  local species="$1"
  local human_ref="$2"
  case "$species" in
    hg38)
      if [[ "$human_ref" == "t2t" ]]; then
        TX2GENE_PATH='/path/to/reference/human_t2t/chm13v2.0.annotation.gtf'
        TE_ANNOTATION_TSV='/path/to/reference/TE_GTF/hs1_te_annotation.tsv'
      else
        TX2GENE_PATH='/path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf'
        TE_ANNOTATION_TSV='/path/to/reference/TE_GTF/hg38_te_annotation.tsv'
      fi
      ;;
    mm10)
      TX2GENE_PATH='/path/to/reference/mouse_mm10/gencode.vM25.primary_assembly.annotation.gtf'
      TE_ANNOTATION_TSV='/path/to/reference/TE_GTF/mm10_te_annotation.tsv'
      ;;
    mm39)
      TX2GENE_PATH='/path/to/reference/mouse_mm39/gencode.vM38.primary_assembly.annotation.gtf'
      TE_ANNOTATION_TSV='/path/to/reference/TE_GTF/mm39_te_annotation.tsv'
      ;;
  esac
}
resolve_refs "${SPECIES}" "${HUMAN_REF}"
if [[ "${RUN_DOWNSTREAM}" == "true" ]]; then
  [[ -f "${TX2GENE_PATH}" ]] || { echo "[run_auto_rnaseq] 错误：未找到 GTF: ${TX2GENE_PATH}" >&2; exit 1; }
  [[ -f "${TE_ANNOTATION_TSV}" ]] || { echo "[run_auto_rnaseq] 错误：未找到 TE 注释: ${TE_ANNOTATION_TSV}" >&2; exit 1; }
fi
if [[ "${HUMAN_REF}" == "t2t" && "${RUN_SALMONTE}" == "true" ]]; then
  log_msg "提示：T2T 模式下 SalmonTE 仍使用其内置 hs 参考；如需完全一致，建议手动关闭 --run-salmonte。" | tee -a "${AUTOMATION_LOG}"
fi

build_metadata_from_cli() {
  python3 - <<'PY' "${SAMPLE_METADATA_AUTO}" "${DEFAULT_CONDITION}" "${DEFAULT_REPLICATE}" "${SAMPLE_METADATA}" "${#SAMPLE_META_ARGS[@]}" "${#CONDITION_MAP_ARGS[@]}" "${SAMPLE_META_ARGS[@]:-}" --CONDSEP-- "${CONDITION_MAP_ARGS[@]:-}"
import csv, sys
from collections import OrderedDict

out_csv = sys.argv[1]
default_condition = sys.argv[2]
default_replicate = sys.argv[3]
base_file = sys.argv[4]
n_sample_meta = int(sys.argv[5])
n_condition_map = int(sys.argv[6])
rest = sys.argv[7:]
sep_idx = rest.index('--CONDSEP--')
sample_meta_args = rest[:sep_idx]
condition_map_args = rest[sep_idx+1:]

meta = OrderedDict()

def add_row(sample, condition, replicate):
    sample = str(sample).strip()
    if not sample:
        return
    meta[sample] = {
        'sample': sample,
        'condition': str(condition).strip() or default_condition,
        'replicate': str(replicate).strip() or default_replicate,
    }

if base_file:
    with open(base_file, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        fields = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if 'sample' not in fields:
            raise SystemExit('sample_metadata 文件缺少 sample 列')
        s_col = fields['sample']
        c_col = fields.get('condition')
        r_col = fields.get('replicate')
        for row in reader:
            add_row(row.get(s_col, ''), row.get(c_col, default_condition) if c_col else default_condition, row.get(r_col, default_replicate) if r_col else default_replicate)

for raw in sample_meta_args[:n_sample_meta]:
    parts = raw.split(':')
    if len(parts) != 3:
        raise SystemExit(f'--sample-meta 格式必须为 sample:condition:replicate，收到: {raw}')
    add_row(parts[0], parts[1], parts[2])

for raw in condition_map_args[:n_condition_map]:
    parts = raw.split(':', 1)
    if len(parts) != 2:
        raise SystemExit(f'--condition-map 格式必须为 condition:sample1,sample2,...，收到: {raw}')
    cond, sample_blob = parts
    samples = [x.strip() for x in sample_blob.split(',') if x.strip()]
    for idx, sample in enumerate(samples, start=1):
        replicate = idx
        if sample in meta and str(meta[sample].get('replicate', '')).strip() not in {'', default_replicate, 'NA'}:
            replicate = meta[sample]['replicate']
        add_row(sample, cond, replicate)

with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['sample', 'condition', 'replicate'])
    writer.writeheader()
    for row in meta.values():
        writer.writerow(row)
print(f'[run_auto_rnaseq] 已生成 sample metadata: {out_csv}')
PY
}

build_contrast_file() {
  if [[ -n "${CONTRAST_FILE}" ]]; then
    [[ -f "${CONTRAST_FILE}" ]] || { echo "[run_auto_rnaseq] 错误：contrast 文件不存在: ${CONTRAST_FILE}" >&2; exit 1; }
    cp "${CONTRAST_FILE}" "${CONTRAST_OUT}"
    return 0
  fi

  python3 - <<'PY' "${CONTRAST_OUT}" "${#CONTRAST_ARGS[@]}" "${CONTRAST_ARGS[@]:-}"
import csv, sys
out_csv = sys.argv[1]
n = int(sys.argv[2])
items = sys.argv[3:3+n]
if n == 0:
    raise SystemExit('需要 --contrast-file 或至少一个 --contrast')
rows = []
for raw in items:
    parts = raw.split(':')
    if len(parts) not in (2, 3):
        raise SystemExit(f'--contrast 格式必须为 CASE:CONTROL 或兼容旧格式 group_col:CASE:CONTROL，收到: {raw}')
    if len(parts) == 2:
        rows.append({'case': parts[0], 'control': parts[1]})
    else:
        rows.append({'case': parts[1], 'control': parts[2]})
with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['case', 'control'])
    writer.writeheader()
    writer.writerows(rows)
print(f'[run_auto_rnaseq] 已生成 contrasts: {out_csv}')
PY
}

backup_design_file_if_needed() {
  local f="$1"
  if [[ -f "${f}" && "${REPLACE_DESIGN}" == "true" ]]; then
    local backup="${f}.bak_${RUN_TS}"
    cp "${f}" "${backup}"
    log_msg "检测到已有设计文件，已备份: ${backup}" | tee -a "${AUTOMATION_LOG}"
  fi
}

write_condition_template_from_samplesheet() {
  [[ -f "${SAMPLESHEET}" ]] || return 0
  if [[ -f "${SAMPLE_METADATA_AUTO}" && "${REPLACE_DESIGN}" != "true" ]]; then
    log_msg "保留已有 condition.csv（使用 --replace-design 才会覆盖）: ${SAMPLE_METADATA_AUTO}"
    return 0
  fi
  backup_design_file_if_needed "${SAMPLE_METADATA_AUTO}"
  python3 - <<'PY' "${SAMPLESHEET}" "${SAMPLE_METADATA_AUTO}" "${LEGACY_SAMPLE_METADATA_AUTO}"
import csv, sys
samplesheet, out_csv, legacy_csv = sys.argv[1:4]
rows = []
seen = set()
with open(samplesheet, newline='', encoding='utf-8-sig') as fh:
    reader = csv.DictReader(fh)
    need = ['sample', 'condition', 'replicate']
    missing = [c for c in need if c not in (reader.fieldnames or [])]
    if missing:
        raise SystemExit(f'samplesheet 缺少列: {missing}')
    for row in reader:
        sample = str(row.get('sample', '')).strip()
        if not sample or sample in seen:
            continue
        seen.add(sample)
        rows.append({
            'sample': sample,
            'condition': str(row.get('condition', '')).strip() or 'NA',
            'replicate': str(row.get('replicate', '')).strip() or 'NA'
        })
for path in [out_csv, legacy_csv]:
    with open(path, 'w', newline='', encoding='utf-8') as wf:
        writer = csv.DictWriter(wf, fieldnames=['sample', 'condition', 'replicate'])
        writer.writeheader()
        writer.writerows(rows)
print(f'[run_auto_rnaseq] 已生成 condition 模板: {out_csv}')
PY
}

write_contrast_template_from_condition() {
  [[ -f "${SAMPLE_METADATA_AUTO}" ]] || return 0
  if [[ -f "${CONTRAST_OUT}" && "${REPLACE_DESIGN}" != "true" ]]; then
    log_msg "保留已有 contrast.csv（使用 --replace-design 才会覆盖）: ${CONTRAST_OUT}"
    return 0
  fi
  backup_design_file_if_needed "${CONTRAST_OUT}"
  python3 - <<'PY' "${SAMPLE_METADATA_AUTO}" "${CONTRAST_OUT}" "${LEGACY_CONTRAST_OUT}"
import csv, itertools, sys
cond_csv, out_csv, legacy_csv = sys.argv[1:4]
conds = []
with open(cond_csv, newline='', encoding='utf-8-sig') as fh:
    reader = csv.DictReader(fh)
    if 'condition' not in (reader.fieldnames or []):
        raise SystemExit('condition.csv 缺少 condition 列')
    for row in reader:
        cond = str(row.get('condition', '')).strip()
        if cond and cond.upper() != 'NA' and cond not in conds:
            conds.append(cond)
rows = [{'case': a, 'control': b} for a, b in itertools.permutations(conds, 2)]
for path in [out_csv, legacy_csv]:
    with open(path, 'w', newline='', encoding='utf-8') as wf:
        writer = csv.DictWriter(wf, fieldnames=['case', 'control'])
        writer.writeheader()
        writer.writerows(rows)
print(f'[run_auto_rnaseq] 已生成 contrast 模板: {out_csv}; rows={len(rows)}')
PY
}

prepare_design_files_for_downstream() {
  if [[ -f "${LEGACY_SAMPLE_METADATA_AUTO}" && ! -f "${SAMPLE_METADATA_AUTO}" ]]; then
    cp "${LEGACY_SAMPLE_METADATA_AUTO}" "${SAMPLE_METADATA_AUTO}"
    log_msg "已从旧路径恢复 condition.csv: ${SAMPLE_METADATA_AUTO}" | tee -a "${AUTOMATION_LOG}"
  fi
  if [[ -f "${LEGACY_CONTRAST_OUT}" && ! -f "${CONTRAST_OUT}" ]]; then
    cp "${LEGACY_CONTRAST_OUT}" "${CONTRAST_OUT}"
    log_msg "已从旧路径恢复 contrast.csv: ${CONTRAST_OUT}" | tee -a "${AUTOMATION_LOG}"
  fi
  if [[ ! -f "${SAMPLE_METADATA_AUTO}" && -f "${SAMPLESHEET}" ]]; then
    write_condition_template_from_samplesheet | tee -a "${AUTOMATION_LOG}"
  fi
  if [[ ! -f "${CONTRAST_OUT}" && -f "${SAMPLE_METADATA_AUTO}" ]]; then
    write_contrast_template_from_condition | tee -a "${AUTOMATION_LOG}"
  fi
}

prepare_fastq_and_metadata() {
  if [[ -n "${MANIFEST}" ]]; then
    [[ -f "${MANIFEST}" ]] || { echo "[run_auto_rnaseq] 错误：manifest 不存在: ${MANIFEST}" >&2; exit 1; }
    [[ -x "${MANIFEST_HELPER}" ]] || { echo "[run_auto_rnaseq] 错误：manifest helper 不存在: ${MANIFEST_HELPER}" >&2; exit 1; }
    mkdir -p "${DOWNLOAD_DIR}"
    log_msg "步骤1：manifest -> FASTQ（下载/整理）" | tee -a "${AUTOMATION_LOG}"
    python3 "${MANIFEST_HELPER}" \
      --manifest "${MANIFEST}" \
      --fastq-dir "${DOWNLOAD_DIR}" \
      --sample-metadata-out "${SAMPLE_METADATA_AUTO}" \
      --resolved-manifest-out "${RESOLVED_MANIFEST_OUT}" \
      --sra-threads "${SRA_THREADS}" | tee -a "${AUTOMATION_LOG}"

    if [[ -n "${SAMPLE_METADATA}" || ${#SAMPLE_META_ARGS[@]} -gt 0 || ${#CONDITION_MAP_ARGS[@]} -gt 0 ]]; then
      local manifest_meta_bak="${INPUTS_DIR}/sample_metadata.from_manifest.csv"
      local merged_meta="${INPUTS_DIR}/sample_metadata.merged_base.csv"
      cp "${SAMPLE_METADATA_AUTO}" "${manifest_meta_bak}"
      if [[ -n "${SAMPLE_METADATA}" && -f "${SAMPLE_METADATA}" ]]; then
        python3 - <<'PY2' "${manifest_meta_bak}" "${SAMPLE_METADATA}" "${merged_meta}"
import csv, sys
base_f, user_f, out_f = sys.argv[1:4]
rows = []
for idx, path in enumerate([base_f, user_f]):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if 'sample' not in cols:
            raise SystemExit(f'metadata 文件缺少 sample 列: {path}')
        s_col = cols['sample']
        c_col = cols.get('condition')
        r_col = cols.get('replicate')
        for row in reader:
            sample = str(row.get(s_col, '')).strip()
            if not sample:
                continue
            rows.append({'sample': sample, 'condition': str(row.get(c_col, '')).strip() if c_col else '', 'replicate': str(row.get(r_col, '')).strip() if r_col else ''})
with open(out_f, 'w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['sample','condition','replicate'])
    writer.writeheader()
    writer.writerows(rows)
PY2
        SAMPLE_METADATA="${merged_meta}"
      else
        SAMPLE_METADATA="${manifest_meta_bak}"
      fi
      build_metadata_from_cli | tee -a "${AUTOMATION_LOG}"
    fi
    FASTQ_DIR="${DOWNLOAD_DIR}"
  else
    [[ -n "${FASTQ_DIR}" ]] || { echo "[run_auto_rnaseq] 错误：运行上游时必须提供 --manifest 或 --fastq-dir" >&2; exit 1; }
    [[ -d "${FASTQ_DIR}" ]] || { echo "[run_auto_rnaseq] 错误：fastq 目录不存在: ${FASTQ_DIR}" >&2; exit 1; }
    if [[ -n "${SAMPLE_METADATA}" || ${#SAMPLE_META_ARGS[@]} -gt 0 || ${#CONDITION_MAP_ARGS[@]} -gt 0 ]]; then
      build_metadata_from_cli | tee -a "${AUTOMATION_LOG}"
    fi
  fi

  log_msg "步骤2：生成 samplesheet" | tee -a "${AUTOMATION_LOG}"
  GEN_CMD=(bash "${SAMPLESHEET_GEN}"
    --input "${FASTQ_DIR}"
    --output "${SAMPLESHEET}"
    --layout "${LAYOUT}"
    --r1-pattern "${R1_PATTERN}"
    --r2-pattern "${R2_PATTERN}"
    --condition "${DEFAULT_CONDITION}"
    --replicate "${DEFAULT_REPLICATE}")
  if [[ -f "${SAMPLE_METADATA_AUTO}" ]]; then
    GEN_CMD+=(--metadata-csv "${SAMPLE_METADATA_AUTO}")
  elif [[ -n "${SAMPLE_METADATA}" && -f "${SAMPLE_METADATA}" ]]; then
    GEN_CMD+=(--metadata-csv "${SAMPLE_METADATA}")
  fi
  "${GEN_CMD[@]}" | tee -a "${AUTOMATION_LOG}"
  write_condition_template_from_samplesheet | tee -a "${AUTOMATION_LOG}"
  if [[ ! -f "${CONTRAST_OUT}" ]]; then
    write_contrast_template_from_condition | tee -a "${AUTOMATION_LOG}"
  fi
}

if [[ "${RUN_UPSTREAM}" == "true" ]]; then
  prepare_fastq_and_metadata
fi

if [[ "${RUN_DOWNSTREAM}" == "true" ]]; then
  if [[ -f "${SAMPLE_METADATA_AUTO}" ]]; then
    :
  elif [[ -n "${SAMPLE_METADATA}" && -f "${SAMPLE_METADATA}" ]]; then
    cp "${SAMPLE_METADATA}" "${SAMPLE_METADATA_AUTO}"
  elif [[ "${RUN_UPSTREAM}" == "false" && -f "${SAMPLESHEET}" ]]; then
    write_condition_template_from_samplesheet | tee -a "${AUTOMATION_LOG}"
  elif [[ ${#SAMPLE_META_ARGS[@]} -gt 0 || ${#CONDITION_MAP_ARGS[@]} -gt 0 ]]; then
    build_metadata_from_cli | tee -a "${AUTOMATION_LOG}"
  fi
  prepare_design_files_for_downstream
  if [[ -n "${CONTRAST_FILE}" || ${#CONTRAST_ARGS[@]} -gt 0 ]]; then
    build_contrast_file | tee -a "${AUTOMATION_LOG}"
  fi
fi

if [[ "${RUN_UPSTREAM}" == "true" ]]; then
  UP_CMD=(bash "${PIPELINE_RUNNER}"
    --samplesheet "${SAMPLESHEET}"
    --species "${SPECIES}"
    --strand "${STRANDEDNESS}"
    --outdir "${RESULTS_DIR}"
    --work-dir "${WORK_DIR}"
    --log-dir "${LOG_DIR}"
    --log-file "${PIPELINE_LOG_FILE}"
    --nxf-home "${NXF_HOME_DIR}"
    --nf-launch-dir "${NF_LAUNCH_DIR}"
    --conda-cache-dir "${CONDA_CACHE_DIR}"
    --tmp-dir "${TMP_DIR}"
    --progress-interval "${PROGRESS_INTERVAL}"
    --failure-policy "${FAILURE_POLICY}"
    --run-fastp "${RUN_FASTP}"
    --run-gene-count-branch "${RUN_STAR_FC}"
    --run-salmon "${RUN_SALMON}"
    --run-stringtie "${RUN_STRINGTIE}"
    --run-tecount "${RUN_TECOUNT}"
    --run-telocal "${RUN_TELOCAL}"
    --run-tetranscripts "${RUN_TETRANSCRIPTS}"
    --run-dedup "${RUN_DEDUP}"
    --run-markdup-qc "${RUN_MARKDUP_QC}"
    --run-rediscoverte "${RUN_REDISCOVERTE}"
    --run-rediscoverte-rollup "${RUN_REDISCOVERTE_ROLLUP}"
    --rediscoverte-rollup-conda-prefix "${ROLLUP_CONDA_PREFIX}"
    --picard-markdup-java-heap "${PICARD_MARKDUP_JAVA_HEAP}"
    --run-salmonte "${RUN_SALMONTE}"
    --run-telescope "${RUN_TELESCOPE}" \
    --run-fastqc "${RUN_FASTQC}" \
    --run-rnaseq-metrics "${RUN_RNASEQ_METRICS}" \
    --run-multiqc "${RUN_MULTIQC}" \
    --multiqc-cmd "${MULTIQC_CMD}")
  [[ -n "${REF_FLAT}" ]] && UP_CMD+=(--ref-flat "${REF_FLAT}")
  [[ -n "${RIBOSOMAL_INTERVALS}" ]] && UP_CMD+=(--ribosomal-intervals "${RIBOSOMAL_INTERVALS}")
  [[ "${RESUME}" == "true" ]] && UP_CMD+=(--resume) || UP_CMD+=(--no-resume)
  [[ -z "${RESUME_SESSION}" ]] || UP_CMD+=(--resume-session "${RESUME_SESSION}")
  [[ -n "${PROFILE}" ]] && UP_CMD+=(--profile "${PROFILE}")
  [[ -n "${QUEUE_SIZE}" ]] && UP_CMD+=(--queue-size "${QUEUE_SIZE}")
  [[ -n "${MAX_CPUS}" ]] && UP_CMD+=(--max-cpus "${MAX_CPUS}")
  [[ -n "${MAX_MEMORY}" ]] && UP_CMD+=(--max-memory "${MAX_MEMORY}")
  [[ "${DRY_RUN}" == "true" ]] && UP_CMD+=(--dry-run)
  if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
    for x in "${EXTRA_ARGS[@]}"; do
      UP_CMD+=(--extra "${x}")
    done
  fi
  log_msg "步骤3：运行上游 Nextflow（开始）" | tee -a "${AUTOMATION_LOG}"
  log_msg "上游 wrapper 日志: ${PIPELINE_LOG_FILE}" | tee -a "${AUTOMATION_LOG}"
  log_msg "Nextflow 详细日志/报告目录: ${RESULTS_DIR}/_automation/logs" | tee -a "${AUTOMATION_LOG}"
  set +e
  conda run -p "${RNASEQ_ENV}" --no-capture-output "${UP_CMD[@]}" 2>&1 | tee -a "${AUTOMATION_LOG}"
  up_rc=${PIPESTATUS[0]}
  set -e
  log_msg "步骤3：上游结束，exit_code=${up_rc}" | tee -a "${AUTOMATION_LOG}"
  if [[ "${up_rc}" -ne 0 ]]; then
    UPSTREAM_STATUS="FAILED"
    log_msg "上游失败，请重点查看 ${PIPELINE_LOG_FILE} 与 ${RESULTS_DIR}/_automation/logs/nextflow_*.log" | tee -a "${AUTOMATION_LOG}"
    exit "${up_rc}"
  fi
  if summarize_nextflow_failures; then
    UPSTREAM_STATUS="PARTIAL_SUCCESS"
  else
    UPSTREAM_STATUS="SUCCESS"
    log_msg "上游成功，未发现被忽略的失败任务。可查看 trace/report/timeline：${RESULTS_DIR}/_automation/logs/nextflow_*.trace.tsv|*.report.html|*.timeline.html" | tee -a "${AUTOMATION_LOG}"
  fi
  if [[ "${DRY_RUN}" == "true" ]]; then
    log_msg "DRY-RUN 完成：未执行任何 Nextflow process。" | tee -a "${AUTOMATION_LOG}"
    log_msg "计划模块：featureCounts=${RUN_STAR_FC}, Salmon=${RUN_SALMON}, REdiscoverTE=${RUN_REDISCOVERTE}, REdiscoverTE-rollup=${RUN_REDISCOVERTE_ROLLUP}, TEtranscripts=${RUN_TETRANSCRIPTS}, Telescope=${RUN_TELESCOPE}, TElocal=${RUN_TELOCAL}, SalmonTE=${RUN_SALMONTE}" | tee -a "${AUTOMATION_LOG}"
    exit 0
  fi
else
  UPSTREAM_STATUS="SKIPPED"
  [[ -d "${RESULTS_DIR}" ]] || { echo "[run_auto_rnaseq] 错误：downstream-only 模式下 results-dir 不存在: ${RESULTS_DIR}" >&2; exit 1; }
  log_msg "跳过上游，仅使用已有结果目录: ${RESULTS_DIR}" | tee -a "${AUTOMATION_LOG}"
fi

if [[ "${RUN_DOWNSTREAM}" == "true" ]]; then
  if [[ "${DRY_RUN}" == "true" ]]; then
    log_msg "DRY-RUN downstream：将读取 ${SAMPLE_METADATA_AUTO} 与 ${CONTRAST_OUT}，输出 ${PLOT_OUTDIR}；未执行 R 分析。"
    exit 0
  fi
  mkdir -p "${PLOT_OUTDIR}"
  if [[ -z "${PLOT_THREADS}" ]]; then
    if [[ -n "${MAX_CPUS}" ]]; then
      PLOT_THREADS="${MAX_CPUS}"
      if [[ "${PLOT_THREADS}" =~ ^[0-9]+$ ]] && (( PLOT_THREADS > 8 )); then
        PLOT_THREADS="8"
      fi
    else
      PLOT_THREADS="1"
    fi
  fi
  if [[ -z "${PLOT_THREADS}" ]]; then
    PLOT_THREADS="1"
  fi
  SAMPLE_TABLE_FOR_DS="${SAMPLE_METADATA_AUTO}"
  [[ -f "${SAMPLE_TABLE_FOR_DS}" ]] || { echo "[run_auto_rnaseq] 错误：下游缺少 sample metadata: ${SAMPLE_TABLE_FOR_DS}" >&2; exit 1; }
  [[ -f "${CONTRAST_OUT}" ]] || { echo "[run_auto_rnaseq] 错误：下游缺少 contrasts: ${CONTRAST_OUT}" >&2; exit 1; }
  log_msg "步骤4：运行下游差异分析与可视化" | tee -a "${AUTOMATION_LOG}"
  DOWNSTREAM_TOOL_ARGS=()
  [[ -z "${ONLY_TOOLS}" ]] || DOWNSTREAM_TOOL_ARGS+=(--only-tools "${ONLY_TOOLS}")
  [[ -z "${SKIP_TOOLS}" ]] || DOWNSTREAM_TOOL_ARGS+=(--skip-tools "${SKIP_TOOLS}")
  DOWNSTREAM_TOOL_ARGS+=(--partial-input-policy "${PARTIAL_INPUT_POLICY}")
  conda run -p "${DOWNSTREAM_ENV}" --no-capture-output Rscript "${DOWNSTREAM_RUNNER}" \
    --results-dir "${RESULTS_DIR}" \
    --outdir "${PLOT_OUTDIR}" \
    --species "${SPECIES}" \
    --te-annotation-tsv "${TE_ANNOTATION_TSV}" \
    --tx2gene-path "${TX2GENE_PATH}" \
    --sample-table "${SAMPLE_TABLE_FOR_DS}" \
    --contrast-file "${CONTRAST_OUT}" \
    --tecount-strand "${TECOUNT_STRAND}" \
    --padj-cutoff "${PADJ_CUTOFF}" \
    --lfc-cutoff "${LFC_CUTOFF}" \
    --baseMean-min "${BASEMEAN_MIN}" \
    --label-top-n "${LABEL_TOP_N}" \
    --heatmap-top-n "${HEATMAP_TOP_N}" \
    --volcano-orientation "${VOLCANO_ORIENTATION}" \
    --gray-nonsig "${GRAY_NONSIG}" \
    --exploratory-method "${EXPLORATORY_METHOD}" \
    --exploratory-fixed-bcv "${EXPLORATORY_FIXED_BCV}" \
    --plot-threads "${PLOT_THREADS}" \
    "${DOWNSTREAM_TOOL_ARGS[@]}" | tee -a "${AUTOMATION_LOG}"
else
  log_msg "按要求跳过下游" | tee -a "${AUTOMATION_LOG}"
fi

log_msg "流程完成，状态=${UPSTREAM_STATUS}。inputs: ${INPUTS_DIR}" | tee -a "${AUTOMATION_LOG}"
log_msg "自动化中间目录: ${RESULTS_DIR}/_automation" | tee -a "${AUTOMATION_LOG}"
log_msg "上游结果目录: ${RESULTS_DIR}" | tee -a "${AUTOMATION_LOG}"
log_msg "下游绘图目录: ${PLOT_OUTDIR}" | tee -a "${AUTOMATION_LOG}"
