#!/usr/bin/env bash
set -euo pipefail

ORIG_ARGS=("$@")
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${SCRIPT_DIR}"
PIPELINE_RUNNER="${PROJECT_ROOT}/run_pipeline.sh"
MANIFEST_GEN="${PROJECT_ROOT}/bin/prepare_manifest_fastq.py"

FASTQ_DIR=""
SPECIES="hg38"
RESOURCE_TIER="standard"
ASSAY="cutrun"
LAYOUT="auto"
OUTDIR=""
WORK_DIR=""
MANIFEST=""
CONTROL_SAMPLE=""
CONTROL_REGEX=""
NO_AUTO_CONTROL="false"
FORCE_MANIFEST="false"
INIT_ONLY="false"
BACKGROUND="false"
BACKGROUND_INTERNAL="false"
RESUME="true"
DRY_RUN="false"
PREVIEW="false"
SKIP_PREFLIGHT="false"
PROFILE="conda"
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
T3E_MAX_BED_READS="1000000"
T3E_ITERATIONS="100"
ALLO_COMMAND=""
REPENTOOLS_COMMAND=""
REPENTOOLS_INDEX_DIR=""
REPENTOOLS_GTF=""
REPENTOOLS_RET=""
REPENTOOLS_TARGET_GROUPS=""
REPENTOOLS_INPUT_SAMPLES=""
MACS_QVALUE=""
MACS_PVALUE=""
MACS_BROAD_CUTOFF=""
EXTRA_ARGS=()

usage() {
  cat <<'USAGE'
ChIP-seq / CUT&RUN / CUT&Tag 自动入口：只指定 FASTQ 目录，自动生成 manifest 并运行 Nextflow。

常用：
  run_auto_chipseq --fastq-dir FASTQ_DIR --species hg38 --outdir results
  run_auto_chipseq --fastq-dir FASTQ_DIR --outdir results --init-only
  run_auto_chipseq --fastq-dir FASTQ_DIR --outdir results --preview

参数：
  --fastq-dir DIR              FASTQ 所在目录，递归扫描
  --species STR                hg38 | mm39，默认 hg38
  --resource-tier STR          downstream resources: small | standard | full
  --assay STR                  chipseq | cutrun | cuttag，默认 cutrun
  --layout STR                 auto | PE | SE，默认 auto
  --outdir DIR                 默认 <fastq-dir同级>/chipseq_auto_results
  --work-dir DIR               默认 <outdir>/_automation/work
  --manifest FILE              默认 <outdir>/_automation/inputs/manifest.csv
  --control-sample NAME        手动指定 IgG/Input/control 样本
  --control-regex REGEX        自动识别 control 的正则
  --no-auto-control            不自动给 target 样本填 igg 列
  --force-manifest             覆盖重建 manifest
  --init-only                  只生成 manifest，不启动流程
  --background                 后台运行
  --resume / --no-resume
  --preview                    nextflow -preview，不真正运行 jobs
  --dry-run                    只打印命令
  --skip-preflight             跳过 manifest/FASTQ/参考/磁盘检查（不建议）
  --profile STR                Nextflow profile，默认 conda
  --trim BOOL                  默认 true
  --run-analysis BOOL          默认 true
  --make-rpgc-track BOOL       默认 true
  --make-te-tracks BOOL        默认 true；为 TE/L1 relaxed BAM 生成 bigWig
  --make-te-locus-best-track BOOL
                               默认 true；生成可复现的单一最佳位置 5 bp CPM 轨道
  --te-k INT                   TE multi-mapping 数，默认 25
  --te-duplicate-policy STR    mark | keep | remove，默认 mark
  --te-remove-blacklist BOOL   默认 false
  --run-downstream             主流程完成后运行下游工具
  --te-methods LIST            可选 TE 适配器：t3e,allo,repentools
  --te-methods-execute         执行已安装的外部 TE 工具
  --te-repeat-bed FILE         T3E/Allo 使用的 RepeatMasker BED
  --te-saf FILE                TE SAF，用于自动生成 RepeatMasker BED
  --te-anno FILE               TE annotation TSV，用于自动生成 RepeatMasker BED
  --t3e-dir DIR                T3E 仓库目录
  --t3e-python FILE            T3E 隔离环境 Python
  --t3e-max-bed-reads INT      T3E BED 行数上限（默认 1000000；0=全部）
  --t3e-iterations INT         T3E 迭代次数（默认 100）
  --allo-command CMD           Allo 命令模板
  --repentools-command CMD     RepEnTools 命令模板
  --repentools-index-dir DIR   CHM13 index directory for FASTQ adapter
  --repentools-gtf FILE        RepeatMasker GTF for FASTQ adapter
  --repentools-ret FILE        RepEnTools ret executable
  --repentools-target-groups SPEC  group=chip1,chip2;... mapping
  --repentools-input-samples SPEC  input1,input2 mapping
  --macs-qvalue FLOAT           MACS -q cutoff；数值越小越严格
  --macs-pvalue FLOAT           MACS -p cutoff；设置后用 p-value 而不是 q-value
  --macs-broad-cutoff FLOAT     MACS --broad-cutoff；数值越小越严格
  --extra "ARGS"               继续传给 run_pipeline/Nextflow
USAGE
}

log_msg() { echo "[run_auto_chipseq] $*"; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
    --species) SPECIES="$2"; shift 2 ;;
    --resource-tier) RESOURCE_TIER="$2"; shift 2 ;;
    --assay) ASSAY="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --work-dir) WORK_DIR="$2"; shift 2 ;;
    --manifest) MANIFEST="$2"; shift 2 ;;
    --control-sample) CONTROL_SAMPLE="$2"; shift 2 ;;
    --control-regex) CONTROL_REGEX="$2"; shift 2 ;;
    --no-auto-control) NO_AUTO_CONTROL="true"; shift 1 ;;
    --force-manifest) FORCE_MANIFEST="true"; shift 1 ;;
    --init-only) INIT_ONLY="true"; shift 1 ;;
    --background) BACKGROUND="true"; shift 1 ;;
    --background-internal) BACKGROUND_INTERNAL="true"; BACKGROUND="false"; shift 1 ;;
    --resume) RESUME="true"; shift 1 ;;
    --no-resume) RESUME="false"; shift 1 ;;
    --preview) PREVIEW="true"; shift 1 ;;
    --dry-run) DRY_RUN="true"; shift 1 ;;
    --skip-preflight) SKIP_PREFLIGHT="true"; shift 1 ;;
    --profile) PROFILE="$2"; shift 2 ;;
    --trim) TRIM="$2"; shift 2 ;;
    --run-analysis) RUN_ANALYSIS="$2"; shift 2 ;;
    --make-rpgc-track) MAKE_RPGC_TRACK="$2"; shift 2 ;;
    --make-te-tracks|--make-te-track) MAKE_TE_TRACK="$2"; shift 2 ;;
    --make-te-locus-best-track|--make-te-locus-best-tracks) MAKE_TE_LOCUS_BEST_TRACK="$2"; shift 2 ;;
    --te-k) TE_K="$2"; shift 2 ;;
    --te-duplicate-policy) TE_DUPLICATE_POLICY="$2"; shift 2 ;;
    --te-remove-blacklist) TE_REMOVE_BLACKLIST="$2"; shift 2 ;;
    --run-downstream) RUN_DOWNSTREAM="true"; shift 1 ;;
    --te-methods) TE_METHODS="$2"; shift 2 ;;
    --te-methods-execute) TE_METHODS_EXECUTE="true"; shift ;;
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
    --macs-qvalue) MACS_QVALUE="$2"; shift 2 ;;
    --macs-pvalue) MACS_PVALUE="$2"; shift 2 ;;
    --macs-broad-cutoff) MACS_BROAD_CUTOFF="$2"; shift 2 ;;
    --extra) EXTRA_ARGS+=("$2"); shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[run_auto_chipseq] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$FASTQ_DIR" ]] || { echo "[run_auto_chipseq] ERROR: --fastq-dir is required" >&2; exit 1; }
[[ -d "$FASTQ_DIR" ]] || { echo "[run_auto_chipseq] ERROR: fastq dir not found: $FASTQ_DIR" >&2; exit 1; }
[[ -x "$PIPELINE_RUNNER" ]] || { echo "[run_auto_chipseq] ERROR: runner not executable: $PIPELINE_RUNNER" >&2; exit 1; }
[[ -f "$MANIFEST_GEN" ]] || { echo "[run_auto_chipseq] ERROR: manifest generator not found: $MANIFEST_GEN" >&2; exit 1; }

BASE_DIR="$(cd "$(dirname "$FASTQ_DIR")" && pwd)"
[[ -n "$OUTDIR" ]] || OUTDIR="${BASE_DIR}/chipseq_auto_results"
[[ -n "$WORK_DIR" ]] || WORK_DIR="${OUTDIR}/_automation/work"
INPUTS_DIR="${OUTDIR}/_automation/inputs"
LOG_DIR="${OUTDIR}/_automation/logs"
mkdir -p "$INPUTS_DIR" "$LOG_DIR"
[[ -n "$MANIFEST" ]] || MANIFEST="${INPUTS_DIR}/manifest.csv"
AUTOMATION_LOG="${LOG_DIR}/automation_$(date +%Y%m%d_%H%M%S).log"
PID_FILE="${LOG_DIR}/automation.pid"

if [[ "$BACKGROUND" == "true" && "$BACKGROUND_INTERNAL" != "true" && "${RUN_AUTO_CHIPSEQ_BG_CHILD:-0}" != "1" ]]; then
  FILTERED_ARGS=()
  for x in "${ORIG_ARGS[@]}"; do
    [[ "$x" == "--background" ]] && continue
    [[ "$x" == "--background-internal" ]] && continue
    FILTERED_ARGS+=("$x")
  done
  log_msg "后台运行，日志: ${AUTOMATION_LOG}"
  nohup env RUN_AUTO_CHIPSEQ_BG_CHILD=1 bash "$0" "${FILTERED_ARGS[@]}" --background-internal > "$AUTOMATION_LOG" 2>&1 &
  echo $! > "$PID_FILE"
  log_msg "PID: $(cat "$PID_FILE")"
  exit 0
fi

if [[ "$FORCE_MANIFEST" == "true" || ! -f "$MANIFEST" ]]; then
  CMD_GEN=(python "$MANIFEST_GEN" --fastq-dir "$FASTQ_DIR" --output "$MANIFEST" --species "$SPECIES" --assay "$ASSAY" --layout "$LAYOUT")
  if [[ -n "$CONTROL_SAMPLE" ]]; then CMD_GEN+=(--control-sample "$CONTROL_SAMPLE"); fi
  if [[ -n "$CONTROL_REGEX" ]]; then CMD_GEN+=(--control-regex "$CONTROL_REGEX"); fi
  if [[ "$NO_AUTO_CONTROL" == "true" ]]; then CMD_GEN+=(--no-auto-control); fi
  log_msg "扫描 FASTQ 并生成 manifest"
  log_msg "CMD=${CMD_GEN[*]}"
  "${CMD_GEN[@]}"
else
  log_msg "检测到已有 manifest，保留不覆盖: $MANIFEST"
fi

log_msg "manifest 预览:"
head -n 12 "$MANIFEST" || true

if [[ "$INIT_ONLY" == "true" ]]; then
  log_msg "初始化完成。你可以检查/编辑 manifest 后再运行。"
  exit 0
fi

CMD=(bash "$PIPELINE_RUNNER" --manifest "$MANIFEST" --outdir "$OUTDIR" --work-dir "$WORK_DIR" --profile "$PROFILE" --species "$SPECIES" --resource-tier "$RESOURCE_TIER" --trim "$TRIM" --run-analysis "$RUN_ANALYSIS" --make-rpgc-track "$MAKE_RPGC_TRACK" --make-te-tracks "$MAKE_TE_TRACK" --make-te-locus-best-track "$MAKE_TE_LOCUS_BEST_TRACK" --te-k "$TE_K" --te-duplicate-policy "$TE_DUPLICATE_POLICY" --te-remove-blacklist "$TE_REMOVE_BLACKLIST")
if [[ "$RUN_DOWNSTREAM" == "true" ]]; then CMD+=(--run-downstream); fi
if [[ -n "$TE_METHODS" ]]; then CMD+=(--te-methods "$TE_METHODS"); fi
if [[ "$TE_METHODS_EXECUTE" == "true" ]]; then CMD+=(--te-methods-execute); fi
[[ -n "$TE_REPEAT_BED" ]] && CMD+=(--te-repeat-bed "$TE_REPEAT_BED")
[[ -n "$TE_SAF" ]] && CMD+=(--te-saf "$TE_SAF")
[[ -n "$TE_ANNO" ]] && CMD+=(--te-anno "$TE_ANNO")
[[ -n "$T3E_DIR" ]] && CMD+=(--t3e-dir "$T3E_DIR")
[[ -n "$T3E_PYTHON" ]] && CMD+=(--t3e-python "$T3E_PYTHON")
CMD+=(--t3e-max-bed-reads "$T3E_MAX_BED_READS" --t3e-iterations "$T3E_ITERATIONS")
[[ -n "$ALLO_COMMAND" ]] && CMD+=(--allo-command "$ALLO_COMMAND")
[[ -n "$REPENTOOLS_COMMAND" ]] && CMD+=(--repentools-command "$REPENTOOLS_COMMAND")
[[ -n "$REPENTOOLS_INDEX_DIR" ]] && CMD+=(--repentools-index-dir "$REPENTOOLS_INDEX_DIR")
[[ -n "$REPENTOOLS_GTF" ]] && CMD+=(--repentools-gtf "$REPENTOOLS_GTF")
[[ -n "$REPENTOOLS_RET" ]] && CMD+=(--repentools-ret "$REPENTOOLS_RET")
[[ -n "$REPENTOOLS_TARGET_GROUPS" ]] && CMD+=(--repentools-target-groups "$REPENTOOLS_TARGET_GROUPS")
[[ -n "$REPENTOOLS_INPUT_SAMPLES" ]] && CMD+=(--repentools-input-samples "$REPENTOOLS_INPUT_SAMPLES")
if [[ -n "$MACS_QVALUE" ]]; then CMD+=(--macs-qvalue "$MACS_QVALUE"); fi
if [[ -n "$MACS_PVALUE" ]]; then CMD+=(--macs-pvalue "$MACS_PVALUE"); fi
if [[ -n "$MACS_BROAD_CUTOFF" ]]; then CMD+=(--macs-broad-cutoff "$MACS_BROAD_CUTOFF"); fi
if [[ "$RESUME" == "true" ]]; then CMD+=(--resume); else CMD+=(--no-resume); fi
if [[ "$PREVIEW" == "true" ]]; then CMD+=(--preview); fi
if [[ "$DRY_RUN" == "true" ]]; then CMD+=(--dry-run); fi
if [[ "$SKIP_PREFLIGHT" == "true" ]]; then CMD+=(--skip-preflight); fi
if [[ "$BACKGROUND_INTERNAL" == "true" ]]; then CMD+=(--background-internal); fi
if [[ ${#EXTRA_ARGS[@]} -gt 0 ]]; then
  for x in "${EXTRA_ARGS[@]}"; do CMD+=(--extra "$x"); done
fi

log_msg "启动流程"
log_msg "CMD=${CMD[*]}"
"${CMD[@]}"
