#!/usr/bin/env bash
set -euo pipefail

# Run the upstream RepEnTools FASTQ workflow on a resolved CUT&RUN manifest.
# RepEnTools is intentionally kept as a FASTQ/T2T adapter: its reference is
# CHM13/hs1 and its statistical model expects two ChIP and two input libraries.
# This wrapper stages symlinks with the names required by `ret`, records the
# exact command, and never silently treats hg38 BAMs as CHM13 data.

MANIFEST=""
OUT_DIR=""
INDEX_DIR="/path/to/CUTnRUN/resources/RepEnTools/chm13v2/indexes"
GTF_FILE="/path/to/CUTnRUN/resources/RepEnTools/chm13v2/annotation/rmsk.gtf"
RET_BIN="/path/to/.local/share/CUTnRUN-tools/RepEnTools/ret"
THREADS=8
# If omitted, both mappings are derived from the manifest.  Hard-coded sample
# names are intentionally not used because they can silently mix projects.
TARGET_GROUPS=""
INPUT_SAMPLES=""
RUN="false"
FORCE="false"

usage() {
  cat <<'USAGE'
Run RepEnTools using FASTQs named by a resolved CUT&RUN manifest.

Required:
  --manifest FILE       resolved_manifest.csv with fastq_1/fastq_2 columns
  --out-dir DIR         output root for staged inputs and result folders

Options:
  --index-dir DIR       directory containing chm13-2.{1..8}.ht2
  --gtf FILE            RepeatMasker-derived GTF (featureCounts exon/gene_id)
  --ret FILE            RepEnTools ret executable
  --target-groups SPEC  semicolon-separated group=sample1,sample2 pairs
                        (default: derive two target replicates per manifest group)
  --input-samples SPEC   comma-separated input samples (default: first two controls)
  --threads INT          threads passed to ret (default: 8)
  --force                rerun groups with completed report sentinels
  --execute              run ret; without this flag only stage/validate inputs
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --manifest) MANIFEST="$2"; shift 2 ;;
    --out-dir) OUT_DIR="$2"; shift 2 ;;
    --index-dir) INDEX_DIR="$2"; shift 2 ;;
    --gtf) GTF_FILE="$2"; shift 2 ;;
    --ret) RET_BIN="$2"; shift 2 ;;
    --target-groups) TARGET_GROUPS="$2"; shift 2 ;;
    --input-samples) INPUT_SAMPLES="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --force) FORCE="true"; shift ;;
    --execute) RUN="true"; shift ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage; exit 1 ;;
  esac
done

[[ -s "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: --out-dir is required" >&2; exit 1; }
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VISUAL_SCRIPT="$SCRIPT_DIR/render_te_method_visuals.py"
PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python || true)}"
[[ -n "$PYTHON_BIN" ]] || { echo "ERROR: python3/python not found" >&2; exit 1; }
mkdir -p "$OUT_DIR"
STATUS_FILE="$OUT_DIR/method_status.tsv"
if [[ -s "$STATUS_FILE" ]]; then
  # Keep one unambiguous current attempt.  Historical rows are retained as a
  # timestamped legacy file and never participate in the exit-code decision.
  mv "$STATUS_FILE" "${STATUS_FILE%.tsv}.legacy_$(date -u +%Y%m%dT%H%M%SZ)_$$.tsv"
fi
printf 'group\tstatus\tdetail\n' > "$STATUS_FILE"
skip_preflight() {
  printf 'preflight\tSKIP\t%s\n' "$1" >> "$STATUS_FILE"
  echo "[repentools] SKIP: $1" >&2
  if [[ -f "$VISUAL_SCRIPT" ]]; then
    "$PYTHON_BIN" "$VISUAL_SCRIPT" --methods-dir "$(dirname "$OUT_DIR")" \
      --output-dir "$(dirname "$OUT_DIR")/visuals" >/dev/null || true
  fi
  exit 42
}
[[ -d "$INDEX_DIR" ]] || skip_preflight "RepEnTools index directory not found: $INDEX_DIR"
for suffix in 1 2 3 4 5 6 7 8; do
  [[ -s "$INDEX_DIR/chm13-2.${suffix}.ht2" ]] || {
    skip_preflight "missing RepEnTools index: $INDEX_DIR/chm13-2.${suffix}.ht2"
  }
done
[[ -s "$GTF_FILE" ]] || skip_preflight "RepEnTools GTF not found: $GTF_FILE"
[[ -x "$RET_BIN" ]] || skip_preflight "RepEnTools ret not executable: $RET_BIN"
fastq_for() {
  # Avoid Bash-4 associative arrays so the wrapper remains usable in local
  # macOS shells as well as the Linux Pod.  The manifest is small and this
  # lookup runs only once per staged FASTQ.
  "$PYTHON_BIN" - "$MANIFEST" "$1" "$2" <<'PY'
import csv, sys
manifest, wanted, column = sys.argv[1], sys.argv[2], sys.argv[3]
with open(manifest, newline='') as handle:
    for row in csv.DictReader(handle):
        if (row.get('sample') or '').strip() == wanted:
            print((row.get(column) or '').strip())
            break
PY
}

if [[ -z "$INPUT_SAMPLES" || -z "$TARGET_GROUPS" ]]; then
  derived=$("$PYTHON_BIN" - "$MANIFEST" <<'PY'
import csv, sys
from collections import defaultdict
rows = list(csv.DictReader(open(sys.argv[1], newline='')))
truthy = {'1', 'true', 'yes', 'y'}
controls = [str(r.get('sample','')).strip() for r in rows if str(r.get('is_igg','')).strip().lower() in truthy]
groups = defaultdict(list)
for r in rows:
    if str(r.get('is_igg','')).strip().lower() in truthy:
        continue
    sample = str(r.get('sample','')).strip()
    group = str(r.get('group','')).strip() or 'all_targets'
    if sample:
        groups[group].append(sample)
if not controls:
    print('ERROR\tmanifest has no is_igg=true controls')
    raise SystemExit(0)
if len(controls) < 2:
    print('ERROR\tRepEnTools requires two input/control samples; found ' + ','.join(controls))
    raise SystemExit(0)
target_pairs = []
for group, samples in sorted(groups.items()):
    if len(samples) >= 2:
        target_pairs.append(group + '=' + ','.join(samples[:2]))
    else:
        print('WARN\tgroup ' + group + ' has fewer than two target replicates', file=sys.stderr)
if not target_pairs:
    print('ERROR\tno group with two target replicates')
    raise SystemExit(0)
print('OK\t' + ';'.join(target_pairs) + '\t' + ','.join(controls[:2]))
PY
)
  if [[ "$derived" == OK$'\t'* ]]; then
    if [[ -z "$TARGET_GROUPS" ]]; then TARGET_GROUPS="${derived#*$'\t'}"; TARGET_GROUPS="${TARGET_GROUPS%%$'\t'*}"; fi
    if [[ -z "$INPUT_SAMPLES" ]]; then INPUT_SAMPLES="${derived##*$'\t'}"; fi
  else
    detail="${derived#*$'\t'}"
    printf 'manifest_validation\tSKIP\t%s\n' "$detail" >> "$STATUS_FILE"
    echo "[repentools] SKIP: $detail" >&2
    exit 42
  fi
fi

printf '# derived_target_groups=%s\n# derived_input_samples=%s\n' "$TARGET_GROUPS" "$INPUT_SAMPLES" > "$OUT_DIR/manifest_mapping.txt"

stage_sample() {
  local src="$1" dst="$2"
  [[ -s "$src" ]] || { echo "ERROR: FASTQ missing or empty: $src" >&2; return 1; }
  ln -sfn "$src" "$dst"
}

report_ready() {
  local report="$1"
  [[ -s "$report" ]] || return 1
  awk 'NR > 1 && $0 !~ /^[[:space:]]*$/ {found=1; exit} END {exit(found ? 0 : 1)}' "$report"
}

run_group() {
  local group="$1" chip1="$2" chip2="$3" input1="$4" input2="$5"
  local group_dir="$OUT_DIR/$group" stage_dir="$OUT_DIR/$group/input"
  mkdir -p "$stage_dir"
  stage_sample "$(fastq_for "$input1" fastq_1)" "$stage_dir/Rep1_input_R1.fq.gz" || return 1
  stage_sample "$(fastq_for "$input1" fastq_2)" "$stage_dir/Rep1_input_R2.fq.gz" || return 1
  stage_sample "$(fastq_for "$input2" fastq_1)" "$stage_dir/Rep2_input_R1.fq.gz" || return 1
  stage_sample "$(fastq_for "$input2" fastq_2)" "$stage_dir/Rep2_input_R2.fq.gz" || return 1
  stage_sample "$(fastq_for "$chip1" fastq_1)" "$stage_dir/Rep1_ChIP_R1.fq.gz" || return 1
  stage_sample "$(fastq_for "$chip1" fastq_2)" "$stage_dir/Rep1_ChIP_R2.fq.gz" || return 1
  stage_sample "$(fastq_for "$chip2" fastq_1)" "$stage_dir/Rep2_ChIP_R1.fq.gz" || return 1
  stage_sample "$(fastq_for "$chip2" fastq_2)" "$stage_dir/Rep2_ChIP_R2.fq.gz" || return 1

  # RepEnTools emits multi-gigabyte reports.  Treat the experiment summary
  # and report as a resumable completion sentinel so an interrupted or
  # repeated pipeline invocation does not remap/count the same group again.
  if [[ "$RUN" == true && "$FORCE" != true ]] && report_ready "$stage_dir/ret_experiment_summary.csv" && report_ready "$stage_dir/ret_report.csv"; then
    # Recreate stable group-scoped copies as well. Older runs may have
    # completed in the staging directory before the copy step was added.
    find "$stage_dir" -maxdepth 1 -type f \( -name 'ret_*' -o -name '*multiple_feature_counts.txt.summary' \) -exec cp -f {} "$group_dir/" \;
    printf '%s\tPASS\tRepEnTools FASTQ workflow already completed (resume)\n' "$group" >> "$STATUS_FILE"
    echo "[repentools] SKIP $group (completed reports found; use --force to rerun)"
    return 0
  fi

  local command_text
  # The upstream ret has an interactive shebang and sources ~/opt/conda.
  # Invoke it explicitly with bash and set HOME to the Pod user's conda home;
  # otherwise a root-launched pipeline silently exits before running qc().
  index_q="$(printf '%q' "$INDEX_DIR")"
  ret_dir="$(dirname "$RET_BIN")"
  ret_dir_q="$(printf '%q' "$ret_dir")"
  stage_q="$(printf '%q' "$stage_dir")"
  ret_q="$(printf '%q' "$RET_BIN")"
  gtf_q="$(printf '%q' "$GTF_FILE")"
  command_text="export HOME=/path/to; export HISAT2_INDEXES=${index_q}; export REPENTOOLS_DIR=${ret_dir_q}; export ADAPTERS=\"\$REPENTOOLS_DIR\"; cd ${stage_q}; bash ${ret_q} -s . -l '500;500;500;500' -g ${gtf_q} -n ${index_q} -p ${THREADS}"
  printf '%s\n' "$command_text" > "$group_dir/COMMAND.txt"
  if [[ "$RUN" != true ]]; then
    printf '%s\tREADY\tstaged FASTQs; rerun with --execute\n' "$group" >> "$STATUS_FILE"
    return 0
  fi
  echo "[repentools] START $group"
  if bash -lc "$command_text" >"$group_dir/repentools.stdout.log" 2>"$group_dir/repentools.stderr.log"; then
    # ret writes reports in the stage directory; copy a stable, group-scoped
    # manifest so downstream inventory can discover the outputs recursively.
    find "$stage_dir" -maxdepth 1 -type f \( -name 'ret_*' -o -name '*multiple_feature_counts.txt.summary' \) -exec cp -f {} "$group_dir/" \;
    # A zero exit code from an upstream shell workflow is not sufficient:
    # ret can finish before qc/counting and leave an empty directory.  Require
    # the two stable report sentinels before exposing PASS to downstream.
    missing=()
    report_ready "$group_dir/ret_experiment_summary.csv" || missing+=(ret_experiment_summary.csv)
    report_ready "$group_dir/ret_report.csv" || missing+=(ret_report.csv)
    if [[ ${#missing[@]} -gt 0 ]]; then
      printf '%s\tFAIL\tcommand returned 0 but required output is missing: %s\n' "$group" "${missing[*]}" >> "$STATUS_FILE"
      echo "[repentools] FAIL $group: required report missing (${missing[*]})" >&2
      overall_status=1
      return 0
    fi
    printf '%s\tPASS\tRepEnTools FASTQ workflow completed; validated ret_experiment_summary.csv and ret_report.csv\n' "$group" >> "$STATUS_FILE"
    echo "[repentools] DONE $group"
  else
    printf '%s\tFAIL\tsee %s/repentools.stderr.log\n' "$group" "$group_dir" >> "$STATUS_FILE"
    echo "[repentools] FAIL $group" >&2
    return 1
  fi
}

overall_status=0
IFS=';' read -r -a pairs <<< "$TARGET_GROUPS"
for pair in "${pairs[@]}"; do
  [[ -n "$pair" ]] || continue
  group="${pair%%=*}"
  samples="${pair#*=}"
  IFS=',' read -r chip1 chip2 <<< "$samples"
  IFS=',' read -r input1 input2 <<< "$INPUT_SAMPLES"
  chip1_fq1="$(fastq_for "$chip1" fastq_1)"
  chip2_fq1="$(fastq_for "$chip2" fastq_1)"
  input1_fq1="$(fastq_for "$input1" fastq_1)"
  input2_fq1="$(fastq_for "$input2" fastq_1)"
  [[ -n "$chip1_fq1" && -n "$chip2_fq1" ]] || {
    printf '%s\tSKIP\tmanifest lacks target pair %s,%s\n' "$group" "$chip1" "$chip2" >> "$STATUS_FILE"; continue;
  }
  [[ -n "$input1_fq1" && -n "$input2_fq1" ]] || {
    printf '%s\tSKIP\tmanifest lacks input pair %s,%s\n' "$group" "$input1" "$input2" >> "$STATUS_FILE"; continue;
  }
  if ! run_group "$group" "$chip1" "$chip2" "$input1" "$input2"; then
    overall_status=1
  fi
done

# Keep the RepEnTools-only entry point consistent with the full downstream
# pipeline: it also creates a validated visualization status and, when report
# columns are parseable, a group-level SVG summary.
if [[ -f "$SCRIPT_DIR/render_te_method_visuals.py" ]]; then
  "$PYTHON_BIN" "$SCRIPT_DIR/render_te_method_visuals.py" \
    --methods-dir "$(dirname "$OUT_DIR")" \
    --output-dir "$(dirname "$OUT_DIR")/visuals" >/dev/null || true
fi

echo "[repentools] status: $STATUS_FILE"
if [[ "$overall_status" -ne 0 ]]; then
  exit "$overall_status"
fi
# 42 is the pipeline-wide optional-module SKIP code.  This prevents a
# manifest with no compatible CHM13 group (or a dry-run with only READY rows)
# from being recorded as PASS by run_downstream.sh.
if awk -F '\t' 'NR > 1 && $2 == "PASS" {found=1} END {exit found ? 0 : 1}' "$STATUS_FILE"; then
  exit 0
fi
if awk -F '\t' 'NR > 1 && $2 == "FAIL" {found=1} END {exit found ? 0 : 1}' "$STATUS_FILE"; then
  exit 1
fi
exit 42
