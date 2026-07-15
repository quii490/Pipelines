#!/usr/bin/env bash
set -euo pipefail

RESULTS_DIR=""
TRACE_FILE=""
ENGINE_LOG=""
PIPELINE_PID=""
OUTPUT_FILE=""
INTERVAL=30
WATCH=false

usage() {
  cat <<'USAGE'
Usage:
  bash rnaseq_progress.sh --results-dir DIR [options]

Options:
  --trace FILE         Trace file for this run. Defaults to the newest trace.
  --engine-log FILE    Nextflow engine log for this run. Defaults to the newest log.
  --pid PID            Pipeline wrapper PID; used to identify a live run.
  --output FILE        Write a progress snapshot to FILE.
  --watch              Refresh until --pid exits (default interval: 30 seconds).
  --interval SEC       Refresh interval used with --watch.
  -h, --help           Show this help.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --trace) TRACE_FILE="$2"; shift 2 ;;
    --engine-log) ENGINE_LOG="$2"; shift 2 ;;
    --pid) PIPELINE_PID="$2"; shift 2 ;;
    --output) OUTPUT_FILE="$2"; shift 2 ;;
    --watch) WATCH=true; shift ;;
    --interval) INTERVAL="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$RESULTS_DIR" ]] || { echo "--results-dir is required" >&2; exit 1; }
[[ "$INTERVAL" =~ ^[0-9]+$ ]] && (( INTERVAL > 0 )) || { echo "--interval must be a positive integer" >&2; exit 1; }

LOG_DIR="${RESULTS_DIR}/_automation/logs"

newest_file() {
  local pattern="$1"
  compgen -G "$pattern" | sort -r | head -n 1 || true
}

resolve_defaults() {
  [[ -n "$TRACE_FILE" ]] || TRACE_FILE="$(newest_file "${LOG_DIR}/nextflow_*.trace.tsv")"
  [[ -n "$ENGINE_LOG" ]] || ENGINE_LOG="$(newest_file "${LOG_DIR}/nextflow_*.log")"
}

run_is_alive() {
  [[ -n "$PIPELINE_PID" ]] && kill -0 "$PIPELINE_PID" 2>/dev/null
}

workflow_state() {
  local failed_count=0
  if [[ -n "$TRACE_FILE" && -s "$TRACE_FILE" ]]; then
    failed_count="$(awk -F '\t' '
      NR == 1 { for (i = 1; i <= NF; i++) if ($i == "status") s=i; next }
      s && $s != "" && $s !~ /^(COMPLETED|CACHED)$/ { n++ }
      END { print n + 0 }
    ' "$TRACE_FILE")"
  fi
  if run_is_alive; then
    printf 'RUNNING'
  elif [[ "$failed_count" -gt 0 ]]; then
    printf 'PARTIAL_SUCCESS'
  elif [[ -n "$ENGINE_LOG" && -f "$ENGINE_LOG" ]] && grep -q 'Session aborted' "$ENGINE_LOG"; then
    printf 'FAILED'
  elif [[ -n "$ENGINE_LOG" && -f "$ENGINE_LOG" ]] && grep -q 'Execution complete -- Goodbye' "$ENGINE_LOG"; then
    printf 'COMPLETED'
  else
    printf 'UNKNOWN'
  fi
}

snapshot() {
  local state now task_summary stage_summary latest
  now="$(date '+%Y-%m-%d %H:%M:%S')"
  state="$(workflow_state)"

  task_summary='no completed tasks recorded yet'
  stage_summary=''
  if [[ -n "$TRACE_FILE" && -s "$TRACE_FILE" ]]; then
    task_summary="$(awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) { if ($i == "status") s=i; if ($i == "process" || $i == "name") p=i }
        next
      }
      s && $s != "" { count[$s]++ }
      END { for (x in count) printf "%s=%d ", x, count[x] }
    ' "$TRACE_FILE")"
    stage_summary="$(awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) { if ($i == "status") s=i; if ($i == "process" || $i == "name") p=i }
        next
      }
      s && p && $p != "" { proc=$p; sub(/ \(.*/, "", proc); key=proc SUBSEP $s; count[key]++ }
      END {
        for (key in count) {
          split(key, parts, SUBSEP)
          proc=parts[1]; status=parts[2]
          item[proc]=item[proc] sprintf(" %s=%d", status, count[key])
        }
        for (proc in item) print proc ":" item[proc]
      }
    ' "$TRACE_FILE" | sort)"
  fi

  latest=''
  if [[ -n "$ENGINE_LOG" && -f "$ENGINE_LOG" ]]; then
    latest="$(grep -E "ERROR|Error executing process|Process .*failed|Submitted process|Completed process|Session aborted" "$ENGINE_LOG" | tail -n 1 || true)"
  fi

  cat <<EOF
RNA-seq pipeline progress
Updated: ${now}
State: ${state}
Results: ${RESULTS_DIR}
Pipeline PID: ${PIPELINE_PID:-not supplied}
Trace: ${TRACE_FILE:-not created yet}
Engine log: ${ENGINE_LOG:-not created yet}

Task status: ${task_summary}
Non-completed trace records: $(awk -F '\t' '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "status") s=i; next }
  s && $s != "" && $s !~ /^(COMPLETED|CACHED)$/ { n++ }
  END { print n + 0 }
' "${TRACE_FILE:-/dev/null}" 2>/dev/null)
Per-process status:
${stage_summary:-  no trace records yet}

Latest notable event:
${latest:-  no error or task event recorded yet}
EOF
}

write_snapshot() {
  local text tmp latest
  text="$(snapshot)"
  if [[ -n "$OUTPUT_FILE" ]]; then
    mkdir -p "$(dirname "$OUTPUT_FILE")"
    tmp="${OUTPUT_FILE}.tmp.$$"
    printf '%s\n' "$text" > "$tmp"
    mv -f "$tmp" "$OUTPUT_FILE"
    latest="${LOG_DIR}/progress_latest.txt"
    if [[ "$OUTPUT_FILE" != "$latest" ]]; then
      cp -f "$OUTPUT_FILE" "$latest"
    fi
  else
    printf '%s\n' "$text"
  fi
}

resolve_defaults
if [[ "$WATCH" == true ]]; then
  while true; do
    write_snapshot
    run_is_alive || break
    sleep "$INTERVAL"
  done
  write_snapshot
else
  write_snapshot
fi
