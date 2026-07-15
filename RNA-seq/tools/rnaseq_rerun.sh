#!/usr/bin/env bash
set -euo pipefail
SCRIPT_PATH="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(cd "$(dirname "${SCRIPT_PATH}")" && pwd)"
RUNNER="${SCRIPT_DIR}/../rnaseq/run_auto_rnaseq.sh"
RESULTS_DIR=""
MODULES=""
EXTRA=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --results-dir) RESULTS_DIR="$2"; shift 2 ;;
    --modules) MODULES="$2"; shift 2 ;;
    *) EXTRA+=("$1"); shift ;;
  esac
done
[[ -n "$RESULTS_DIR" && -n "$MODULES" ]] || { echo "Usage: rnaseq-rerun --results-dir DIR --modules MODULE1,MODULE2 [downstream options]" >&2; exit 2; }
exec bash "$RUNNER" --results-dir "$RESULTS_DIR" --downstream --only-tools "$MODULES" "${EXTRA[@]}"
