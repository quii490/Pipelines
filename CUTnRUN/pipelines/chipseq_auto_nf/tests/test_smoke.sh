#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

python3 "$ROOT/tools/cutrun_cli/cutrun_status.py" init \
  --status-file "$TMP/status.tsv" --run-id smoke
python3 "$ROOT/tools/cutrun_cli/cutrun_status.py" record \
  --status-file "$TMP/status.tsv" --run-id smoke --module optional:test \
  --status PASS --outputs-ok
python3 "$ROOT/tools/cutrun_cli/cutrun_status.py" finalize \
  --status-file "$TMP/status.tsv" --run-id smoke --output "$TMP/summary.json"
grep -q $'optional:test\t1\tPASS' "$TMP/status_latest.tsv"

mkdir -p "$TMP/results/manifest"
cp "$ROOT/pipelines/chipseq_auto_nf/tests/fixtures/manifest.control_map.csv" \
  "$TMP/results/manifest/resolved_manifest.csv"
bash "$ROOT/pipelines/chipseq_auto_nf/run_downstream.sh" \
  --results-dir "$TMP/results" --species hg38 --skip-bw-cor --skip-count-draw
test -s "$TMP/results/09_downstream/run_manifest.json"
test -s "$TMP/results/09_downstream/module_status_latest.tsv"
test -s "$TMP/results/09_downstream/run_report.html"

# Optional TE tools have heterogeneous upstream formats.  Verify that the
# adapter creates real method-specific SVGs from a tiny fixture and records
# their status for the report.
cp -R "$ROOT/pipelines/chipseq_auto_nf/tests/fixtures/te_methods" "$TMP/te_methods"
python3 "$ROOT/tools/te_methods/render_te_method_visuals.py" \
  --methods-dir "$TMP/te_methods" --output-dir "$TMP/te_methods/visuals"
test -s "$TMP/te_methods/visualization_status.tsv"
test -s "$TMP/te_methods/visuals/allo/allo_mapped_fraction.svg"
test -s "$TMP/te_methods/visuals/t3e/t3e_top_features.svg"
test -s "$TMP/te_methods/visuals/repentools_fastq/G1/repentools/repentools_top_features.svg"
printf '%s\n' 'smoke test PASS'
