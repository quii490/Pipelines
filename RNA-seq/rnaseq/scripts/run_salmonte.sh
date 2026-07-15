#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Run SalmonTE with a patched salmon binary from a conda env.

Required:
  --salmonte-dir DIR
  --reference REF
  --input-dir DIR
  --conda-prefix DIR

Optional:
  --exprtype STR
  --outdir DIR
EOF
}

SALMONTE_DIR=""
REFERENCE=""
INPUT_DIR=""
CONDA_PREFIX_DIR=""
EXPRTYPE="count"
OUTDIR="SalmonTE_output"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --salmonte-dir) SALMONTE_DIR="$2"; shift 2 ;;
    --reference) REFERENCE="$2"; shift 2 ;;
    --input-dir) INPUT_DIR="$2"; shift 2 ;;
    --conda-prefix) CONDA_PREFIX_DIR="$2"; shift 2 ;;
    --exprtype) EXPRTYPE="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown arg: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "$SALMONTE_DIR" ]] || { echo "Missing --salmonte-dir" >&2; exit 1; }
[[ -n "$REFERENCE" ]] || { echo "Missing --reference" >&2; exit 1; }
[[ -n "$INPUT_DIR" ]] || { echo "Missing --input-dir" >&2; exit 1; }
[[ -n "$CONDA_PREFIX_DIR" ]] || { echo "Missing --conda-prefix" >&2; exit 1; }

SALMON_BIN="$CONDA_PREFIX_DIR/bin/salmon"
TARGET_BIN="$SALMONTE_DIR/salmon/linux/bin/salmon"
TARGET_BAK="$SALMONTE_DIR/salmon/linux/bin/salmon.bak_chatgpt"

[[ -x "$SALMON_BIN" ]] || { echo "salmon not found in env: $SALMON_BIN" >&2; exit 1; }
[[ -e "$TARGET_BIN" ]] || { echo "Bundled salmon not found: $TARGET_BIN" >&2; exit 1; }

restore() {
  if [[ -L "$TARGET_BIN" && -e "$TARGET_BAK" ]]; then
    rm -f "$TARGET_BIN"
    mv "$TARGET_BAK" "$TARGET_BIN"
  fi
}
trap restore EXIT

if [[ ! -e "$TARGET_BAK" ]]; then
  mv "$TARGET_BIN" "$TARGET_BAK"
fi
ln -sf "$SALMON_BIN" "$TARGET_BIN"

export PATH="$CONDA_PREFIX_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX_DIR/lib:${LD_LIBRARY_PATH:-}"

python "$SALMONTE_DIR/SalmonTE.py" quant \
  --reference="$REFERENCE" \
  --exprtype="$EXPRTYPE" \
  "$INPUT_DIR"

if [[ "$OUTDIR" != "SalmonTE_output" && -d "SalmonTE_output" ]]; then
  rm -rf "$OUTDIR"
  mv SalmonTE_output "$OUTDIR"
fi
