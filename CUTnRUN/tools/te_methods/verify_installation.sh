#!/usr/bin/env bash
set -euo pipefail

OUT="${1:-tool_installation.tsv}"
mkdir -p "$(dirname "$OUT")"
printf 'tool\tstatus\tversion\tpath\tnotes\n' > "$OUT"

record() {
  printf '%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" >> "$OUT"
}

t3e_py="/path/to/.conda/envs/cutrun-t3e/bin/python"
if [[ -x "$t3e_py" ]]; then
  version="$($t3e_py --version 2>&1)"
  commit="unknown"
  if [[ -d /path/to/.local/share/CUTnRUN-tools/T3E/.git ]]; then
    commit="$(git -C /path/to/.local/share/CUTnRUN-tools/T3E rev-parse --short HEAD)"
  fi
  record t3e installed "$version" "$t3e_py" "repo_commit=$commit"
else
  record t3e missing unknown "$t3e_py" "install environment-t3e.yml"
fi

allo="/path/to/.conda/envs/cutrun-allo/bin/allo"
if [[ -x "$allo" ]]; then
  version="$(sed -n 's/^version = //p' /path/to/.conda/envs/cutrun-allo/lib/python*/site-packages/Allo/allo | head -1 | tr -d '"')"
  tensorflow_version="$(/path/to/.conda/envs/cutrun-allo/bin/python -c 'from importlib.metadata import version; print(version("tensorflow"))')"
  record allo installed "${version:-unknown}" "$allo" "tensorflow=${tensorflow_version}"
else
  record allo missing unknown "$allo" "install environment-allo.yml"
fi

ret="/path/to/.local/share/CUTnRUN-tools/RepEnTools/ret"
if [[ -x "$ret" ]]; then
  commit="unknown"
  if [[ -d /path/to/.local/share/CUTnRUN-tools/RepEnTools/.git ]]; then
    commit="$(git -C /path/to/.local/share/CUTnRUN-tools/RepEnTools rev-parse --short HEAD)"
  fi
  record repentools installed repository "$ret" "repo_commit=$commit;reference=chm13v2_required"
else
  record repentools missing unknown "$ret" "run RepEnTools/installation"
fi

echo "[te_methods] installation status: $OUT"
