#!/usr/bin/env bash
set -euo pipefail

SNAKEMAKE_ROOT="/path/to/softwares/snakemake-9.22.0"
SNAKEMAKE_VENV="/path/to/softwares/snakemake-9.22.0-venv"
MODE="user-bin"
BIN_DIR="${HOME}/.local/bin"
COMMAND_NAME="snakemake"
PYTHON_BIN="python3"

usage() {
  cat <<'USAGE'
Install or expose /path/to/softwares/snakemake-9.22.0

Recommended:
  bash install_snakemake_922.sh --mode user-bin

Modes:
  --mode create-venv        Create a dedicated venv and pip-install the local source tree.
  --mode user-bin          Create ~/.local/bin/snakemake wrapper. Safe and reversible.
  --mode current-conda-env Create $CONDA_PREFIX/bin/snakemake wrapper for current env.
  --mode all-conda-envs    Create wrappers in every writable conda env from `conda env list`.

Options:
  --snakemake-root DIR     Default: /path/to/softwares/snakemake-9.22.0
  --snakemake-venv DIR     Default: /path/to/softwares/snakemake-9.22.0-venv
  --bin-dir DIR            Used by user-bin mode. Default: ~/.local/bin
  --command-name NAME      Default: snakemake
  --python PATH            Python used for create-venv mode. Default: python3
  -h, --help               Show help

Notes:
  - /path/to/softwares/snakemake-9.22.0 is usually a source tree.
    Run --mode create-venv once before installing wrappers if no executable
    snakemake exists under that directory.
  - user-bin mode works in ordinary shells if ~/.local/bin is before older snakemake commands.
  - if an activated conda env has its own snakemake earlier in PATH, use current-conda-env
    or all-conda-envs mode to make version 9.22.0 win inside those envs.
USAGE
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --snakemake-root) SNAKEMAKE_ROOT="$2"; shift 2 ;;
    --snakemake-venv) SNAKEMAKE_VENV="$2"; shift 2 ;;
    --mode) MODE="$2"; shift 2 ;;
    --bin-dir) BIN_DIR="$2"; shift 2 ;;
    --command-name) COMMAND_NAME="$2"; shift 2 ;;
    --python) PYTHON_BIN="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "[install_snakemake_922] Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

find_snakemake_bin() {
  local root="$1"
  local cand
  for cand in \
    "${SNAKEMAKE_VENV}/bin/snakemake" \
    "${root}/bin/snakemake" \
    "${root}/snakemake" \
    "${root}/.venv/bin/snakemake" \
    "${root}/venv/bin/snakemake"; do
    if [[ -x "${cand}" ]]; then
      echo "${cand}"
      return 0
    fi
  done
  cand="$(find "${root}" -maxdepth 4 -type f -name snakemake -perm -111 2>/dev/null | head -1 || true)"
  if [[ -n "${cand}" ]]; then
    echo "${cand}"
    return 0
  fi
  return 1
}

if [[ ! -d "${SNAKEMAKE_ROOT}" ]]; then
  echo "[install_snakemake_922] ERROR: snakemake root not found: ${SNAKEMAKE_ROOT}" >&2
  exit 1
fi

write_wrapper() {
  local dest="$1"
  mkdir -p "$(dirname "${dest}")"
  cat > "${dest}" <<EOF
#!/usr/bin/env bash
exec "${SNAKEMAKE_BIN}" "\$@"
EOF
  chmod +x "${dest}"
  echo "[install_snakemake_922] installed ${dest} -> ${SNAKEMAKE_BIN}"
}

ensure_snakemake_bin() {
  SNAKEMAKE_BIN="$(find_snakemake_bin "${SNAKEMAKE_ROOT}")" || {
    echo "[install_snakemake_922] ERROR: cannot find executable snakemake." >&2
    echo "If ${SNAKEMAKE_ROOT} is a source tree, run:" >&2
    echo "  bash install_snakemake_922.sh --mode create-venv" >&2
    exit 1
  }
}

case "${MODE}" in
  create-venv)
    command -v "${PYTHON_BIN}" >/dev/null 2>&1 || {
      echo "[install_snakemake_922] ERROR: python not found: ${PYTHON_BIN}" >&2
      exit 1
    }
    "${PYTHON_BIN}" -m venv "${SNAKEMAKE_VENV}"
    "${SNAKEMAKE_VENV}/bin/python" -m pip install --upgrade pip setuptools wheel
    "${SNAKEMAKE_VENV}/bin/python" -m pip install "${SNAKEMAKE_ROOT}"
    "${SNAKEMAKE_VENV}/bin/snakemake" --version
    echo "[install_snakemake_922] venv ready: ${SNAKEMAKE_VENV}"
    ;;
  user-bin)
    ensure_snakemake_bin
    write_wrapper "${BIN_DIR}/${COMMAND_NAME}"
    cat <<EOF

Add this to ~/.bashrc or ~/.zshrc if it is not already there:
  export PATH="${BIN_DIR}:\$PATH"

Then reopen shell or run:
  source ~/.bashrc
  ${COMMAND_NAME} --version
EOF
    ;;
  current-conda-env)
    ensure_snakemake_bin
    if [[ -z "${CONDA_PREFIX:-}" ]]; then
      echo "[install_snakemake_922] ERROR: CONDA_PREFIX is empty; activate a conda env first" >&2
      exit 1
    fi
    write_wrapper "${CONDA_PREFIX}/bin/${COMMAND_NAME}"
    "${CONDA_PREFIX}/bin/${COMMAND_NAME}" --version || true
    ;;
  all-conda-envs)
    ensure_snakemake_bin
    command -v conda >/dev/null 2>&1 || {
      echo "[install_snakemake_922] ERROR: conda not found in PATH" >&2
      exit 1
    }
    while IFS= read -r env_path; do
      [[ -n "${env_path}" && -d "${env_path}/bin" && -w "${env_path}/bin" ]] || continue
      write_wrapper "${env_path}/bin/${COMMAND_NAME}"
    done < <(conda env list | awk 'NF && $1 !~ /^#/ {print $NF}' | grep '^/')
    ;;
  *)
    echo "[install_snakemake_922] ERROR: unknown mode: ${MODE}" >&2
    usage
    exit 1
    ;;
esac
