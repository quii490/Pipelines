#!/usr/bin/env bash
set -euo pipefail

resolve_script_dir() {
  local src="${BASH_SOURCE[0]}"
  while [ -L "$src" ]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    [[ "$src" != /* ]] && src="$dir/$src"
  done
  cd -P "$(dirname "$src")" && pwd
}

get_project_root() {
  local script_dir
  script_dir="$(resolve_script_dir)"
  cd "$script_dir/.." && pwd
}

get_species_param() {
  local species="$1"
  local key="$2"
  local config_file="$3"
  local value
  value="$(awk -v sp="$species" -v key="$key" '
    $0 ~ "if \\(params\\.species == " && $0 ~ sp {in_block=1; next}
    in_block && $0 ~ /^}/ {exit}
    in_block && $0 ~ "params\\." key "[[:space:]]*=" {
      if (match($0, /\x27[^\x27]*\x27/)) {
        print substr($0, RSTART + 1, RLENGTH - 2)
        exit
      }
      if (match($0, /"[^"]*"/)) {
        print substr($0, RSTART + 1, RLENGTH - 2)
        exit
      }
    }
  ' "$config_file")"
  local base_dir_token='${baseDir}'
  if [[ "$value" == "${base_dir_token}"/* ]]; then
    local project_root
    project_root="$(cd "$(dirname "$config_file")/.." && pwd)"
    value="${value/"$base_dir_token"/$project_root}"
  fi
  printf '%s\n' "$value"
}
