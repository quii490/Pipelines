#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

fail=0
private_pattern='(/home/[A-Za-z0-9._-]+|/Users/[A-Za-z0-9._-]+|/miniconda3|192\.168\.|10\.[0-9]+\.[0-9]+\.[0-9]+|BEGIN (RSA |OPENSSH )?PRIVATE KEY|github_pat_|ghp_)'

while IFS= read -r -d '' file; do
  # `git ls-files` also reports tracked files deleted in the current change.
  [[ -f "$file" ]] || continue
  [[ "$file" == "scripts/validate_public_repo.sh" ]] && continue

  if grep -I -n -E "$private_pattern" "$file"; then
    echo "ERROR: sensitive content detected in: $file" >&2
    fail=1
  fi

  size=$(wc -c < "$file")
  if (( size > 10000000 )); then
    echo "ERROR: file larger than 10 MB: $file ($size bytes)" >&2
    fail=1
  fi

  case "$file" in
    */tests/fixtures/*) ;;
    *.fastq|*.fastq.gz|*.fq|*.fq.gz|*.bam|*.cram|*.bw|*.bigWig)
      echo "ERROR: sequencing data or generated artifact: $file" >&2
      fail=1
      ;;
  esac
done < <(git ls-files --cached --others --exclude-standard -z)

if (( fail != 0 )); then
  exit 1
fi

echo "Public repository validation passed."
