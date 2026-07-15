#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

fail=0
scan_args=(--hidden -g '!**/.git/**' -g '!.venv/**' -g '!scripts/validate_public_repo.sh' -g '!docs/site/**')

if rg -n "(/home/[A-Za-z0-9._-]+|/Users/[A-Za-z0-9._-]+|/miniconda3|192\\.168\\.|10\\.[0-9]+\\.[0-9]+\\.[0-9]+|BEGIN (RSA |OPENSSH )?PRIVATE KEY|github_pat_|ghp_)" "${scan_args[@]}" .; then
  echo "ERROR: detected a private path, private-network address, or secret-shaped value." >&2
  fail=1
fi

while IFS= read -r -d '' file; do
  size=$(wc -c < "$file")
  if (( size > 10000000 )); then
    echo "ERROR: file larger than 10 MB: $file ($size bytes)" >&2
    fail=1
  fi
done < <(find . -type f -not -path './.git/*' -not -path './.venv/*' -not -path './docs/site/*' -print0)

if find . -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq.gz' -o -name '*.bam' -o -name '*.cram' -o -name '*.bw' -o -name '*.bigWig' \) -not -path '*/tests/fixtures/*' -print | grep -q .; then
  echo "ERROR: detected sequencing data or generated analysis artifacts." >&2
  fail=1
fi

if (( fail != 0 )); then
  exit 1
fi

echo "Public repository validation passed."
