#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
自动扫描 FASTQ 并生成可编辑的 ATAC-seq 设计文件。

用法：
  bash generate_samplesheet.sh \
    --input /path/to/fastq_dir \
    --samplesheet samplesheet.csv \
    --contrasts contrasts.csv

参数：
  --input         FASTQ 所在目录（递归扫描）
  --samplesheet   输出 samplesheet.csv，默认 ./samplesheet.csv
  --contrasts     输出 contrasts.csv，默认 ./contrasts.csv
  --layout        auto / PE / SE，默认 auto
  --condition     默认 condition，默认 NA
  --replicate     默认 replicate，默认 NA
  --metadata-csv  可选；若提供 sample/condition/replicate，可自动回填

说明：
1. 自动识别常见命名：_R1/_R2、_r1/_r2、_1/_2、.R1/.R2、.r1/.r2、-R1/-R2 等。
2. 自动识别扩展名：.fastq.gz / .fq.gz / .fastq / .fq。
3. 会输出两个文件：
   - samplesheet.csv：你只需要手动补 condition 和 replicate
   - contrasts.csv：你手动写 case / control 比较
4. 如果 metadata-csv 已提供 condition / replicate，则会自动回填。
USAGE
}

log_msg() {
  echo "[generate_samplesheet] $*"
}

INPUT=""
SAMPLESHEET="./samplesheet.csv"
CONTRASTS="./contrasts.csv"
LAYOUT="auto"
DEFAULT_CONDITION="NA"
DEFAULT_REPLICATE="NA"
METADATA_CSV=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --samplesheet|--output) SAMPLESHEET="$2"; shift 2 ;;
    --contrasts) CONTRASTS="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --condition) DEFAULT_CONDITION="$2"; shift 2 ;;
    --replicate) DEFAULT_REPLICATE="$2"; shift 2 ;;
    --metadata-csv) METADATA_CSV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "未知参数: $1" >&2; usage; exit 1 ;;
  esac
done

[[ -n "${INPUT}" ]] || { echo "错误：必须提供 --input" >&2; usage; exit 1; }
[[ -d "${INPUT}" ]] || { echo "错误：输入目录不存在: ${INPUT}" >&2; exit 1; }
if [[ -n "${METADATA_CSV}" && ! -f "${METADATA_CSV}" ]]; then
  echo "错误：metadata 文件不存在: ${METADATA_CSV}" >&2
  exit 1
fi

mkdir -p "$(dirname "${SAMPLESHEET}")" "$(dirname "${CONTRASTS}")"

python3 - <<'PY' "${INPUT}" "${SAMPLESHEET}" "${CONTRASTS}" "${LAYOUT}" "${DEFAULT_CONDITION}" "${DEFAULT_REPLICATE}" "${METADATA_CSV}"
import csv
import os
import re
import sys
from pathlib import Path

input_dir, samplesheet, contrasts, layout, default_condition, default_replicate, metadata_csv = sys.argv[1:8]
layout = layout.strip().upper()
if layout not in {"AUTO", "PE", "SE"}:
    raise SystemExit("--layout 只能是 auto / PE / SE")

exts = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
read_patterns = [
    re.compile(r"^(?P<prefix>.+?)[._-]R(?P<read>[12])(?:[._-]L\d{3})?$", re.IGNORECASE),
    re.compile(r"^(?P<prefix>.+?)(?:[._-]L\d{3})[._-]R(?P<read>[12])$", re.IGNORECASE),
    re.compile(r"^(?P<prefix>.+?)[._-](?P<read>[12])(?:[._-]L\d{3})?$", re.IGNORECASE),
]

meta_map = {}
if metadata_csv:
    with open(metadata_csv, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if "sample" not in cols:
            raise SystemExit("metadata-csv 缺少 sample 列")
        sample_col = cols["sample"]
        cond_col = cols.get("condition")
        rep_col = cols.get("replicate")
        for row in reader:
            sample = str(row.get(sample_col, "")).strip()
            if not sample:
                continue
            meta_map[sample] = {
                "condition": str(row.get(cond_col, "")).strip() if cond_col else "",
                "replicate": str(row.get(rep_col, "")).strip() if rep_col else "",
            }


def strip_ext(name: str) -> str:
    for ext in exts:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def normalize_sample_name(prefix: str) -> str:
    prefix = re.sub(r"(?:[._-]S\d+)?(?:[._-]L\d{3})$", "", prefix, flags=re.IGNORECASE)
    prefix = re.sub(r"[._-]+$", "", prefix)
    return prefix


def parse_read(path: Path):
    stem = strip_ext(path.name)
    stem = re.sub(r"[._-]001$", "", stem, flags=re.IGNORECASE)
    m_simple = re.match(r"^(?:R|read)?(?P<read>[12])$", stem, flags=re.IGNORECASE)
    if m_simple:
        sample_from_dir = normalize_sample_name(path.parent.name)
        return sample_from_dir, m_simple.group("read")
    for pattern in read_patterns:
        m = pattern.match(stem)
        if m:
            prefix = normalize_sample_name(m.group("prefix"))
            read = m.group("read")
            return prefix, read
    return stem, None

files = sorted(
    p for p in Path(input_dir).rglob("*")
    if p.is_file() and any(str(p).endswith(ext) for ext in exts)
)
if not files:
    raise SystemExit(f"未找到 FASTQ 文件: {input_dir}")

sample_map = {}
used = set()
for path in files:
    sample, read = parse_read(path)
    rec = sample_map.setdefault(sample, {"R1": [], "R2": [], "SE": []})
    if read == "1":
        rec["R1"].append(str(path.resolve()))
    elif read == "2":
        rec["R2"].append(str(path.resolve()))
    else:
        rec["SE"].append(str(path.resolve()))

rows = []
for sample in sorted(sample_map):
    rec = sample_map[sample]
    r1s = sorted(rec["R1"])
    r2s = sorted(rec["R2"])
    ses = sorted(rec["SE"])

    if layout == "PE":
        n = min(len(r1s), len(r2s))
        if n == 0:
            continue
        if len(r1s) != len(r2s):
            print(f"[generate_samplesheet] 警告: {sample} 的 R1/R2 数量不一致，按最小配对数 {n} 处理", file=sys.stderr)
        rows.append([sample, "PE", default_condition, default_replicate, ",".join(r1s[:n]), ",".join(r2s[:n])])
        used.update(r1s[:n]); used.update(r2s[:n])
        continue

    if layout == "SE":
        se_files = ses if ses else r1s
        if not se_files:
            continue
        rows.append([sample, "SE", default_condition, default_replicate, ",".join(se_files), ""])
        used.update(se_files)
        continue

    if r1s and r2s:
        n = min(len(r1s), len(r2s))
        if len(r1s) != len(r2s):
            print(f"[generate_samplesheet] 警告: {sample} 的 R1/R2 数量不一致，按最小配对数 {n} 处理", file=sys.stderr)
        rows.append([sample, "PE", default_condition, default_replicate, ",".join(r1s[:n]), ",".join(r2s[:n])])
        used.update(r1s[:n]); used.update(r2s[:n])
    elif ses:
        rows.append([sample, "SE", default_condition, default_replicate, ",".join(ses), ""])
        used.update(ses)
    elif r1s:
        rows.append([sample, "SE", default_condition, default_replicate, ",".join(r1s), ""])
        used.update(r1s)

if not rows:
    raise SystemExit("未能根据当前 layout 检测到有效样本")

rows.sort(key=lambda x: x[0])

with open(samplesheet, "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["sample", "layout", "condition", "replicate", "r1", "r2"])
    for row in rows:
        sample = row[0]
        meta = meta_map.get(sample, {})
        row[2] = meta.get("condition") or row[2]
        row[3] = meta.get("replicate") or row[3]
        writer.writerow(row)

conditions = sorted({r[2] for r in rows if r[2] and r[2] != "NA"})
with open(contrasts, "w", newline="", encoding="utf-8") as fh:
    writer = csv.writer(fh)
    writer.writerow(["case", "control"])
    if len(conditions) >= 2:
        for case in conditions[1:]:
            writer.writerow([case, conditions[0]])

print(f"[generate_samplesheet] 样本数: {len(rows)}", file=sys.stderr)
print(f"[generate_samplesheet] 已写出 samplesheet: {samplesheet}", file=sys.stderr)
print(f"[generate_samplesheet] 已写出 contrasts: {contrasts}", file=sys.stderr)
PY

log_msg "已生成: ${SAMPLESHEET}"
log_msg "已生成: ${CONTRASTS}"
log_msg "samplesheet 预览："
head -n 8 "${SAMPLESHEET}" || true
log_msg "contrasts 预览："
head -n 8 "${CONTRASTS}" || true
