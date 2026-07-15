#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
自动扫描 FASTQ 并生成 samplesheet.csv

用法：
  bash generate_samplesheet.sh \
    --input /path/to/project_or_fastq_dir \
    --output samplesheet.csv \
    --layout auto \
    --r1-pattern "_R1" \
    --r2-pattern "_R2" \
    --condition NA \
    --metadata-csv sample_metadata.csv

参数：
  --input         FASTQ 所在目录（递归扫描；支持所有 FASTQ 在同一目录，或项目目录下多个子目录）
  --output        输出 csv 路径，默认 ./samplesheet.csv
  --layout        PE / SE / auto，默认 auto
  --r1-pattern    可选人工指定 R1 标记，默认 _R1
  --r2-pattern    可选人工指定 R2 标记，默认 _R2
  --condition     condition 默认值，默认 NA
  --replicate     replicate 默认值，默认 NA
  --metadata-csv  可选，样本元信息表；至少包含 sample，可选 condition/replicate

说明：
1. 递归支持常见后缀：.fastq.gz / .fq.gz / .fastq / .fq
2. 自动识别常见双端命名：_1/_2、_R1/_R2、_r1/_r2、.1/.2、-R1/-R2 等
3. 未成功配对且不是明显 R2 的文件，在 auto / SE 模式下按单端处理
4. 若提供 --metadata-csv，则按 sample 列回填 condition / replicate；未命中时使用默认值
USAGE
}

log_msg() {
  echo "[generate_samplesheet] $*"
}

INPUT=""
OUTPUT="./samplesheet.csv"
LAYOUT="auto"
R1_PATTERN="_R1"
R2_PATTERN="_R2"
DEFAULT_CONDITION="NA"
DEFAULT_REPLICATE="NA"
METADATA_CSV=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input) INPUT="$2"; shift 2 ;;
    --output) OUTPUT="$2"; shift 2 ;;
    --layout) LAYOUT="$2"; shift 2 ;;
    --r1-pattern) R1_PATTERN="$2"; shift 2 ;;
    --r2-pattern) R2_PATTERN="$2"; shift 2 ;;
    --condition) DEFAULT_CONDITION="$2"; shift 2 ;;
    --replicate) DEFAULT_REPLICATE="$2"; shift 2 ;;
    --metadata-csv) METADATA_CSV="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "未知参数: $1" >&2; usage; exit 1 ;;
  esac
done

if [[ -z "${INPUT}" ]]; then
  echo "错误：必须提供 --input" >&2
  usage
  exit 1
fi
if [[ ! -d "${INPUT}" ]]; then
  echo "错误：输入目录不存在: ${INPUT}" >&2
  exit 1
fi
if [[ -n "${METADATA_CSV}" && ! -f "${METADATA_CSV}" ]]; then
  echo "错误：metadata 文件不存在: ${METADATA_CSV}" >&2
  exit 1
fi

LAYOUT_UPPER="$(echo "${LAYOUT}" | tr '[:lower:]' '[:upper:]')"
if [[ "${LAYOUT_UPPER}" != "AUTO" && "${LAYOUT_UPPER}" != "PE" && "${LAYOUT_UPPER}" != "SE" ]]; then
  echo "错误：--layout 只能是 auto / PE / SE" >&2
  exit 1
fi

mkdir -p "$(dirname "${OUTPUT}")"
INPUT_REAL="$(python3 - <<'PY' "${INPUT}"
import os, sys
print(os.path.realpath(sys.argv[1]))
PY
)"

log_msg "开始扫描 FASTQ: ${INPUT_REAL}"
python3 - <<'PY' "${INPUT_REAL}" "${OUTPUT}" "${LAYOUT_UPPER}" "${R1_PATTERN}" "${R2_PATTERN}" "${DEFAULT_CONDITION}" "${DEFAULT_REPLICATE}" "${METADATA_CSV}"
import csv, re, sys
from pathlib import Path

input_dir, out_csv, layout_mode, r1_pattern, r2_pattern, default_condition, default_replicate, meta_csv = sys.argv[1:9]
FASTQ_SUFFIXES = ('.fastq.gz', '.fq.gz', '.fastq', '.fq')
PAIR_RE = re.compile(r'(?i)^(?P<sample>.+?)(?P<sep>[._-])(?P<read>r?[12])(?:[._-]?00[12])?$')


def log(*args):
    print('[generate_samplesheet]', *args, flush=True)


def strip_fastq_suffix(name: str):
    for suf in FASTQ_SUFFIXES:
        if name.endswith(suf):
            return name[:-len(suf)], suf
    return name, ''


def infer_sample_and_read(path: Path, r1_pat: str, r2_pat: str):
    stem, _ = strip_fastq_suffix(path.name)
    if r1_pat and r1_pat in stem:
        sample = stem.split(r1_pat, 1)[0].rstrip('._-') or path.parent.name
        mate_stem = stem.replace(r1_pat, r2_pat, 1)
        return sample, 'R1', mate_stem
    if r2_pat and r2_pat in stem:
        sample = stem.split(r2_pat, 1)[0].rstrip('._-') or path.parent.name
        mate_stem = stem.replace(r2_pat, r1_pat, 1)
        return sample, 'R2', mate_stem
    m = PAIR_RE.match(stem)
    if m:
        sample = m.group('sample').rstrip('._-') or path.parent.name
        read_token = m.group('read')
        read_upper = read_token.upper()
        read = 'R1' if read_upper in {'1', 'R1'} else 'R2'
        mate_token_map = {'1': '2', '2': '1', 'R1': 'R2', 'R2': 'R1', 'r1': 'r2', 'r2': 'r1'}
        mate_token = mate_token_map.get(read_token, 'R2' if read == 'R1' else 'R1')
        mate_stem = f"{m.group('sample')}{m.group('sep')}{mate_token}"
        return sample, read, mate_stem
    sample = stem.rstrip('._-') or path.parent.name
    return sample, 'SE', None


def build_path_lookup(files):
    lookup = {}
    for fp in files:
        stem, _ = strip_fastq_suffix(fp.name)
        lookup[(str(fp.parent), stem)] = fp
    return lookup


root = Path(input_dir)
files = sorted([p for p in root.rglob('*') if p.is_file() and p.name.endswith(FASTQ_SUFFIXES)])
if not files:
    raise SystemExit(f'错误：未在目录中找到 FASTQ 文件: {input_dir}')

path_lookup = build_path_lookup(files)
used = set()
rows = []
seen_samples = set()

for fp in files:
    if str(fp) in used:
        continue
    sample, read_type, mate_stem = infer_sample_and_read(fp, r1_pattern, r2_pattern)
    mate = path_lookup.get((str(fp.parent), mate_stem)) if mate_stem else None

    if layout_mode == 'PE':
        if read_type != 'R1':
            continue
        if mate is None:
            log(f'警告：缺少配对 R2，跳过: {fp}')
            continue
        rows.append({'sample': sample, 'layout': 'PE', 'condition': default_condition, 'replicate': default_replicate, 'r1': str(fp), 'r2': str(mate)})
        used.add(str(fp)); used.add(str(mate))
    elif layout_mode == 'SE':
        if read_type == 'R2':
            continue
        rows.append({'sample': sample, 'layout': 'SE', 'condition': default_condition, 'replicate': default_replicate, 'r1': str(fp), 'r2': ''})
        used.add(str(fp))
    else:
        if read_type == 'R1' and mate is not None:
            rows.append({'sample': sample, 'layout': 'PE', 'condition': default_condition, 'replicate': default_replicate, 'r1': str(fp), 'r2': str(mate)})
            used.add(str(fp)); used.add(str(mate))
        elif read_type == 'R2' and mate is not None:
            continue
        elif read_type == 'R2':
            log(f'警告：检测到疑似 R2 但未找到 R1，跳过: {fp}')
            used.add(str(fp))
        else:
            rows.append({'sample': sample, 'layout': 'SE', 'condition': default_condition, 'replicate': default_replicate, 'r1': str(fp), 'r2': ''})
            used.add(str(fp))

    if sample in seen_samples:
        log(f'警告：检测到重复 sample 名，保留多行: {sample}')
    seen_samples.add(sample)

meta_map = {}
if meta_csv:
    with open(meta_csv, newline='', encoding='utf-8-sig') as fh:
        reader = csv.DictReader(fh)
        cols = {c.lower().strip(): c for c in (reader.fieldnames or [])}
        if 'sample' not in cols:
            raise SystemExit('metadata-csv 缺少 sample 列')
        s_col = cols['sample']
        c_col = cols.get('condition')
        r_col = cols.get('replicate')
        for row in reader:
            sample = str(row.get(s_col, '')).strip()
            if not sample:
                continue
            meta_map[sample] = {
                'condition': str(row.get(c_col, '')).strip() if c_col else '',
                'replicate': str(row.get(r_col, '')).strip() if r_col else '',
            }

with open(out_csv, 'w', newline='', encoding='utf-8') as fh:
    writer = csv.DictWriter(fh, fieldnames=['sample','layout','condition','replicate','r1','r2'])
    writer.writeheader()
    for row in rows:
        meta = meta_map.get(row['sample'], {})
        row['condition'] = meta.get('condition') or row['condition']
        row['replicate'] = meta.get('replicate') or row['replicate']
        writer.writerow(row)

log(f'已生成: {out_csv}; rows={len(rows)}')
PY

log_msg "前几行预览："
head -n 6 "${OUTPUT}"
