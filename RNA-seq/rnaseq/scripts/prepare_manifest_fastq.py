#!/usr/bin/env python3
import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
from collections import OrderedDict
from pathlib import Path


def log(*args):
    print('[prepare_manifest_fastq]', *args, flush=True)


def read_table(path):
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in {'.csv', '.tsv', '.txt'}:
        sep = '\t' if suffix in {'.tsv', '.txt'} else ','
        with path.open('r', encoding='utf-8-sig', newline='') as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            rows = [dict(r) for r in reader]
            fields = reader.fieldnames or []
        return fields, rows
    if suffix in {'.xlsx', '.xls'}:
        try:
            import pandas as pd
        except Exception as e:
            raise SystemExit(f'无法读取 Excel manifest，请安装 pandas/openpyxl: {e}')
        df = pd.read_excel(path)
        rows = df.fillna('').to_dict(orient='records')
        return list(df.columns), rows
    raise SystemExit(f'不支持的 manifest 格式: {path}')


def norm_key(s):
    return str(s).strip().lower().replace(' ', '_').replace('-', '_')


def norm_val(x):
    return str(x).strip() if x is not None else ''


def split_multi(val):
    val = norm_val(val)
    if not val:
        return []
    parts = []
    for chunk in val.replace('|', ';').split(';'):
        for x in chunk.split(','):
            x = x.strip()
            if x:
                parts.append(x)
    return parts


def pick(row, *aliases):
    for a in aliases:
        if a in row and norm_val(row[a]):
            return norm_val(row[a])
    return ''


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def append_fastq(src, dest):
    ensure_parent(dest)
    src = str(src)
    dest = str(dest)
    if src.endswith('.gz'):
        with open(src, 'rb') as rf, open(dest, 'ab') as wf:
            shutil.copyfileobj(rf, wf)
    else:
        with open(dest, 'ab') as wf:
            proc = subprocess.run(['gzip', '-c', src], stdout=wf)
            if proc.returncode != 0:
                raise RuntimeError(f'gzip 压缩失败: {src}')


def download_url(url, tmp_dir):
    name = os.path.basename(urllib.parse.urlparse(url).path) or 'download.fastq.gz'
    target = os.path.join(tmp_dir, name)
    cmd = ['curl', '-L', '--fail', '-o', target, url]
    if shutil.which('curl') is None:
        cmd = ['wget', '-O', target, url]
        if shutil.which('wget') is None:
            raise RuntimeError('未找到 curl 或 wget，无法下载 URL FASTQ')
    log('下载 URL:', url)
    subprocess.run(cmd, check=True)
    return target


def fetch_text(url):
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode('utf-8', errors='replace')


def resolve_runs_from_gsm(gsm):
    term = urllib.parse.quote(f'{gsm}[All Fields]')
    esearch = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&retmode=json&term={term}'
    data = json.loads(fetch_text(esearch))
    ids = data.get('esearchresult', {}).get('idlist', [])
    if not ids:
        return []
    efetch = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id=' + ','.join(ids) + '&rettype=runinfo&retmode=text'
    txt = fetch_text(efetch)
    rows = list(csv.DictReader(txt.splitlines()))
    runs = []
    for row in rows:
        run = norm_val(row.get('Run', ''))
        if run:
          runs.append(run)
    return sorted(set(runs))


def sra_download_to_fastq(run_acc, out_r1, out_r2, threads=8):
    if shutil.which('fasterq-dump') is None:
        raise RuntimeError('未找到 fasterq-dump，请先安装 SRA Toolkit')
    with tempfile.TemporaryDirectory(prefix=f'sra_{run_acc}_') as tmpdir:
        tmpdir_p = Path(tmpdir)
        if shutil.which('prefetch') is not None:
            log('prefetch:', run_acc)
            subprocess.run(['prefetch', run_acc, '-O', tmpdir], check=True)
            run_input = str(next(tmpdir_p.glob(f'**/{run_acc}.sra'), run_acc))
        else:
            run_input = run_acc
        log('fasterq-dump:', run_acc)
        subprocess.run(['fasterq-dump', '--split-files', '--threads', str(threads), '-O', tmpdir, run_input], check=True)
        single = tmpdir_p / f'{run_acc}.fastq'
        r1 = tmpdir_p / f'{run_acc}_1.fastq'
        r2 = tmpdir_p / f'{run_acc}_2.fastq'
        if r1.exists() and r2.exists():
            append_fastq(r1, out_r1)
            append_fastq(r2, out_r2)
            return 'PE'
        if single.exists():
            append_fastq(single, out_r1)
            return 'SE'
        raise RuntimeError(f'未找到 fasterq-dump 输出: {run_acc}')


def main():
    ap = argparse.ArgumentParser(description='从 manifest 准备 FASTQ，并输出样本元信息')
    ap.add_argument('--manifest', required=True)
    ap.add_argument('--fastq-dir', required=True)
    ap.add_argument('--sample-metadata-out', required=True)
    ap.add_argument('--resolved-manifest-out', required=True)
    ap.add_argument('--sra-threads', type=int, default=8)
    args = ap.parse_args()

    fields, rows = read_table(args.manifest)
    if not rows:
        raise SystemExit('manifest 为空')

    lower_alias = {norm_key(c): c for c in fields}
    norm_rows = []
    for row in rows:
        norm_rows.append({norm_key(k): norm_val(v) for k, v in row.items()})

    fastq_dir = Path(args.fastq_dir)
    fastq_dir.mkdir(parents=True, exist_ok=True)

    meta = OrderedDict()
    resolved = []

    for idx, row in enumerate(norm_rows, start=1):
        sample = pick(row, 'sample', 'sample_name', 'sampleid', 'sample_id', 'name')
        gsm = pick(row, 'gsm', 'geo_accession', 'geo_sample', 'sample_accession')
        run_val = pick(row, 'run', 'srr', 'run_accession', 'accession')
        condition = pick(row, 'condition', 'group', 'treatment', 'case_control', 'status')
        replicate = pick(row, 'replicate', 'rep', 'repeat')
        layout = pick(row, 'layout', 'library_layout', 'paired_end')
        if not sample:
            sample = gsm or run_val or f'sample_{idx}'
        sample = sample.replace(' ', '_').replace('/', '_')
        if not condition:
            condition = 'NA'
        if not replicate:
            replicate = 'NA'

        local_r1 = split_multi(pick(row, 'r1', 'read1', 'fastq_1', 'fq1', 'path1', 'file1'))
        local_r2 = split_multi(pick(row, 'r2', 'read2', 'fastq_2', 'fq2', 'path2', 'file2'))
        local_single = split_multi(pick(row, 'fastq', 'fq', 'path', 'file', 'fastq_path'))
        url_r1 = split_multi(pick(row, 'url1', 'fastq_url_1', 'r1_url', 'ftp_1', 'fastq_ftp_1'))
        url_r2 = split_multi(pick(row, 'url2', 'fastq_url_2', 'r2_url', 'ftp_2', 'fastq_ftp_2'))
        url_single = split_multi(pick(row, 'url', 'fastq_url', 'ftp', 'fastq_ftp'))
        runs = split_multi(run_val)
        if not runs and gsm:
            try:
                runs = resolve_runs_from_gsm(gsm)
                if runs:
                    log(f'GSM 解析到 runs: {gsm} -> {",".join(runs)}')
            except Exception as e:
                log(f'警告：GSM 解析失败 {gsm}: {e}')

        out_r1 = fastq_dir / f'{sample}_R1.fastq.gz'
        out_r2 = fastq_dir / f'{sample}_R2.fastq.gz'
        out_se = fastq_dir / f'{sample}.fastq.gz'

        # 清理旧的目标，避免重复附加旧文件
        for tgt in (out_r1, out_r2, out_se):
            if tgt.exists():
                tgt.unlink()

        final_layout = ''
        with tempfile.TemporaryDirectory(prefix=f'manifest_{sample}_') as tmpdir:
            if local_r1 or local_r2:
                if len(local_r1) != len(local_r2):
                    raise RuntimeError(f'{sample}: 本地 R1/R2 数量不一致')
                for p1, p2 in zip(local_r1, local_r2):
                    append_fastq(Path(args.manifest).parent / p1 if not os.path.isabs(p1) else p1, out_r1)
                    append_fastq(Path(args.manifest).parent / p2 if not os.path.isabs(p2) else p2, out_r2)
                final_layout = 'PE'
            elif local_single:
                for p in local_single:
                    append_fastq(Path(args.manifest).parent / p if not os.path.isabs(p) else p, out_se)
                final_layout = 'SE'
            elif url_r1 or url_r2:
                if len(url_r1) != len(url_r2):
                    raise RuntimeError(f'{sample}: URL R1/R2 数量不一致')
                for u1, u2 in zip(url_r1, url_r2):
                    f1 = download_url(u1, tmpdir)
                    f2 = download_url(u2, tmpdir)
                    append_fastq(f1, out_r1)
                    append_fastq(f2, out_r2)
                final_layout = 'PE'
            elif url_single:
                for u in url_single:
                    f = download_url(u, tmpdir)
                    append_fastq(f, out_se)
                final_layout = 'SE'
            elif runs:
                seen_layout = set()
                for run_acc in runs:
                    lyt = sra_download_to_fastq(run_acc, out_r1, out_r2, threads=args.sra_threads)
                    seen_layout.add(lyt)
                if len(seen_layout) != 1:
                    raise RuntimeError(f'{sample}: 多个 run 的 layout 不一致: {seen_layout}')
                final_layout = list(seen_layout)[0]
                if final_layout == 'SE' and out_r1.exists() and not out_se.exists():
                    out_r1.rename(out_se)
            else:
                raise RuntimeError(f'{sample}: manifest 行无法解析到 fastq/path/url/run/GSM')

        if layout:
            layout_norm = layout.upper()
            if layout_norm.startswith('P'):
                layout_norm = 'PE'
            elif layout_norm.startswith('S'):
                layout_norm = 'SE'
            elif layout_norm in {'TRUE', 'T', 'YES', 'Y', '1'}:
                layout_norm = 'PE'
            elif layout_norm in {'FALSE', 'F', 'NO', 'N', '0'}:
                layout_norm = 'SE'
            if layout_norm in {'PE', 'SE'} and final_layout and layout_norm != final_layout:
                raise RuntimeError(f'{sample}: manifest layout={layout_norm} 与下载结果 {final_layout} 不一致')
            final_layout = layout_norm

        meta[sample] = {'sample': sample, 'condition': condition, 'replicate': replicate}
        resolved.append({
            'sample': sample,
            'condition': condition,
            'replicate': replicate,
            'layout': final_layout,
            'r1': str(out_r1) if final_layout == 'PE' else str(out_se),
            'r2': str(out_r2) if final_layout == 'PE' else '',
            'gsm': gsm,
            'runs': ';'.join(runs)
        })
        log(f'完成: sample={sample} layout={final_layout}')

    with open(args.sample_metadata_out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['sample', 'condition', 'replicate'])
        writer.writeheader()
        for row in meta.values():
            writer.writerow(row)

    with open(args.resolved_manifest_out, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=['sample', 'condition', 'replicate', 'layout', 'r1', 'r2', 'gsm', 'runs'])
        writer.writeheader()
        writer.writerows(resolved)

    log('metadata 输出:', args.sample_metadata_out)
    log('resolved manifest 输出:', args.resolved_manifest_out)


if __name__ == '__main__':
    main()
