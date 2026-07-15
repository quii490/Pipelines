#!/usr/bin/env python3
import csv
import os
import sys
from pathlib import Path

print('[prepare_manifest_fastq] start', file=sys.stderr)
if len(sys.argv) < 3:
    print('Usage: prepare_manifest_fastq.py <input_samplesheet.csv> <output_metadata.csv>', file=sys.stderr)
    sys.exit(1)

inp, outp = sys.argv[1], sys.argv[2]
rows = []
with open(inp, newline='', encoding='utf-8-sig') as f:
    reader = csv.DictReader(f)
    for r in reader:
        rows.append({
            'sample': r.get('sample', ''),
            'condition': r.get('condition', 'NA') or 'NA',
            'replicate': r.get('replicate', 'NA') or 'NA',
            'layout': r.get('layout', 'PE') or 'PE',
            'note': ''
        })
with open(outp, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=['sample', 'condition', 'replicate', 'layout', 'note'])
    writer.writeheader()
    writer.writerows(rows)
print(f'[prepare_manifest_fastq] wrote {len(rows)} rows to {outp}', file=sys.stderr)
