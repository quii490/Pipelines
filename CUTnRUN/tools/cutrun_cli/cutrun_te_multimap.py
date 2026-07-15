#!/usr/bin/env python3
"""Summarize the effective TE multi-mapping policy for every TE BAM."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--max-records", type=int, default=0, help="0 scans all alignments")
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    samtools = shutil.which("samtools")
    rows: list[dict[str, object]] = []
    if not samtools:
        args.output.write_text("sample\tstatus\treason\n\tSKIP\tsamtools missing\n")
        return 0
    for bam in sorted(args.bam_dir.glob("*_te.bam")):
        total = secondary = nh_tagged = nh1 = nh_multi = 0
        proc = subprocess.Popen([samtools, "view", "-F", "4", str(bam)], stdout=subprocess.PIPE, text=True)
        assert proc.stdout is not None
        for line in proc.stdout:
            if args.max_records and total >= args.max_records:
                break
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            total += 1
            flag = int(fields[1])
            secondary += int(bool(flag & 256))
            nh = None
            for tag in fields[11:]:
                if tag.startswith("NH:i:"):
                    try:
                        nh = int(tag.split(":", 2)[2])
                    except ValueError:
                        nh = None
                    break
            if nh is not None:
                nh_tagged += 1
                if nh == 1:
                    nh1 += 1
                elif nh > 1:
                    nh_multi += 1
        proc.stdout.close()
        proc.wait()
        rows.append({
            "sample": bam.name[:-7] if bam.name.endswith("_te.bam") else bam.stem,
            "status": "PASS" if proc.returncode == 0 else "FAIL",
            "total_mapped_alignments": total,
            "secondary_alignments": secondary,
            "secondary_fraction": secondary / total if total else "",
            "nh_tagged": nh_tagged,
            "nh_tagged_fraction": nh_tagged / total if total else "",
            "nh1_alignments": nh1,
            "nh_gt1_alignments": nh_multi,
            "nh_gt1_fraction_of_tagged": nh_multi / nh_tagged if nh_tagged else "",
            "records_scanned": total,
            "scan_scope": "all" if not args.max_records else f"first_{args.max_records}",
        })
    fields = list(rows[0]) if rows else ["sample", "status", "reason"]
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    summary = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "bam_dir": str(args.bam_dir.resolve()),
        "sample_count": len(rows),
        "policy": "retain TE multimappers; report NH/secondary sensitivity; no cross-tool mixing",
        "output": str(args.output.resolve()),
    }
    args.output.with_suffix(".json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
