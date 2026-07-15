#!/usr/bin/env python3
"""Generate compact classical and TE QC metrics from an existing run.

The script is deliberately resumable and reports unavailable metrics as SKIP
instead of silently dropping them.  FRiP is calculated with mapped alignment
counts and BEDTools read-overlap counts when both tools are available.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def run(cmd: list[str]) -> tuple[int, str]:
    try:
        proc = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        return proc.returncode, proc.stdout
    except OSError:
        return 127, ""


def flagstat(path: Path) -> dict[str, float | int]:
    if not path.is_file():
        return {}
    text = path.read_text(errors="replace")
    result: dict[str, float | int] = {}
    patterns = {
        "total_alignments": r"^(\d+) \+ \d+ in total",
        "mapped_alignments": r"^(\d+) \+ \d+ mapped \(([^)]+)%",
        "properly_paired": r"^(\d+) \+ \d+ properly paired \(([^)]+)%",
        "duplicates": r"^(\d+) \+ \d+ duplicates",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            result[key] = int(match.group(1))
            if len(match.groups()) > 1:
                try:
                    result[f"{key}_percent"] = float(match.group(2))
                except ValueError:
                    pass
    return result


def bed_stats(path: Path) -> tuple[int, int]:
    n = 0
    bases = 0
    if not path.is_file():
        return 0, 0
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                bases += max(0, int(fields[2]) - int(fields[1]))
                n += 1
            except ValueError:
                continue
    return n, bases


def frip(bam: Path, peak: Path, mapped: int) -> tuple[int | None, float | None, str]:
    if mapped <= 0 or not bam.is_file() or not peak.is_file():
        return None, None, "SKIP:missing_bam_peak_or_mapped_count"
    if not shutil.which("bedtools") or not shutil.which("samtools"):
        return None, None, "SKIP:bedtools_or_samtools_missing"
    code, text = run(["bash", "-lc", f"bedtools intersect -abam {shlex_quote(str(bam))} -b {shlex_quote(str(peak))} -u | samtools view -c -"])
    if code != 0:
        return None, None, f"SKIP:intersection_failed_rc{code}"
    try:
        count = int(text.strip())
    except ValueError:
        return None, None, "SKIP:invalid_intersection_count"
    return count, count / mapped, "PASS"


def shlex_quote(value: str) -> str:
    return "'" + value.replace("'", "'\\''") + "'"


def main() -> int:
    args = parse_args()
    results = args.results_dir.resolve()
    manifest = (args.manifest or results / "manifest" / "resolved_manifest.csv").resolve()
    if not manifest.is_file():
        fallback = results / "_automation" / "inputs" / "manifest.csv"
        if fallback.is_file():
            manifest = fallback
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")

    with manifest.open(newline="") as handle:
        samples = list(csv.DictReader(handle))
    output = (args.output or results / "09_downstream" / "qc_metrics.tsv").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample", "is_igg", "branch", "total_alignments", "mapped_alignments",
        "mapped_percent", "properly_paired", "properly_paired_percent", "duplicates",
        "duplicate_fraction", "nonduplicate_fraction",
        "peak_type", "peak_count", "peak_bases", "frip_reads", "frip", "frip_status",
        "complexity_status", "nsc_status", "rsc_status",
    ]
    rows: list[dict[str, object]] = []
    for sample_row in samples:
        sample = str(sample_row.get("sample", "")).strip()
        is_igg = str(sample_row.get("is_igg", "")).lower() in {"1", "true", "yes", "y"}
        for branch, bam_name in (("standard", f"{sample}_clean.bam"), ("te", f"{sample}_te.bam")):
            bam_dir = results / ("04_clean_bam" if branch == "standard" else "04_te_bam")
            bam = bam_dir / bam_name
            metrics = flagstat(bam.with_suffix(".flagstat.txt"))
            mapped = int(metrics.get("mapped_alignments", 0) or 0)
            peak_root = results / "06_peaks" / sample
            peak_specs = [
                ("broad", peak_root / "broad" / f"{sample}_peaks.broadPeak"),
                ("narrow", peak_root / "narrow" / f"{sample}_peaks.narrowPeak"),
            ]
            for peak_type, peak in peak_specs:
                count, bases = bed_stats(peak)
                hit, fraction, status = frip(bam, peak, mapped)
                row = {
                    "sample": sample,
                    "is_igg": str(is_igg).lower(),
                    "branch": branch,
                    **metrics,
                    "peak_type": peak_type,
                    "peak_count": count if peak.is_file() else "",
                    "peak_bases": bases if peak.is_file() else "",
                    "frip_reads": hit if hit is not None else "",
                    "frip": fraction if fraction is not None else "",
                    "frip_status": status if peak.is_file() else "SKIP:peak_missing",
                    "duplicate_fraction": (metrics.get("duplicates", 0) / metrics.get("total_alignments", 1)) if metrics.get("total_alignments") else "",
                    "nonduplicate_fraction": (1 - metrics.get("duplicates", 0) / metrics.get("total_alignments", 1)) if metrics.get("total_alignments") else "",
                    "complexity_status": "PASS" if mapped else "SKIP:no_mapped_reads",
                    "nsc_status": "SKIP:phantompeakqualtools_not_integrated",
                    "rsc_status": "SKIP:phantompeakqualtools_not_integrated",
                }
                rows.append(row)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "results_dir": str(results),
        "manifest": str(manifest),
        "rows": len(rows),
        "frip_pass": sum(row.get("frip_status") == "PASS" for row in rows),
        "frip_skipped": sum(str(row.get("frip_status", "")).startswith("SKIP") for row in rows),
        "output": str(output),
    }
    (output.with_suffix(".json")).write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
