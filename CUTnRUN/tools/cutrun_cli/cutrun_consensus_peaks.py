#!/usr/bin/env python3
"""Build replicate-supported consensus peaks for broad and narrow calls."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def truthy(value: str) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if end > start:
                result[f[0]].append((start, end))
    return result


def union_support(peaks: dict[str, dict[str, list[tuple[int, int]]]], min_support: int) -> list[tuple[str, int, int, int]]:
    events: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for sample, by_chrom in peaks.items():
        for chrom, intervals in by_chrom.items():
            for start, end in intervals:
                events[chrom].append((start, 1))
                events[chrom].append((end, -1))
    result = []
    for chrom, chrom_events in events.items():
        chrom_events.sort(key=lambda x: (x[0], -x[1]))
        depth = 0
        prev = None
        segment_start = None
        for coordinate, delta in chrom_events:
            if prev is not None and coordinate > prev and depth >= min_support:
                if segment_start is None:
                    segment_start = prev
            elif segment_start is not None and coordinate > prev:
                result.append((chrom, segment_start, prev, depth))
                segment_start = None
            depth += delta
            prev = coordinate
        if segment_start is not None and prev is not None:
            result.append((chrom, segment_start, prev, depth))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--min-support", type=int, default=2)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with args.manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    groups: dict[str, list[str]] = defaultdict(list)
    for row in rows:
        if truthy(row.get("is_igg", "")):
            continue
        sample = str(row.get("sample", "")).strip()
        group = str(row.get("group", "")).strip() or "all_targets"
        if sample:
            groups[group].append(sample)
    summary = []
    for peak_type, suffix in (("broad", ".broadPeak"), ("narrow", ".narrowPeak")):
        for group, samples in sorted(groups.items()):
            sample_peaks: dict[str, dict[str, list[tuple[int, int]]]] = {}
            for sample in samples:
                path = args.results_dir / "06_peaks" / sample / peak_type / f"{sample}_peaks{suffix}"
                if path.is_file() and path.stat().st_size:
                    sample_peaks[sample] = read_bed(path)
            if len(sample_peaks) < args.min_support:
                summary.append({"group": group, "peak_type": peak_type, "replicates": len(sample_peaks), "consensus_peaks": 0, "status": "SKIP"})
                continue
            consensus = union_support(sample_peaks, args.min_support)
            output = args.out_dir / f"{group}.{peak_type}.consensus.bed"
            with output.open("w") as handle:
                for i, (chrom, start, end, support) in enumerate(consensus, 1):
                    handle.write(f"{chrom}\t{start}\t{end}\t{group}_{peak_type}_consensus_{i}\t{support}\n")
            summary.append({"group": group, "peak_type": peak_type, "replicates": len(sample_peaks), "consensus_peaks": len(consensus), "status": "PASS"})
    summary_path = args.out_dir / "consensus_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["group", "peak_type", "replicates", "consensus_peaks", "status"], delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
