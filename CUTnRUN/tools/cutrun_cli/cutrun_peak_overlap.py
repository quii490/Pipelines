#!/usr/bin/env python3
"""Calculate transparent pairwise peak overlap statistics without pybedtools."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_peak(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--peak expects LABEL=PATH")
    label, path = value.split("=", 1)
    if not label or not path:
        raise argparse.ArgumentTypeError("--peak expects non-empty LABEL=PATH")
    return label, Path(path)


def intervals(path: Path) -> dict[str, list[tuple[int, int]]]:
    result: dict[str, list[tuple[int, int]]] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if end > start:
                result.setdefault(fields[0], []).append((start, end))
    for chrom in result:
        result[chrom].sort()
    return result


def overlap(a: dict[str, list[tuple[int, int]]], b: dict[str, list[tuple[int, int]]]) -> tuple[int, int, int]:
    hit_a = 0
    hit_b = 0
    bases = 0
    for chrom, a_items in a.items():
        b_items = b.get(chrom, [])
        j = 0
        for start, end in a_items:
            while j < len(b_items) and b_items[j][1] <= start:
                j += 1
            k = j
            has = False
            while k < len(b_items) and b_items[k][0] < end:
                left, right = max(start, b_items[k][0]), min(end, b_items[k][1])
                if right > left:
                    has = True
                    bases += right - left
                k += 1
            hit_a += int(has)
    for chrom, b_items in b.items():
        a_items = a.get(chrom, [])
        j = 0
        for start, end in b_items:
            while j < len(a_items) and a_items[j][1] <= start:
                j += 1
            k = j
            if any(a_items[k2][0] < end and a_items[k2][1] > start for k2 in range(k, len(a_items))):
                hit_b += 1
    return hit_a, hit_b, bases


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--peak", action="append", required=True, type=parse_peak)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    parsed: dict[str, dict[str, list[tuple[int, int]]]] = {}
    for label, path in args.peak:
        if not path.is_file():
            raise SystemExit(f"peak file not found: {path}")
        parsed[label] = intervals(path)
    labels = list(parsed)
    rows: list[dict[str, object]] = []
    for i, first in enumerate(labels):
        first_count = sum(len(x) for x in parsed[first].values())
        for j, second in enumerate(labels):
            second_count = sum(len(x) for x in parsed[second].values())
            if i == j:
                hit_a, hit_b, bases = first_count, second_count, sum(e - s for xs in parsed[first].values() for s, e in xs)
            else:
                hit_a, hit_b, bases = overlap(parsed[first], parsed[second])
            rows.append({
                "first": first, "second": second,
                "first_peaks": first_count, "second_peaks": second_count,
                "first_overlapping": hit_a, "second_overlapping": hit_b,
                "first_fraction": hit_a / first_count if first_count else "",
                "second_fraction": hit_b / second_count if second_count else "",
                "overlap_bases": bases,
            })
    output = args.out_dir / "pairwise_overlap.tsv"
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    (args.out_dir / "labels.tsv").write_text("label\tpath\n" + "\n".join(f"{label}\t{path}" for label, path in args.peak) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
