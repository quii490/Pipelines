#!/usr/bin/env python3
"""Create a compact, machine-readable inventory of a completed run."""

from __future__ import annotations

import argparse
import csv
import json
import re
from datetime import datetime, timezone
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-dir", required=True, type=Path)
    p.add_argument("--manifest", type=Path, default=None)
    return p.parse_args()


def first_match(directory: Path, pattern: str) -> Path | None:
    matches = sorted(directory.glob(pattern)) if directory.is_dir() else []
    return matches[0] if matches else None


def first_recursive_file(directory: Path, pattern: str = "*") -> Path | None:
    """Return the first regular file below a directory (including subdirs)."""
    if not directory.is_dir():
        return None
    matches = sorted(path for path in directory.rglob(pattern) if path.is_file())
    return matches[0] if matches else None


def count_lines(path: Path | None) -> int | None:
    if not path or not path.is_file():
        return None
    with path.open(errors="replace") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def flagstat_metrics(path: Path | None) -> dict[str, object]:
    if not path or not path.is_file():
        return {}
    text = path.read_text(errors="replace")
    result: dict[str, object] = {}
    patterns = {
        "total": r"^(\d+) \+ \d+ in total",
        "mapped": r"^(\d+) \+ \d+ mapped \(([^)]+)%",
        "properly_paired": r"^(\d+) \+ \d+ properly paired \(([^)]+)%",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text, re.MULTILINE)
        if match:
            result[key] = int(match.group(1))
            if len(match.groups()) > 1:
                result[f"{key}_percent"] = float(match.group(2))
    return result


def te_count_metrics(path: Path | None) -> dict[str, object]:
    if not path or not path.is_file():
        return {}
    result: dict[str, object] = {}
    with path.open() as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if len(row) == 2 and row[0] not in {"metric", "count"}:
                try:
                    result[row[0]] = int(row[1])
                except ValueError:
                    result[row[0]] = row[1]
    return result


def status_for(path: Path | None) -> str:
    return "PASS" if path and path.is_file() and path.stat().st_size > 0 else "MISSING"


def main() -> int:
    args = parse_args()
    results = args.results_dir.resolve()
    manifest = args.manifest or (results / "manifest" / "resolved_manifest.csv")
    if not manifest.is_file():
        manifest = results / "_automation" / "inputs" / "manifest.csv"
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")

    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    outdir = results / "09_downstream"
    outdir.mkdir(parents=True, exist_ok=True)
    inventory_path = outdir / "results_manifest.tsv"
    qc_path = outdir / "qc_metrics.tsv"

    inventory_fields = [
        "sample", "is_igg", "igg", "layout",
        "standard_bam", "te_bam", "standard_bam_bytes", "te_bam_bytes",
        "standard_bai", "te_bai",
        "narrow_peak", "broad_peak", "narrow_peaks", "broad_peaks",
        "standard_track", "te_track", "locus_track",
        "annotation", "homer", "status", "required_reason",
    ]
    qc_fields = ["sample", "branch", "total", "mapped", "mapped_percent", "properly_paired", "properly_paired_percent", "final_secondary", "final_nh_tagged"]
    inventory: list[dict[str, object]] = []
    qc_rows: list[dict[str, object]] = []
    expected_targets = 0
    missing_required = 0

    for row in rows:
        sample = str(row.get("sample", "")).strip()
        is_igg = str(row.get("is_igg", "")).strip().lower() in {"1", "true", "yes", "y"}
        if not is_igg:
            expected_targets += 1
        clean_bam = results / "04_clean_bam" / f"{sample}_clean.bam"
        clean_bai = results / "04_clean_bam" / f"{sample}_clean.bam.bai"
        te_bam = results / "04_te_bam" / f"{sample}_te.bam"
        te_bai = results / "04_te_bam" / f"{sample}_te.bam.bai"
        narrow = first_match(results / "06_peaks" / sample, "narrow/*_peaks.narrowPeak")
        broad = first_match(results / "06_peaks" / sample, "broad/*_peaks.broadPeak")
        standard_track = first_match(results / "05_tracks", f"{sample}_*_rpgc.bw") or first_match(results / "05_tracks", f"{sample}_*_rpkm.bw")
        te_track = first_match(results / "05_tracks_te", f"{sample}_te_*.bw")
        locus_track = first_match(results / "05_tracks_te_locus_best", f"{sample}_te_locus_best_*.bw")
        annotation = first_recursive_file(outdir / "annotation" / sample, "*.gene_structure.tsv")
        homer = first_recursive_file(outdir / "homer" / sample, "STATUS.tsv") or first_recursive_file(outdir / "homer" / sample, "knownResults.html")
        required_reasons = []
        required_ok = clean_bam.is_file() and clean_bai.is_file() and te_bam.is_file() and te_bai.is_file() and standard_track is not None and te_track is not None
        if not clean_bam.is_file(): required_reasons.append("standard_bam_missing")
        if not clean_bai.is_file(): required_reasons.append("standard_bai_missing")
        if not te_bam.is_file(): required_reasons.append("te_bam_missing")
        if not te_bai.is_file(): required_reasons.append("te_bai_missing")
        if standard_track is None: required_reasons.append("standard_track_missing")
        if te_track is None: required_reasons.append("te_track_missing")
        if not is_igg:
            required_ok = required_ok and narrow is not None and broad is not None
            if narrow is None: required_reasons.append("narrow_peak_missing")
            if broad is None: required_reasons.append("broad_peak_missing")
        if not required_ok:
            missing_required += 1
        entry = {
            "sample": sample,
            "is_igg": str(is_igg).lower(),
            "igg": str(row.get("igg", "")).strip(),
            "layout": str(row.get("layout", "")).strip(),
            "standard_bam": str(clean_bam) if clean_bam.is_file() else "",
            "te_bam": str(te_bam) if te_bam.is_file() else "",
            "standard_bam_bytes": clean_bam.stat().st_size if clean_bam.is_file() else "",
            "te_bam_bytes": te_bam.stat().st_size if te_bam.is_file() else "",
            "standard_bai": str(clean_bai) if clean_bai.is_file() else "",
            "te_bai": str(te_bai) if te_bai.is_file() else "",
            "narrow_peak": str(narrow) if narrow else "",
            "broad_peak": str(broad) if broad else "",
            "narrow_peaks": count_lines(narrow),
            "broad_peaks": count_lines(broad),
            "standard_track": str(standard_track) if standard_track else "",
            "te_track": str(te_track) if te_track else "",
            "locus_track": str(locus_track) if locus_track else "",
            "annotation": status_for(annotation),
            "homer": status_for(homer),
            "status": "PASS" if required_ok else "MISSING_REQUIRED",
            "required_reason": ";".join(required_reasons),
        }
        inventory.append(entry)

        for branch, flagstat, counts in (
            ("standard", results / "04_clean_bam" / f"{sample}_clean.flagstat.txt", None),
            ("te", results / "04_te_bam" / f"{sample}_te.flagstat.txt", results / "04_te_bam" / f"{sample}_te.clean_counts.tsv"),
        ):
            metrics = flagstat_metrics(flagstat)
            metrics.update(te_count_metrics(counts))
            qc_rows.append({"sample": sample, "branch": branch, **metrics})

    with inventory_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=inventory_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(inventory)
    with qc_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=qc_fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(qc_rows)

    summary = {
        "schema_version": 2,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "results_dir": str(results),
        "manifest": str(manifest.resolve()),
        "sample_count": len(rows),
        "target_count": expected_targets,
        "control_count": len(rows) - expected_targets,
        "samples_missing_required_outputs": missing_required,
        "inventory": str(inventory_path),
        "qc_metrics": str(qc_path),
        "status": "PASS" if missing_required == 0 else "INCOMPLETE",
    }
    run_manifest = outdir / "run_manifest.json"
    if run_manifest.is_file():
        try:
            run_payload = json.loads(run_manifest.read_text())
            summary["run_id"] = run_payload.get("run_id", "")
            summary["run_manifest"] = str(run_manifest)
        except json.JSONDecodeError:
            summary["run_manifest_status"] = "INVALID_JSON"
    status_latest = outdir / "module_status_latest.tsv"
    if status_latest.is_file():
        summary["module_status_latest"] = str(status_latest)
    method_statuses = sorted(outdir.glob("te_methods/**/method_status.tsv"))
    summary["te_method_status_files"] = [str(path) for path in method_statuses]
    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 0 if missing_required == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
