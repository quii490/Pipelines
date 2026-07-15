#!/usr/bin/env python3
"""Fail-fast validation for a CUTnRUN/chipseq manifest and its references.

The check is intentionally independent of Nextflow so a user can run it before
spending hours on trimming/alignment.  Structural problems are errors; a target
without a control is reported as a warning because MACS3 can still run without
one, but the choice is made visible in the run metadata.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample",
    "species",
    "is_igg",
    "layout",
    "fastq_1",
    "fastq_2",
    "igg",
}
REFERENCE_KEYS = (
    "bowtie2_index",
    "gene_saf",
    "te_saf",
    "gene_anno",
    "te_anno",
    "blacklist",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--species", choices=("hg38", "mm39"), default=None)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--repentools-index-dir", type=Path, default=None)
    parser.add_argument("--repentools-gtf", type=Path, default=None)
    parser.add_argument(
        "--skip-gzip",
        action="store_true",
        help="skip full gzip integrity checks (not recommended for first runs)",
    )
    parser.add_argument(
        "--min-free-gb",
        type=float,
        default=50.0,
        help="minimum free space required before starting (default: 50 GB)",
    )
    return parser.parse_args()


def bool_value(value: str | None) -> bool:
    return str(value or "").strip().lower() in {"1", "true", "yes", "y"}


def file_digest(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def reference_block(config_text: str, species: str) -> str:
    match = re.search(rf"(?m)^\s*{re.escape(species)}\s*:\s*\[(.*?)(?=^\s*\]\s*$)", config_text, re.S)
    return match.group(1) if match else ""


def reference_paths(config_text: str, species: str) -> dict[str, Path]:
    block = reference_block(config_text, species)
    result: dict[str, Path] = {}
    for key in REFERENCE_KEYS:
        match = re.search(rf"{re.escape(key)}\s*:\s*'([^']+)'", block)
        if match:
            result[key] = Path(match.group(1))
    return result


def check_gzip(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as handle:
            while handle.read(1024 * 1024):
                pass
    except (OSError, EOFError) as exc:
        return False, str(exc)
    return True, "ok"


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    warnings: list[str] = []
    checks: list[dict[str, object]] = []
    started = datetime.now(timezone.utc).isoformat()

    def record(name: str, status: str, detail: str) -> None:
        checks.append({"name": name, "status": status, "detail": detail})
        stream = sys.stderr if status == "ERROR" else sys.stdout
        print(f"[preflight] {status:7s} {name}: {detail}", file=stream)

    if not args.manifest.is_file():
        errors.append(f"manifest not found: {args.manifest}")
        record("manifest", "ERROR", str(args.manifest))
    if not args.config.is_file():
        errors.append(f"Nextflow config not found: {args.config}")
        record("config", "ERROR", str(args.config))

    rows: list[dict[str, str]] = []
    if not errors:
        try:
            with args.manifest.open(newline="") as handle:
                reader = csv.DictReader(handle)
                columns = set(reader.fieldnames or [])
                missing = sorted(REQUIRED_COLUMNS - columns)
                if missing:
                    errors.append("manifest missing columns: " + ", ".join(missing))
                    record("manifest columns", "ERROR", ", ".join(missing))
                rows = list(reader)
        except (OSError, csv.Error) as exc:
            errors.append(f"cannot read manifest: {exc}")
            record("manifest", "ERROR", str(exc))

    samples = [str(row.get("sample", "")).strip() for row in rows]
    duplicates = sorted({sample for sample in samples if sample and samples.count(sample) > 1})
    if duplicates:
        errors.append("duplicate sample names: " + ", ".join(duplicates))
        record("sample names", "ERROR", ", ".join(duplicates))
    elif rows:
        record("sample names", "PASS", f"{len(samples)} unique samples")

    observed_species = sorted({str(row.get("species", "")).strip() for row in rows if row.get("species")})
    if any(species not in {"hg38", "mm39"} for species in observed_species):
        errors.append("unsupported species in manifest: " + ", ".join(observed_species))
        record("species", "ERROR", ", ".join(observed_species))
    elif args.species and observed_species and observed_species != [args.species]:
        errors.append(f"manifest species {observed_species} does not match --species {args.species}")
        record("species", "ERROR", f"manifest={observed_species}, requested={args.species}")
    elif observed_species:
        record("species", "PASS", ", ".join(observed_species))

    controls = {
        str(row.get("sample", "")).strip()
        for row in rows
        if bool_value(row.get("is_igg"))
    }
    sample_set = set(samples)
    if not controls:
        warnings.append("no sample is marked is_igg=true; peak calling will run without control BAM")
        record("controls", "WARN", "none marked is_igg=true")
    else:
        record("controls", "PASS", ", ".join(sorted(controls)))
        mapped_controls = {
            str(row.get("igg", "")).strip()
            for row in rows
            if not bool_value(row.get("is_igg")) and str(row.get("igg", "")).strip()
        }
        if len(controls) > 1 and len(mapped_controls) == 1:
            only_control = next(iter(mapped_controls))
            warning = (
                f"multiple controls detected ({', '.join(sorted(controls))}) but all targets "
                f"map to {only_control}; verify the experimental design"
            )
            warnings.append(warning)
            record("control design", "WARN", warning)

    fastq_paths: list[Path] = []
    for row in rows:
        sample = str(row.get("sample", "")).strip() or "<empty sample>"
        layout = str(row.get("layout", "")).strip().upper()
        if layout not in {"PE", "SE"}:
            errors.append(f"{sample}: layout must be PE or SE, got {layout!r}")
            record(f"layout {sample}", "ERROR", layout or "empty")
        f1 = Path(str(row.get("fastq_1", "")).strip())
        f2 = Path(str(row.get("fastq_2", "")).strip())
        expected = [f1] if layout == "SE" else [f1, f2]
        for path in expected:
            if not str(path) or not path.is_file() or path.stat().st_size == 0:
                errors.append(f"{sample}: FASTQ missing or empty: {path}")
                record(f"FASTQ {sample}", "ERROR", str(path))
            else:
                fastq_paths.append(path)
        if layout == "PE" and (not str(f1) or not str(f2)):
            errors.append(f"{sample}: PE sample requires both fastq_1 and fastq_2")
        if layout == "PE" and f1.is_file() and f2.is_file():
            record(f"FASTQ pair {sample}", "PASS", f"{f1.name}, {f2.name}")
        elif layout == "SE" and f1.is_file():
            record(f"FASTQ {sample}", "PASS", f1.name)

        is_control = bool_value(row.get("is_igg"))
        igg = str(row.get("igg", "")).strip()
        if not is_control and igg:
            if igg not in sample_set:
                errors.append(f"{sample}: igg={igg} is not present in manifest")
                record(f"control link {sample}", "ERROR", igg)
            elif igg not in controls:
                errors.append(f"{sample}: igg={igg} is not marked is_igg=true")
                record(f"control link {sample}", "ERROR", igg)
            else:
                record(f"control link {sample}", "PASS", igg)
        elif not is_control and not igg:
            warnings.append(f"{sample}: no igg control assigned")
            record(f"control link {sample}", "WARN", "none")

    if not args.skip_gzip:
        for path in sorted(set(fastq_paths)):
            ok, detail = check_gzip(path)
            if not ok:
                errors.append(f"gzip integrity failed: {path}: {detail}")
                record(f"gzip {path.name}", "ERROR", detail)
            else:
                record(f"gzip {path.name}", "PASS", "ok")
    else:
        warnings.append("gzip integrity checks were skipped")
        record("gzip integrity", "WARN", "skipped by --skip-gzip")

    species = args.species or (observed_species[0] if len(observed_species) == 1 else "")
    refs: dict[str, Path] = {}
    if species and args.config.is_file():
        refs = reference_paths(args.config.read_text(), species)
        missing_keys = sorted(set(REFERENCE_KEYS) - set(refs))
        if missing_keys:
            errors.append(f"reference keys missing from config for {species}: {', '.join(missing_keys)}")
            record("reference config", "ERROR", ", ".join(missing_keys))
        for key, path in refs.items():
            if key == "bowtie2_index":
                index_files = sorted(path.parent.glob(path.name + "*.bt2")) + sorted(path.parent.glob(path.name + "*.bt2l"))
                if len(index_files) < 6:
                    errors.append(f"bowtie2 index incomplete for {species}: {path}")
                    record(f"reference {key}", "ERROR", str(path))
                else:
                    record(f"reference {key}", "PASS", f"{len(index_files)} index files")
            elif not path.is_file() or path.stat().st_size == 0:
                errors.append(f"reference missing or empty: {path}")
                record(f"reference {key}", "ERROR", str(path))
            else:
                record(f"reference {key}", "PASS", str(path))

    if args.repentools_index_dir or args.repentools_gtf:
        if args.repentools_index_dir:
            shards = [args.repentools_index_dir / f"chm13-2.{suffix}.ht2" for suffix in range(1, 9)]
            missing = [str(path) for path in shards if not path.is_file() or path.stat().st_size == 0]
            if missing:
                errors.append("RepEnTools index incomplete: " + ", ".join(missing))
                record("RepEnTools index", "ERROR", ", ".join(missing))
            else:
                record("RepEnTools index", "PASS", f"8 shards in {args.repentools_index_dir}")
        if args.repentools_gtf:
            if not args.repentools_gtf.is_file() or args.repentools_gtf.stat().st_size == 0:
                errors.append(f"RepEnTools GTF missing or empty: {args.repentools_gtf}")
                record("RepEnTools GTF", "ERROR", str(args.repentools_gtf))
            else:
                record("RepEnTools GTF", "PASS", str(args.repentools_gtf))

    required_commands = ("samtools", "bowtie2", "bedtools", "macs3", "featureCounts", "Rscript")
    for command in required_commands:
        path = shutil.which(command)
        if path:
            record(f"tool {command}", "PASS", path)
        else:
            errors.append(f"required tool not found in PATH: {command}")
            record(f"tool {command}", "ERROR", "not found in PATH")
    for command in ("bamCoverage", "computeMatrix", "plotHeatmap"):
        path = shutil.which(command)
        record(f"optional tool {command}", "PASS" if path else "WARN", path or "not found; related downstream module will SKIP")

    usage = shutil.disk_usage(args.outdir.parent if args.outdir.parent.exists() else Path.cwd())
    free_gb = usage.free / (1024**3)
    fastq_gb = sum(path.stat().st_size for path in set(fastq_paths) if path.is_file()) / (1024**3)
    required_gb = max(args.min_free_gb, fastq_gb * 2.0)
    if free_gb < required_gb:
        errors.append(f"free space {free_gb:.1f} GB is below estimated minimum {required_gb:.1f} GB")
        record("disk space", "ERROR", f"free={free_gb:.1f} GB, required>={required_gb:.1f} GB")
    else:
        record("disk space", "PASS", f"free={free_gb:.1f} GB, estimated minimum={required_gb:.1f} GB")

    metadata = {
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(args.manifest.resolve()),
        "config": str(args.config.resolve()),
        "species": species,
        "sample_count": len(rows),
        "target_count": sum(not bool_value(row.get("is_igg")) for row in rows),
        "control_count": len(controls),
        "fastq_count": len(set(fastq_paths)),
        "fastq_gb": round(fastq_gb, 3),
        "references": {key: str(path) for key, path in refs.items()},
        "repentools": {
            "index_dir": str(args.repentools_index_dir) if args.repentools_index_dir else "",
            "gtf": str(args.repentools_gtf) if args.repentools_gtf else "",
        },
        "checks": checks,
        "warnings": warnings,
        "errors": errors,
        "status": "PASS" if not errors else "FAIL",
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(f"[preflight] {'PASS' if not errors else 'FAIL'}: {len(rows)} samples, {len(controls)} controls; report={args.json_out}")
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
