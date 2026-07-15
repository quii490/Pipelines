#!/usr/bin/env python3
"""Create a reproducible run manifest from a CUT&RUN results directory.

Large FASTQ/BAM files are represented by size, mtime and a streaming content
fingerprint unless ``--hash-large`` is requested.  Reference and configuration
files are always SHA256 hashed.  The result is intentionally independent of
Nextflow so it can also describe resumed or partially completed runs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path, max_bytes: int | None = None) -> str:
    digest = hashlib.sha256()
    remaining = max_bytes
    with path.open("rb") as handle:
        while True:
            size = 1024 * 1024 if remaining is None else min(1024 * 1024, remaining)
            if size <= 0:
                break
            block = handle.read(size)
            if not block:
                break
            digest.update(block)
            if remaining is not None:
                remaining -= len(block)
    return digest.hexdigest()


def file_info(value: str | Path, full: bool = False, always_hash: bool = False) -> dict[str, object]:
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        return {"path": str(path), "exists": False}
    result: dict[str, object] = {
        "path": str(path),
        "exists": True,
        "bytes": path.stat().st_size,
        "mtime_utc": datetime.fromtimestamp(path.stat().st_mtime, timezone.utc).isoformat(),
    }
    if always_hash or full or path.stat().st_size <= 128 * 1024 * 1024:
        result["sha256"] = sha256(path)
    else:
        result["sha256_prefix_128MiB"] = sha256(path, 128 * 1024 * 1024)
        result["hash_scope"] = "first_128MiB"
    return result


def version(command: list[str]) -> list[str]:
    try:
        proc = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False, timeout=10)
        return proc.stdout.strip().splitlines()[:3]
    except (OSError, subprocess.TimeoutExpired) as exc:
        return [f"UNAVAILABLE: {exc}"]


def read_manifest(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def reference_candidates(config_text: str) -> list[Path]:
    values: list[Path] = []
    for line in config_text.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip()
        value = value.split("//", 1)[0].strip().strip("',\"")
        if value.startswith("/") and key in {"bowtie2_index", "gene_saf", "te_saf", "gene_anno", "te_anno", "blacklist"}:
            values.append(Path(value))
    return values


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--run-id")
    parser.add_argument("--hash-large", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    results = args.results_dir.resolve()
    manifest = (args.manifest or results / "manifest" / "resolved_manifest.csv").resolve()
    if not manifest.is_file():
        fallback = results / "_automation" / "inputs" / "manifest.csv"
        if fallback.is_file():
            manifest = fallback
    if not manifest.is_file():
        raise SystemExit(f"manifest not found: {manifest}")

    run_id = args.run_id or os.environ.get("CUTRUN_RUN_ID") or datetime.now(timezone.utc).strftime("run_%Y%m%dT%H%M%SZ")
    config = args.config.resolve() if args.config else None
    rows = read_manifest(manifest)
    inputs: dict[str, dict[str, object]] = {"manifest": file_info(manifest, full=True)}
    for row in rows:
        sample = str(row.get("sample", "")).strip()
        for key in ("fastq_1", "fastq_2"):
            value = str(row.get(key, "")).strip()
            if value:
                inputs[f"{sample}:{key}"] = file_info(value, full=args.hash_large)

    references: dict[str, dict[str, object]] = {}
    if config and config.is_file():
        inputs["nextflow_config"] = file_info(config, full=True)
        for path in reference_candidates(config.read_text(errors="replace")):
            if not path.is_file() and list(path.parent.glob(path.name + "*.bt2*")):
                references[str(path)] = {
                    "path": str(path),
                    "exists": True,
                    "kind": "bowtie2_prefix",
                    "shards": [file_info(shard, full=True, always_hash=True) for shard in sorted(path.parent.glob(path.name + "*.bt2*"))],
                }
            else:
                references[str(path)] = file_info(path, full=True, always_hash=True)

    tools = {
        "python": version(["python3", "--version"]),
        "nextflow": version(["nextflow", "-version"]),
        "samtools": version(["samtools", "--version"]),
        "macs3": version(["macs3", "--version"]),
        "bowtie2": version(["bowtie2", "--version"]),
        "bedtools": version(["bedtools", "--version"]),
        "computeMatrix": version(["computeMatrix", "--version"]),
        "Rscript": version(["Rscript", "--version"]),
    }
    external_tools = {
        "t3e_python": Path(os.environ.get("T3E_PYTHON", "/path/to/.conda/envs/cutrun-t3e/bin/python")),
        "allo": Path(os.environ.get("ALLO_BIN", "/path/to/.conda/envs/cutrun-allo/bin/allo")),
        "repentools_ret": Path(os.environ.get("REPENTOOLS_BIN", "/path/to/.local/share/CUTnRUN-tools/RepEnTools/ret")),
    }
    for name, path in external_tools.items():
        tools[name] = {"path": str(path), "exists": path.is_file(), "version": version([str(path), "--version"]) if path.is_file() else []}

    payload = {
        "schema_version": 2,
        "run_id": run_id,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "results_dir": str(results),
        "host": platform.node(),
        "platform": platform.platform(),
        "inputs": inputs,
        "references": references,
        "sample_count": len(rows),
        "target_count": sum(str(row.get("is_igg", "")).lower() not in {"1", "true", "yes", "y"} for row in rows),
        "control_count": sum(str(row.get("is_igg", "")).lower() in {"1", "true", "yes", "y"} for row in rows),
        "parameters": {
            "macs3": {"always_broad_and_narrow": True},
            "peak_annotator": os.environ.get("PEAK_ANNOTATOR", "chipseeker"),
            "te_multimap_policy": os.environ.get("TE_MULTIMAP_POLICY", "retain_and_report_sensitivity"),
            "te_locus_max_regions": os.environ.get("TE_LOCUS_MAX_REGIONS", "100000"),
            "te_locus_bin_size": os.environ.get("TE_LOCUS_BIN_SIZE", "25"),
            "peak_heatmap_max_regions": os.environ.get("PEAK_HEATMAP_MAX_REGIONS", "100000"),
            "resource_tier": os.environ.get("CUTRUN_RESOURCE_TIER", "standard"),
            "te_heatmap_write_values": os.environ.get("TE_HEATMAP_WRITE_VALUES", "false"),
        },
        "software": tools,
    }
    output = (args.output or results / "09_downstream" / "run_manifest.json").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(output)
    print(json.dumps({"run_id": run_id, "output": str(output), "status": "PASS"}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
