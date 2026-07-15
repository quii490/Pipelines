#!/usr/bin/env python3
import argparse
import hashlib
import json
import shutil
from datetime import datetime
from pathlib import Path


def copy_if_exists(src, dest):
    if not src.exists():
        return
    if src.is_dir():
        shutil.copytree(src, dest, dirs_exist_ok=True, ignore=shutil.ignore_patterns("*.bam", "*.bai", "*.fastq*", "*.fq*"))
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)


def main():
    ap = argparse.ArgumentParser(description="Publish compact RNA-seq deliverables")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--outdir")
    args = ap.parse_args()
    results = Path(args.results_dir).resolve()
    out = Path(args.outdir).resolve() if args.outdir else results / "deliverables" / datetime.now().strftime("%Y%m%d_%H%M%S")
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"publish directory is not empty: {out}")
    out.mkdir(parents=True, exist_ok=True)
    for name in ["plots", "multiqc", "condition.csv", "contrast.csv", "pipeline_report.html", "pipeline_timeline.html"]:
        copy_if_exists(results / name, out / name)
    copy_if_exists(results / "_automation" / "run_manifest_latest.json", out / "provenance" / "run_manifest.json")
    logs = results / "_automation" / "logs"
    if logs.exists():
        for p in logs.glob("progress_latest.txt"):
            copy_if_exists(p, out / "provenance" / p.name)
        for p in logs.glob("*.trace.tsv"):
            copy_if_exists(p, out / "provenance" / p.name)
    manifest = []
    for p in sorted(x for x in out.rglob("*") if x.is_file()):
        digest = hashlib.sha256(p.read_bytes()).hexdigest()
        manifest.append({"path": str(p.relative_to(out)), "size": p.stat().st_size, "sha256": digest})
    (out / "publish_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"published={out} files={len(manifest)}")


if __name__ == "__main__":
    main()
