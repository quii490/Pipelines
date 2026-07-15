#!/usr/bin/env python3
import argparse
import getpass
import json
import os
import platform
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def git_info(repo):
    try:
        commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True, stderr=subprocess.DEVNULL).strip())
        return {"commit": commit, "dirty": dirty}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--repo", required=True)
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--command", required=True)
    ap.add_argument("--param", action="append", default=[])
    args = ap.parse_args()
    params = {}
    for item in args.param:
        key, _, value = item.partition("=")
        params[key] = value
    sample_path = Path(args.samplesheet).resolve()
    payload = {
        "created_at": datetime.now(timezone.utc).isoformat(),
        "user": getpass.getuser(),
        "uid": os.getuid(),
        "host": platform.node(),
        "platform": platform.platform(),
        "command": args.command,
        "params": params,
        "samplesheet": {"path": str(sample_path), "size": sample_path.stat().st_size if sample_path.exists() else None},
        "executables": {name: shutil.which(name) for name in ["nextflow", "java", "conda", "python3", "Rscript"]},
        "git": git_info(Path(args.repo).resolve()),
    }
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
