#!/usr/bin/env python3
"""Write an immutable, machine-readable description of a pipeline launch."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import subprocess
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def capture(command: list[str], cwd: Path | None = None) -> str:
    try:
        result = subprocess.run(
            command,
            cwd=cwd,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        return result.stdout.strip()
    except OSError as error:
        return f"UNAVAILABLE: {error}"


def capture_first_lines(command: list[str], n: int = 3) -> list[str]:
    """Capture a short version banner without making optional tools fatal."""
    return capture(command).splitlines()[:n]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--main-nf", required=True)
    parser.add_argument("--nextflow-config", required=True)
    parser.add_argument("--run-id", default=os.environ.get("CUTRUN_RUN_ID", ""))
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    output = Path(args.output).resolve()
    manifest = Path(args.manifest).resolve()
    main_nf = Path(args.main_nf).resolve()
    config = Path(args.nextflow_config).resolve()
    project = main_nf.parent.parent.parent
    command = args.command[1:] if args.command[:1] == ["--"] else args.command

    git_status = capture(
        ["git", "-c", f"safe.directory={project}", "status", "--short"], cwd=project
    )
    git_diff = capture(
        ["git", "-c", f"safe.directory={project}", "diff", "--binary"], cwd=project
    )
    payload = {
        "schema_version": 1,
        "run_id": args.run_id,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "host": platform.node(),
        "platform": platform.platform(),
        "git": {
            "commit": capture(
                ["git", "-c", f"safe.directory={project}", "rev-parse", "HEAD"], cwd=project
            ),
            "branch": capture(
                ["git", "-c", f"safe.directory={project}", "branch", "--show-current"], cwd=project
            ),
            "dirty": bool(git_status),
            "status": git_status.splitlines(),
            "tracked_diff_sha256": hashlib.sha256(git_diff.encode()).hexdigest(),
        },
        "inputs": {
            "manifest": {"path": str(manifest), "sha256": sha256(manifest)},
            "main_nf": {"path": str(main_nf), "sha256": sha256(main_nf)},
            "nextflow_config": {"path": str(config), "sha256": sha256(config)},
        },
        "software": {
            "java": capture_first_lines(["java", "-version"], 2),
            "nextflow": capture_first_lines(["nextflow", "-version"], 3),
            "python": platform.python_version(),
            "tools": {
                "samtools": capture_first_lines(["samtools", "--version"]),
                "macs3": capture_first_lines(["macs3", "--version"]),
                "bowtie2": capture_first_lines(["bowtie2", "--version"]),
                "bedtools": capture_first_lines(["bedtools", "--version"]),
                "deepTools": capture_first_lines(["computeMatrix", "--version"]),
                "R": capture_first_lines(["Rscript", "--version"]),
            },
        },
        "environment": {
            key: os.environ.get(key, "")
            for key in ("JAVA_HOME", "NXF_HOME", "NXF_CONDA_CACHEDIR", "CHIPSEQ_ENV")
        },
        "wrapper_args_shell_escaped": os.environ.get("CUTRUN_WRAPPER_ARGS", ""),
        "command": command,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(output.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(output)


if __name__ == "__main__":
    main()
