#!/usr/bin/env python3
import argparse
import os
import shutil
import signal
import subprocess
from pathlib import Path


def alive(pid_file):
    try:
        pid = int(pid_file.read_text().strip())
        os.kill(pid, 0)
        return pid
    except (OSError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description="Safely preview or clean Nextflow work cache")
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--confirm", action="store_true")
    args = ap.parse_args()
    results = Path(args.results_dir).resolve()
    automation = results / "_automation"
    work = automation / "work"
    launch = automation / "nextflow_launch"
    if results == Path("/") or automation not in work.parents:
        raise SystemExit("unsafe results/work path")
    live = [(p, alive(p)) for p in (automation / "logs").glob("*.pid")]
    live = [(p, pid) for p, pid in live if pid]
    if live:
        raise SystemExit("refusing to clean: active PID(s): " + ", ".join(f"{pid} ({p.name})" for p, pid in live))
    size = "0"
    if work.exists():
        size = subprocess.run(["du", "-sh", str(work)], text=True, capture_output=True, check=False).stdout.split()[0]
    print(f"results={results}\nwork={work}\nwork_size={size}\nlaunch={launch}")
    if not args.confirm:
        print("preview only; add --confirm to clean completed Nextflow sessions")
        return
    if not launch.exists():
        raise SystemExit("persistent Nextflow launch directory not found; no cache was removed")
    subprocess.run(["nextflow", "clean", "-f"], cwd=launch, check=True)
    print("Nextflow cache cleanup completed")


if __name__ == "__main__":
    main()
