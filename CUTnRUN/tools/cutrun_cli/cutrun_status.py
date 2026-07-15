#!/usr/bin/env python3
"""Attempt-aware status ledger for CUT&RUN downstream modules.

The old status file was append-only but had no run identifier and was
recreated on every invocation.  This utility keeps every attempt auditable and
also writes a latest-status view for humans and automation.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path


FIELDS = [
    "run_id",
    "module",
    "attempt",
    "status",
    "started_utc",
    "finished_utc",
    "reason",
    "outputs_ok",
    "output_paths",
]
VALID = {"PASS", "FAIL", "SKIP", "RUNNING"}


def now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        return []
    with path.open(newline="", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def ensure_ledger(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.is_file() and path.stat().st_size:
        with path.open(errors="replace") as handle:
            header = handle.readline().rstrip("\n").split("\t")
        if header == FIELDS:
            return
        legacy = path.with_name(f"{path.stem}.legacy_{datetime.now().strftime('%Y%m%d_%H%M%S')}{path.suffix}")
        shutil.move(str(path), str(legacy))
    with path.open("w", newline="") as handle:
        csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n").writeheader()


def write_latest(path: Path, rows: list[dict[str, str]], run_id: str | None = None) -> Path:
    latest_path = path.with_name(path.stem + "_latest" + path.suffix)
    latest: dict[tuple[str, str], dict[str, str]] = {}
    for row in rows:
        if run_id and row.get("run_id") != run_id:
            continue
        key = (row.get("run_id", ""), row.get("module", ""))
        latest[key] = row
    with latest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(latest.values())
    return latest_path


def cmd_init(args: argparse.Namespace) -> int:
    ensure_ledger(args.status_file)
    rows = read_rows(args.status_file)
    latest = write_latest(args.status_file, rows, args.run_id)
    meta = {
        "schema_version": 2,
        "run_id": args.run_id,
        "status_file": str(args.status_file.resolve()),
        "latest_file": str(latest.resolve()),
        "created_utc": now(),
    }
    if args.meta:
        args.meta.parent.mkdir(parents=True, exist_ok=True)
        args.meta.write_text(json.dumps(meta, indent=2, ensure_ascii=False) + "\n")
    return 0


def cmd_record(args: argparse.Namespace) -> int:
    if args.status not in VALID:
        raise SystemExit(f"invalid status: {args.status}")
    ensure_ledger(args.status_file)
    rows = read_rows(args.status_file)
    attempt = 1 + sum(
        1 for row in rows if row.get("run_id") == args.run_id and row.get("module") == args.module
    )
    row = {
        "run_id": args.run_id,
        "module": args.module,
        "attempt": str(attempt),
        "status": args.status,
        "started_utc": args.started or now(),
        "finished_utc": args.finished or now(),
        "reason": args.reason or "",
        "outputs_ok": "true" if args.outputs_ok else "false",
        "output_paths": args.output_paths or "",
    }
    with args.status_file.open("a", newline="") as handle:
        csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n").writerow(row)
    write_latest(args.status_file, rows + [row], args.run_id)
    return 0


def cmd_finalize(args: argparse.Namespace) -> int:
    rows = read_rows(args.status_file)
    latest_path = write_latest(args.status_file, rows, args.run_id)
    relevant = [row for row in rows if not args.run_id or row.get("run_id") == args.run_id]
    latest: dict[str, dict[str, str]] = {}
    for row in relevant:
        latest[row.get("module", "")] = row
    counts = {status: sum(row.get("status") == status for row in latest.values()) for status in sorted(VALID)}
    required_fail = [module for module, row in latest.items() if row.get("status") == "FAIL" and not module.startswith("optional:")]
    summary = {
        "schema_version": 2,
        "run_id": args.run_id or "",
        "generated_utc": now(),
        "latest_status_file": str(latest_path.resolve()),
        "module_count": len(latest),
        "counts": counts,
        "required_failures": sorted(required_fail),
        "status": "FAIL" if required_fail else ("INCOMPLETE" if counts.get("RUNNING", 0) else "PASS"),
    }
    inventory_summary = args.status_file.parent / "run_summary.json"
    if inventory_summary.is_file():
        try:
            inventory = json.loads(inventory_summary.read_text())
            if inventory.get("status") and inventory.get("status") != "PASS" and summary["status"] == "PASS":
                summary["status"] = "INCOMPLETE"
                summary["incomplete_reason"] = "required outputs missing according to results summary"
        except json.JSONDecodeError:
            summary["incomplete_reason"] = "run_summary.json is invalid"
    output = args.output or args.status_file.with_name("module_status_summary.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, ensure_ascii=False))
    return 1 if summary["status"] == "FAIL" else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    init = sub.add_parser("init")
    init.add_argument("--status-file", required=True, type=Path)
    init.add_argument("--run-id", required=True)
    init.add_argument("--meta", type=Path)
    init.set_defaults(func=cmd_init)

    record = sub.add_parser("record")
    record.add_argument("--status-file", required=True, type=Path)
    record.add_argument("--run-id", required=True)
    record.add_argument("--module", required=True)
    record.add_argument("--status", required=True)
    record.add_argument("--started")
    record.add_argument("--finished")
    record.add_argument("--reason", default="")
    record.add_argument("--outputs-ok", action="store_true")
    record.add_argument("--output-paths", default="")
    record.set_defaults(func=cmd_record)

    finalize = sub.add_parser("finalize")
    finalize.add_argument("--status-file", required=True, type=Path)
    finalize.add_argument("--run-id")
    finalize.add_argument("--output", type=Path)
    finalize.set_defaults(func=cmd_finalize)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
