#!/usr/bin/env python3
"""Create a normalized target/control plan shared by TE method adapters."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path


def truthy(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    manifest = Path(args.manifest).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError("manifest is empty")
    required = {"sample", "species", "is_igg"}
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"manifest is missing columns: {sorted(missing)}")

    species = {row["species"].strip() for row in rows if row["species"].strip()}
    if len(species) != 1:
        raise ValueError(f"TE method adapters require one species per run, got {sorted(species)}")
    controls = {row["sample"].strip() for row in rows if truthy(row["is_igg"])}
    samples = {row["sample"].strip() for row in rows}
    plan = []
    for row in rows:
        sample = row["sample"].strip()
        if truthy(row["is_igg"]):
            plan.append({"sample": sample, "is_igg": True, "control": ""})
            continue
        control = row.get("igg", "").strip()
        if not control and len(controls) == 1:
            control = next(iter(controls))
        if control and control not in samples:
            raise ValueError(f"sample {sample} refers to unknown control {control}")
        if not control and len(controls) > 1:
            raise ValueError(f"sample {sample} has no unambiguous control")
        plan.append({"sample": sample, "is_igg": False, "control": control})

    with (output_dir / "method_plan.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["sample", "is_igg", "control"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(plan)
    payload = {
        "schema_version": 1,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "manifest": str(manifest),
        "species": next(iter(species)),
        "controls": sorted(controls),
        "targets": [item["sample"] for item in plan if not item["is_igg"]],
        "plan": plan,
    }
    (output_dir / "method_plan.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(f"species={next(iter(species))}")
    print(f"controls={','.join(sorted(controls))}")
    print(f"targets={','.join(payload['targets'])}")


if __name__ == "__main__":
    main()
