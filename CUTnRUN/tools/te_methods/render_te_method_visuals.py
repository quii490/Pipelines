#!/usr/bin/env python3
"""Render small, auditable visual summaries for the optional TE methods.

The upstream tools do not share an output format.  This adapter deliberately
does not invent values: a plot is written only when a real, parseable result
file is found.  Every method receives a row in ``visualization_status.tsv`` so
that a missing/unsupported tool is visible in the final report rather than
being mistaken for a successful analysis.

The renderer uses hand-written SVG so it works in the lightweight chipseq
environment without requiring matplotlib.  SVG is directly viewable in the
HTML report and remains editable for publication figures.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable, Iterator, Sequence


STATUS_HEADER = "method\tstatus\tdetail\toutputs\n"


def safe_text(value: object) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value)).strip("_") or "item"


def numeric(value: str) -> float | None:
    value = value.strip().replace(",", "")
    if not value or value.lower() in {"na", "nan", "none", "null", "-"}:
        return None
    value = value.rstrip("%")
    try:
        result = float(value)
    except ValueError:
        return None
    return result if math.isfinite(result) else None


def read_status(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    try:
        with path.open(newline="", errors="replace") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    except OSError:
        return []


def last_status(rows: Sequence[dict[str, str]], names: Iterable[str]) -> str:
    wanted = {name.lower() for name in names}
    for row in reversed(rows):
        method = (row.get("method") or row.get("group") or "").strip().lower()
        if method in wanted or any(method.startswith(name + ":") for name in wanted):
            return (row.get("status") or "UNKNOWN").strip().upper()
    return "MISSING"


def iter_table(path: Path, limit: int = 50000) -> tuple[list[str], list[list[str]], int]:
    """Read at most *limit* data rows from TSV/CSV/whitespace tables."""
    rows: list[list[str]] = []
    skipped = 0
    try:
        with path.open(newline="", errors="replace") as handle:
            for raw in handle:
                if not raw.strip() or raw.lstrip().startswith("#"):
                    continue
                line = raw.rstrip("\r\n")
                if "\t" in line:
                    parts = line.split("\t")
                elif "," in line:
                    parts = next(csv.reader([line]))
                else:
                    parts = line.split()
                if parts:
                    rows.append([part.strip() for part in parts])
                if len(rows) >= limit + 1:
                    break
    except (OSError, UnicodeError):
        return [], [], 0
    if not rows:
        return [], [], 0
    header = rows[0]
    data = rows[1:]
    # A table without a header is still useful when it has a clear numeric
    # second column; give the columns stable names instead of guessing labels.
    header_has_numeric = sum(numeric(value) is not None for value in header)
    if header_has_numeric > max(0, len(header) // 2):
        width = len(header)
        data = rows
        header = [f"column_{index + 1}" for index in range(width)]
    width = len(header)
    normalized = [row[:width] + [""] * max(0, width - len(row)) for row in data]
    truncated = max(0, len(rows) - 1 - limit)
    return header, normalized[:limit], truncated


def svg_bar(path: Path, title: str, labels: Sequence[str], values: Sequence[float],
            x_label: str = "value") -> None:
    width, height = 1000, max(260, 90 + len(labels) * 26)
    margin_left, margin_right = 230, 45
    plot_width = width - margin_left - margin_right
    vmax = max([0.0] + [abs(value) for value in values])
    vmax = vmax if vmax > 0 else 1.0
    rows = []
    for index, (label, value) in enumerate(zip(labels, values)):
        y = 68 + index * 26
        bar_width = max(1.0, plot_width * abs(value) / vmax)
        color = "#2563eb" if value >= 0 else "#dc2626"
        rows.append(
            f'<text x="{margin_left - 10}" y="{y + 5}" text-anchor="end" '
            f'font-size="12">{xml_escape(label[:34])}</text>'
            f'<rect x="{margin_left}" y="{y - 12}" width="{bar_width:.2f}" height="16" '
            f'fill="{color}" rx="3"/>'
            f'<text x="{min(width - 5, margin_left + bar_width + 6):.2f}" y="{y + 1}" '
            f'font-size="11">{value:.4g}</text>'
        )
    svg = (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">'
        '<rect width="100%" height="100%" fill="white"/>'
        f'<text x="24" y="28" font-size="18" font-family="Arial" font-weight="bold">'
        f'{xml_escape(title)}</text>'
        f'<text x="{margin_left}" y="{height - 14}" font-size="11" fill="#64748b">'
        f'{xml_escape(x_label)}</text>'
        + "".join(rows) + "</svg>"
    )
    path.write_text(svg)


def xml_escape(value: str) -> str:
    return (value.replace("&", "&amp;").replace("<", "&lt;")
            .replace(">", "&gt;").replace('"', "&quot;"))


def parse_flagstat(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    patterns = {
        "total": r"^(\d+) \+ \d+ in total",
        "mapped": r"^(\d+) \+ \d+ mapped",
        "properly_paired": r"^(\d+) \+ \d+ properly paired",
        "duplicates": r"^(\d+) \+ \d+ duplicates",
    }
    try:
        text = path.read_text(errors="replace")
    except OSError:
        return values
    for key, pattern in patterns.items():
        match = re.search(pattern, text, flags=re.MULTILINE)
        if match:
            values[key] = float(match.group(1))
    return values


def allo_visuals(methods_dir: Path, output_dir: Path, status_rows: list[dict[str, str]]) -> tuple[str, str]:
    allo_dir = methods_dir / "allo"
    flagstats = sorted(allo_dir.glob("*_allo.flagstat.txt"))
    out_dir = output_dir / "allo"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for flagstat in flagstats:
        values = parse_flagstat(flagstat)
        if "total" not in values:
            continue
        sample = flagstat.name[: -len("_allo.flagstat.txt")]
        total = values["total"]
        mapped = values.get("mapped", 0.0)
        rows.append({
            "sample": sample,
            "total": str(int(total)),
            "mapped": str(int(mapped)),
            "mapped_fraction": f"{mapped / total:.6f}" if total else "0",
            "properly_paired": str(int(values.get("properly_paired", 0))),
            "duplicates": str(int(values.get("duplicates", 0))),
        })
    table_path = out_dir / "allo_qc.tsv"
    with table_path.open("w") as handle:
        handle.write("sample\ttotal\tmapped\tmapped_fraction\tproperly_paired\tduplicates\n")
        for row in rows:
            handle.write("\t".join(row[key] for key in ("sample", "total", "mapped", "mapped_fraction", "properly_paired", "duplicates")) + "\n")
    if not rows:
        method_state = last_status(status_rows, {"allo"})
        if method_state == "PASS":
            return "FAIL", "Allo is PASS but no validated *_allo.flagstat.txt was found"
        return "SKIP", f"no validated Allo flagstat output (adapter status: {method_state})"
    labels = [row["sample"] for row in rows]
    values = [float(row["mapped_fraction"]) * 100.0 for row in rows]
    figure = out_dir / "allo_mapped_fraction.svg"
    svg_bar(figure, "Allo mapped fraction", labels, values, "mapped reads (%)")
    return "PASS", f"validated {len(rows)} Allo flagstat files; mapped-fraction plot written"


def candidate_tables(root: Path, method: str) -> Iterator[Path]:
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".txt", ".tsv", ".csv"}:
            continue
        name = path.name.lower()
        if any(token in name for token in ("log", "command", "stdout", "stderr", "flagstat")):
            continue
        if method == "t3e" and any(token in name for token in ("_counts", "_background", "probability")):
            continue
        yield path


def table_visuals(method: str, root: Path, output_dir: Path, status_rows: list[dict[str, str]]) -> tuple[str, str]:
    out_dir = output_dir / method
    out_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, str]] = []
    plot_candidates: list[tuple[float, Path, list[str], list[list[str]], int, int]] = []
    for path in candidate_tables(root, method):
        header, data, truncated = iter_table(path)
        if not header or not data:
            continue
        numeric_columns: list[tuple[int, int]] = []
        for index, column in enumerate(header):
            count = sum(numeric(row[index]) is not None for row in data if index < len(row))
            if count >= 2:
                numeric_columns.append((index, count))
        summaries.append({
            "file": str(path.relative_to(root)),
            "rows_read": str(len(data)),
            "truncated": str(truncated),
            "columns": str(len(header)),
            "numeric_columns": ",".join(header[index] for index, _ in numeric_columns),
        })
        if not numeric_columns:
            continue
        # Prefer a result/enrichment table and a named numeric effect column;
        # otherwise use the first numeric column with at least two values.
        name_score = 0
        lowered_name = path.name.lower()
        if any(token in lowered_name for token in ("enrich", "report", "summary", "result")):
            name_score += 10
        chosen_index, _ = max(numeric_columns, key=lambda item: (name_score, item[1]))
        label_index = 0 if chosen_index != 0 else (1 if len(header) > 1 else 0)
        valid = []
        for row in data:
            value = numeric(row[chosen_index]) if chosen_index < len(row) else None
            label = row[label_index] if label_index < len(row) else ""
            if value is not None and label:
                valid.append((abs(value), label, value))
        if len(valid) >= 2:
            valid.sort(reverse=True)
            plot_candidates.append((float(name_score), path, header, [[x[1], str(x[2])] for x in valid[:20]], chosen_index, label_index))
    summary_path = out_dir / f"{method}_table_summary.tsv"
    with summary_path.open("w") as handle:
        handle.write("file\trows_read\ttruncated\tcolumns\tnumeric_columns\n")
        for row in summaries:
            handle.write("\t".join(row[key] for key in ("file", "rows_read", "truncated", "columns", "numeric_columns")) + "\n")
    if not plot_candidates:
        method_state = last_status(status_rows, {method})
        if method_state == "PASS":
            return "FAIL", f"{method} is PASS but no parseable numeric result table was found"
        return "SKIP", f"no parseable {method} result table (adapter status: {method_state})"
    _, source, header, values, chosen_index, _ = max(plot_candidates, key=lambda item: (item[0], len(item[3])))
    labels = [value[0] for value in values]
    numbers = [float(value[1]) for value in values]
    figure = out_dir / f"{method}_top_features.svg"
    svg_bar(figure, f"{method} top features ({source.name})", labels, numbers, header[chosen_index])
    return "PASS", f"parsed {len(summaries)} result table(s); top-feature plot written"


def repentools_visuals(methods_dir: Path, output_dir: Path, status_rows: list[dict[str, str]]) -> tuple[str, str]:
    root = methods_dir / "repentools_fastq"
    out_root = output_dir / "repentools_fastq"
    out_root.mkdir(parents=True, exist_ok=True)
    groups = sorted(path for path in root.iterdir() if path.is_dir()) if root.is_dir() else []
    if not groups:
        state = last_status(status_rows, {"repentools", "manifest_validation"})
        return "SKIP", f"no RepEnTools group directory found (adapter status: {state})"
    group_rows: list[dict[str, str]] = []
    statuses: list[str] = []
    for group in groups:
        group_out = out_root / safe_text(group.name)
        group_out.mkdir(parents=True, exist_ok=True)
        group_status_rows = read_status(group / "method_status.tsv")
        if not group_status_rows:
            group_status_rows = status_rows
        # Report files are copied by the wrapper.  Requiring the two primary
        # sentinels prevents an empty ret directory from being labelled PASS.
        required = [group / "ret_experiment_summary.csv", group / "ret_report.csv"]
        missing = [path.name for path in required if not path.is_file() or path.stat().st_size == 0]
        if missing:
            state = last_status(group_status_rows, {group.name})
            statuses.append("FAIL" if state == "PASS" else "SKIP")
            group_rows.append({"group": group.name, "status": statuses[-1], "detail": "missing " + ",".join(missing)})
            continue
        state, detail = table_visuals("repentools", group, group_out, group_status_rows)
        statuses.append(state)
        group_rows.append({"group": group.name, "status": state, "detail": detail})
    summary_path = out_root / "visualization_summary.tsv"
    with summary_path.open("w") as handle:
        handle.write("group\tstatus\tdetail\n")
        for row in group_rows:
            handle.write(f"{row['group']}\t{row['status']}\t{row['detail']}\n")
    if any(state == "PASS" for state in statuses):
        return "PASS", f"rendered visual summaries for {sum(state == 'PASS' for state in statuses)} RepEnTools group(s)"
    if any(state == "FAIL" for state in statuses):
        return "FAIL", "RepEnTools reports exist but one or more groups lack parseable visual input"
    return "SKIP", "no validated RepEnTools report could be visualized"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methods-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    methods_dir = args.methods_dir.resolve()
    output_dir = (args.output_dir or methods_dir / "visuals").resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    adapter_status = read_status(methods_dir / "method_status.tsv")
    rows: list[tuple[str, str, str, str]] = []
    allo_state, allo_detail = allo_visuals(methods_dir, output_dir, adapter_status)
    rows.append(("allo", allo_state, allo_detail, "allo"))
    t3e_state, t3e_detail = table_visuals("t3e", methods_dir / "t3e" / "results", output_dir, adapter_status)
    rows.append(("t3e", t3e_state, t3e_detail, "t3e"))
    rep_status = read_status(methods_dir / "repentools_fastq" / "method_status.tsv")
    rep_state, rep_detail = repentools_visuals(methods_dir, output_dir, rep_status or adapter_status)
    rows.append(("repentools", rep_state, rep_detail, "repentools_fastq"))
    status_path = methods_dir / "visualization_status.tsv"
    with status_path.open("w") as handle:
        handle.write(STATUS_HEADER)
        for method, state, detail, rel in rows:
            output_path = output_dir / rel
            if output_path.exists():
                try:
                    outputs = str(output_path.relative_to(methods_dir))
                except ValueError:
                    outputs = str(output_path)
            else:
                outputs = ""
            handle.write(f"{method}\t{state}\t{detail}\t{outputs}\n")
    with (output_dir / "index.tsv").open("w") as handle:
        handle.write("method\tstatus\tdetail\toutputs\n")
        for method, state, detail, rel in rows:
            handle.write(f"{method}\t{state}\t{detail}\t{rel}\n")
    print(f"[te_visuals] status: {status_path}")
    states = [state for _, state, _, _ in rows]
    # Unsupported methods are a normal optional outcome (42 => SKIP).  A
    # method claiming PASS while lacking parseable output is a real contract
    # violation and should be visible as FAIL in module_status_latest.tsv.
    if any(state == "FAIL" for state in states):
        return 1
    if not any(state == "PASS" for state in states):
        return 42
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
