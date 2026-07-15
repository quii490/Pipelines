#!/usr/bin/env python3
"""Build a self-contained HTML/PDF run report from pipeline outputs."""

from __future__ import annotations

import argparse
import csv
import html
import json
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path


def read_tsv(path: Path | None) -> list[dict[str, str]]:
    if not path or not path.is_file():
        return []
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def json_file(path: Path | None) -> dict[str, object]:
    if not path or not path.is_file():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def table(rows: list[dict[str, str]], limit: int = 200) -> str:
    if not rows:
        return "<p class='muted'>No records available.</p>"
    keys = list(rows[0].keys())
    header = "".join(f"<th>{html.escape(str(key))}</th>" for key in keys)
    body = []
    for row in rows[:limit]:
        cells = []
        for key in keys:
            value = str(row.get(key, ""))
            cls = ""
            if key == "status" or key.endswith("_status"):
                if value == "PASS":
                    cls = " class='pass'"
                elif value == "FAIL":
                    cls = " class='fail'"
                elif value.startswith("SKIP"):
                    cls = " class='skip'"
            cells.append(f"<td{cls}>{html.escape(value)}</td>")
        body.append("<tr>" + "".join(cells) + "</tr>")
    suffix = f"<p class='muted'>Showing {min(limit, len(rows))} of {len(rows)} rows.</p>" if len(rows) > limit else ""
    return f"<div class='table-wrap'><table><thead><tr>{header}</tr></thead><tbody>{''.join(body)}</tbody></table></div>{suffix}"


def link(path: Path, root: Path) -> str:
    try:
        rel = path.resolve().relative_to(root.resolve())
    except ValueError:
        rel = path.name
    return html.escape(str(rel))


def make_pdf(html_path: Path, pdf_path: Path) -> str:
    candidates = [
        ["weasyprint", str(html_path), str(pdf_path)],
        ["wkhtmltopdf", "--enable-local-file-access", str(html_path), str(pdf_path)],
        ["chromium", "--headless", "--no-sandbox", f"--print-to-pdf={pdf_path}", str(html_path)],
        ["chromium-browser", "--headless", "--no-sandbox", f"--print-to-pdf={pdf_path}", str(html_path)],
    ]
    for command in candidates:
        if not shutil.which(command[0]):
            continue
        try:
            proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False, timeout=180)
        except (OSError, subprocess.TimeoutExpired):
            continue
        if proc.returncode == 0 and pdf_path.is_file() and pdf_path.stat().st_size > 0:
            return f"PASS:{command[0]}"
    return "SKIP:no HTML-to-PDF converter installed"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    results = args.results_dir.resolve()
    downstream = results / "09_downstream"
    output = (args.output or downstream / "run_report.html").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    status_rows = read_tsv(downstream / "module_status_latest.tsv")
    if not status_rows:
        status_rows = read_tsv(downstream / "module_status.tsv")
    inventory_rows = read_tsv(downstream / "results_manifest.tsv")
    qc_rows = read_tsv(downstream / "qc_metrics.tsv")
    multimap_rows = read_tsv(downstream / "te_multimap_sensitivity.tsv")
    consensus_rows = read_tsv(downstream / "consensus_peaks" / "consensus_summary.tsv")
    te_method_rows = read_tsv(downstream / "te_methods" / "method_status.tsv")
    repentools_rows = read_tsv(downstream / "te_methods" / "repentools_fastq" / "method_status.tsv")
    te_visual_rows = read_tsv(downstream / "te_methods" / "visualization_status.tsv")
    run_manifest = json_file(downstream / "run_manifest.json")
    status_summary = json_file(downstream / "module_status_summary.json")
    inventory_summary = json_file(downstream / "run_summary.json")
    generated = datetime.now(timezone.utc).isoformat()

    counts: dict[str, int] = {}
    for row in status_rows:
        counts[row.get("status", "UNKNOWN")] = counts.get(row.get("status", "UNKNOWN"), 0) + 1
    cards = "".join(
        f"<div class='card'><div class='card-value'>{value}</div><div class='card-label'>{html.escape(key)}</div></div>"
        for key, value in sorted(counts.items())
    )
    if not cards:
        cards = "<div class='card'><div class='card-value'>0</div><div class='card-label'>modules</div></div>"

    # Do not recursively walk raw BAM/FASTQ/RepEnTools report trees on NFS.
    # Figure-producing modules are confined to these directories; limiting
    # discovery keeps report generation responsive for multi-GB runs.
    figure_roots = [
        downstream / "figures",
        downstream / "heatmaps",
        downstream / "correlation",
        downstream / "te_methods" / "visuals",
        downstream / "te_methods" / "repentools_fastq",
    ]
    figure_paths = sorted(
        p for root in figure_roots if root.is_dir()
        for p in root.rglob("*")
        if p.is_file() and p.suffix.lower() in {".pdf", ".png", ".svg", ".html"}
    )
    figures = "".join(
        f"<li><a href='{html.escape(link(path, results))}'>{html.escape(str(path.relative_to(results)))}</a></li>"
        for path in figure_paths[:1000]
    ) or "<li class='muted'>No figures found.</li>"
    manifest_summary = {
        "run_id": run_manifest.get("run_id", ""),
        "created_utc": run_manifest.get("created_utc", ""),
        "sample_count": run_manifest.get("sample_count", ""),
        "target_count": run_manifest.get("target_count", ""),
        "control_count": run_manifest.get("control_count", ""),
        "results_dir": str(results),
    }
    manifest_table = table([{k: str(v) for k, v in manifest_summary.items()}], limit=1)
    report_status = status_summary.get("status", "UNKNOWN")
    if inventory_summary.get("status") and inventory_summary.get("status") != "PASS":
        report_status = "INCOMPLETE"
    te_rows = te_visual_rows + te_method_rows + repentools_rows
    if any((row.get("status") or "").upper() == "FAIL" for row in te_rows):
        report_status = "INCOMPLETE"
    elif report_status in {"", "UNKNOWN"} and any((row.get("status") or "").upper() == "SKIP" for row in te_rows):
        report_status = "INCOMPLETE"
    css = """
    :root { color-scheme: light; font-family: Inter, Arial, sans-serif; }
    body { margin: 0; color: #1f2937; background: #f5f7fb; }
    header { background: linear-gradient(120deg,#172554,#2563eb); color: white; padding: 32px 6vw; }
    main { max-width: 1500px; margin: 24px auto; padding: 0 24px 64px; }
    h1 { margin: 0 0 8px; font-size: 30px; } h2 { margin-top: 34px; color: #172554; }
    .muted { color: #64748b; font-size: 13px; } .cards { display: flex; flex-wrap: wrap; gap: 12px; }
    .card { background: white; border-radius: 12px; padding: 16px 22px; min-width: 120px; box-shadow: 0 2px 12px #0f172a12; }
    .card-value { font-size: 26px; font-weight: 700; } .card-label { color: #64748b; font-size: 12px; text-transform: uppercase; }
    .pass { color: #047857; font-weight: 700; } .fail { color: #b91c1c; font-weight: 700; } .skip { color: #b45309; font-weight: 700; }
    .table-wrap { overflow: auto; background: white; border-radius: 10px; box-shadow: 0 2px 12px #0f172a12; }
    table { width: 100%; border-collapse: collapse; font-size: 12px; } th, td { padding: 8px 10px; border-bottom: 1px solid #e2e8f0; text-align: left; white-space: nowrap; }
    th { background: #eef2ff; position: sticky; top: 0; } a { color: #1d4ed8; } ul { line-height: 1.7; }
    .notice { padding: 12px 16px; border-left: 4px solid #2563eb; background: #eff6ff; border-radius: 6px; }
    """
    document = f"""<!doctype html>
<html lang='en'><head><meta charset='utf-8'><title>CUT&RUN/ChIP-seq run report</title><style>{css}</style></head>
<body><header><h1>CUT&amp;RUN / ChIP-seq + TE run report</h1><div>Generated {html.escape(generated)} · status: <strong>{html.escape(report_status)}</strong></div></header>
<main>
<div class='notice'>This report keeps hg38/mm39 core outputs separate from CHM13/RepEnTools outputs. SKIP means a method was not scientifically valid or was not available; it is not treated as a successful result.</div>
<h2>Module status</h2><div class='cards'>{cards}</div>
<h2>Run manifest</h2>{manifest_table}
<h2>Attempt-aware status table</h2>{table(status_rows)}
<h2>Sample inventory</h2>{table(inventory_rows)}
<h2>QC metrics</h2>{table(qc_rows)}
<h2>TE multi-mapping sensitivity</h2>{table(multimap_rows)}
<h2>Replicate consensus peaks</h2>{table(consensus_rows)}
<h2>TE method adapter status</h2>{table(te_method_rows)}{table(repentools_rows)}
<h2>TE method visualization status</h2>{table(te_visual_rows)}
<h2>Figures and tracks</h2><ul>{figures}</ul>
</main></body></html>"""
    output.write_text(document)
    pdf = output.with_suffix(".pdf")
    pdf_status = make_pdf(output, pdf)
    (output.with_suffix(".pdf.status.txt")).write_text(pdf_status + "\n")
    print(json.dumps({"status": "PASS", "html": str(output), "pdf": str(pdf), "pdf_status": pdf_status}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
