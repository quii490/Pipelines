#!/usr/bin/env python3
"""Build compact DEG/TE summaries and an offline HTML index from downstream outputs."""

import argparse
import csv
import html
import json
import math
import os
import re
from pathlib import Path
from statistics import median


def number(value):
    try:
        x = float(value)
        return x if math.isfinite(x) else None
    except (TypeError, ValueError):
        return None


def choose(fields, names):
    lower = {x.lower(): x for x in fields}
    for name in names:
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def summarize_de(path, root, padj_cutoff, lfc_cutoff, base_mean_min):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        rows = list(reader)
    lfc_col = choose(fields, ["log2FoldChange", "log2FC", "logFC"])
    padj_col = choose(fields, ["padj", "FDR", "adj.P.Val"])
    mean_col = choose(fields, ["baseMean", "mean_expression", "logCPM"])
    if not lfc_col:
        return None
    lfc_values, eligible, significant = [], [], []
    has_padj = False
    for row in rows:
        lfc = number(row.get(lfc_col))
        if lfc is None:
            continue
        lfc_values.append(lfc)
        mean = number(row.get(mean_col)) if mean_col else None
        keep = mean is None or mean >= base_mean_min
        if not keep:
            continue
        eligible.append(lfc)
        padj = number(row.get(padj_col)) if padj_col else None
        if padj is not None:
            has_padj = True
            if padj <= padj_cutoff and abs(lfc) >= lfc_cutoff:
                significant.append(lfc)
    selected = significant if has_padj else [x for x in eligible if abs(x) >= lfc_cutoff]
    rel = path.relative_to(root)
    parts = rel.parts
    module = parts[0] if len(parts) > 1 else "unknown"
    name = path.name[:-7] if path.name.endswith("_DE.csv") else path.stem
    contrast = name[len(module) + 1 :] if name.startswith(module + "_") else name
    return {
        "module": module,
        "contrast": contrast,
        "analysis_mode": "deseq2" if has_padj else "exploratory",
        "tested_features": len(lfc_values),
        "eligible_features": len(eligible),
        "up": sum(x >= lfc_cutoff for x in selected),
        "down": sum(x <= -lfc_cutoff for x in selected),
        "median_lfc": f"{median(eligible):.6g}" if eligible else "",
        "positive_fraction": f"{sum(x > 0 for x in eligible) / len(eligible):.6g}" if eligible else "",
        "padj_cutoff": padj_cutoff if has_padj else "NA",
        "lfc_cutoff": lfc_cutoff,
        "base_mean_min": base_mean_min,
        "de_matrix": str(rel),
    }


def read_lfc(path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        id_col = choose(fields, ["feature_id", "gene_id", "gene_name", "Geneid", fields[0] if fields else ""])
        lfc_col = choose(fields, ["log2FoldChange", "log2FC", "logFC"])
        if not id_col or not lfc_col:
            return {}
        out = {}
        for row in reader:
            feature = str(row.get(id_col, "")).strip()
            if re.match(r"^ENS[A-Z0-9]+\.\d+$", feature):
                feature = re.sub(r"\.\d+$", "", feature)
            value = number(row.get(lfc_col))
            if feature and value is not None:
                out[feature] = value
        return out


def pearson(x, y):
    if len(x) < 3:
        return None
    mx, my = sum(x) / len(x), sum(y) / len(y)
    xx = sum((v - mx) ** 2 for v in x)
    yy = sum((v - my) ** 2 for v in y)
    if xx == 0 or yy == 0:
        return None
    return sum((a - mx) * (b - my) for a, b in zip(x, y)) / math.sqrt(xx * yy)


def module_level(module):
    for level in ["repClass", "repFamily", "repName"]:
        if module.endswith("_" + level):
            return level
    return "gene" if module.startswith("gene_") else None


def concordance_rows(de_files, summaries, root):
    summary_by_path = {x["de_matrix"]: x for x in summaries}
    grouped = {}
    for path in de_files:
        rel = str(path.relative_to(root))
        summary = summary_by_path.get(rel)
        if not summary:
            continue
        level = module_level(summary["module"])
        if not level:
            continue
        grouped.setdefault((summary["contrast"], level), []).append((summary["module"], path))
    output = []
    for (contrast, level), entries in sorted(grouped.items()):
        cache = {module: read_lfc(path) for module, path in entries}
        modules = sorted(cache)
        for i, left in enumerate(modules):
            for right in modules[i + 1 :]:
                common = sorted(set(cache[left]) & set(cache[right]))
                value = pearson([cache[left][x] for x in common], [cache[right][x] for x in common])
                output.append({
                    "contrast": contrast, "feature_level": level, "module_a": left, "module_b": right,
                    "matched_features": len(common), "pearson_lfc": f"{value:.6g}" if value is not None else "NA"
                })
    return output


def write_csv(path, rows, fields):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_status(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def build_html(root, summaries, statuses, pdfs, tables):
    modules = sorted({x["module"] for x in summaries} | {p.relative_to(root).parts[0] for p in pdfs if p.relative_to(root).parts})
    status_rows = "".join(
        f"<tr><td>{html.escape(x.get('module',''))}</td><td><span class='status {html.escape(x.get('status',''))}'>{html.escape(x.get('status',''))}</span></td><td>{html.escape(x.get('reason',''))}</td></tr>"
        for x in statuses[-200:]
    ) or "<tr><td colspan='3'>No status records</td></tr>"
    deg_rows = "".join(
        "<tr>" + "".join(f"<td>{html.escape(str(x[k]))}</td>" for k in ["module", "contrast", "analysis_mode", "tested_features", "up", "down", "median_lfc", "positive_fraction"]) +
        f"<td><a href='{html.escape(x['de_matrix'])}'>CSV</a></td></tr>" for x in summaries
    ) or "<tr><td colspan='9'>No DE matrices found</td></tr>"
    cards = []
    for module in modules:
        items = [p for p in pdfs if p.relative_to(root).parts[0] == module]
        links = "".join(f"<li><a href='{html.escape(str(p.relative_to(root)))}'>{html.escape(p.name)}</a></li>" for p in items)
        cards.append(f"<section class='module' data-module='{html.escape(module.lower())}'><h2>{html.escape(module)} <small>{len(items)} PDFs</small></h2><ul>{links or '<li>No PDF</li>'}</ul></section>")
    payload = {
        "root": str(root), "modules": len(modules), "pdfs": len(pdfs), "tables": len(tables), "comparisons": len(summaries)
    }
    return f"""<!doctype html><html lang='zh-CN'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>
<title>RNA-seq results</title><style>
:root{{--ink:#25313b;--muted:#66727c;--line:#d9dfe3;--bg:#f6f7f7;--paper:#fff;--accent:#3c6e71;--ok:#4f7b63;--warn:#a2763b;--bad:#a45454}}*{{box-sizing:border-box}}body{{margin:0;font:14px/1.5 Arial,sans-serif;color:var(--ink);background:var(--bg)}}header{{padding:28px 4vw 20px;background:var(--paper);border-bottom:1px solid var(--line)}}h1{{margin:0 0 8px;font-size:28px}}main{{padding:20px 4vw 50px}}.metrics{{display:flex;gap:24px;flex-wrap:wrap;color:var(--muted)}}.toolbar{{position:sticky;top:0;background:var(--bg);padding:12px 0;z-index:2}}input{{width:min(520px,100%);padding:10px;border:1px solid var(--line);border-radius:4px}}table{{width:100%;border-collapse:collapse;background:var(--paper);margin:10px 0 24px}}th,td{{padding:8px;border-bottom:1px solid var(--line);text-align:left}}th{{background:#eef1f1}}.module{{background:var(--paper);border:1px solid var(--line);padding:14px 18px;margin:12px 0;border-radius:6px}}.module h2{{font-size:18px;margin:0 0 8px}}small{{color:var(--muted);font-weight:normal}}ul{{columns:3;column-gap:28px;margin:0;padding-left:18px}}li{{break-inside:avoid;margin:4px 0}}a{{color:var(--accent)}}.status{{font-weight:bold}}.success{{color:var(--ok)}}.failed,.incomplete{{color:var(--bad)}}.skipped{{color:var(--warn)}}@media(max-width:900px){{ul{{columns:1}}.wide{{overflow-x:auto}}}}
</style></head><body><header><h1>RNA-seq results</h1><div class='metrics'><span>Modules: {payload['modules']}</span><span>Comparisons: {payload['comparisons']}</span><span>PDF: {payload['pdfs']}</span><span>Tables: {payload['tables']}</span></div></header><main>
<h2>Module status</h2><div class='wide'><table><tr><th>Module/step</th><th>Status</th><th>Reason</th></tr>{status_rows}</table></div>
<h2>DE/TE summary</h2><div class='wide'><table><tr><th>Module</th><th>Contrast</th><th>Mode</th><th>Tested</th><th>Up</th><th>Down</th><th>Median LFC</th><th>Positive fraction</th><th>Matrix</th></tr>{deg_rows}</table></div>
<div class='toolbar'><input id='q' placeholder='Filter module or filename'></div>{''.join(cards)}
</main><script>const q=document.getElementById('q');q.addEventListener('input',()=>{{const s=q.value.toLowerCase();document.querySelectorAll('.module').forEach(x=>x.hidden=!x.textContent.toLowerCase().includes(s));}});</script></body></html>"""


def write_module_indexes(root, pdfs):
    modules = sorted({p.relative_to(root).parts[0] for p in pdfs if len(p.relative_to(root).parts) > 1})
    for module in modules:
        module_dir = root / module
        links = []
        for path in pdfs:
            rel = path.relative_to(module_dir) if module_dir in path.parents else None
            if rel is not None:
                links.append(f"<li><a href='{html.escape(str(rel))}'>{html.escape(str(rel))}</a></li>")
        page = f"<!doctype html><html lang='zh-CN'><meta charset='utf-8'><title>{html.escape(module)}</title><style>body{{font:14px Arial;max-width:1100px;margin:30px auto;padding:0 20px}}li{{margin:6px 0}}a{{color:#3c6e71}}</style><h1>{html.escape(module)}</h1><p><a href='../index.html'>返回总报告</a></p><ul>{''.join(links)}</ul></html>"
        (module_dir / "index.html").write_text(page, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plot-dir", required=True)
    ap.add_argument("--padj-cutoff", type=float, default=0.05)
    ap.add_argument("--lfc-cutoff", type=float, default=0.58)
    ap.add_argument("--base-mean-min", type=float, default=5.0)
    args = ap.parse_args()
    root = Path(args.plot_dir).resolve()
    summary_dir = root / "_summary"
    de_files = sorted(root.glob("*/de_matrices/*_DE.csv"))
    summaries = [x for x in (summarize_de(p, root, args.padj_cutoff, args.lfc_cutoff, args.base_mean_min) for p in de_files) if x]
    fields = ["module", "contrast", "analysis_mode", "tested_features", "eligible_features", "up", "down", "median_lfc", "positive_fraction", "padj_cutoff", "lfc_cutoff", "base_mean_min", "de_matrix"]
    write_csv(summary_dir / "deg_summary.csv", summaries, fields)
    write_csv(summary_dir / "lfc_symmetry_summary.csv", summaries, ["module", "contrast", "analysis_mode", "eligible_features", "median_lfc", "positive_fraction", "de_matrix"])
    concordance = concordance_rows(de_files, summaries, root)
    write_csv(summary_dir / "lfc_method_concordance.csv", concordance, ["contrast", "feature_level", "module_a", "module_b", "matched_features", "pearson_lfc"])
    statuses = read_status(summary_dir / "module_status.csv")
    pdfs = sorted(root.rglob("*.pdf"))
    tables = sorted([*root.rglob("*.csv"), *root.rglob("*.tsv"), *root.rglob("*.txt")])
    (root / "index.html").write_text(build_html(root, summaries, statuses, pdfs, tables), encoding="utf-8")
    write_module_indexes(root, pdfs)
    (summary_dir / "report_manifest.json").write_text(json.dumps({"pdf": len(pdfs), "tables": len(tables), "de_matrices": len(summaries)}, indent=2), encoding="utf-8")
    print(f"[rnaseq-report] index={root / 'index.html'} comparisons={len(summaries)} pdf={len(pdfs)}")


if __name__ == "__main__":
    main()
