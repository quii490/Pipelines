#!/usr/bin/env python3
"""Organize raw ATAC downstream outputs into a user-facing v3 layout."""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


FIG_EXT = {".pdf", ".png"}
TABLE_EXT = {".csv", ".tsv", ".txt", ".xlsx"}


def mkdir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def copy_file(src: Path, dst: Path, records: list[tuple[str, str, str]]) -> None:
    if not src.exists() or not src.is_file():
        return
    mkdir(dst.parent)
    shutil.copy2(src, dst)
    records.append((str(dst), str(src), classify_dst(dst)))


def classify_dst(path: Path) -> str:
    parts = path.parts
    if "figures" in parts:
        return "figure"
    if "results" in parts:
        return "result"
    if "reports" in parts:
        return "report"
    return "other"


def strip_contrast_suffix(stem: str, contrast: str) -> str:
    return re.sub(rf"_{re.escape(contrast)}$", "", stem)


def te_distribution_subcategory(name: str) -> str:
    low = name.lower()
    if "status" in low:
        return "status"
    if "enrichment_vs_genome" in low:
        return "genome_enrichment"
    if "distribution" in low:
        return "composition"
    return "misc"


def figure_category(rel: Path, contrast: str) -> tuple[str, str]:
    name = rel.name
    stem = strip_contrast_suffix(rel.stem, contrast)
    low = "/".join(rel.parts).lower()

    if "annotation" in rel.parts:
        if "te_distribution" in low:
            direction = next((x for x in rel.parts if x in {"sig", "up", "down"}), "all")
            subcat = te_distribution_subcategory(name)
            return (f"differential/annotation/te_distribution/{subcat}", f"{contrast}.{direction}.{stem}{rel.suffix}")
        direction = next((x for x in rel.parts if x in {"sig", "up", "down"}), "all")
        return ("differential/annotation/gene_structure", f"{contrast}.{direction}.{stem}{rel.suffix}")

    if "TE" in rel.parts or name.startswith("TE_"):
        if "violin" in name.lower():
            return ("differential/te/violin", f"{contrast}.{stem}{rel.suffix}")
        if "class" in name.lower():
            return ("differential/te/class", f"{contrast}.{stem}{rel.suffix}")
        if "family" in name.lower() and "subfamily" not in name.lower():
            return ("differential/te/family", f"{contrast}.{stem}{rel.suffix}")
        if "subfamily" in name.lower():
            return ("differential/te/subfamily", f"{contrast}.{stem}{rel.suffix}")
        return ("differential/te/peak", f"{contrast}.{stem}{rel.suffix}")

    if "gene" in rel.parts or name.startswith("Gene_") or name.startswith("GO_") or name.startswith("GSEA_"):
        if "volcano" in name.lower():
            return ("differential/gene/volcano", f"{contrast}.{stem}{rel.suffix}")
        if "effect" in name.lower():
            return ("differential/gene/effect_size", f"{contrast}.{stem}{rel.suffix}")
        return ("differential/gene", f"{contrast}.{stem}{rel.suffix}")

    if name.startswith("MA_"):
        return ("differential/regions/MA", f"{contrast}.{stem}{rel.suffix}")
    if "volcano" in name.lower():
        return ("differential/regions/volcano", f"{contrast}.{stem}{rel.suffix}")
    if "effect" in name.lower() or name.startswith("Peak_class_"):
        return ("differential/regions/effect_size", f"{contrast}.{stem}{rel.suffix}")
    if "heatmap" in name.lower():
        return ("differential/regions/heatmap", f"{contrast}.{stem}{rel.suffix}")
    return ("misc", f"{contrast}.{stem}{rel.suffix}")


def result_category(src: Path, contrast: str) -> tuple[str, str]:
    name = src.name
    stem = src.stem
    if name in {"up.bed", "down.bed"}:
        return ("differential/beds", f"{contrast}.{name}")
    if name.startswith("differential_peaks_"):
        return ("differential/regions", f"{contrast}.differential_regions{src.suffix}")
    if name.startswith("peak_annotation_") or name.startswith("peak_to_gene_annotation"):
        return ("differential/annotation", f"{contrast}.{stem}{src.suffix}")
    if name.startswith("TE_"):
        direction = next((x for x in src.parts if x in {"sig", "up", "down"}), "all")
        if "TE_distribution" in src.parts or "TE_distribution" in str(src):
            subcat = te_distribution_subcategory(name)
            return (f"differential/te_distribution/{subcat}", f"{contrast}.{direction}.{stem}{src.suffix}")
        return ("differential/te", f"{contrast}.{direction}.{stem}{src.suffix}")
    if name.startswith("gene_") or name.startswith("GO_") or name.startswith("GSEA_"):
        return ("differential/gene", f"{contrast}.{stem}{src.suffix}")
    return ("misc", f"{contrast}.{name}")


def organize_contrast(raw: Path, out: Path, contrast_dir: Path, records: list[tuple[str, str, str]]) -> None:
    contrast = contrast_dir.name

    for src in contrast_dir.rglob("*"):
        if not src.is_file():
            continue
        rel = src.relative_to(contrast_dir)
        if src.suffix.lower() in FIG_EXT:
            cat, fname = figure_category(rel, contrast)
            copy_file(src, out / "figures" / cat / fname, records)
        elif src.suffix.lower() in TABLE_EXT or src.suffix.lower() == ".bed":
            cat, fname = result_category(src, contrast)
            copy_file(src, out / "results" / cat / fname, records)


def organize_top_level(raw: Path, out: Path, records: list[tuple[str, str, str]]) -> None:
    for src in (raw / "figures" / "qc").glob("*"):
        if src.is_file() and src.suffix.lower() in FIG_EXT:
            copy_file(src, out / "figures" / "qc" / src.name, records)
    for src in (raw / "figures" / "overview").glob("*"):
        if src.is_file() and src.suffix.lower() in FIG_EXT:
            copy_file(src, out / "figures" / "overview" / src.name, records)
    for src in (raw / "results" / "overview").glob("*"):
        if src.is_file() and src.suffix.lower() in TABLE_EXT:
            copy_file(src, out / "results" / "overview" / src.name, records)

    for cbed in (raw / "contrast_beds").glob("*"):
        if not cbed.is_dir():
            continue
        for src in cbed.glob("*.bed"):
            copy_file(src, out / "results" / "differential" / "beds" / f"{cbed.name}.{src.name}", records)


def write_reports(out: Path, raw: Path, records: list[tuple[str, str, str]]) -> None:
    reports = mkdir(out / "reports")
    index = reports / "downstream_file_index.tsv"
    with index.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["kind", "organized_path", "raw_path"])
        for dst, src, kind in sorted(records):
            writer.writerow([kind, dst, src])

    fig_count = sum(1 for _, _, kind in records if kind == "figure")
    res_count = sum(1 for _, _, kind in records if kind == "result")
    readme = reports / "README_downstream_outputs.md"
    readme.write_text(
        "\n".join(
            [
                "# ATAC 下游结果说明",
                "",
                "这里是面向使用者整理后的下游结果。R 原始输出保存在 `legacy/raw_r_output`，用于追溯。",
                "",
                "## 目录结构",
                "",
                "- `results/differential/regions`：peak/bin/region 层面的核心差异表。",
                "- `results/differential/annotation`：差异 region 的基因结构和 TE overlap 注释。",
                "- `results/differential/gene`：由差异 region 派生的基因聚合、GO 和 GSEA。",
                "- `results/differential/te`：由差异 region 派生的 TE class/family/subfamily 聚合。",
                "- `results/differential/te_distribution/{status,composition,genome_enrichment}`：TE overlap 状态、组成比例，以及相对全基因组 TE 背景的富集/缺失。",
                "- `results/differential/beds`：Up/Down region BED。",
                "- `results/overview`：跨 contrast 汇总表。",
                "- `figures/qc`：library size 和 logCPM 分布。",
                "- `figures/overview`：PCA、样本相关性和跨 contrast 总览。",
                "- `figures/differential/{regions,gene,te,annotation}`：各类差异图；其中 `annotation/te_distribution/{status,composition,genome_enrichment}` 已按 TE distribution 图类型拆分。",
                "- 全局 TSS/gene body/TE/L1 signal profile 位于上一级 `profile_heatmaps`。",
                "",
                f"R 原始输出：`{raw}`",
                f"整理后的图：{fig_count}",
                f"整理后的结果文件：{res_count}",
                "",
                "## differential、gene 和 TE 的关系",
                "",
                "`differential` 是父层，每行代表一个被检验或评分的 peak/bin/region。`gene` 和 `TE` 是由这些 region 派生出的解释层。一个 region 可能同时进入 gene 和 TE 汇总，也可能只进入其中一个或都不进入，所以不能把三者的数量相加。",
                "",
                "## TE 多重比对边界",
                "",
                "TE 差异结果来自标准 clean BAM 的 peak/bin counts 与 TE 注释 overlap，代表可定位的 TE-associated regulatory regions。`04_bw_te` relaxed track 会保留 MAPQ=0 reads，适合观察重复区域附近的信号；但它没有使用 EM/fractional allocation，不能当作严格的 TE family/subfamily 总量定量。",
                "",
                "原始文件到整理文件的对应关系见 `downstream_file_index.tsv`。",
                "",
            ]
        ),
        encoding="utf-8",
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-dir", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--clean", action="store_true", help="Remove organized results/figures/reports before rebuilding.")
    args = parser.parse_args()

    raw = args.raw_dir.resolve()
    out = args.outdir.resolve()
    if not raw.exists():
        raise SystemExit(f"raw dir not found: {raw}")

    if args.clean:
        for name in ("results", "figures", "reports", "contrasts", "contrast_beds", "differential"):
            path = out / name
            if path.exists():
                shutil.rmtree(path)

    records: list[tuple[str, str, str]] = []
    organize_top_level(raw, out, records)
    contrasts = raw / "contrasts"
    if contrasts.exists():
        for contrast_dir in sorted(p for p in contrasts.iterdir() if p.is_dir()):
            organize_contrast(raw, out, contrast_dir, records)
    write_reports(out, raw, records)
    print(f"[organize_atac_downstream_outputs] organized {len(records)} files into {out}")


if __name__ == "__main__":
    main()
