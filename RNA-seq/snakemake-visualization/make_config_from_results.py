#!/usr/bin/env python3
"""Create a Snakemake visualization config from a standard RNA-seq results dir."""

from __future__ import annotations

import argparse
import os
import re
from pathlib import Path


DEFAULT_GTF = {
    "hg38": "/path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf",
    "t2t": "/path/to/reference/human_t2t/chm13v2.0.annotation.gtf",
    "mm10": "/path/to/reference/mouse_mm10/gencode.vM25.primary_assembly.annotation.gtf",
    "mm39": "/path/to/reference/mouse_mm39/gencode.vM38.primary_assembly.annotation.gtf",
}


def yaml_scalar(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    text = str(value)
    if text == "":
        return '""'
    if re.fullmatch(r"[A-Za-z0-9_./:+-]+", text):
        return text
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def choose_plot_dir(results_dir: Path, requested: str | None) -> Path:
    if requested:
        path = Path(requested).expanduser()
        if not path.is_absolute():
            path = results_dir / path
        if not path.is_dir():
            raise SystemExit(f"plot dir not found: {path}")
        return path.resolve()

    for name in ("plots", "plot_with_replicates", "plot_4.9", "plot_with_replicates_4.28"):
        path = results_dir / name
        if path.is_dir() and list(path.rglob("*_matrix.csv")):
            return path.resolve()
    raise SystemExit(
        "Cannot find plot directory with *_matrix.csv. Tried: plots, "
        "plot_with_replicates, plot_4.9, plot_with_replicates_4.28"
    )


def choose_annotation(results_dir: Path, plot_dir: Path, requested: str | None) -> str:
    if requested:
        return str(Path(requested).expanduser().resolve())
    candidates = [
        plot_dir / "te_annotation.preview.csv",
        results_dir / "plots" / "te_annotation.preview.csv",
        results_dir / "plot_with_replicates" / "te_annotation.preview.csv",
        results_dir / "plot_4.9" / "te_annotation.preview.csv",
    ]
    for path in candidates:
        if path.is_file():
            return str(path.resolve())
    hits = sorted(results_dir.rglob("te_annotation.preview.csv"))
    return str(hits[0].resolve()) if hits else ""


def infer_matrix_type(name: str) -> str:
    if name.startswith("gene_"):
        return "gene"
    if name.startswith("TE_"):
        return "te"
    return "generic"


def infer_te_levels(name: str) -> tuple[str, str]:
    if name.endswith("_repName") or name.endswith("_subfamily"):
        return "repName", "repFamily"
    if name.endswith("_repFamily") or name.endswith("_family"):
        return "repFamily", "repClass"
    if name.endswith("_repClass") or name.endswith("_class"):
        return "repClass", "repClass"
    if name.endswith("_locus"):
        return "locus_id", "repFamily"
    return "repName", "repFamily"


def matrix_name(path: Path) -> str:
    return path.name.removesuffix("_matrix.csv")


def prefer_matrix_path(existing: Path, candidate: Path, name: str, plot_dir: Path) -> Path:
    existing_parent_match = existing.parent.name == name
    candidate_parent_match = candidate.parent.name == name
    if candidate_parent_match and not existing_parent_match:
        return candidate
    if existing_parent_match and not candidate_parent_match:
        return existing
    # Prefer the shallower file after parent-name preference, then stable lexical order.
    existing_depth = len(existing.relative_to(plot_dir).parts)
    candidate_depth = len(candidate.relative_to(plot_dir).parts)
    if candidate_depth < existing_depth:
        return candidate
    if candidate_depth == existing_depth and str(candidate) < str(existing):
        return candidate
    return existing


def collect_matrices(plot_dir: Path, include: str | None, exclude: str | None) -> dict[str, Path]:
    include_re = re.compile(include) if include else None
    exclude_re = re.compile(exclude) if exclude else None
    found: dict[str, Path] = {}
    for path in sorted(plot_dir.rglob("*_matrix.csv")):
        rel = path.relative_to(plot_dir)
        if "_panels" in rel.parts:
            continue
        name = matrix_name(path)
        if include_re and not include_re.search(name):
            continue
        if exclude_re and exclude_re.search(name):
            continue
        if name in found:
            found[name] = prefer_matrix_path(found[name], path, name, plot_dir)
        else:
            found[name] = path
    if not found:
        raise SystemExit(f"No *_matrix.csv found in {plot_dir}")
    return found


def build_config(args: argparse.Namespace) -> str:
    results_dir = Path(args.results_dir).expanduser().resolve()
    if not results_dir.is_dir():
        raise SystemExit(f"results dir not found: {results_dir}")

    sample_table = Path(args.sample_table).expanduser().resolve() if args.sample_table else results_dir / "condition.csv"
    contrast_file = Path(args.contrast_file).expanduser().resolve() if args.contrast_file else results_dir / "contrast.csv"
    if not sample_table.is_file():
        raise SystemExit(f"sample table not found: {sample_table}")
    if not contrast_file.is_file():
        raise SystemExit(f"contrast file not found: {contrast_file}")

    plot_dir = choose_plot_dir(results_dir, args.plot_dir)
    te_annotation = choose_annotation(results_dir, plot_dir, args.te_annotation_tsv)
    matrices = collect_matrices(plot_dir, args.include_matrix, args.exclude_matrix)
    tx2gene = args.tx2gene_path or DEFAULT_GTF.get(args.species, DEFAULT_GTF["hg38"])

    lines = [
        f"outdir: {yaml_scalar(args.outdir or str(results_dir / 'snakemake_visualization'))}",
        f"species: {yaml_scalar(args.species)}",
        f"rscript: {yaml_scalar(args.rscript)}",
        "",
        f"sample_table: {yaml_scalar(sample_table)}",
        f"contrast_file: {yaml_scalar(contrast_file)}",
        f"tx2gene_path: {yaml_scalar(tx2gene)}",
        f"te_annotation_tsv: {yaml_scalar(te_annotation)}",
        "",
        f"padj_cutoff: {yaml_scalar(args.padj_cutoff)}",
        f"lfc_cutoff: {yaml_scalar(args.lfc_cutoff)}",
        f"baseMean_min: {yaml_scalar(args.base_mean_min)}",
        f"label_top_n: {yaml_scalar(args.label_top_n)}",
        f"volcano_orientation: {yaml_scalar(args.volcano_orientation)}",
        f"gray_nonsig: {yaml_scalar(not args.keep_nonsig_color)}",
        f"diff_threads: {yaml_scalar(args.diff_threads)}",
        "",
        f"run_pathway: {yaml_scalar(not args.no_pathway)}",
        f"run_go: {yaml_scalar(not args.no_go)}",
        f"run_gsea: {yaml_scalar(not args.no_gsea)}",
        f"disable_gseaplot2: {yaml_scalar(args.disable_gseaplot2)}",
        f"run_te_analysis: {yaml_scalar(not args.no_te_analysis)}",
        "",
        "matrices:",
    ]
    for name, path in matrices.items():
        typ = infer_matrix_type(name)
        lines.append(f"  {name}:")
        lines.append(f"    type: {typ}")
        lines.append(f"    path: {yaml_scalar(path.resolve())}")
        if typ == "gene":
            lines.append(f"    tx2gene_path: {yaml_scalar(tx2gene)}")
        elif typ == "te":
            label_level, color_level = infer_te_levels(name)
            if te_annotation:
                lines.append(f"    te_annotation_tsv: {yaml_scalar(te_annotation)}")
            lines.append(f"    te_label_level: {yaml_scalar(label_level)}")
            lines.append(f"    te_color_level: {yaml_scalar(color_level)}")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate snakemake-visualization config from RNA-seq pipeline results."
    )
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--plot-dir", help="Default auto-detects plots/plot_with_replicates/plot_4.9.")
    parser.add_argument("--out", required=True, help="Output YAML path.")
    parser.add_argument("--outdir", help="Snakemake output directory.")
    parser.add_argument("--species", default="hg38", choices=sorted(DEFAULT_GTF))
    parser.add_argument("--rscript", default="conda run -n downstream Rscript")
    parser.add_argument("--sample-table")
    parser.add_argument("--contrast-file")
    parser.add_argument("--tx2gene-path")
    parser.add_argument("--te-annotation-tsv")
    parser.add_argument("--include-matrix", help="Regex over matrix names.")
    parser.add_argument("--exclude-matrix", help="Regex over matrix names.")
    parser.add_argument("--padj-cutoff", type=float, default=0.05)
    parser.add_argument("--lfc-cutoff", type=float, default=0.58)
    parser.add_argument("--base-mean-min", type=float, default=5)
    parser.add_argument("--label-top-n", type=int, default=40)
    parser.add_argument("--volcano-orientation", choices=("classic", "horizontal"), default="classic")
    parser.add_argument("--keep-nonsig-color", action="store_true", help="Do not force non-significant points to gray.")
    parser.add_argument("--diff-threads", type=int, default=1)
    parser.add_argument("--no-pathway", action="store_true")
    parser.add_argument("--no-go", action="store_true")
    parser.add_argument("--no-gsea", action="store_true")
    parser.add_argument("--disable-gseaplot2", action="store_true")
    parser.add_argument("--no-te-analysis", action="store_true")
    args = parser.parse_args()

    out = Path(args.out).expanduser()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(build_config(args), encoding="utf-8")
    print(f"[make_config_from_results] wrote {out}")


if __name__ == "__main__":
    main()
