#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path


def read_tsv(path):
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value):
    try:
        if value in (None, "", "NA", "nan", "NaN"):
            return None
        x = float(value)
        if math.isnan(x):
            return None
        return x
    except Exception:
        return None


def fmt(value, digits=3):
    x = as_float(value)
    if x is None:
        return "NA"
    if abs(x) >= 1000:
        return f"{x:,.0f}"
    return f"{x:.{digits}f}"


def status_counts(rows):
    counts = {"PASS": 0, "WARN": 0, "FAIL": 0, "NA": 0}
    for row in rows:
        counts[row.get("overall_qc", "NA")] = counts.get(row.get("overall_qc", "NA"), 0) + 1
    return counts


def write_markdown(outdir, species, qc_rows, pass_rows):
    report = outdir / "QC_REPORT.md"
    counts = status_counts(pass_rows)
    lines = [
        "# ATAC-seq QC 中文报告",
        "",
        f"- 物种/基因组：`{species}`",
        f"- 结果目录：`{outdir}`",
        "",
        "## 1. 总览",
        "",
        f"- 样本数：{len(qc_rows)}",
        f"- PASS：{counts.get('PASS', 0)}",
        f"- WARN：{counts.get('WARN', 0)}",
        f"- FAIL：{counts.get('FAIL', 0)}",
        f"- NA：{counts.get('NA', 0)}",
        "",
        "优先查看：",
        "",
        "- `03_qc/qc_pass_fail.tsv`：每个样本的判读表。",
        "- `03_qc/atac_qc_summary_with_metrics.tsv`：完整 QC 数值。",
        "- `03_qc/FRiP_barplot.pdf/png`：peak 信噪比。",
        "- `03_qc/Mito_fraction_barplot.pdf/png`：线粒体比例。",
        "- `03_qc/Library_complexity_barplot.pdf/png`：文库复杂度。",
        "- `03_qc/Fragment_class_fraction_barplot.pdf/png`：PE fragment 分类。",
        "",
        "## 2. 样本级判读",
        "",
        "| sample | layout | overall | FRiP | mt_fraction | NRF | PBC1 | PBC2 |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    pass_by_sample = {row.get("sample"): row for row in pass_rows}
    for row in qc_rows:
        sample = row.get("sample", "")
        prow = pass_by_sample.get(sample, {})
        lines.append(
            "| {sample} | {layout} | {overall} | {frip} | {mt} | {nrf} | {pbc1} | {pbc2} |".format(
                sample=sample,
                layout=row.get("layout", prow.get("layout", "NA")),
                overall=prow.get("overall_qc", "NA"),
                frip=fmt(row.get("frip")),
                mt=fmt(row.get("mt_fraction")),
                nrf=fmt(row.get("NRF")),
                pbc1=fmt(row.get("PBC1")),
                pbc2=fmt(row.get("PBC2")),
            )
        )

    lines.extend(
        [
            "",
            "## 3. 指标解释",
            "",
            "- **FRiP**：落在 peak 中的 reads/fragments 比例。PE 按 fragment 统计，SE 按 read 统计。通常 `>0.2` 可接受，`>0.3` 更理想，但强依赖细胞类型和实验设计。",
            "- **mt_fraction**：线粒体比例。PE 默认使用 `mitochondrial_fragments / (mitochondrial_fragments + nuclear_fragments_before_blacklist)`。过高常提示细胞状态、裂解或建库问题。",
            "- **NRF**：Non-Redundant Fraction，唯一 fragment/read position 数量除以总 fragment/read position 数量。越接近 1，文库复杂度越好。",
            "- **PBC1**：只出现一次的位置数除以所有非重复位置数。越高越好。",
            "- **PBC2**：只出现一次的位置数除以出现两次的位置数。越高越好，低值提示 bottleneck 或 PCR duplicate 问题。",
            "- **Fragment classes**：PE ATAC 中常看 NFR、mono-nucleosome、di-nucleosome 和更长 fragment 的比例。好的 ATAC 通常有短片段开放区信号，同时保留核小体周期性。",
            "",
            "## 4. 重要概念",
            "",
            "### Reproducible peak 是什么",
            "",
            "单个样本 call 出来的 peak 可能包含随机噪声。**reproducible peak** 指在同一 condition 的多个 biological replicates 中重复出现、或者通过 IDR/pseudo-replicate 分析支持的 peak。它比简单把所有样本 peak 取并集更严格，更适合做正式 count matrix 和差异分析。",
            "",
            "现有 `consensus_peaks.bed` 是所有样本 summit 扩展后合并得到的共同计数区域。它适合作为通用分析骨架，但如果有生物学重复，后续应升级为 `condition reproducible peaks -> union`，减少噪声 peak 进入差异分析。",
            "",
            "### TSS profile 和 TSS enrichment score 的区别",
            "",
            "当前 deeTools 图展示的是 TSS 周围 bigWig signal profile，能看 promoter 附近是否有平均开放峰。严格 TSS enrichment score 需要 Tn5-shifted insertion sites，并用远端背景归一化得到每个样本一个分数。后续可以把它作为独立 QC 模块加入。",
            "",
            "### TE branch 应如何理解",
            "",
            "标准 ATAC 分支使用 MAPQ 过滤和 clean BAM，适合 high-confidence peak、FRiP 和 gene/TE overlap。TE 区域有大量 multi-mapping reads，因此 TE-aware branch 需要单独处理：",
            "",
            "- `unique-mappable TE loci`：只解释能唯一定位的 TE locus，可信但不完整。",
            "- `fractional TE subfamily counts`：用 multi-mapper fractional counting 或 EM 分配，适合看 TE family/subfamily 总体变化。",
            "- `relaxed TE bigWig`：用于可视化探索，不能当作严格 locus-level 定量。",
            "",
            "## 5. 目前仍建议补充的 gold-standard 模块",
            "",
            "- Tn5-shifted TSS enrichment score。",
            "- replicate-supported reproducible peaks 或 IDR。",
            "- sample peak、consensus peak、reproducible peak 三套 FRiP。",
            "- fragment CPM、fragment RPGC、Tn5 insertion CPM 三类 bigWig。",
            "- TE-aware fractional/EM subfamily count branch。",
            "",
            "## 6. 从中间结果补跑或改参数",
            "",
            "这些命令不重新跑 FASTQ/比对，适合只改图形参数、TE filter、contrast 或重新注释：",
            "",
            "```bash",
            f"# 只重跑 downstream 差异分析、基因注释、TE volcano/MA/violin",
            f"atacseq downstream --result-dir {outdir} --species {species} \\",
            f"  --contrast-file {outdir}/_automation/inputs/contrasts.csv \\",
            f"  --outdir {outdir}/08_downstream_manual",
            "",
            "# 对任意 BED 区域画 ATAC heatmap/profile",
            f"atacseq heatmap --regions target_regions.bed \\",
            f"  --bw-glob '{outdir}/04_bw/*.bw' \\",
            f"  --outdir {outdir}/manual_heatmap --name target_regions \\",
            "  --before 3000 --after 3000 --threads 8",
            "",
            "# 重新画 TE heatmap，并可按 TE class/family/name 过滤",
            f"atacseq te --bw-glob '{outdir}/04_bw/*.bw' \\",
            f"  --global-peak-bed {outdir}/06_consensus_peaks/consensus_peaks.bed \\",
            f"  --outdir {outdir}/manual_te_heatmap \\",
            f"  --species {species} --te-class-filter LINE --cores 8",
            "",
            "# 重新生成 TSS accessibility profile",
            f"atacseq tss --bw-glob '{outdir}/04_bw/*.bw' \\",
            f"  --species {species} --outdir {outdir}/manual_tss --cores 8",
            "",
            "# 从已有 clean BAM 重新 call peak 或改 MACS3 参数",
            f"atacseq callpeak --bam-dir {outdir}/02_align \\",
            f"  --species {species} --outdir {outdir}/manual_callpeak --qvalue 0.01",
            "```",
            "",
        ]
    )
    report.write_text("\n".join(lines), encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description="Generate a Chinese ATAC-seq QC report from pipeline outputs.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--species", default="unknown")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    qc_dir = outdir / "03_qc"
    qc_rows = read_tsv(qc_dir / "atac_qc_summary_with_metrics.tsv")
    if not qc_rows:
        qc_rows = read_tsv(qc_dir / "atac_qc_summary.tsv")
    pass_rows = read_tsv(qc_dir / "qc_pass_fail.tsv")
    report = write_markdown(outdir, args.species, qc_rows, pass_rows)
    print(f"[generate_qc_report] wrote {report}")


if __name__ == "__main__":
    main()
