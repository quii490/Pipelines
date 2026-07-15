# ATAC-seq 区域、注释与 motif 工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## 任意 BED 的 heatmap/profile

```bash
bash ATAC-seq/scripts/run_region_heatmap.sh \
  --regions /path/to/regions.bed \
  --bw-glob "/path/to/04_bw/*.bw" \
  --outdir /path/to/region_heatmap \
  --name target_regions \
  --mode reference-point \
  --reference-point center \
  --before 3000 \
  --after 3000 \
  --cores 8
```

用于 promoter、enhancer、motif sites、differential peaks 或自定义 TE subset。若 regions 已按某个样本信号筛选并排序，再展示同一样本会产生选择偏倚；应保留 anchor 和排序规则。

## 差异 peak heatmap

```bash
bash ATAC-seq/scripts/run_diff_peak_heatmap.sh \
  --contrast-bed-dir /path/to/08_downstream/peak_level/contrast_beds \
  --bw-glob "/path/to/04_bw/*.bw" \
  --outdir /path/to/KO_vs_WT_heatmap \
  --before 3000 --after 3000 \
  --cores 8
```

该图展示已选 differential regions 的信号，不重新计算差异显著性。行选择应基于预先声明的 padj/logFC/top-N 规则。

## Peak annotation 与 overlap

```bash
bash ATAC-seq/scripts/run_peak_annotation.sh \
  --input /path/to/peaks.bed \
  --species hg38 \
  --out-prefix /path/to/annotation/sample

bash ATAC-seq/scripts/run_peak_overlap.sh \
  --peaks "A=/path/to/A.bed B=/path/to/B.bed" \
  --outdir /path/to/overlap \
  --distance 0
```

Peak annotation 还可设置 `--gtf`、`--te-bed`、`--promoter-up`、`--promoter-down`、`--txdb-cache` 和 `--no-plots`。最近基因注释不等于该基因是功能靶点。

Overlap 也可用 `--peak-dir /path/to/peaks` 自动扫描 BED/narrowPeak/broadPeak；`--distance` 会先合并相距不超过该距离的 intervals，`--top-n` 只保留每个输入前 N 行。两者都会改变结果，必须记录。Overlap 受 peak width、caller、threshold 和 universe 影响；同时报告交集数量、各自比例和采用的 overlap 定义。

## HOMER motif

```bash
bash ATAC-seq/scripts/run_motif_homer.sh \
  --bed-dir /path/to/contrast_beds \
  --genome hg38 \
  --outdir /path/to/homer \
  --size 200 \
  --run-annotation true
```

motif enrichment 是序列富集证据，不证明该 TF 在本实验中直接结合。background regions 的长度、GC 和可及性会影响结果。

`--bed-dir` 支持组织后的 `<contrast>.up.bed/.down.bed` 和 legacy `<contrast>/up.bed` 结构。`--size 200` 表示以 peak center 为中心的窗口；`--run-annotation false` 可跳过 `annotatePeaks.pl`，但不会改变 motif discovery 输入。

## TOBIAS footprinting

```bash
bash ATAC-seq/scripts/run_tobias.sh \
  --bam-dir /path/to/clean_bam \
  --peaks /path/to/consensus_peaks.bed \
  --sample-meta /path/to/metadata.csv \
  --contrast-file /path/to/contrasts.csv \
  --genome /path/to/genome.fa \
  --motif-meme /path/to/motifs.meme \
  --outdir /path/to/tobias \
  --cores 8
```

Footprinting 对 Tn5 bias、深度、motif quality 和 replicate 设计敏感。当前 wrapper 的 BINDetect 会为每个 condition 选择 metadata 中第一个样本，并写入 `bindetect_sample_choice.tsv`；解释前必须审查该表。它不是 replicate-aware differential model，单样本 footprint score 不应替代重复设计分析。
