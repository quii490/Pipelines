# ATAC-seq 实战教程

| 状态 | 数据 | 最后审查 |
|---|---|---|
| Draft review | 使用公开/已授权 PE ATAC 数据 | 2026-07-16 |

选择两组、每组至少 2 个 biological replicates 的小型 PE 数据。本教程练习完整运行、从 clean BAM 重建、region heatmap 与 TE relaxed track。

## A. 初始化

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode init \
  --fastq-dir /path/to/tutorial/atac/fastq \
  --species hg38 \
  --outdir /path/to/tutorial/atac/results
```

准备 metadata：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

确认 samplesheet 与 metadata sample ID 完全一致。

## B. 正式运行

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode auto \
  --fastq-dir /path/to/tutorial/atac/fastq \
  --metadata-csv /path/to/tutorial/atac/metadata.csv \
  --contrast KO,WT \
  --species hg38 \
  --outdir /path/to/tutorial/atac/results \
  --preset standard \
  --resume
```

## C. 先读 QC

打开 `QC_REPORT.md`，检查 mapping、mitochondrial fraction、duplicates/library complexity、FRiP、fragment-size pattern、样本相关性/PCA。回答：哪个指标的分母是 read，哪个是 fragment？是否所有样本使用相同过滤阶段？

当前 TSS 图是 bigWig signal profile，不应在练习报告中写成严格 ENCODE TSS enrichment score。

## D. Peak 与 fixed-bin

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode downstream \
  --outdir /path/to/tutorial/atac/results \
  --metadata-csv /path/to/tutorial/atac/metadata.csv \
  --contrast KO,WT \
  --species hg38 \
  --levels both
```

比较 peak-level 与 bin-level：feature universe、稀疏度、显著 feature 数和生物学解释是否一致。不要只选择更显著的一套。

## E. 从 clean BAM 重建

```bash
bash ATAC-seq/scripts/run_atac_from_bam.sh \
  --bam-dir /path/to/tutorial/atac/results/02_align \
  --sample-meta /path/to/tutorial/atac/metadata.csv \
  --contrast-file /path/to/tutorial/atac/contrasts.csv \
  --species hg38 \
  --outdir /path/to/tutorial/atac/from_bam \
  --run-tss true \
  --run-nuc-phasing true \
  --cores 16
```

确认 clean BAM 定义与主流程一致。nucleosome phasing 只对 PE 解释。

## F. 自定义 regions

```bash
bash ATAC-seq/scripts/run_region_heatmap.sh \
  --regions /path/to/tutorial/atac/regions.bed \
  --bw-glob "/path/to/tutorial/atac/results/04_bw/*.bw" \
  --outdir /path/to/tutorial/atac/region_heatmap \
  --name tutorial_regions \
  --mode reference-point --reference-point center \
  --before 3000 --after 3000 --cores 8
```

检查 BED assembly、坐标、排序规则和 bigWig normalization。

## G. TE relaxed 练习

从 pre-clean sorted BAM 而不是 strict clean BAM 生成 relaxed track，再画 TE heatmap。比较每一步 clean counts。解释 relaxed track 为何可能有更高背景，以及为什么它不是 locus-specific significance。

## 教程验收

- [ ] 能解释 PE fragment 与 SE read 计数差别；
- [ ] 能区分 peaks、fixed bins 和 bigWig；
- [ ] 能说明 `skipZeros` 如何改变 heatmap；
- [ ] 能判断何时必须从 pre-clean BAM/FASTQ 开始；
- [ ] 能用多个 QC 指标而非单一 FRiP 决定是否调查样本。
