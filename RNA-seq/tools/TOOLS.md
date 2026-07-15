# RNA-seq 常用小工具

这些工具是主流程之外的补充入口，适合从已有 BAM、count matrix、DE 表或 bigWig 快速补分析和补图。它们不会改变 `run_auto_rnaseq.sh` / Nextflow 主流程。

结果汇总与维护：

```bash
rnaseq-report --plot-dir /path/to/results/plots
rnaseq-clean-work --results-dir /path/to/results
rnaseq-publish --results-dir /path/to/results --outdir /path/to/publish
rnaseq-rerun --results-dir /path/to/results --modules TE_REdiscoverTE,TE_Telescope
```

## 0. 已有 BAM 补做 QC metrics

只统计重复率及 RNA-seq 区域分布，不去重复，也不改变任何 count：

```bash
rnaseq-qc-metrics-bam \
  --bam-dir /path/to/02_gene_star \
  --gtf /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --strand reverse \
  --threads 4 \
  --outdir /path/to/qc_metrics
```

输出：

- `markdup/*.markdup.metrics.txt`：重复率 QC；原 BAM 不会被替换。
- `rnaseq_metrics/*.rnaseq.metrics.txt`：exonic、intronic、intergenic 比例及 5'/3' bias。
- `reference/generated.refFlat.txt`：由 GTF 自动生成，可在后续样本中复用。
- `logs/*.qc_metrics.log`：每个样本的完整 Picard 日志。

单个 BAM 用 `--bam FILE`；已有 refFlat 时可用 `--ref-flat FILE` 代替 `--gtf`。普通 reverse-stranded paired-end RNA-seq 使用 `--strand reverse`。

## 1. BAM 得到 counts

```bash
rnaseq-bam-to-counts \
  --bam-dir /path/to/bam \
  --gene-gtf /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --te-gtf /path/to/reference/TE_GTF/hg38_rmsk_TE.gtf \
  --outdir counts_from_bam \
  --prefix project \
  --layout PE \
  --strandedness 2 \
  --threads 16
```

输出 `project.gene.count_matrix.csv` 和 `project.TE.count_matrix.csv`。`--strandedness 0/1/2` 分别对应 unstranded/forward/reverse。若 TE 或 gene 注释是 SAF，用 `--gene-saf` / `--te-saf`。

## 2. counts 得到通用差异矩阵

只做差异分析，不画图、不做注释：

```bash
rnaseq-diff-counts \
  --matrix counts_from_bam/project.gene.count_matrix.csv \
  --sample-table condition.csv \
  --contrast-file contrast.csv \
  --outdir diff_only \
  --matrix-name gene_featureCounts
```

输出 `*.DE_matrix.csv`、`*.normalized_counts.csv`、`*.vst_or_log_matrix.csv`。

集合工具仍然保留，适合从 counts 一步得到 DE 矩阵、基础图和可选 pathway：

```bash
rnaseq-counts-to-de \
  --matrix counts_from_bam/project.gene.count_matrix.csv \
  --sample-table condition.csv \
  --contrast-file contrast.csv \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --outdir de_from_counts \
  --matrix-name gene_featureCounts \
  --run-go true \
  --run-gsea true
```

TE 多层级矩阵同样可用：

```bash
rnaseq-counts-to-de \
  --matrix plots/TE_TEtranscripts_repFamily_matrix.csv \
  --sample-table condition.csv \
  --contrast-file contrast.csv \
  --te-annotation-tsv plots/te_annotation.preview.csv \
  --te-label-level repFamily \
  --te-color-level repClass \
  --outdir de_from_counts \
  --matrix-name TE_TEtranscripts_repFamily
```

每个 contrast 会输出：

- `*.DE_matrix.csv`
- `*.normalized_counts.csv`
- `*.vst_or_log_matrix.csv`
- volcano / MA / top heatmap
- gene 模式可选 GO 和 GSEA Hallmark/KEGG/Reactome/GO

## 3. DE 矩阵注释、补图和通路分析

只注释，不画图：

```bash
rnaseq-annotate-de \
  --de-matrix diff_only/gene_featureCounts/gene_featureCounts_KO_vs_NC/gene_featureCounts_KO_vs_NC.DE_matrix.csv \
  --annotation-mode gene \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --outdir annotate/gene_KO_vs_NC
```

只做通路分析：

```bash
rnaseq-pathway-de \
  --de-matrix diff_only/gene_featureCounts/gene_featureCounts_KO_vs_NC/gene_featureCounts_KO_vs_NC.DE_matrix.csv \
  --species hg38 \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --case KO \
  --control NC \
  --outdir pathway/gene_KO_vs_NC
```

从 DE 矩阵补 volcano/MA 和组成图：

```bash
rnaseq-de-visuals \
  --de-matrix de_from_counts/gene_featureCounts/gene_featureCounts_KO_vs_NC/gene_featureCounts_KO_vs_NC.DE_matrix.csv \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --padj-cutoff 0.05 \
  --lfc-cutoff 0.58 \
  --label-top-n 30 \
  --volcano-orientation classic \
  --gray-nonsig true \
  --outdir de_visuals/gene_KO_vs_NC \
  --annotation-mode gene
```

TE 结果：

```bash
rnaseq-de-visuals \
  --de-matrix de_from_counts/TE_TEtranscripts_repName/TE_TEtranscripts_repName_KO_vs_NC/TE_TEtranscripts_repName_KO_vs_NC.DE_matrix.csv \
  --te-annotation-tsv plots/te_annotation.preview.csv \
  --annotation-mode te \
  --te-label-level repName \
  --te-color-level repFamily \
  --label-top-n 30 \
  --volcano-orientation classic \
  --gray-nonsig true \
  --outdir de_visuals/TE_repName_KO_vs_NC
```

TE 图包括 repClass/repFamily/repName 的上下调数量与比例，以及参考 CUT&RUN 逻辑的 TE age：`milliDiv < 50` 为 young，`50-150` 为 middle，`>=150` 为 old。

只做 TE 多层级和年龄分析：

```bash
rnaseq-te-analysis \
  --de-matrix diff_only/TE_TEtranscripts_repName/TE_TEtranscripts_repName_KO_vs_NC/TE_TEtranscripts_repName_KO_vs_NC.DE_matrix.csv \
  --te-annotation-tsv plots/te_annotation.preview.csv \
  --te-label-level repName \
  --te-color-level repFamily \
  --outdir te_analysis/TE_repName_KO_vs_NC
```

## 4. Snakemake 下游可视化流程

如果一个项目里有多个 gene/TE 矩阵和多个 contrast，推荐用 Snakemake 组织下游：

```bash
cd Pipeline/RNA-seq/snakemake-visualization
cp config.example.yaml my_project.yaml

snakemake \
  -s visualization.smk \
  --configfile my_project.yaml \
  --cores 8 \
  -p \
  --rerun-incomplete
```

说明见：

```text
snakemake-visualization/README.md
```

## 5. bigWig 相关性分析

```bash
rnaseq-bw-cor \
  --bw-dir /path/to/bw \
  --out-prefix qc/rnaseq_bw_cor \
  --binsize 10000 \
  --threads 8 \
  --methods pearson,spearman
```

输出 deepTools 的 `multiBigwigSummary` 矩阵、Pearson/Spearman heatmap、scatterplot 和 PCA，PDF/PNG 都会生成。

## 6. 两个样本相关性 scatter

```bash
rnaseq-two-sample-scatter \
  --matrix counts_from_bam/project.gene.count_matrix.csv \
  --sample-x KO_1 \
  --sample-y NC_1 \
  --out-prefix qc/KO_1_vs_NC_1 \
  --transform log2cpm \
  --label-top-n 20
```

输入可以是 gene count、TE count、归一化矩阵或其它第一列为 feature ID 的 CSV。输出 `*.scatter.pdf/png` 和用于复查的 `*.scatter_values.csv`。

## 7. 已有主流程下游

如果已经跑完上游，优先使用主入口：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream
```

只跑某个模块：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream \
  --only-tools TE_TEtranscripts
```

主下游会自动生成 gene/TE DE、PCA、Pearson、volcano、MA、heatmap、GO、GSEA，以及 TE_TElocal / TE_TEtranscripts / TEcount / REdiscoverTE / SalmonTE 的多层级结果。
