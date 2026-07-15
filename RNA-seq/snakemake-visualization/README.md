# RNA-seq Snakemake 下游可视化流程

这个 Snakemake 流程只负责下游分析和可视化，不改变 FASTQ 到 BAM/counts 的主 Nextflow 流程。它适合从已有 count matrix 批量跑：

```text
counts -> differential matrix -> annotation -> volcano/MA/composition plots
       -> gene pathway analysis or TE class/family/age analysis
```

如果你手里是 BAM，先用独立工具得到 counts：

```bash
rnaseq-bam-to-counts \
  --bam-dir /path/to/bam \
  --gene-gtf /path/to/gencode.gtf \
  --te-gtf /path/to/rmsk_TE.gtf \
  --outdir counts_from_bam \
  --prefix project \
  --layout PE \
  --strandedness 2 \
  --threads 16
```

## 1. 从主流程结果目录直接运行

主流程和传统下游会在结果目录中生成标准文件，例如：

```text
results/
  condition.csv
  contrast.csv
  plots/
    gene_featureCounts/gene_featureCounts_matrix.csv
    gene_Salmon/gene_Salmon_matrix.csv
    TE_TEcount_class/TE_TEcount_class_matrix.csv
    TE_TEcount_family/TE_TEcount_family_matrix.csv
    TE_TEcount_subfamily/TE_TEcount_subfamily_matrix.csv
    TE_TElocal_repName/TE_TElocal_repName_matrix.csv
    TE_TEtranscripts_repName/TE_TEtranscripts_repName_matrix.csv
    te_annotation.preview.csv
```

推荐直接从 `results-dir` 生成配置：

```bash
cd Pipeline/RNA-seq/snakemake-visualization

bash run_snakemake_visualization.sh \
  --results-dir /path/to/lab-data/RNAseq/20260402_C174nSTE_clone/SEQ2603854BJ/results_4.7 \
  --plot-dir plots \
  --outdir /path/to/lab-data/RNAseq/20260402_C174nSTE_clone/SEQ2603854BJ/results_4.7/snakemake_visualization \
  --cores 8 \
  --dry-run
```

如果不指定 `--plot-dir`，会按顺序寻找 `plots`、`plot_with_replicates`、`plot_4.9`、`plot_with_replicates_4.28`。自动配置工具会识别 `*_matrix.csv`，按名称推断 gene/TE 类型，并为 `repName`、`repFamily`、`repClass`、`locus` 等 TE 层次设置默认注释参数。

可以用正则限制矩阵：

```bash
bash run_snakemake_visualization.sh \
  --results-dir /path/to/results \
  --include-matrix '^(gene_featureCounts|TE_TEcount_|TE_TEtranscripts_)' \
  --cores 8
```

## 2. 手动准备配置

复制模板：

```bash
cd Pipeline/RNA-seq/snakemake-visualization
cp config.example.yaml my_project.yaml
```

编辑 `my_project.yaml`：

- `sample_table`: `condition.csv`
- `contrast_file`: `contrast.csv`
- `rscript`: 可选，指定运行 R 脚本的命令
- `tx2gene_path`: gene GTF 或 tx2gene
- `te_annotation_tsv`: TE 注释 TSV，推荐主流程输出的 `te_annotation.preview.csv`
- `matrices`: 一个或多个 count matrix
- `padj_cutoff` / `lfc_cutoff` / `baseMean_min`: volcano、MA、上下调统计共用阈值
- `label_top_n`: volcano/MA 标注 feature 数；设为 `0` 不标注
- `volcano_orientation`: `classic` 为默认经典竖版；`horizontal` 为旧版横放
- `gray_nonsig`: `true` 时非显著点统一灰色；`false` 时保留分组颜色

矩阵类型：

```yaml
matrices:
  gene_featureCounts:
    type: gene
    path: /path/gene_featureCounts_matrix.csv

  TE_TEtranscripts_repName:
    type: te
    path: /path/TE_TEtranscripts_repName_matrix.csv
    te_label_level: repName
    te_color_level: repFamily
```

如果 Snakemake 和 R 包不在同一个环境，可以在配置里指定 R 入口。例如当前 pod 里 `rnaseq` 环境有 Snakemake，`downstream` 环境有 DESeq2/ggplot2/pheatmap/readr：

```yaml
rscript: conda run -n downstream Rscript
```

## 3. 运行

先 dry-run：

```bash
snakemake \
  -s visualization.smk \
  --configfile my_project.yaml \
  -n -p
```

正式运行：

```bash
snakemake \
  -s visualization.smk \
  --configfile my_project.yaml \
  --cores 8 \
  -p \
  --rerun-incomplete
```

也可以用包装脚本：

```bash
bash run_snakemake_visualization.sh --config my_project.yaml --cores 8
```

只重跑某个目标也可以直接指定输出文件，例如：

```bash
snakemake \
  -s visualization.smk \
  --configfile my_project.yaml \
  --cores 4 \
  rnaseq_snakemake_downstream/visuals/gene_featureCounts/gene_featureCounts_KO_vs_NC/gene_featureCounts_KO_vs_NC.visuals.summary.txt
```

## 4. 输出结构

```text
<outdir>/
  diff/<matrix>/<matrix>_<case>_vs_<control>/
    *.DE_matrix.csv
    *.normalized_counts.csv
    *.vst_or_log_matrix.csv

  annotate/<matrix>/<contrast>/
    *.annotated_DE_matrix.csv
    *.direction_summary.csv

  visuals/<matrix>/<contrast>/
    *.volcano.pdf/png
    *.MA.pdf/png
    *.up_down_bar.pdf/png
    *.annotated_DE_matrix.csv

  pathway/<gene_matrix>/<contrast>/
    *_GO_up.csv
    *_GO_down.csv
    *_GSEA/...

  te/<te_matrix>/<contrast>/
    *.TE_annotated_DE_matrix.csv
    *.repClass.direction_counts.csv
    *.repFamily.direction_counts.csv
    *.TE_age.direction_counts.csv
    *.TE_family_class_mean_log2FC_heatmap.pdf
```

## 4. 为什么保留独立工具

Snakemake 适合项目级批量流程；独立工具适合临时补图、单个矩阵排查、或在 notebook/服务器交互环境中快速调用。因此当前设计同时保留：

- 拆分工具：`rnaseq-diff-counts`、`rnaseq-annotate-de`、`rnaseq-pathway-de`、`rnaseq-te-analysis`
- 集合工具：`rnaseq-counts-to-de`、`rnaseq-de-visuals`
- QC 工具：`rnaseq-bw-cor`、`rnaseq-two-sample-scatter`

## 5. 用 `/path/to/lab-data/RNAseq` 测试

如果服务器挂载了 `/path/to/lab-data/RNAseq`，可以先挑一个已有 BAM 目录，从 BAM 跳过 alignment：

```bash
find /path/to/lab-data/RNAseq -name "*.bam" | head
```

假设 BAM 在 `/path/to/lab-data/RNAseq/project1/bam`，先计数：

```bash
rnaseq-bam-to-counts \
  --bam-dir /path/to/lab-data/RNAseq/project1/bam \
  --gene-gtf /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --te-gtf /path/to/reference/TE_GTF/hg38_rmsk_TE.gtf \
  --outdir /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/counts \
  --prefix project1 \
  --layout auto \
  --strandedness 0 \
  --threads 8
```

准备 `condition.csv` 和 `contrast.csv` 后，写 Snakemake 配置：

```yaml
outdir: /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/smk
species: hg38
sample_table: /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/condition.csv
contrast_file: /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/contrast.csv
tx2gene_path: /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf
te_annotation_tsv: /path/to/te_annotation.preview.csv
matrices:
  project1_gene:
    type: gene
    path: /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/counts/project1.gene.count_matrix.csv
  project1_TE:
    type: te
    path: /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/counts/project1.TE.count_matrix.csv
    te_label_level: repName
    te_color_level: repFamily
```

运行：

```bash
bash run_snakemake_visualization.sh \
  --config /path/to/lab-data/RNAseq/project1/rnaseq_downstream_test/config.yaml \
  --cores 8
```
