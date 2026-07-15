# RNA-seq 实战教程

| 状态 | 数据 | 最后审查 |
|---|---|---|
| Draft review | 使用公开/已授权数据自行替换路径 | 2026-07-16 |

本教程不绑定可能失效的下载链接，重点练习一次完整分析和三种中间文件入口。建议选择两组、每组至少 2 个 biological replicates 的小型 bulk RNA-seq 数据。

## 目标

完成后应能：识别 PE/SE 和 strandedness；定义 CASE vs CONTROL；从 FASTQ 运行；从 BAM/counts/DE 表续做；找到核心输出；区分 raw counts、normalized counts 和 TPM。

## A. 从 FASTQ 完整运行

```bash
cd /path/to/Pipelines/RNA-seq/rnaseq

bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/tutorial/rnaseq/fastq \
  --results-dir /path/to/tutorial/rnaseq/results \
  --species hg38 \
  --strand reverse \
  --layout auto \
  --dry-run
```

检查样本数、R1/R2、reference 和计划模块。通过后去掉 `--dry-run`，加入资源和后台参数：

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/tutorial/rnaseq/fastq \
  --results-dir /path/to/tutorial/rnaseq/results \
  --species hg38 \
  --strand reverse \
  --background --resume \
  --max-cpus 16 --max-memory "64 GB"
```

普通 bulk RNA-seq 保持 `--run-dedup false`；duplicate rate 由默认 `--run-markdup-qc true` 作为 QC 记录。

## B. 审查设计并运行下游

编辑 `results/condition.csv`：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

编辑 `results/contrast.csv`：

```csv
case,control
KO,WT
```

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/tutorial/rnaseq/results \
  --species hg38 \
  --downstream \
  --plot-outdir /path/to/tutorial/rnaseq/results/plots_tutorial
```

## C. 阅读结果

按顺序回答：

1. MultiQC 中是否缺样本？
2. gene STAR mapping、featureCounts assigned 和 strand 是否支持当前输入？
3. TE STAR 为什么 multi-mapping 高于 gene branch？
4. PCA 是否按 condition 聚类，是否存在 batch/离群？
5. `KO_vs_WT` 正 log2FC 的方向是什么？
6. volcano、MA 和 heatmap 是否由同一个 DE table 生成？
7. pathway 是 ORA 还是 GSEA，使用什么 universe/ranking？

## D. 从 BAM 重新得到 counts

```bash
rnaseq-bam-to-counts \
  --bam-dir /path/to/tutorial/rnaseq/bam \
  --gene-gtf /path/to/genes.gtf \
  --te-gtf /path/to/repeats.gtf \
  --outdir /path/to/tutorial/rnaseq/counts_from_bam \
  --prefix tutorial \
  --layout PE --strandedness 2 --threads 8
```

比较其 gene assigned rate 和主流程结果。若差异大，检查 BAM 来源、strand、GTF、PE flags 和 filtering，而不是直接选择数值更大的版本。

## E. 从 raw count matrix 重新做 DE

```bash
rnaseq-counts-to-de \
  --matrix /path/to/tutorial/rnaseq/gene_count_matrix.csv \
  --sample-table /path/to/tutorial/rnaseq/condition.csv \
  --contrast-file /path/to/tutorial/rnaseq/contrast.csv \
  --tx2gene-path /path/to/genes.gtf \
  --species hg38 \
  --outdir /path/to/tutorial/rnaseq/counts_downstream \
  --matrix-name gene_featureCounts \
  --make-plots true --run-go true --run-gsea true --threads 4
```

确认输入是 raw counts，不是 TPM/VST。记录实际纳入样本与过滤 feature 数。

## 教程验收

- [ ] 能解释 strandedness 写错为何降低 assigned rate；
- [ ] 能说明 BAM 不能直接画 volcano；
- [ ] 能区分 normalized counts、VST 和 TPM；
- [ ] 能解释无重复结果为什么只能 exploratory；
- [ ] 能从完整 DE 表重建图，而不只保留筛选后的显著表。
