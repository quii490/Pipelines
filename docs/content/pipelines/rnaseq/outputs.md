# RNA-seq 输出指南
实际目录随启用模块变化；不存在的可选目录不代表核心流程失败。先看状态和 QC，再看差异图。
## 结果目录概览

```text
results/
├── 01_fastp*/                 # clean FASTQ 和 fastp 报告（名称以实际流程为准）
├── 02_gene_star/              # gene alignment、BAM、STAR 日志
├── 03_gene_featurecounts/     # gene counts 和 featureCounts 日志
├── 04_salmon*/                # transcript abundance
├── 06_te_star/                # TE alignment/logs
├── 07_te*/                    # 启用的 TE quantification 输出
├── condition.csv
├── contrast.csv
├── plots/                     # 下游矩阵、差异、图和 pathway
└── _automation/
    ├── inputs/                # resolved input/运行输入快照
    ├── logs/                  # 自动化和 Nextflow 日志
    └── work/                  # Nextflow cache，不是最终结果
```

!!! note
    目录编号和可选模块可能随版本变化。发布结果时以本次 `run manifest` 和实际目录为准，不要仅根据旧截图判断。

## 第一次完成后先看什么

1. 自动化/Nextflow 最终状态和第一个失败模块。
2. MultiQC、fastp/FastQC、STAR/HISAT2 和 featureCounts summaries。
3. gene count matrix 的样本数、列名和非零值。
4. TE 模块状态及其 annotation/层级。
5. condition/contrast 与 matrix 样本的一致性。
6. PCA/correlation，然后才是 DEG、TE、GSEA 和展示图。

## Gene counts

featureCounts 输出通常包含 feature identifier、坐标/长度信息和每个样本的 raw counts。用于统计的矩阵应明确：

- 行 ID 是 Ensembl gene ID、symbol 还是其他 attribute。
- 是否去除 annotation version suffix。
- PE 是 fragment 还是 read counting。
- 使用的 strandedness 和 GTF build。

## Differential expression

常见字段：

| 字段 | 含义 |
|---|---|
| `baseMean` | normalization 后的平均表达量指标 |
| `log2FoldChange` | case 相对 control 的 effect size |
| `lfcSE` | log2FC standard error（若方法提供） |
| `pvalue` | raw P value |
| `padj` | multiple-testing adjusted P value/FDR |

不要只发布筛选后的 DEG 表；同时保留完整结果、过滤阈值和 contrast 定义。

## TE 输出

TEtranscripts、Telescope、REdiscoverTE 等输出不可只凭文件名互换。每个输出应记录：工具、版本、annotation、family/subfamily/locus 层级、raw/normalized 单位和 multi-mapping 策略。

## `plots/`

可能包括 PCA、sample correlation、volcano、MA、heatmap、composition、GSEA/pathway、TE hierarchy/age 图和汇总报告。图用于解释和展示；统计事实仍以对应完整表和参数为准。

## 发布和清理

发布工具只应复制必要的图、报告、状态、版本和小型表，不包含 FASTQ、BAM、work/cache。删除 `_automation/work/` 前确认流程不再需要 resume，并保存 trace/report/log。
