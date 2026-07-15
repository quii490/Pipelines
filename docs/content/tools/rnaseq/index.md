# RNA-seq 独立工具

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | RNA-seq maintainers | 2026-07-16 | `RNA-seq/tools/` |

独立工具用于从已有中间结果补分析。若仍有 FASTQ、完整设计和可用缓存，优先用主入口，以保留一致的 reference、状态和报告。

| 已有输入 | 推荐入口 | 主要结果 |
|---|---|---|
| BAM | [BAM 到 counts](bam-to-counts.md) | gene/TE raw counts |
| Raw count matrix | [Counts 到差异分析](counts-to-de.md) | normalized matrix、DE、PCA、图和 pathway |
| DE matrix | [DE、通路与绘图工具](de-tools.md) | annotation、volcano/MA、pathway、TE 图 |
| 多个 results/matrix | [Snakemake 批量可视化](snakemake.md) | 多矩阵、多 contrast 的一致输出 |
| bigWig | `rnaseq-bw-cor` | correlation/PCA |
| 两个样本 count | `rnaseq-two-sample-scatter` | exploratory scatter |

安装统一命令前先查看 `RNA-seq/tools/README.md`；每个命令的 `--help` 是参数事实来源。
