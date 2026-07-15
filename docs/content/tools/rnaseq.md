# RNA-seq 工具

| 工具/入口 | 输入 | 用途 |
|---|---|---|
| `rnaseq-bam-to-counts` | BAM + GTF | 从已有 BAM 补 gene/TE counts |
| `rnaseq-counts-to-de` | count matrix + design | 差异、PCA、heatmap、pathway |
| `run_qc_metrics_bam.sh` | BAM | 补充 BAM/Picard QC |
| `run_telescope_bam.sh` | BAM + TE annotation | locus-level TE quantification |
| Snakemake visualization | results/count matrices | 批量组织下游可视化 |

工具参数以 `RNA-seq/tools/README.md`、`TOOLS.md` 和每个命令的 `--help` 为准。输入不完整时不得把工具输出描述为完整 pipeline 结果。
