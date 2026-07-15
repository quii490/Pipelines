# 输出文件参考

| 内容 | RNA-seq | ATAC-seq | CUT&RUN |
|---|---|---|---|
| Alignment | gene/TE BAM | clean BAM | standard/TE/locus BAM |
| Signal track | 可选 bigWig | standard/TE bigWig | standard/TE/locus bigWig |
| Quantification | gene/TE counts | peak/bin counts | gene/TE counts |
| Regions | — | peaks/consensus | narrow/broad/consensus peaks |
| Report/QC | MultiQC/plots | `QC_REPORT.md` | run report + status ledger |

最终结果与 Nextflow work/cache 分开归档。删除缓存前确认无需 `--resume`，并保留输入表、参数、日志、版本和核心结果。
