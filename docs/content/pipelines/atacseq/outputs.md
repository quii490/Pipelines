# ATAC-seq 输出

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

```text
results/
├── 02_align/                 # clean/final BAM 与比对日志
├── 04_bw/                    # 标准 bigWig
├── 05_peaks/                 # sample/consensus peaks
├── 08_downstream*/           # 差异、annotation、enrichment、TE
├── QC_REPORT.md              # 非生信用户优先查看
└── _automation/              # 输入快照、日志和缓存
```

优先检查 `QC_REPORT.md`、样本是否完整、mapping/mitochondrial/FRiP/TSS、BAM 和 peaks，再查看 differential accessibility、annotation、motif 或 heatmap。
