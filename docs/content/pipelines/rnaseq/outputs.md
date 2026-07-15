# RNA-seq 输出

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

典型结果包含：

```text
results/
├── 02_gene_star/             # gene alignment 与日志
├── 03_gene_featurecounts/    # gene count matrix
├── 06_te_star/               # TE alignment 与日志
├── condition.csv
├── contrast.csv
├── plots/                    # PCA、差异、热图、GSEA 等
└── _automation/              # 输入快照、日志和 Nextflow work/cache
```

优先检查：模块完成状态与日志、MultiQC、样本完整性、gene count matrix、PCA，再解释 DEG/TE/GSEA。`_automation/work/` 是缓存，不是最终结果。
