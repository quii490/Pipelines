# 从 raw counts 到差异分析

| 状态 | 维护人 | 最后验证 | 适用版本 |
|---|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 | `rnaseq-counts-to-de` |

输入必须是 raw count matrix，不是 TPM/FPKM/CPM/VST。矩阵列名必须与 sample table 完全一致。

```bash
rnaseq-counts-to-de \
  --matrix /path/to/count_matrix.csv \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --outdir /path/to/downstream \
  --matrix-name gene_featureCounts \
  --make-plots true \
  --threads 4
```

gene 分析可增加 `--tx2gene-path /path/to/genes.gtf --species hg38 --run-go true --run-gsea true`。TE 分析使用 `--te-annotation-tsv`、`--te-label-level` 和 `--te-color-level`。

```csv
gene_id,KO_1,KO_2,WT_1,WT_2
ENSG000001,120,133,80,90
ENSG000002,5,6,30,28
```

运行后先检查纳入样本、condition、contrast 方向和过滤数量，再看：

- `DE_matrix.csv`：完整差异结果；
- `normalized_counts.csv`：展示用，不作为 raw count model 输入；
- `vst_or_log_matrix.csv`：PCA/heatmap；
- volcano/MA/heatmap：不同视角，不是三个独立统计检验；
- pathway status：区分无显著结果、输入不足与执行失败。

没有 biological replicates 时只能使用明确标注的探索性方法。不得把 fixed BCV 或 logCPM difference 描述为标准重复设计检验。
