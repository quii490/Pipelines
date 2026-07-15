# Snakemake 批量下游可视化

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 |

用于从已有主流程 results 或多个 count matrix 批量组织一致输出，不替代 FASTQ→BAM/counts 的主 Nextflow 流程。

```bash
cd RNA-seq/snakemake-visualization

bash run_snakemake_visualization.sh \
  --results-dir /path/to/results \
  --outdir /path/to/results/snakemake_visualization \
  --cores 8 \
  --dry-run
```

从 results 生成配置后，人工确认 matrices、condition、contrasts 和启用模块。参数以当前 `--help` 为准。只纳入部分矩阵可使用 `--include-matrix`；关闭 pathway/GO/GSEA/TE 使用 wrapper 的 `--no-*` 选项。

```text
diff/<matrix>/<contrast>/
annotate/<matrix>/<contrast>/
visuals/<matrix>/<contrast>/
pathway/<gene_matrix>/<contrast>/
te/<te_matrix>/<contrast>/
```

输出与主流程 `plots/` 分开，避免覆盖。保留配置、dry-run、DAG/日志和软件版本；配置改变后先重新 dry-run。
