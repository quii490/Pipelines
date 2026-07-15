# DE、通路和绘图工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 |

| 命令 | 输入 | 用途 |
|---|---|---|
| `rnaseq-annotate-de` | DE table | 补 gene/feature annotation |
| `rnaseq-de-visuals` | DE table | volcano、MA、组成图 |
| `rnaseq-pathway-de` | DE table/ranked list | GO/ORA/GSEA |
| `rnaseq-te-analysis` | TE DE table | class/family/name/age 图 |
| `rnaseq-gsea` | ranked genes | 独立 GSEA |
| `rnaseq-bw-cor` | bigWig | signal correlation/PCA |
| `rnaseq-two-sample-scatter` | 两列 counts | exploratory scatter |

## DE matrix 补图

```bash
rnaseq-de-visuals \
  --de-matrix /path/to/KO_vs_WT.DE_matrix.csv \
  --tx2gene-path /path/to/genes.gtf \
  --annotation-mode gene \
  --padj-cutoff 0.05 \
  --lfc-cutoff 0.58 \
  --label-top-n 30 \
  --outdir /path/to/de_visuals \
  --prefix KO_vs_WT
```

改变 padj/LFC cutoff 只改变展示/筛选，不重新拟合模型。保留未筛选 DE 表。

## Pathway

```bash
rnaseq-pathway-de \
  --de-matrix /path/to/KO_vs_WT.DE_matrix.csv \
  --species hg38 \
  --tx2gene-path /path/to/genes.gtf \
  --case KO --control WT \
  --outdir /path/to/pathway/KO_vs_WT \
  --run-go true --run-gsea true
```

ORA/GO 需要显著 gene list 和正确 background universe；GSEA 使用完整 ranked list。两者问题不同，不应只挑有显著结果的一种报告。

## bigWig correlation

```bash
rnaseq-bw-cor \
  --bw-dir /path/to/bw \
  --out-prefix /path/to/qc/rnaseq_bw \
  --binsize 10000 \
  --threads 8 \
  --methods pearson,spearman
```

bigWig 必须同一 assembly、normalization、bin 和 blacklist 口径。

## 两样本 scatter

```bash
rnaseq-two-sample-scatter \
  --matrix /path/to/count_matrix.csv \
  --sample-x KO_1 --sample-y WT_1 \
  --out-prefix /path/to/qc/KO_1_vs_WT_1 \
  --transform log2cpm \
  --label-top-n 20
```

Pearson 强调线性，Spearman 强调秩关系。相关性高不等于没有系统性偏移，也不能替代 biological replicates。
