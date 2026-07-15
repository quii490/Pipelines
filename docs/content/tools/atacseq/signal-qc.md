# ATAC-seq 信号与 QC 工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## TSS profile

```bash
bash ATAC-seq/scripts/run_tss_qc.sh \
  --bw-glob "/path/to/04_bw/*.bw" \
  --species hg38 \
  --outdir /path/to/qc/tss \
  --cores 8 \
  --upstream 3000 \
  --downstream 3000
```

输出 matrix、profile 和 heatmap。这里是 bigWig-based TSS signal profile，不等同于严格基于 Tn5 insertions 的 ENCODE TSS enrichment score。

## Gene-body profile

```bash
bash ATAC-seq/scripts/run_gene_body_profile.sh \
  --bw-glob "/path/to/04_bw/*.bw" \
  --species hg38 \
  --outdir /path/to/qc/gene_body \
  --body-length 5000 \
  --cores 8
```

用于观察转录单元周围 aggregate signal。scale-regions 会把不同长度 gene body 缩放到相同长度，不能解释为真实碱基距离。

## bigWig correlation 与 PCA

```bash
bash ATAC-seq/scripts/run_bw_correlation.sh \
  --bw-glob "/path/to/04_bw/*.bw" \
  --outdir /path/to/qc/bw_correlation \
  --cores 8
```

只比较同一 assembly、normalization、bin size 和过滤策略的 tracks。相关性高表示全局形状相似，不证明 peak-level 变化不存在，也不替代 biological replicates。

## Nucleosome phasing

```bash
bash ATAC-seq/scripts/run_nuc_phasing.sh \
  --bam-glob "/path/to/clean_bam/*.clean.bam" \
  --outdir /path/to/qc/nucleosome \
  --cores 8
```

主要用于 PE fragment lengths。SE 数据或深度不足时应报告不适用/NA，而不是用 read length 代替 fragment distribution。

## 常见空图原因

- BED 与 bigWig chromosome naming 不一致；
- bigWig 不是目标 assembly；
- regions 超出 contig 范围；
- `--skipZeros` 删除大量无信号区域；
- track normalization 或过滤策略不同；
- 输入 glob 没有匹配文件。
