# 从 BAM 生成 gene/TE counts

| 状态 | 维护人 | 最后验证 | 适用版本 |
|---|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 | `rnaseq-bam-to-counts` |

适用于已有 coordinate-sorted BAM，但缺少标准 raw count matrix 的情况。它不会重新评估 FASTQ trimming、污染或 alignment 参数。

```bash
rnaseq-bam-to-counts \
  --bam-dir /path/to/bam \
  --gene-gtf /path/to/reference/genes.gtf \
  --te-gtf /path/to/reference/repeats.gtf \
  --outdir /path/to/counts \
  --prefix project \
  --layout auto \
  --strandedness 0 \
  --threads 8
```

`--strandedness` 使用 featureCounts 数字：`0=unstranded`、`1=forward`、`2=reverse`；与主入口的文字值不同。

运行前检查：

```bash
samtools quickcheck /path/to/bam/*.bam
samtools view -H /path/to/bam/sample.bam | head
```

确认 BAM/GTF 同一 assembly，BAM 已 coordinate-sort，并知道 PE/SE、strand、duplicate 与 secondary alignment filtering。

典型输出：

```text
project.gene.featureCounts.txt
project.gene.count_matrix.csv
project.TE.featureCounts.txt
project.TE.count_matrix.csv
project.bam_to_counts.summary.txt
```

查看 `Assigned`、`Unassigned_NoFeatures`、`Unassigned_MultiMapping` 和各样本总 counts。gene assigned rate 低优先检查 strand、GTF 和 assembly。

!!! warning

    从 unique-only BAM 无法恢复 TE multi-mapping 信息。能输出 TE counts 不代表输入适合 locus-aware TE 解释。
