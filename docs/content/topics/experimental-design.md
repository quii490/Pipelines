# 实验设计、重复与 contrast

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | Analysis maintainers | 2026-07-16 |

## Biological 与 technical replicate

biological replicates 是独立生物样本，用于估计生物变异；technical replicates/lane 主要反映技术重复，通常不能替代 biological replicates。先按 pipeline 规则合并 lanes，再保留 biological sample identity。

## 最小 metadata

```csv
sample,condition,replicate,batch
WT_1,WT,1,B1
WT_2,WT,2,B2
KO_1,KO,1,B1
KO_2,KO,2,B2
```

复杂设计还需要 donor、sex、time、treatment、pair/block 等。sample ID 必须匹配 count matrix/BAM/track。

## Contrast 方向

本文档统一使用 `CASE vs CONTROL`。`KO,WT` 或 `KO:WT` 表示正 log2FC 为 KO 高于 WT。把方向写入分析计划和文件名。

## 配对和 batch

同一 donor 前后处理应使用 paired/block design；仅把 donor 写在表里但不加入模型没有作用。batch 可以加入模型的前提是它与 condition 不完全混杂。

## 无重复数据

没有 biological replicates 时不能可靠估计组内变异。log fold change、scatter 或 fixed BCV 只能作为探索性结果，必须明确降级措辞，不能报告成标准 FDR-supported differential result。

## Control 的 assay 差异

- RNA-seq：对照通常是 biological condition；
- ATAC-seq：通常没有 IgG，但需批次和实验质量对照；
- CUT&RUN/CUT&Tag/ChIP-seq：IgG/input/matched control 影响 peak calling 与背景解释，manifest 必须逐 target 配对。

## 避免先看结果再改设计

先声明 primary contrasts、过滤、协变量和排除规则。查看 PCA 后发现标签错误可以纠正；因为结果“不显著”而反复重定义 group/cutoff 会增加假阳性风险。
