# ATAC-seq Quality Control

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

QC 阈值是调查触发器，不是跨物种、细胞类型和实验条件通用的硬门槛。

| 指标 | 看什么 | 异常时优先检查 |
|---|---|---|
| FASTQ quality | base quality、adapter | 测序批次、read trimming |
| Mapping rate | 与同批样本比较 | 物种/index、污染、read 长度 |
| Mitochondrial fraction | 线粒体 reads 占比定义一致 | 细胞状态、建库与过滤口径 |
| Duplicate/library complexity | 可用独立 fragments | 低起始量、PCR、测序过深 |
| FRiP | fragments in peaks | peak set、背景、peak caller |
| Fragment-size pattern | nucleosome-free/mono/di | PE 数据、文库质量、深度 |
| TSS profile | TSS 周围 signal shape | 注释版本、track normalization |
| Replicate correlation/PCA | 同组是否接近 | 标签、批次、异常样本 |

## 评估原则

1. 先确认分母、过滤阶段和 PE/SE 计数单位；
2. 同一批次、相同 assay 和相似深度下比较；
3. 结合多个指标，不能用单一 FRiP 或 mapping rate 决定去留；
4. 排除样本必须记录原因，并分别评估含/不含该样本的结论稳健性。

!!! warning "TSS 指标命名"

    当前 bigWig-based TSS profile 不能直接写成 ENCODE TSS enrichment score。若项目需要该指标，应采用明确的 Tn5 insertion 定义和相应实现。
