# ATAC-seq Quality Control

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

| 指标 | 判断重点 | 异常时排查 |
|---|---|---|
| Mapping rate | 同批次一致性 | build、污染、read 质量 |
| Mitochondrial fraction | 是否异常偏高 | 建库质量、细胞状态 |
| FRiP | signal 是否集中在 peaks | peak 参数、深度、背景 |
| TSS profile/enrichment | TSS 附近富集形态 | TSS build、Tn5 处理、深度 |
| Fragment distribution | PE 是否有核小体周期 | layout、insert size、建库 |
| Replicate correlation | 重复是否一致 | 样本标签、批次、离群样本 |

不要把 profile 曲线数值误写成标准 insertion-site TSS enrichment score。QC 阈值是调查触发器，不是自动删除规则。
