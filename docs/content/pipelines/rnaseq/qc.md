# RNA-seq Quality Control

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

| 指标 | 经验性检查 | 异常时优先排查 |
|---|---|---|
| Q30/adapter | 多数 bases 高质量、adapter 可控 | 测序质量、read length、trim 设置 |
| Mapping rate | 同批次应接近；明显偏低需调查 | 物种/build、污染、read 质量 |
| Assigned reads | 与 mapping rate 和 annotation 一致 | GTF build、strand、feature type |
| rRNA fraction | 异常升高提示建库问题 | depletion、污染、annotation |
| PCA/correlation | 生物重复应比跨组更接近 | 样本标签、批次、离群样本 |

阈值依赖物种、组织、建库和实验目标，不能仅凭单一百分比删除样本。所有排除决定必须记录。
