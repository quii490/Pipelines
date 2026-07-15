# RNA-seq 下游与 TE 分析

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

下游读取 `condition.csv` 和 `contrast.csv`，生成 gene/TE differential expression、PCA、volcano、heatmap、GSEA 和汇总图。默认遇到某模块缺样本时跳过该模块，避免生成残缺 panel；只有明确接受部分样本时才启用 allow 策略。

TE family/subfamily 与 locus-level 方法回答不同问题。优先在多个方法间比较方向和稳健性，不以单一工具的显著结果作为唯一证据。参见 [TE 分析专题](../../topics/te-analysis.md)。
