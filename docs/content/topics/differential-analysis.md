# 差异分析与可视化

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft | Analysis maintainers | 2026-07-16 |

## 分析前

1. 确认每行样本的 condition、replicate、batch 和配对关系；
2. 预先写明 `CASE vs CONTROL` 方向；
3. 使用 raw counts 进入 count-based model，不把 TPM/CPM 当作输入 counts；
4. 先检查 library size、sample correlation、PCA 和异常样本；
5. 有批次或配对设计时使用相应 design，而不是事后“校正图”。

## 解释结果

同时报告 effect size（如 log2FC）、uncertainty/p-value 与 adjusted p-value。FDR 不显著不等同于“完全无效应”，尤其在低重复或低计数时；大 fold change 也可能来自不稳定的低计数。

| 图 | 用途 | 常见误用 |
|---|---|---|
| PCA | 发现全局结构、批次和标签问题 | 把分离当作因果证明 |
| MA plot | 检查 abundance 与 effect | 忽略低计数不稳定性 |
| Volcano | 同时展示 effect 与 significance | 只标最显著而不看表达量 |
| Heatmap | 展示选定 feature 的模式 | 用全数据先筛再夸大聚类 |

保存用于绘图的表、筛选规则、随机种子和脚本版本，保证图可以重建。
