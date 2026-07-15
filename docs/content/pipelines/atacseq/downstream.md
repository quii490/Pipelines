# ATAC-seq 下游与 TE 分析

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

## Peak 与 fixed-bin

- peak-level：聚焦明确开放区域，结果依赖 peak universe 和 calling 策略；
- fixed-bin：无须预先定义 peaks，适合全局扫描，但多重检验更多且解释粒度不同；
- `--levels both` 可同时生成，但不能挑选更显著的一套而忽略分析计划。

## 推荐分析顺序

1. 检查 metadata、contrast 方向和生物学重复；
2. 查看 count distribution、sample correlation 与 PCA；
3. 运行差异可及性，报告 log2FC 与 FDR；
4. 将差异区域关联到基因、motif 或已知调控元件；
5. 在 bigWig 和 heatmap 中验证区域级信号；
6. 使用独立证据支持功能解释。

## TE 分析

TE overlap、metaprofile 和 relaxed tracks 可回答家族/亚家族层面的信号富集问题。TE 区域高度重复，低 MAPQ reads 和 locus assignment 存在歧义。

```text
问题：某 TE family 周围是否整体更开放？ → metaprofile / aggregate heatmap
问题：哪些 TE loci 发生变化？           → locus-aware quantification + 严格多重比对策略
问题：某个浏览器轨迹是否更高？          → 先核对 normalization、MAPQ 和可比性
```

任何 TE 结论都应记录 annotation 版本、overlap 规则、multi-mapping 策略、MAPQ 与 normalization。
