# TE 分析

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft | TE analysis maintainers | 2026-07-16 |

先定义问题层级，再选择工具。family/subfamily abundance、单 locus 信号和 TE 周围 aggregate profile 是三类不同问题。

| 问题 | 典型策略 | 主要限制 |
|---|---|---|
| 哪些 TE subfamilies 表达变化？ | family/subfamily count model | annotation、multi-mapping 分配 |
| 哪些 TE loci 表达变化？ | locus-aware quantification | 唯一定位能力、低计数 |
| 某 mark 是否富集于 TE？ | overlap / metaprofile | region background、mappability |
| 某 locus 是否有 binding/accessibility？ | locus track + orthogonal evidence | locus ambiguity、track normalization |

## Multi-mapping

`unique-only`、best-hit、fractional/EM allocation 和保留多重比对 reads 会产生不同估计。没有一种策略适用于所有问题。报告必须写明 aligner、允许的多重比对数、MAPQ、分配方法和 annotation release。

## 最小报告清单

- genome assembly 与 TE annotation 来源/版本；
- family/subfamily/locus 层级；
- read/fragment 计数单位和 strandedness；
- multi-mapping 与 duplicate 策略；
- normalization、contrast 和统计模型；
- control/background 与 multiple-testing correction；
- 在浏览器或独立实验中的验证证据。

!!! warning

    relaxed bigWig 或 locus-best track 适合发现候选，不能单独证明 locus-specific biological effect。
