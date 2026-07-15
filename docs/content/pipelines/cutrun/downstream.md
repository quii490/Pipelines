# CUT&RUN 下游与 TE 分析

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | CUT&RUN maintainers | 2026-07-16 | `main` |

## 常规下游

1. 用 standard BAM/track 完成基础 QC；
2. 采用与 target 匹配的 narrow/broad peak 策略；
3. 建立一致的 region universe 与 count matrix；
4. 检查 correlation/PCA 和 metadata；
5. 进行 differential binding、annotation、motif 或 metaprofile；
6. 在浏览器中核验代表性区域。

## TE 问题与方法选择

| 问题 | 优先证据 |
|---|---|
| 某 TE subfamily 是否整体富集？ | TE overlap、aggregate profile、family-level statistics |
| 某 locus 是否有候选 binding？ | locus-best track + mappability + orthogonal evidence |
| 处理组是否变化？ | replicate-aware count model + predeclared contrast |
| 外部 TE 方法是否成功？ | module status、真实输出与方法专属 QC |

T3E、Allo、RepEnTools 等方法的输入定义和支持 assembly 不同。只比较可比的统计量；不把某方法 skipped 当作阴性结果。

TE 结果必须记录 annotation release、multi-mapping 策略、BAM 类型、control、peak/overlap 规则和多重检验方法。
