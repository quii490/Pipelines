# ATAC-seq 下游与 TE 分析

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

下游支持 peak-level、fixed-bin differential accessibility、peak annotation、pathway、motif、TOBIAS footprinting、region/TE/gene-body heatmap 与 nucleosome phasing。

`--run-te-relaxed-tracks` 生成探索性 TE/L1 track，但不是 EM/fractional multi-mapper allocation。只有问题确实需要 multi-mapping signal 时才启用，并在报告中说明方法限制。
