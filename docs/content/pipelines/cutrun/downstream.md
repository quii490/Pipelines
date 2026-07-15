# CUT&RUN 下游与 TE 分析

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

下游包括 BAM/bigWig correlation、consensus peaks、ChIPseeker annotation、pathway、HOMER、signal heatmap、TE family/subfamily 与 locus heatmap。

T3E、Allo 和 RepEnTools 是可选外部工具，只有安装与输入验证通过时才执行。RepEnTools FASTQ 工作流使用 CHM13/T2T 资源，不能把 hg38/mm39 BAM 静默当作 CHM13 数据。

高分辨率 `full` 档位可能显著增加 CPU、内存和磁盘，应先用 `standard` 检查结果。
