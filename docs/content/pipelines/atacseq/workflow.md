# ATAC-seq Workflow

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

主流程执行 trimming、alignment、mitochondrial/blacklist/duplicate 处理、BAM 与 bigWig、MACS peak calling、feature counting、peak/bin differential analysis 和 QC。

PE 的 blacklist 过滤按 read name 删除整对 mates；PE counting/FRiP 以 fragment 为单位，SE 以 read 为单位。TSS profile 不等同于 ENCODE insertion-site TSS enrichment，解释时必须区分。
