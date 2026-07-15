# CUT&RUN 常用工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

工具位于 `CUTnRUN/tools/cutrun_cli/` 和 `CUTnRUN/tools/te_methods/`。主 pipeline 能完成时优先使用主入口；独立工具适合恢复中断结果、只补报告/QC，或从已有 BAM/bigWig/peaks 做透明的小分析。

| 任务 | 页面 |
|---|---|
| 恢复 peaks、下游、QC 和报告 | [恢复、状态与报告](recovery-report.md) |
| BAM/bigWig correlation、heatmap、consensus、annotation、HOMER | [信号、peaks 与 annotation](signal-peaks.md) |
| TE multi-mapping、locus heatmap、T3E/Allo/RepEnTools | [TE 工具与外部方法](te-methods.md) |

## 三类输入先分清

- `04_clean_bam`：标准高可信分支，常规 peak/gene/track；
- `04_te_bam`：保留更多 multi-mapping，用于 TE counts/signal；
- `04_te_locus_best_bam`：确定性 best-locus 候选定位，不替代 fractional/EM counts。

工具返回 0 也要验证输出。`PASS`、`FAIL`、`SKIP`、`INCOMPLETE` 的含义见[输出指南](../../pipelines/cutrun/outputs.md)。
