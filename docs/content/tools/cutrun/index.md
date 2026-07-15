# CUT&RUN 常用工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

工具位于 `CUTnRUN/tools/cutrun_cli/` 和 `CUTnRUN/tools/te_methods/`。主 pipeline 能完成时优先使用主入口；独立工具适合恢复中断结果、只补报告/QC，或从已有 BAM/bigWig/peaks 做透明的小分析。CUTnRUN 当前仍在开发，运行前以实际 `--help` 和小数据验证为准。

| 任务 | 页面 |
|---|---|
| 运行前检查、状态 ledger、完整性 inventory | [Preflight、状态与完整性](validation-status.md) |
| 恢复 peaks、下游、QC 和报告 | [恢复、状态与报告](recovery-report.md) |
| BAM/bigWig correlation、heatmap、consensus、annotation、HOMER | [信号、peaks 与 annotation](signal-peaks.md) |
| Gene/TE annotation、peak overlap、GO/Reactome | [注释、overlap 与通路](annotation-enrichment.md) |
| TE multi-mapping、locus heatmap、T3E/Allo/RepEnTools | [TE 工具与外部方法](te-methods.md) |

## 完整脚本清单

| 工具 | 类型 | 作用 |
|---|---|---|
| `cutrun_preflight.py` | 用户/流程入口 | manifest、FASTQ、reference、空间的 fail-fast 检查 |
| `cutrun_call_peaks_published.sh` | 恢复工具 | 从已发布 clean BAM 补 narrow/broad peaks |
| `cutrun_bam_cor.sh` | 用户工具 | clean BAM bins correlation |
| `cutrun_bw_cor.sh` | 用户工具 | 同类 bigWig bins correlation |
| `cutrun_dt_heatmap.sh` | 用户工具 | 任意 region 的 deepTools heatmap/profile |
| `cutrun_consensus_peaks.py` | 用户/流程工具 | replicate-supported narrow/broad consensus |
| `cutrun_peak_overlap.py` | 用户工具 | 多 peak 集合成对 overlap 统计 |
| `cutrun_annotate_peaks` | 用户工具 | ChIPseeker/ChIPpeakAnno gene 与 TE annotation |
| `cutrun_homer.sh` | 用户工具 | HOMER annotation 和 de novo motif |
| `cutrun_pathway_enrichment.R` | 用户/流程工具 | annotation gene IDs 的 GO BP/Reactome ORA |
| `cutrun_te_qc.py` | 用户/报告工具 | classical 与 TE QC 汇总 |
| `cutrun_te_multimap.py` | 用户工具 | TE BAM secondary/NH 等策略核查 |
| `cutrun_te_locus_heatmap.sh` | 用户工具 | strand-aware TE 5′/body heatmap |
| `cutrun_results_summary.py` | 报告工具 | 核心输出 inventory 和完整性状态 |
| `cutrun_run_manifest.py` | 归档工具 | 输入、reference、版本和 fingerprints |
| `cutrun_report.py` | 报告工具 | 生成 self-contained HTML/PDF report |
| `cutrun_status.py` | 内部/高级 | downstream attempt-aware status ledger |
| `run_te_methods.sh` | 高级适配器 | T3E、Allo、RepEnTools plan/execute/status |
| `run_repentools_fastq.sh` | 高级适配器 | manifest FASTQ → RepEnTools CHM13 workflow |
| `prepare_repentools_reference.sh` | 参考构建 | CHM13 FASTA/RepeatMasker → index/GTF profile |
| `build_repeat_bed.py` | 内部辅助 | TE SAF + annotation → RepeatMasker-like BED |
| `te_method_plan.py` | 内部辅助 | 从 manifest 生成外部方法 plan |
| `render_te_method_visuals.py` | 报告辅助 | 只对可解析的真实结果生成 SVG 汇总 |
| `verify_installation.sh` | 维护辅助 | 记录维护者本地外部工具安装状态；当前非通用安装器 |
| `*_chipseeker.R`、`*_chippeakanno.R` | 内部后端 | `cutrun_annotate_peaks` 调用的 R backend |

## 三类输入先分清

- `04_clean_bam`：标准高可信分支，常规 peak/gene/track；
- `04_te_bam`：保留更多 multi-mapping，用于 TE counts/signal；
- `04_te_locus_best_bam`：确定性 best-locus 候选定位，不替代 fractional/EM counts。

工具返回 0 也要验证输出。`PASS`、`FAIL`、`SKIP`、`INCOMPLETE` 的含义见[输出指南](../../pipelines/cutrun/outputs.md)。

## 使用原则

1. 优先运行 wrapper，而不是直接调用内部 R/Python backend；
2. 同一个 correlation/heatmap 中只放相同 assembly、过滤和 normalization 的输入；
3. recovery 输出写入新目录或保留 run ID，不覆盖最后一个已验证结果；
4. 外部 TE 工具先 plan，再 `--execute`；
5. CUTnRUN 稳定前，所有命令都应在小型公开/授权数据上核对输出和失败语义。
