# CUT&RUN 工具

| 工具 | 用途 |
|---|---|
| `cutrun_status.py` | 查看模块运行状态 |
| `cutrun_results_summary.py` | 汇总结果 |
| `cutrun_report.py` | 生成运行报告 |
| `cutrun_consensus_peaks.py` | 生成重复支持的 consensus peaks |
| `cutrun_peak_overlap.py` | Peak overlap |
| `cutrun_te_qc.py` | TE 分支 QC |
| `cutrun_te_locus_heatmap.sh` | TE locus heatmap |

工具位于 `CUTnRUN/tools/cutrun_cli/`。优先通过主入口的 `--run-downstream` 调度；手工调用时必须保存输入、参考、参数和软件版本。
