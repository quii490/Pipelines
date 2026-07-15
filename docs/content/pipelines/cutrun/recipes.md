# CUT&RUN：按现有文件选择入口

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

## 我有 FASTQ

按 [Quick Start](quick-start.md) 生成 manifest、preview 并运行 `run_auto_chipseq.sh --run-downstream`。

## 核心流程完成，只想运行/恢复下游

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_downstream.sh \
  --results-dir /path/to/results \
  --species hg38 \
  --resource-tier standard \
  --resume
```

downstream 使用 fingerprint/sentinel 复用已通过模块；失败模块会重跑。只有确实要强制全部重算才用 `--no-resume`。

## 已有 published `04_clean_bam`，peak 缺失

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_call_peaks_published.sh \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --species hg38
```

这是中断恢复工具；新项目应让 Nextflow 主流程生成 peaks。运行前确认 matched control 与 narrow/broad 策略。

## 已有 counts，只想重画 gene/TE 图

使用 `cutrun_count_draw --skip-featurecounts --force-r-reprocess`，输出到新的 analysis 目录。counts 必须来自相同样本集合和 counting 口径。

## 已有 BAM/bigWig/peaks，只做独立分析

- BAM correlation：`cutrun_bam_cor.sh`；
- bigWig correlation：`cutrun_bw_cor.sh`；
- region heatmap：`cutrun_dt_heatmap.sh`；
- consensus peaks：`cutrun_consensus_peaks.py`；
- annotation/HOMER/overlap：相应 `cutrun_*` 工具；
- 详细用法见[CUT&RUN 工具手册](../../tools/cutrun/index.md)。

## 只补报告、QC 或可复现清单

- `cutrun_te_qc.py`：从现有 run 汇总 classical/TE QC；
- `cutrun_results_summary.py`：机器可读结果 inventory；
- `cutrun_run_manifest.py`：记录输入、reference、配置和 fingerprint；
- `cutrun_report.py`：重建 self-contained HTML/PDF report。

## 只跑外部 TE 方法

必须从核心流程的 resolved manifest 和 TE BAM 开始。T3E/Allo 可由 `run_te_methods.sh` 调度；RepEnTools 是独立 CHM13/T2T FASTQ 工作流，不能把 hg38/mm39 BAM 当作其输入。失败、输入不足或不支持必须写为 `FAIL/SKIP`。
