# CUT&RUN 恢复、状态与报告

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

## 只恢复 downstream

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_downstream.sh \
  --results-dir /path/to/results \
  --species hg38 \
  --resource-tier standard \
  --resume
```

可设置 `--threads`、`--peak-jobs`。资源不足时保持 `peak-jobs=1`；`PEAK_JOBS=2` 只适合 CPU/内存和 I/O 充足时。`--no-resume` 会强制更多重算，不是普通排错第一步。

## 从 published clean BAM 恢复 peaks

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_call_peaks_published.sh \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --species hg38
```

它从 `04_clean_bam/*.bam` 重新调用 MACS3 narrow+broad peaks，适合 Nextflow 在发布 BAM 后中断。新 run 应让主流程产生 peaks。先核对 manifest 的 control mapping。

## 补 QC

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_te_qc.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --output /path/to/results/09_downstream/qc_metrics.tsv
```

不可计算指标应记录 `SKIP` 而不是 0。FRiP 的 mapped alignment/read/fragment 口径必须与项目报告一致。

## 结果 inventory

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_results_summary.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv
```

用于找缺失/空文件，不替代科学 QC。

需要在提交作业前检查 manifest/reference/空间时，使用[Preflight、状态与完整性](validation-status.md)，不要等到 alignment 后才发现 FASTQ 或 control mapping 错误。

## 可复现 run manifest

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_run_manifest.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --config /path/to/run.config \
  --output /path/to/results/09_downstream/run_manifest.json
```

可用 `--run-id` 绑定一次 downstream attempt。大型文件默认记录 size、mtime 和 streaming fingerprint；`--hash-large` 会显著增加 I/O，仅在归档政策要求完整 hash 时使用。reference/config 文件仍应完整 SHA256。

## 重建报告

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_report.py \
  --results-dir /path/to/results \
  --output /path/to/results/09_downstream/run_report.html
```

HTML 是主要交付；PDF 缺失可能只是转换器未安装。报告不会把缺失的核心结果变成成功，先看 module status。
