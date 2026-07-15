# CUT&RUN 输出与结果解释

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

```text
results/
├── manifest/resolved_manifest.csv
├── 01_prepared_reads/
├── 02_align/                    # standard alignment
├── 02_align_te/                 # TE alignment
├── 03_sorted_bam/
├── 03_sorted_bam_te/
├── 04_clean_bam/                # standard final BAM
├── 04_te_bam/                   # TE-aware final BAM
├── 04_te_locus_best_bam/        # deterministic best-locus BAM
├── 05_tracks/                   # standard RPKM/RPGC bigWig
├── 05_tracks_te/                # TE RPGC bigWig
├── 05_tracks_te_locus_best/     # 5 bp CPM best-locus bigWig
├── 06_peaks/<sample>/broad|narrow/
├── 07_featurecounts/
├── 08_analysis/results/
├── 09_downstream/
└── _automation/
```

## 先看四个文件

```text
09_downstream/run_report.html
09_downstream/module_status_latest.tsv
09_downstream/module_status_summary.json
09_downstream/run_manifest.json
```

`module_status.tsv` 保存 attempt 历史；`latest` 是当前状态。`PASS` 是命令和基本输出检查通过；`FAIL` 是失败/输出缺失；`SKIP` 是不适用或条件不足；`INCOMPLETE` 表示核心文件缺失，不能交付为完成结果。

## 三类 BAM/track

| 分支 | 主要用途 | 不能直接声称 |
|---|---|---|
| standard clean | 常规 peaks、gene counts、RPGC/RPKM track | 完整保留重复区域 reads |
| TE-aware | TE counts、TE aggregate signal | 唯一定位到某 locus |
| locus-best | TE locus 候选 heatmap/IGV | 等价于 fractional/EM counts |

三类 signal 不放进同一个 correlation matrix。图注注明 BAM/track 类型、normalization、bin size、MAPQ/duplicate/blacklist 策略。

## Peaks

每个 target 可有 narrow 和 broad。TF/尖锐信号通常重点看 narrow；broad histone marks 重点看 broad。不能把两类 peaks 合并成一个“更完整”的集合。检查 matched control 和 MACS 日志。

## FeatureCounts 与 analysis

`featurecounts_gene.txt` 和 `featurecounts_te.txt` 是统一分析输入之一。TE 表需结合 SAF/annotation、multi-mapping 和 control normalization 解释。`08_analysis` 中的 lfc-over-IgG/normalized visualizations 不等于 replicate-aware differential test。

## Downstream

常见内容包括 BAM/bigWig correlation、consensus peaks、annotation、HOMER、heatmaps、QC、TE methods 和报告。Consensus 至少需要两个 target replicates；不足时 SKIP 是正确行为。

外部 TE 方法没有输出先看 `method_status.tsv` 和 `visualization_status.tsv`。工具执行成功但必需文件缺失应为 FAIL；空图不代表阴性。

## 归档

保存 resolved manifest、run metadata、日志、status history/latest、run manifest、HTML report、QC、核心 counts/peaks 和关键图。`latest` 方便查看但可被覆盖，正式结论关联 run ID 的不可变产物。
