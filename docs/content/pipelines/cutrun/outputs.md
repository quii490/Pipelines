# CUT&RUN 输出

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

```text
results/
├── 04_clean_bam/
├── 04_te_bam/
├── 04_te_locus_best_bam/
├── 05_tracks/
├── 05_tracks_te/
├── 05_tracks_te_locus_best/
├── 06_peaks/
├── 07_featurecounts/
├── 08_analysis/results/
├── 09_downstream/
│   ├── module_status_latest.tsv
│   ├── qc_metrics.tsv
│   └── run_report.html
└── _automation/
```

先看模块最新状态、run report、QC 和 paired control，再解释 peaks、annotation、heatmap 或 TE 方法输出。
