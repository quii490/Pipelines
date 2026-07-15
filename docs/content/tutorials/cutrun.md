# CUT&RUN 实战教程

| 状态 | 数据 | 最后审查 |
|---|---|---|
| Draft review | 使用含 matched control 的公开/已授权数据 | 2026-07-16 |

选择一个 target、至少 2 个 biological replicates，并包含明确 IgG/input control。教程练习 manifest、preview、三类 BAM/track、恢复下游和 TE 方法状态。

## A. 生成并审查 manifest

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/tutorial/cutrun/fastq \
  --species hg38 --assay cutrun \
  --outdir /path/to/tutorial/cutrun/results \
  --init-only
```

确认列：`sample,species,assay,group,replicate,igg,is_igg,layout,fastq_1,fastq_2`。每个 target 的 `igg` 必须指向 `is_igg=true` 的实际 sample。

## B. Preview

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_pipeline.sh \
  --manifest /path/to/tutorial/cutrun/results/_automation/inputs/manifest.csv \
  --outdir /path/to/tutorial/cutrun/results \
  --preview
```

修复 FASTQ、control、reference 或磁盘错误，不使用 `--skip-preflight` 掩盖问题。

## C. 正式运行

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/tutorial/cutrun/fastq \
  --manifest /path/to/tutorial/cutrun/results/_automation/inputs/manifest.csv \
  --species hg38 --assay cutrun \
  --outdir /path/to/tutorial/cutrun/results \
  --resource-tier standard \
  --run-downstream --background --resume
```

第一次保持 `te-k=25`、TE duplicate policy `mark` 和默认 tracks，不同时测试多个高级策略。

## D. 判断完成状态

先看：

```text
09_downstream/run_report.html
09_downstream/module_status_latest.tsv
09_downstream/module_status_summary.json
09_downstream/run_manifest.json
```

区分核心结果与可选模块。HOMER/外部 TE 方法失败不会删除 BAM/counts，但必须在报告中解释。

## E. 比较三类 signal

在 IGV 或统一 heatmap 中分别查看：

- standard track：常规高可信 signal/peak；
- TE-aware track：保留更多重复区域证据；
- locus-best track：候选 locus 形状。

不要把三类 track 放进同一个 correlation matrix；不要把 locus-best 当作 fractional counts。

## F. 下游恢复

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_downstream.sh \
  --results-dir /path/to/tutorial/cutrun/results \
  --species hg38 \
  --resource-tier standard \
  --resume
```

观察 attempt-aware status：已通过模块应复用，失败模块重跑。确认 `latest` 与本次 run ID 一致。

## G. TE 方法 plan

先以不执行模式检查 T3E/Allo/RepEnTools 计划和依赖，再决定是否 `--execute`。RepEnTools 若没有 CHM13 reference、2 ChIP + 2 input，应为 SKIP，不应生成空图。

## 教程验收

- [ ] 能修正多个 control 的 `igg` mapping；
- [ ] 能区分 narrow 与 broad peaks；
- [ ] 能解释 PASS/FAIL/SKIP/INCOMPLETE；
- [ ] 能解释 TE secondary/NH tags 为什么重要；
- [ ] 能说明 motif enrichment、peak number 和 correlation 各自不能证明什么。
