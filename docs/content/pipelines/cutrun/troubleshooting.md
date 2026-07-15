# CUT&RUN 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 | `main` |

## Preflight 失败

查看生成的 preflight JSON，修复 manifest 列、FASTQ gzip、reference、control 或磁盘。第一次正式运行不建议 `--skip-preflight`；resume 时已通过检查可复用，只有明确诊断场景才跳过。

## Control 配对错误

先检查 manifest，再检查自动 control 推断日志。使用 `--control-sample` 或 `--control-regex` 显式指定；必要时启用 `--no-auto-control`。不要在结果生成后才修改标签而不重算。

## Preview 通过但正式任务失败

Preview 验证的是计划和部分输入，不保证计算资源、环境或所有文件可读。找到第一个失败 task 的日志，修复根因后从相同 `--work-dir` 使用 `--resume`。

若提示 Nextflow/Java，先运行 `java -version` 与 `nextflow -version`，使用 pipeline 配置的 Java 17+ 环境。不要在多个 shell 中混用不同 Java/Nextflow。

## Peaks 为空或过多

确认 assay、target 属于 narrow 还是 broad、matched control、MACS cutoff、BAM 过滤与有效 reads。不要只通过放宽 p/q cutoff 追求 peak 数量。

## Track 有信号但 counts 为空

核对所用 BAM 类型、region coordinates、chrom names、featureCounts format 和文件路径。track 与 matrix 可能来自不同过滤策略。

## TE 外部方法失败

查看 `module_status_latest.tsv` 和方法日志；确认依赖、annotation、assembly 与输入格式。unsupported genome 应记录 skipped，而不是生成占位“成功”结果。

## TE BAM secondary=0 或没有 NH tag

```bash
samtools view -c -f 256 /path/to/04_te_bam/sample_te.bam
samtools view /path/to/04_te_bam/sample_te.bam | head
```

若没有 secondary/NH，检查 TE alignment 和清理阶段，不能可靠做 fractional/multi-mapping 分析。从 standard clean BAM 无法恢复这些信息，通常需从 FASTQ 重跑 TE 分支。

## Locus-best 结果看似非常特异

检查 read tie、mappability、MAPQ 与 standard/TE-aware tracks。locus-best 是候选定位策略，不能自动消除重复序列歧义。

## 报告是旧版本

核对文件时间、run ID 和 `latest` 指向。正式解释应关联本次 run 的不可变产物，而非仅依赖覆盖式 latest 文件。

## Report HTML 有但 PDF 缺失

HTML 是主要有效报告。PDF 转换依赖 weasyprint/wkhtmltopdf/Chromium；查看 `run_report.pdf.status.txt`。不要因 PDF 缺失判定核心分析失败。

## Heatmap 很慢或磁盘增长

使用 `--resource-tier small/standard`，或调高 `TE_LOCUS_BIN_SIZE`、降低 `TE_LOCUS_MAX_REGIONS`。默认不写大型 raw matrix；只有明确需要时才开启。`te-k` 增大会增加 TE BAM 和下游计算。

## 远端 Pod/节点中断

先确认 results 与 `_automation/work` 在持久化卷、PID 已不存在，再以同一 commit/环境/命令 `--resume`。不要把不同 genome build 或 `te-k` 的恢复结果写入旧目录。
