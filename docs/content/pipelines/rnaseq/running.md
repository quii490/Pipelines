# 运行与恢复 RNA-seq

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | RNA-seq maintainers | 2026-07-16 | `main` |

## 推荐入口与低层入口

日常使用 `run_auto_rnaseq.sh`，它负责扫描/下载 FASTQ、生成设计模板、组织日志和调用 Nextflow。只有已经准备好 samplesheet 且需要直接控制 Nextflow 时才使用 `run_pipeline.sh`。

注意参数名称不同：自动入口使用 `--results-dir`，低层入口使用 `--outdir`。

## 运行模式

| 选项 | 用途 |
|---|---|
| `--dry-run` | 只做预检，不提交任务 |
| `--upstream-only` | 只跑上游；默认行为 |
| `--downstream` | 读取已有 results，只跑差异与绘图 |
| `--background` | 后台运行并写入自动化日志 |
| `--resume` | 复用 Nextflow 已完成任务；默认 |
| `--no-resume` | 不复用缓存，谨慎使用 |

## 日志和状态

优先查看：

```text
results/_automation/logs/automation_*.log
results/_automation/logs/nextflow_*.log
results/02_gene_star/logs/
results/03_gene_featurecounts/logs/
results/06_te_star/logs/
```

`_automation/work/` 是 Nextflow task 工作目录，不是用户结果目录。排障时先看自动化日志和模块日志，不要从 work 目录随机寻找文件。

## 中断后恢复

相同输入、reference 和关键参数下，在同一结果目录重新运行并保留 `--resume`。若修改了 samplesheet、reference build、alignment 或关键计数参数，应先判断旧缓存是否仍有效；必要时使用新结果目录，避免新旧结果混合。

`--resume-session ID` 只用于明确知道目标 Nextflow session 的高级恢复场景。

## 只补部分下游

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/results \
  --species hg38 \
  --downstream \
  --only-tools Gene_featureCounts,TE_Telescope
```

跳过模块使用 `--skip-tools`。缺少样本时默认 `--partial-input-policy skip`，整模块跳过以避免生成残缺 panel；只有在分析设计明确允许部分样本时才使用 `allow`。

## 资源与并发

`--max-cpus` 和 `--max-memory` 是流程上限，不代表每个 task 都占用全部资源。`--plot-threads` 控制下游 contrast 并发，内存不足时优先降低并发，而不是删除输入或关闭 QC。

## 完成后归档

保留设计文件、resolved input、日志、版本信息、核心 counts、QC、差异表和最终图。确认不再需要 `--resume` 后再用清理工具检查 work 占用；不要手工删除正在运行任务的 work/cache。
