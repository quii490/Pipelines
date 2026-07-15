# 运行、监控与恢复 CUT&RUN

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

## 三层入口

| 入口 | 用途 |
|---|---|
| `run_auto_chipseq.sh` | 日常从 FASTQ 自动生成 manifest、preflight、运行和下游 |
| `run_pipeline.sh` | 已审查 manifest，直接控制 Nextflow/preview |
| `run_downstream.sh` | 核心 results 已存在，只恢复/重跑下游 |

日常优先第一层。复杂 control mapping 时，先 `--init-only`，编辑 manifest 后可用第二层 preview。

## 后台运行和监控

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/fastq \
  --species hg38 --assay cutrun \
  --outdir /path/to/results \
  --run-downstream --background --resume
```

```bash
cat /path/to/results/_automation/logs/automation.pid
tail -f /path/to/results/_automation/logs/automation_*.log
tail -f /path/to/results/_automation/logs/nextflow_*.log
```

## 中断恢复

输入、manifest、reference、关键参数和代码版本相同时，保留 `_automation/work` 并重跑原命令 `--resume`。若 Pod 被删除，先确认 results/work 位于持久化卷，再在相同代码/环境中恢复。

修改 genome、manifest control mapping、`te-k`、duplicate/blacklist 或主过滤策略时使用新的结果目录。不要混合旧 BAM、新 peaks 和旧报告。

## 只恢复 downstream

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_downstream.sh \
  --results-dir /path/to/results \
  --species hg38 \
  --resource-tier standard \
  --resume
```

下游是 attempt-aware：查看 `module_status.tsv` 历史和 `module_status_latest.tsv` 当前状态。失败模块可重跑，已通过模块复用。`--no-resume` 只用于确认 cache/sentinel 不应使用时。

## 停止任务

先确认 PID 属于本次流程，再发送 TERM：

```bash
kill "$(cat /path/to/results/_automation/logs/automation.pid)"
```

不要删除 work 目录来“停止”任务。停止后先查看首个失败/中断 task 的 `.command.err`、`.command.out` 和 `.command.sh`。
