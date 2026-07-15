# CUT&RUN Preflight、状态与完整性工具

| 状态 | 事实来源 | 最后验证 |
|---|---|---|
| Draft / under development | `CUTnRUN/tools/cutrun_cli/` | 2026-07-16 |

## 运行前检查：`cutrun_preflight.py`

低层 runner 会调用它；高级用户也可在提交长任务前单独运行：

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_preflight.py \
  --manifest /path/to/resolved_manifest.csv \
  --species hg38 \
  --config /path/to/run.config \
  --outdir /path/to/results \
  --json-out /path/to/results/preflight.json \
  --min-free-gb 100
```

| 参数 | 必需 | 说明 |
|---|---:|---|
| `--manifest` | 是 | 已解析 manifest；检查列、FASTQ 和 control mapping |
| `--species` | 否 | `hg38` 或 `mm39` |
| `--config` | 是 | 本次 reference/runner 配置 |
| `--outdir` | 是 | 用于磁盘空间和输出可写性检查 |
| `--json-out` | 是 | machine-readable 检查结果 |
| `--min-free-gb` | 否 | 最小剩余空间，默认 50 GB |
| `--skip-gzip` | 否 | 跳过完整 gzip integrity；首次运行不推荐 |
| `--repentools-index-dir` | 否 | 只有计划 RepEnTools 时提供 |
| `--repentools-gtf` | 否 | CHM13 RepeatMasker-derived GTF |

target 缺 control 可能被记录为 warning，因为 MACS3 技术上仍可运行；这不表示设计合理。没有合理 control 时必须在结果解释中说明限制。

## 结果完整性：`cutrun_results_summary.py`

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_results_summary.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv
```

它检查所需文件是否存在、非空并生成 machine-readable inventory。它回答“产物是否齐”，不回答“实验质量是否好”。一个完整但低信噪比的 run 仍可能科学上失败。

## Attempt-aware 状态：`cutrun_status.py`

该工具主要供 `run_downstream.sh` 和模块 wrapper 使用。普通用户读取 `module_status_latest.tsv` 即可，不需要手工改 ledger。

```bash
# 初始化一次 attempt ledger
python3 CUTnRUN/tools/cutrun_cli/cutrun_status.py init \
  --status-file /path/to/results/09_downstream/module_status.tsv \
  --run-id 20260716_120000 \
  --meta /path/to/results/09_downstream/status_meta.json

# wrapper 记录一个模块结果
python3 CUTnRUN/tools/cutrun_cli/cutrun_status.py record \
  --status-file /path/to/results/09_downstream/module_status.tsv \
  --run-id 20260716_120000 \
  --module peak_annotation --status PASS \
  --outputs-ok --output-paths "annotation/sample.tsv"

# 汇总本次 run
python3 CUTnRUN/tools/cutrun_cli/cutrun_status.py finalize \
  --status-file /path/to/results/09_downstream/module_status.tsv \
  --run-id 20260716_120000
```

有效模块状态为 `RUNNING/PASS/FAIL/SKIP`。`finalize` 会把必要模块 FAIL 视为整体 FAIL；存在未完成核心输出时可能是 INCOMPLETE。不要手工把 FAIL 改成 PASS，应该修复输入或依赖后记录新的 attempt。

## QC 汇总：`cutrun_te_qc.py`

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_te_qc.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --output /path/to/results/09_downstream/qc_metrics.tsv
```

不可计算指标应为 `SKIP/NA`，不能填 0。该工具的 FRiP 依赖 samtools/bedtools 可用性和现有 peaks；比较不同项目时核对 read、alignment 或 fragment 分母。
