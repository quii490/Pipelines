# CUT&RUN 参数选择

| 状态 | 维护人 | 最后审查 | 事实来源 |
|---|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 | `run_auto_chipseq.sh --help` |

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

## 输入与项目

| 参数 | 默认 | 什么时候改 |
|---|---|---|
| `--fastq-dir` | 必需 | FASTQ 根目录，递归扫描 |
| `--species` | `hg38` | 当前支持 `hg38/mm39`；必须匹配全部参考 |
| `--assay` | `cutrun` | `chipseq/cutrun/cuttag` |
| `--layout` | `auto` | 自动判断错误时明确设 `PE/SE` |
| `--outdir` | FASTQ 同级结果目录 | 建议显式指定可写持久化目录 |
| `--work-dir` | `<outdir>/_automation/work` | 需要独立高速持久化 work 时改 |
| `--manifest` | 自动生成路径 | 使用人工审查的 manifest |

## Manifest 与 control

| 参数 | 用途 |
|---|---|
| `--control-sample NAME` | 全局明确 control，适合单 control 设计 |
| `--control-regex REGEX` | 从样本名识别 control |
| `--no-auto-control` | 多 control/复杂设计时禁止自动填 `igg` |
| `--force-manifest` | 重新生成并覆盖；使用前备份人工 manifest |

## 执行控制

| 参数 | 说明 |
|---|---|
| `--init-only` | 只生成 manifest |
| `--preview` | Nextflow preview，不运行 jobs |
| `--dry-run` | 只打印命令 |
| `--skip-preflight` | 跳过 FASTQ/reference/disk 检查，不建议正式分析 |
| `--background` | 后台运行并保存 PID/log |
| `--resume/--no-resume` | 复用/不复用 work；日常使用 resume |
| `--profile` | Nextflow profile，默认 conda |
| `--resource-tier` | `small/standard/full`，控制下游 heatmap 规模与分辨率 |

## 核心与 tracks

| 参数 | 默认 | 解释 |
|---|---:|---|
| `--trim` | `true` | 输入已可靠 trimming 才关闭 |
| `--run-analysis` | `true` | 关闭后不生成完整 count_draw 分析 |
| `--make-rpgc-track` | `true` | standard RPGC track |
| `--make-te-tracks` | `true` | TE-aware track |
| `--make-te-locus-best-track` | `true` | 5 bp CPM best-locus track |
| `--run-downstream` | 关闭 | 需要 correlation/annotation/motif/heatmap/report 时显式开启 |

## TE 策略

| 参数 | 默认 | 风险与选择 |
|---|---:|---|
| `--te-k` | `25` | 保留的多重定位数；增大 BAM/计算量，改后用新目录 |
| `--te-duplicate-policy` | `mark` | `mark/keep/remove`；改变有效信号，需实验理由 |
| `--te-remove-blacklist` | `false` | 开启会改变 TE counts/track，与旧 run 不可直接混比 |
| `--te-methods` | 空 | `t3e,allo,repentools` 列表 |
| `--te-methods-execute` | 关闭 | 只在外部环境、输入和资源验证后执行 |

T3E 参数包括 `--t3e-dir`、`--t3e-python`、`--t3e-max-bed-reads`、`--t3e-iterations`；Allo 使用 `--allo-command`；RepEnTools 使用 index/GTF/ret 和 target/input mapping 参数。外部路径放在本地配置，不提交仓库。

TE annotation 覆盖包括 `--te-repeat-bed`、`--te-saf` 和 `--te-anno`。RepEnTools 高级参数为 `--repentools-command`、`--repentools-index-dir`、`--repentools-gtf`、`--repentools-ret`、`--repentools-target-groups` 和 `--repentools-input-samples`。

## MACS

| 参数 | 说明 |
|---|---|
| `--macs-qvalue` | q-value cutoff，越小越严格 |
| `--macs-pvalue` | 设置后使用 p-value 而非 q-value |
| `--macs-broad-cutoff` | broad peak cutoff，越小越严格 |

MACS 内部对应 `--broad-cutoff`；用户入口应使用 `--macs-broad-cutoff`。

不要为了增加 peak 数量而反复放宽 cutoff。narrow 和 broad 回答的区域定义不同，应分别保存和报告。

## 低层 runner 的过滤参数

`run_pipeline.sh` 还提供：

| 参数 | 默认 | 选择原则 |
|---|---:|---|
| `--min-mapq` | `30` | standard 分支；降低会增加 ambiguous/background signal |
| `--min-frag` | `30` | 过短 fragment 下限；需结合 assay/layout |
| `--max-frag` | `1200` | 文库片段明显更长且 QC 支持时调整 |
| `--extend-reads-se` | `250` | SE read extension；影响 signal/peak |
| `--ignore-for-normalization` | `chrX` | track normalization 排除 contig；需按实验性别/物种审查 |
| `--download-threads` | `8` | 环境下载并发 |

自动入口没有单独暴露时，可经 `--extra` 传递，但必须确认名称属于 Nextflow/低层 runner 并保存实际命令。任何过滤参数改变都应使用新结果目录。

`--extra "ARGS"` 会继续传给低层 runner/Nextflow，只适合明确知道参数归属的维护者。
