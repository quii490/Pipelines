# ATAC-seq 参数

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

CLI 的事实来源始终是：

```bash
bash ATAC-seq/run_auto_atacseq.sh --help
```

## 输入与运行

| 参数 | 选择/默认 | 说明 |
|---|---|---|
| `--mode` | `auto/init/upstream/downstream` | 执行阶段 |
| `--fastq-dir` | 上游必需 | FASTQ 根目录 |
| `--species` | `hg38/mm10/mm39` | 参考基因组 |
| `--layout` | `auto/PE/SE` | 测序布局 |
| `--outdir` | 路径 | 结果目录 |
| `--work-dir` | 路径 | Nextflow 工作目录 |
| `--samplesheet` | 路径 | 使用已有 samplesheet |
| `--init-only` | flag | 只生成/检查输入 |
| `--background` | flag | 后台运行 |
| `--resume` / `--no-resume` | 推荐 resume | 是否复用缓存 |
| `--extra` | 字符串 | 透传高级参数；必须记录 |

## 资源与预设

| 参数 | 说明 |
|---|---|
| `--preset quick` | 关闭 TSS/gene-body/TE heatmap、fixed-bin、motif、footprinting、nucleosome phasing；用于输入与核心流程测试 |
| `--preset standard` | 保留当前标准模块组合；正式项目的起点 |
| `--preset full` | 在 standard 上开启 fixed-bin、HOMER motif、TOBIAS footprinting 和 nucleosome phasing；依赖和计算更多 |
| `--profile-preset off/quick/standard/full` | 信号 profile 预设 |
| `--profile-cores` | profile 计算线程 |
| `--queue-size` | 最大排队任务数 |
| `--max-cpus` / `--max-memory` | 全局资源上限 |
| `--fastp-max-forks` / `--fastp-timeout` | fastp 并发与超时 |

## Downstream

| 参数 | 说明 |
|---|---|
| `--metadata-csv` | condition/replicate 表 |
| `--contrast CASE,CONTROL` | 可重复指定的比较 |
| `--contrast-file` | 批量比较文件 |
| `--overwrite-inputs` | 覆盖生成的下游输入；使用前备份人工修改 |
| `--levels peak/bin/both` | peak、fixed-bin 或两者 |
| `--downstream-only` | 仅下游兼容开关 |

## TE track

| 参数 | 说明 |
|---|---|
| `--run-te-relaxed-tracks BOOL` | 是否生成 relaxed TE tracks |
| `--te-mapq` | relaxed BAM/track 的 MAPQ 阈值 |
| `--te-track-normalization CPM/RPGC/RPKM/BPM` | track 归一化 |
| `--te-bw-binsize` | bigWig bin size |

## 低层参数：何时才需要改

自动入口可用 `--extra` 传给 `run_pipeline.sh`。常见低层参数：

| 参数 | 当前默认 | 选择原则 |
|---|---:|---|
| `--mapq` | `30` | standard BAM 的可信 mapping 阈值；降低会增加背景/重复区域 reads |
| `--track-normalization` | `RPGC` | standard track；跨样本必须一致 |
| `--consensus-half-width` | `250` | summit ±250 bp；改变 region universe 和 counts |
| `--call-broad` | false | 常规 ATAC peaks 多为 narrow；特定问题才开启 |
| `--fixedbin-size` | `100000` | 更小分辨率高但 feature/计算量大 |
| `--run-markdup` | true | duplicate 处理属于标准清理的一部分 |
| `--remove-mito` / `--remove-blacklist` | true | 改变 clean BAM 定义，跨 run 必须一致 |

示例：

```bash
--extra "--fixedbin_size 50000"
```

透传的实际 Nextflow 参数名以低层 `run_pipeline.sh --help`/配置为准；不确定时不要使用 `--extra`。

!!! danger "高风险参数"

    降低 MAPQ、启用 relaxed tracks 或改变 normalization 会改变信号解释。必须在方法和图注中写明；relaxed track 不代表 EM/fractional multi-mapper 分配。
