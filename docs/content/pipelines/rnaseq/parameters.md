# RNA-seq 参数参考
本页解释推荐值和风险边界；当前参数事实来源始终是：
```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
```
## 输入、路径和设计

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--manifest FILE` | — | accession、URL 或本地 FASTQ manifest |
| `--fastq-dir DIR` | — | 本地 FASTQ 根目录 |
| `--results-dir DIR` | 自动推断 | 上游主结果目录 |
| `--work-root DIR` | 自动推断 | 下载和自动化工作根目录 |
| `--download-dir DIR` | `<work-root>/fastq` | manifest 下载 FASTQ 位置 |
| `--plot-outdir DIR` | `<results>/plots` | 下游绘图目录；`--plot-dir` 为别名 |
| `--sample-metadata FILE` | — | `sample,condition,replicate` |
| `--sample-meta` | — | `sample:condition:replicate`，可重复 |
| `--condition-map` | — | `condition:sample1,sample2`，可重复 |
| `--contrast-file FILE` | — | 推荐列 `case,control` |
| `--contrast CASE:CONTROL` | — | 命令行 contrast，可重复 |
| `--replace-design` | false | 覆盖设计前自动备份；谨慎使用 |
| `--default-condition` / `--default-replicate` | `NA` | metadata 未命中时的占位；正式设计不应依赖 NA |

## Reads 与 reference

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--species` | `hg38` | `hg38/mm10/mm39` |
| `--human-ref` | `hg38` | `hg38/t2t`，仅人类分支 |
| `--aligner` | `star` | `star/hisat2` |
| `--strand` | `reverse` | `unstranded/forward/reverse` |
| `--strandedness` / `--tecount-strand` | — | `--strand` 的兼容别名；新命令统一用 `--strand` |
| `--layout` | `auto` | `auto/PE/SE` |
| `--r1-pattern` | `_R1` | R1 识别模式 |
| `--r2-pattern` | `_R2` | R2 识别模式 |
| `--sra-threads` | `8` | accession 下载转换线程 |

## 运行控制

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--background` | false | 后台运行整个自动化流程 |
| `--resume` | true | 复用任务缓存 |
| `--no-resume` | false | 禁止缓存复用 |
| `--resume-session` | — | 指定 Nextflow session；高级恢复 |
| `--upstream-only` | true | 只跑上游 |
| `--downstream` | false | 只跑下游 |
| `--downstream-only` | false | `--downstream` 的兼容别名 |
| `--profile NAME` | 配置决定 | Nextflow profile |
| `--queue-size` | 配置决定 | Nextflow queue size |
| `--max-cpus` | 配置决定 | 最大 CPU |
| `--max-memory` | 配置决定 | 最大 memory，如 `64 GB` |
| `--dry-run` | false | 预检和计划，不提交任务 |
| `--failure-policy` | `core` | `core/strict`；可选模块失败策略 |
| `--progress-interval` | `30` | 自动化进度文件刷新秒数 |

## 上游模块

| 参数 | 当前默认 | 说明 |
|---|---:|---|
| `--run-fastp` | true | trimming/QC |
| `--run-fastqc` | true | 原始 FASTQ FastQC |
| `--run-gene-count-branch` | true | gene alignment/counts |
| `--run-salmon` | true | transcript quantification |
| `--run-stringtie` | false | StringTie 分支 |
| `--run-tecount` | false | TEcount |
| `--run-telocal` | false | TElocal，按需开启 |
| `--run-tetranscripts` | true | TEtranscripts |
| `--run-telescope` | true | locus-level EM allocation |
| `--run-rediscoverte` | `auto` | hg38 自动运行，其他情况跳过 |
| `--run-rediscoverte-rollup` | true | 生成分层 TE 汇总 |
| `--run-salmonte` | false | 实验性兼容模块 |
| `--run-dedup` | false | 用 dedup BAM counts；普通 RNA-seq 不推荐 |
| `--run-markdup-qc` | true | 只输出 duplicate metrics |
| `--run-multiqc` | true | 汇总 QC |
| `--run-rnaseq-metrics` | true | Picard RNA-seq metrics |

## 下游阈值和绘图

| 参数 | 默认值 | 说明 |
|---|---:|---|
| `--padj-cutoff` | `0.05` | FDR cutoff |
| `--lfc-cutoff` | `0.58` | absolute log2FC cutoff |
| `--baseMean-min` | `5` | 低表达过滤 |
| `--label-top-n` | `40` | 标注数量 |
| `--heatmap-top-n` | `40` | 热图特征数 |
| `--volcano-orientation` | `classic` | `classic/horizontal` |
| `--gray-nonsig` | true | 非显著点统一灰色 |
| `--plot-threads` | 1 或资源上限 | contrast 并发数 |
| `--only-tools` | — | 只运行列出的下游模块 |
| `--skip-tools` | — | 跳过列出的模块 |
| `--partial-input-policy` | `skip` | `skip/error/allow` |
| `--allow-partial-inputs` | false | 旧兼容参数；true 等价于 policy=allow |

这些阈值是当前软件默认而非普遍生物学标准。报告应同时保留连续的 effect size、raw P value/FDR 和过滤前结果。

## 无重复探索性分析

`--exploratory-method` 支持 `logCPM_diff` 或 `edgeR_fixedBCV`，后者由 `--exploratory-fixed-bcv` 提供固定 BCV。此类结果必须标记为探索性，不能声称得到了可靠的组内变异估计。

## 环境和高级覆盖

`--rnaseq-env`、`--downstream-env`、`--rediscoverte-rollup-conda-prefix`、`--multiqc-cmd` 用于提供本地环境位置；公开仓库不保存真实内部路径。

Picard 可用 `--ref-flat`（兼容别名 `--refFlat`）、`--ribosomal-intervals` 和 `--picard-markdup-java-heap` 覆盖自动资源。文件必须与 GTF/reference 匹配，JVM heap 要小于 task memory。

`--extra "ARGS"` 会将原始参数传给 Nextflow，只适合了解低层配置的维护者。
