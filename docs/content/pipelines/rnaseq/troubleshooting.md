# RNA-seq 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 | `main` |

先保存完整命令、日志、设计文件和 Git commit。定位第一个失败模块，不要只复制最后一行 Nextflow 汇总错误。

## 无法写 automation log 或 PID

**现象：** `Permission denied`，流程在提交任务前退出。

**诊断：** 检查 `--results-dir`、其父目录和已有 `_automation/` 的 owner/permissions。

**处理：** 使用当前用户可写的新结果目录，或由管理员修复权限。不要用 `sudo` 运行整个 pipeline。

## STAR index 无法打开或 mapping 很低

检查 index 前缀文件是否完整、物种/build 是否正确、FASTQ 是否污染或损坏。对比同批次样本的 STAR summary。不要混用 hg38、T2T、mm10 和 mm39。

## featureCounts assigned rate 很低

依次核对 BAM/GTF chromosome names、build、`--strand`、feature type 和 gene attribute。若 mapping 正常但 assigned 极低，strandedness 或 annotation 不一致是优先嫌疑。

## TE counts 异常低或为空

检查 TE annotation build、工具输入 BAM 是否保留所需 multi-mapping 信息、strandedness 和模块日志。不要用标准唯一比对 BAM 静默替代 TE 分支输入。

TE STAR 应允许足够的 multi-mapping；把 `outSAMmultNmax` 限成 1 会破坏下游 TE 工具所需信息。Telescope、TEtranscripts、REdiscoverTE、TElocal 的统计层级不同，不能仅按 counts 大小判断哪个“正确”。

## duplicate rate 高

普通 bulk RNA-seq 默认只用 `--run-markdup-qc true` 记录 Picard metrics，不用 dedup BAM counts。高 duplicate 可能来自高表达基因、低复杂度或 PCR；结合 library complexity、表达集中度和建库信息判断。不要为降低 duplicate 指标默认开启 `--run-dedup true`。

## 下游没有差异结果

检查：

- `condition.csv` 样本是否与 matrix 列名完全一致。
- `contrast.csv` 的 case/control 是否存在。
- 是否至少有合理 replicates。
- 模块是否因 partial input policy 被 SKIP。
- 过滤后是否没有足够 feature，而不是脚本失败。

无显著 feature 是合法结果，不应通过放宽阈值制造显著性。

## DE 表没有 p-value/padj

常见原因是每组没有 biological replicates，流程进入探索模式。确认 sample table 的 replicate/condition，而不是把 technical lanes 当重复。探索性 logCPM difference/fixed BCV 只能用于候选排序。

## `--only-tools` 没有运行预期模块

模块名区分大小写并随当前结果命名。先查看下游状态/已有 matrix 名称，再使用当前文档列出的 `gene_featureCounts`、`TE_TEtranscripts` 等。被 `partial-input-policy=skip` 跳过时先修复缺失样本。

## MultiQC 没有生成

确认 `--run-multiqc true`、`--multiqc-cmd` 可执行，并检查上游是否产生 MultiQC 能识别的日志。必要时从已有结果目录单独补跑，但记录使用的 MultiQC 版本和搜索目录。

## 样本不完整

当 condition 有 16 个样本但某 matrix 只有 9 个时，默认整模块跳过。先找缺失发生在哪个上游步骤；只有实验设计明确只分析子集时才允许 partial input，并在结果中列出实际样本。

## Nextflow task 失败后如何恢复

修复输入、资源或环境后，在相同结果/work 位置使用 `--resume`。若修改 reference 或关键分析参数，使用新结果目录或明确的新 session，避免复用不兼容缓存。

## 内存不足或任务被杀

检查 task `.command.err`、scheduler/容器 OOM 状态和 Nextflow trace。先降低 `--queue-size`/`--plot-threads` 或提高对应 process memory；不要同时重复启动多个后台 run。TE 和 MultiQC 扫描大量文件时也会增加 I/O。

## 报告问题时提供

```text
入口命令（删除敏感路径）
Git commit
失败模块和首个错误
对应 automation/Nextflow/module 日志片段
输入表头和匿名示例行
物种、reference build、layout、strand
```
