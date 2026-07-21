# 通用运行流程

```mermaid
flowchart LR
  A[确认 assay 和设计] --> B[检查 FASTQ/参考/磁盘]
  B --> C[生成输入模板]
  C --> D[人工检查分组和 control]
  D --> E[dry-run/preview/preflight]
  E --> F[正式运行]
  F --> G[状态与 QC]
  G --> H[统计和专题分析]
  H --> I[归档和发布]
```

## 三种常见起点

### 从 FASTQ 开始

使用自动入口生成输入表。人工核对后再正式运行，这是信息最完整、最推荐的方式。

### 从 pipeline results 开始

保留原结果和 work/cache；仅修改 metadata、contrast 或绘图参数时，使用 downstream-only 入口。输出到新目录，避免覆盖原始结果。

### 从中间文件开始

先确认文件是哪个 genome build、怎样过滤、是否 deduplicate、PE/SE 单位、是否包含 multi-mapping、使用何种 normalization。然后在工具目录中选择 BAM→counts、BAM→track/peak、matrix→DE 或 bigWig→profile。中间文件缺少 provenance 时，宁可重跑上游。

## 何时使用 resume

输入、关键参数、代码版本和 work 目录不变时使用 `--resume`。如果改变 genome build、manifest、aligner、多重比对策略或主过滤规则，应使用新的结果目录。不同策略不要混写到同一个目录。
