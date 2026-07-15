# CUT&RUN / CUT&Tag / ChIP-seq Quick Start

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 | `run_auto_chipseq.sh --help` |

以下命令从仓库根目录执行。最安全顺序是：**生成 manifest → 人工修正 control → preview → 正式运行 → 检查状态和报告。**

## 1. 查看入口

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

## 2. 检查 FASTQ 和磁盘

```text
/path/to/project/fastq/
├── H3K27ac_WT_1_R1.fastq.gz
├── H3K27ac_WT_1_R2.fastq.gz
├── IgG_WT_1_R1.fastq.gz
└── IgG_WT_1_R2.fastq.gz
```

每个 PE 样本必须有 R1/R2。建议预留明显高于 FASTQ 总量的空间；`--te-k 25` 的 TE BAM 可能很大。

## 3. 生成 manifest

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --assay cutrun \
  --outdir /path/to/project/results \
  --init-only
```

默认 manifest 位于：

```text
results/_automation/inputs/manifest.csv
```

逐行核对 sample、FASTQ、group/target、replicate、`is_igg`/control 配对。多个 IgG 或复杂设计不要依赖文件名自动推断。

## 4. Preview

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --manifest /path/to/project/results/_automation/inputs/manifest.csv \
  --species hg38 \
  --assay cutrun \
  --outdir /path/to/project/results \
  --preview
```

不要日常使用 `--skip-preflight`。preflight 对 FASTQ、manifest、reference 和磁盘的错误应在长任务前修复。

## 5. 正式运行

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --manifest /path/to/project/results/_automation/inputs/manifest.csv \
  --species hg38 \
  --assay cutrun \
  --outdir /path/to/project/results \
  --resource-tier standard \
  --run-downstream \
  --background \
  --resume
```

`--run-downstream` 必须显式开启，才会在主流程后运行 correlation、annotation、motif、heatmap 和总报告。外部 T3E/Allo/RepEnTools 不属于默认必需模块。

## 6. 监控

```bash
cat /path/to/project/results/_automation/logs/automation.pid
tail -f /path/to/project/results/_automation/logs/automation_*.log
```

## 7. 完成后先看四个文件

```text
09_downstream/run_report.html
09_downstream/module_status_latest.tsv
09_downstream/module_status_summary.json
09_downstream/run_manifest.json
```

- 核心 BAM、counts 或 peaks 缺失时属于 `INCOMPLETE`；
- `FAIL` 要查日志；
- `SKIP` 表示不适用或条件不满足，不是“没有信号”；
- standard、TE-aware 和 locus-best track 不能混作同一种定量证据。

下一步读[输出指南](outputs.md)和[QC](qc.md)。
