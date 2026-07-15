# RNA-seq QC、报告、清理与交付

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 |

## 已有 BAM 补 QC

使用 `RNA-seq/tools/run_qc_metrics_bam.sh --help`。BAM、refFlat、rRNA intervals、GTF/FASTA 必须匹配；QC 输出不改变 counts。

## 重建报告

```bash
rnaseq-report \
  --plot-dir /path/to/results/plots \
  --lfc-cutoff 0.58 \
  --padj-cutoff 0.05
```

可重建 `index.html`、DEG summary、LFC symmetry 与 module status。报告是索引，不替代原始差异表。

## 检查和清理 work

```bash
rnaseq-clean-work --results-dir /path/to/results
rnaseq-clean-work --results-dir /path/to/results --confirm
```

第一条只预览。只有确认没有任务运行、也不再需要 resume 时才 `--confirm`。清理后部分任务只能重跑。

## 生成交付目录

```bash
rnaseq-publish \
  --results-dir /path/to/results \
  --outdir /path/to/deliverables
```

交付应包含图、MultiQC、状态、版本和关键表，不包含 FASTQ、BAM、work、内部路径或真实身份信息。发布前人工打开报告并再次脱敏。

## 安装 wrapper

```bash
bash RNA-seq/tools/install_rnaseq_cli_tools.sh --help
```

wrapper 可能从一个 Conda 环境调用另一个 R 环境；归档时同时记录 wrapper 和实际 `Rscript` 环境。
