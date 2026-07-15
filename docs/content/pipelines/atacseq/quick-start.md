# ATAC-seq Quick Start

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

先生成输入模板：

```bash
cd ATAC-seq
bash run_auto_atacseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --outdir /path/to/project/results \
  --init-only
```

检查 `_automation/inputs/samplesheet.csv` 和 `contrasts.csv`，再运行：

```bash
bash run_auto_atacseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --outdir /path/to/project/results \
  --preset standard \
  --background --resume
```

成功后先看 `QC_REPORT.md`、样本级 QC、BAM/bigWig 和 peak，再解释差异结果。
