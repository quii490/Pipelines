# RNA-seq Quick Start

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

```bash
cd RNA-seq/rnaseq

bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/project/fastq \
  --results-dir /path/to/project/results \
  --species hg38 \
  --strand reverse \
  --dry-run
```

确认预检后移除 `--dry-run`，可加 `--background --resume`。上游完成后检查并编辑：

```text
results/condition.csv
results/contrast.csv
```

再运行下游：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/project/results \
  --species hg38 \
  --downstream
```

成功标志：核心 BAM/count matrix 非空、MultiQC/模块日志无核心失败、下游输出与设计表样本一致。
