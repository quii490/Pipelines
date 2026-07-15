# CUT&RUN Quick Start

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

```bash
cd CUTnRUN/pipelines/chipseq_auto_nf

./run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --assay cutrun \
  --outdir /path/to/project/results \
  --init-only
```

检查 manifest 的 group、replicate、control/IgG 后预览：

```bash
./run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --outdir /path/to/project/results \
  --preview
```

正式运行：

```bash
./run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --outdir /path/to/project/results \
  --resource-tier standard \
  --background --resume --run-downstream
```
