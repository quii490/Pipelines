# RNA-seq 主流程

日常使用 `run_auto_rnaseq.sh`，只有已经准备好 samplesheet 且需要直接控制 Nextflow 时才使用 `run_pipeline.sh`。

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/project/fastq \
  --results-dir /path/to/project/results \
  --species hg38 \
  --strand reverse \
  --dry-run
```

参数事实来源：

```bash
bash run_auto_rnaseq.sh --help
```

完整使用说明：<https://quii490.github.io/lab-pipelines/pipelines/rnaseq/>
