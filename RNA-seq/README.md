# RNA-seq Pipeline

bulk RNA-seq 的 FASTQ/manifest 自动化主流程、独立补跑工具和下游可视化。

```bash
bash rnaseq/run_auto_rnaseq.sh \
  --fastq-dir /path/to/project/fastq \
  --results-dir /path/to/project/results \
  --species hg38 \
  --strand reverse \
  --dry-run
```

完整文档：<https://quii490.github.io/Pipelines/pipelines/rnaseq/>

- 推荐主入口：`rnaseq/run_auto_rnaseq.sh`
- 独立补跑工具：`tools/`
- 批量下游可视化：`snakemake-visualization/`
- 低层 Nextflow 入口：`rnaseq/run_pipeline.sh`

不要向仓库提交 FASTQ、BAM、结果目录、日志或本地 reference 配置。
