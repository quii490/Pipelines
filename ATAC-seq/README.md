# ATAC-seq Pipeline

bulk ATAC-seq 从 FASTQ 到 clean BAM、bigWig、peaks、QC、差异可及性和专题分析。

```bash
bash run_auto_atacseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --outdir /path/to/project/results \
  --init-only
```

完整文档：<https://quii490.github.io/Pipelines/pipelines/atacseq/>

- 推荐入口：`run_auto_atacseq.sh`
- 低层 Nextflow 入口：`run_pipeline.sh`
- 独立补跑工具：`scripts/`
- 共享参考模板：`conf/`

开始正式运行前必须检查生成的 samplesheet 和 contrasts。不要提交测序数据或分析结果。
