# CUT&RUN / CUT&Tag / ChIP-seq Pipeline

从 FASTQ 生成标准与 TE 分支的 BAM、bigWig、peaks、counts、QC 和下游报告。

```bash
bash pipelines/chipseq_auto_nf/run_auto_chipseq.sh \
  --fastq-dir /path/to/project/fastq \
  --species hg38 \
  --assay cutrun \
  --outdir /path/to/project/results \
  --init-only
```

完整文档：<https://quii490.github.io/Pipelines/pipelines/cutrun/>

> 当前 CUT&RUN pipeline 和文档仍在开发与验证中；请将其视为 Draft，正式项目使用前审查实际代码、manifest、参数和输出。

- 推荐入口：`pipelines/chipseq_auto_nf/run_auto_chipseq.sh`
- 主流程：`pipelines/chipseq_auto_nf/`
- 独立工具：`tools/cutrun_cli/`
- 可选 TE 方法：`tools/te_methods/`

正式运行前检查 manifest 中的 target、replicate 和 control/IgG。不要提交内部 reference 路径或真实数据。
