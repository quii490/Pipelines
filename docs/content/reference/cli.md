# CLI 参数参考

三个推荐入口：

```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
bash ATAC-seq/run_auto_atacseq.sh --help
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

共有概念包括 `--fastq-dir`、`--species`、结果目录、`--background` 和 `--resume`，但结果目录参数并不统一：RNA-seq 推荐入口使用 `--results-dir`，ATAC-seq/CUT&RUN 使用 `--outdir`。不要把低层 Nextflow 启动器参数直接套到自动化入口。

CLI `--help` 是参数事实来源；站点参数表提供解释和推荐值。
