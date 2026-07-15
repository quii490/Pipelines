# 输出与归档参考

三个 pipeline 的具体目录分别见各自 [RNA-seq 输出](../pipelines/rnaseq/outputs.md)、[ATAC-seq 输出](../pipelines/atacseq/outputs.md) 与 [CUT&RUN 输出](../pipelines/cutrun/outputs.md)。

## 最小归档单元

- 输入 samplesheet/manifest、contrast 和脱敏后的 metadata；
- 完整命令、Git commit/release、环境与 reference 版本；
- run log、module status 和 QC；
- 可重用的 count/region/result 表；
- 关键图与生成图的脚本/参数；
- README，说明哪些模块被跳过或失败。

不要提交 FASTQ、BAM/CRAM、bigWig、Nextflow work、完整结果目录或含真实样本信息的报告到公开 Git 仓库。大型正式结果应存放在受控数据存储中。
