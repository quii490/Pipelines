# 通用故障排查

1. 记录完整命令、Git commit、输入表和失败时间。
2. 找到第一个失败模块，而不是最后一条汇总错误。
3. 检查文件存在性、权限、磁盘、reference build 和 Conda/Nextflow 版本。
4. 修复原因后在原结果目录使用 `--resume`。
5. 若核心结果成功但可选模块失败，明确记录 SKIP/FAIL，不伪造空结果。

专题问题见 [RNA-seq](../pipelines/rnaseq/troubleshooting.md)、[ATAC-seq](../pipelines/atacseq/troubleshooting.md)和 [CUT&RUN](../pipelines/cutrun/troubleshooting.md)。
