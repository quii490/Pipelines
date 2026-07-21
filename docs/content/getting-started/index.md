# 入门指南

## 第一次运行前检查

- [ ] 知道数据属于 RNA-seq、ATAC-seq、CUT&RUN/CUT&Tag/ChIP-seq；
- [ ] 知道物种和选择基因组（人类选hg38，小鼠推荐选mm10，也可选mm39）；
- [ ] 知道单双端测序、建库链特异性 strandedness（RNA-seq）和实验condition、contrast、replicate设计；
- [ ] 检查数据命名是否合理、数据目录结构合理；
- [ ] FASTQ 完整、结果目录可写、磁盘空间足够、运行资源足够；
- [ ] `--help` 可执行，参考配置属于同一 genome build；
- [ ] 已运行 dry-run/preview/preflight；
- [ ] 长任务使用可恢复的 work 目录，完整命令已保存。

如果其中任何一项不清楚，先解决设计或输入问题，不要靠放宽参数“让流程跑起来”。
## 推荐阅读顺序
**建议看一眼激活环境、样本命名即可去看目标pipeline的Quick Start**

1. [激活环境](installation.md)
2. [命令行最小知识](command-line.md)
3. [样本命名与目录](naming.md)
4. [参考基因组](references.md)
5. [通用运行流程](common-workflow.md)
6. 目标 pipeline 的 Quick Start
