# 入门指南

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-16 | `main` |

## 第一次运行前检查

- [ ] 知道数据属于 RNA-seq、ATAC-seq、CUT&RUN、CUT&Tag 还是 ChIP-seq；
- [ ] 知道物种和 genome build；
- [ ] 知道 single-end/paired-end、建库 strandedness（RNA-seq）和 control 设计；
- [ ] biological replicate、condition、batch 和 contrast 已写入 metadata；
- [ ] FASTQ 完整、结果目录可写、磁盘空间足够；
- [ ] `--help` 可执行，参考配置属于同一 genome build；
- [ ] 已运行 dry-run/preview/preflight；
- [ ] 长任务使用可恢复的 work 目录，完整命令已保存。

如果其中任何一项不清楚，先解决设计或输入问题，不要靠放宽参数“让流程跑起来”。

## 推荐阅读顺序

1. [系统与安装](installation.md)
2. [命令行最小知识](command-line.md)
3. [样本命名与目录](naming.md)
4. [参考基因组](references.md)
5. [通用运行流程](common-workflow.md)
6. 目标 pipeline 的 Quick Start
