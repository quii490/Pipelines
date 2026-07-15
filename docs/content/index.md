# Lab Bioinformatics Pipelines

| 项目 | 内容 |
|---|---|
| 状态 | Active |
| 维护人 | Lab pipeline maintainers |
| 最后审查 | 2026-07-15 |
| 文档版本 | `main` |

可复现的 RNA-seq、ATAC-seq 与 CUT&RUN/CUT&Tag/ChIP-seq 分析流程。本网站回答四个问题：应选择哪个流程、怎样准备输入、如何运行、结果是否可信。

## 选择流程

| 数据或问题 | 推荐流程 | 首先阅读 |
|---|---|---|
| bulk RNA 表达、gene/TE differential expression | RNA-seq | [RNA-seq Quick Start](pipelines/rnaseq/quick-start.md) |
| 染色质开放性、peaks、motif、footprinting | ATAC-seq | [ATAC-seq Quick Start](pipelines/atacseq/quick-start.md) |
| 蛋白/组蛋白结合、CUT&RUN、CUT&Tag、ChIP-seq | CUT&RUN | [CUT&RUN Quick Start](pipelines/cutrun/quick-start.md) |

!!! warning "公开仓库"
    不要在命令、截图、日志或配置中提交真实样本信息、密钥、内部服务器地址或私有资源位置。参见[数据安全](getting-started/security.md)。

## 推荐学习顺序

1. 阅读[入门指南](getting-started/index.md)并确认计算资源。
2. 进入目标 pipeline 的 Quick Start。
3. 先初始化输入模板，再核对实验设计。
4. 先看 QC 和汇总报告，再解释差异结果。
5. 遇到问题先查 pipeline 专属页，再查[通用故障排查](troubleshooting/index.md)。
