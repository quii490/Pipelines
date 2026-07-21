# ZhaoLab Bioinformatics Pipelines

| 项目   | 内容         |
| ---- | ---------- |
| 状态   | Active     |
| 维护人  | CC         |
| 最后编辑 | 2026-07-16 |

这里是 ZhaoLab RNA-seq、ATAC-seq、CUT&RUN/CUT&Tag/ChIP-seq  pipelines的中文使用手册。命令和参数保留英文。
## 选择服务器
1. 推荐使用医学院集群，计算、存储资源充足，pipelines部署完整且经过检查。
	如何登录并使用医学院集群：
2. 如果数据量较小、只是做简单分析、医学院集群坏了😓，可以用实验室服务器。
## 先按手里的文件选择

| 你现在有             | RNA-seq                 | ATAC-seq                       | CUT&RUN/CUT&Tag/ChIP-seq      |
| ---------------- | ----------------------- | ------------------------------ | ----------------------------- |
| FASTQ            | 完整主流程                   | 完整主流程                          | 完整主流程                         |
| 已完成 results      | 只重跑 downstream/绘图       | 只重跑 downstream                 | `run_downstream.sh --resume`  |
| BAM              | gene/TE counts、QC       | 重建 peak/count/track/downstream | 恢复 peak、QC、track 或独立工具        |
| raw count matrix | DE、PCA、heatmap、pathway  | 差异可及性（需 regions + metadata）    | count_draw/差异与可视化             |
| DE/差异 region 表   | 火山图/MA图、注释、pathway/TE 图 | heatmap、annotation、motif、火山图   | annotation、motif、heatmap      |
| bigWig + BED     | correlation/profile     | TSS、region/TE heatmap          | correlation、region/TE heatmap |
## 选择 pipeline
看你做了啥实验🧪

| 生物学问题                                       | 推荐入口     | 从这里开始                                           |
| ------------------------------------------- | -------- | ----------------------------------------------- |
| bulk RNA 表达、gene/TE differential expression | RNA-seq  | [Quick Start](pipelines/rnaseq/quick-start.md)  |
| 染色质开放性、peaks、motif、footprinting             | ATAC-seq | [Quick Start](pipelines/atacseq/quick-start.md) |
| 蛋白/组蛋白结合、CUT&RUN、CUT&Tag、ChIP-seq           | CUT&RUN  | [Quick Start](pipelines/cutrun/quick-start.md)  |

## 新成员最短学习路径

1. 阅读[开始前检查](getting-started/index.md)和[命令行最小知识](getting-started/command-line.md)。
2. 进入目标 pipeline 的 Quick Start，先运行 `--help` 和 dry-run/preview。
3. 让入口生成 samplesheet/manifest，再人工核对 condition、replicate、contrast 和 control。
4. 正式运行后先看状态、QC 和样本关系，再看差异结果。
5. 需要从 BAM、matrix、bigWig 或已有 results 开始时，使用对应 pipeline 的 Running 页面和工具目录。

!!! warning "公开仓库"

    不要提交真实样本标识、密钥、内部服务器地址、私有资源位置或未发表结果。参见[数据安全](getting-started/security.md)。
