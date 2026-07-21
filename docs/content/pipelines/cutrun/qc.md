# CUT&RUN Quality Control

| 指标 | 主要问题 | 异常时检查 |
|---|---|---|
| FASTQ/trim | 原始 reads 是否可靠 | adapter、read quality、批次 |
| Alignment | 是否匹配正确 assembly | index、污染、read length |
| Duplicate/complexity | 独立分子是否足够 | 起始量、PCR、测序深度 |
| Fragment size | 是否符合 assay 特征 | layout、实验与过滤 |
| FRiP / signal-to-noise | 信号是否集中于 peaks | control、peak type、背景 |
| Peak number/width | 与 target 类型是否一致 | narrow/broad、cutoff、control |
| Replicate correlation | 生物学重复是否一致 | 标签、批次、异常样本 |
| Control background | IgG/input 是否异常 | 配对、实验背景、归一化 |

TF、narrow histone mark、broad histone mark 不能共用一个机械阈值。先确认指标定义和单位，再结合 track、peak 和重复一致性判断。

TE-aware 信号增加可能同时来自重复序列保留与背景增加；必须与 standard track、control 和 mappability 一起解释。
