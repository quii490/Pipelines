# Counts 与信号标准化

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | Analysis maintainers | 2026-07-16 |

“标准化”不是一个统一步骤。不同方法解决不同偏差，不能因为单位看起来相似就互换。

## RNA expression

| 表示 | 主要用途 | 不适合 |
|---|---|---|
| raw counts | DESeq2/edgeR 等差异模型 | 直接跨样本比较未校正 library size |
| DESeq2 normalized counts / edgeR CPM | 样本间展示、QC | 替代模型的 raw input |
| VST/rlog/logCPM | PCA、correlation、heatmap | 作为 raw count DE input |
| TPM | 同样本内相对 abundance、部分展示 | 直接做标准 count-based DE |
| FPKM/RPKM | 长度和深度校正的旧式表达单位 | 跨项目无条件比较、count model |

TPM 的和在每个样本近似固定；它校正 feature length，但 composition bias 仍可能影响解释。TE length、multi-mapping 和 annotation 定义使 TE TPM 更复杂。

## Genomic signal tracks

| 方法 | 概念 | 使用注意 |
|---|---|---|
| CPM | 每百万 mapped reads/fragments | 不校正有效基因组大小或 input |
| BPM | bins per million，deepTools 定义 | 先确认工具实现 |
| RPKM | CPM 再按 bin/region kb 调整 | 依赖 bin length |
| RPGC | 1× genome coverage 风格 | 需要 effective genome size |
| fold over control / log2 ratio | target 相对 input/control | pseudocount、scale 和负值处理影响图形 |

只比较相同 assembly、过滤、read/fragment 单位、normalization、effective genome size 和 bin size 的 bigWig。

## Peak/count analysis

ATAC/CUT&RUN 的 region count matrix应保留 raw fragment/read counts进入模型；library-size/composition normalization由模型完成。FRiP、peak counts 和 track intensity不是彼此替代的 normalization。

## Batch correction

统计设计中加入 batch 与在可视化 matrix 上去 batch 是两件事。不要把经过 ComBat/`removeBatchEffect` 的数值作为 DESeq2 raw counts。若 batch 与 condition 完全混杂，计算方法无法恢复缺失的实验设计信息。

## 选择流程

```text
差异统计？ → raw counts + replicate-aware model
PCA/heatmap？ → VST/logCPM 或一致的 signal summary
浏览器轨迹？ → 同一 CPM/RPGC 等 track normalization
跨方法/跨项目？ → 先验证 feature、reference、过滤和单位完全可比
```
