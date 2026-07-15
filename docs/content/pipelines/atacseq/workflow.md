# ATAC-seq Workflow

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

| 阶段 | 主要输入 | 主要产物 | 核心检查 |
|---|---|---|---|
| FASTQ QC/trim | FASTQ | clean FASTQ、QC | reads、adapter、质量 |
| Alignment | clean FASTQ | BAM | mapping、proper pairs |
| Filtering | BAM | analysis-ready BAM | mito、blacklist、MAPQ、duplicates |
| Signal | filtered BAM | bigWig | normalization、track 连续性 |
| Peak calling | filtered BAM | peak files | peak 数量与信噪比 |
| Quantification | peaks/bins + BAM | count matrix | PE fragment/SE read 单位 |
| Downstream | matrix + metadata | PCA、correlation、DA | 分组、批次、重复一致性 |

PE blacklist filtering 必须按 read name 同步移除 mates，避免产生人为 orphan reads。PE featureCounts 和 FRiP 使用 fragment 单位；SE 使用 read 单位。

## 特殊模块

- nucleosome phasing 主要适用于 PE；低深度或不可计算时应报告 `NA`，不能伪造数值；
- TSS profile 描述的是 bigWig 信号形状，不等同于严格 Tn5 insertion TSS enrichment；
- relaxed TE track 用于浏览和 metaprofile，不能单独证明某个 TE locus 的可及性变化。

每个结论应能回溯到输入 BAM、过滤策略、归一化方式和生成脚本版本。
