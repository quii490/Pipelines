# CUT&RUN Workflow

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | CUT&RUN maintainers | 2026-07-16 | `main` |

## 主流程

| 阶段 | 目的 | 关键决定 |
|---|---|---|
| Input/preflight | 验证 FASTQ、manifest 与资源 | assay、control、assembly |
| Trim/alignment | 生成比对结果 | layout、aligner 参数 |
| Standard BAM | 高可信常规信号 | MAPQ、duplicates、blacklist |
| TE BAM | 保留更多重复区域证据 | multimapping 策略 |
| Locus-best BAM | 为每条 read 选择最佳 locus | tie/ambiguity 处理 |
| Tracks/peaks | 生成浏览与区域集合 | normalization、narrow/broad、control |
| Quantification | 生成 peak/feature matrices | counting unit 和 region set |
| Analysis/report | QC、比较、图形、模块状态 | metadata、contrast、成功定义 |

## 三类信号不能混称

- standard：适合常规 peak calling 和可复核的高可信轨迹；
- TE-aware：保留更多重复区域 reads，用于 TE 汇总分析；
- locus-best：提供 locus-level 可视化候选，但唯一位置解释仍受多重比对歧义限制。

报告图和表必须注明使用了哪一类 BAM/track，不能把三者作为等价技术重复合并。
