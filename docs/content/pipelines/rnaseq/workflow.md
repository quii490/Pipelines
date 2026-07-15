# RNA-seq Workflow

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

核心阶段包括：FASTQ 完整性和质量控制、clean reads、gene alignment、gene counting、转录本定量、TE 多重比对/定量、MultiQC、差异分析和图形报告。

普通 gene expression 与 TE analysis 对 multi-mapping reads 的处理不同。不要用只保留唯一比对的 BAM 替代 TE 专用分支，也不要把不同 TE 工具的计量层级当作完全等价结果。
