# 独立工具与脚本

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Tool maintainers | 2026-07-15 | `main` |

独立工具用于从已有 BAM、count matrix、bigWig 或 peaks 补跑特定分析。若原始 FASTQ 和完整实验设计仍可用，优先运行主 pipeline；独立工具不应绕过上游 QC 或制造与主流程不一致的结果。

选择：[RNA-seq 工具](rnaseq/index.md)、[ATAC-seq 工具](atacseq/index.md)、[CUT&RUN 工具](cutrun/index.md)。新增稳定工具时使用仓库中的 `docs/content/page-templates/tool.md`；维护模板不发布到网站导航。
