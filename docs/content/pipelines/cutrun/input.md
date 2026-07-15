# CUT&RUN 输入准备

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

`--init-only` 递归扫描 FASTQ 并创建 manifest。必须人工确认：

- 每个 sample 的 layout 和 R1/R2。
- target/group、biological replicate。
- IgG、Input 或其他 control 是否正确关联。
- assay 为 `cutrun`、`cuttag` 或 `chipseq`。

自动 control 正则只提供初稿。识别不可靠时用 `--control-sample` 或关闭自动填充并手工编辑 manifest。没有合理对照时不得把其他样本冒充 control。
