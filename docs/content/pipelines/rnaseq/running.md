# 运行 RNA-seq

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

- 正式运行：在 Quick Start 命令中移除 `--dry-run`。
- 后台运行：加 `--background`，日志写入结果目录的自动化日志目录。
- 恢复：默认启用 `--resume`；需要全新执行时才使用 `--no-resume`。
- 只跑上游：`--upstream-only`。
- 只跑下游：`--downstream`。
- 仅运行部分下游模块：`--only-tools`；跳过模块用 `--skip-tools`。

不要混用 `--results-dir` 与低层 `run_pipeline.sh` 的 `--outdir`。日常使用统一入口。
