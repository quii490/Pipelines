# 运行 CUT&RUN

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

推荐顺序：`--init-only` → 编辑 manifest → `--preview` → 正式运行。正式入口默认进行 FASTQ、control、参考资源、磁盘和 gzip preflight。

| 选项 | 用途 |
|---|---|
| `--resource-tier small` | 测试和受限资源 |
| `--resource-tier standard` | 日常推荐 |
| `--resource-tier full` | 全量高分辨率专题分析 |
| `--run-downstream` | 主流程完成后运行下游 |
| `--resume` | 修复后复用缓存 |
| `--dry-run` | 只打印命令 |

不要常规使用 `--skip-preflight`。
