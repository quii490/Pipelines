# 运行 ATAC-seq

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

| 模式 | 用途 |
|---|---|
| `--init-only` | 只生成可编辑输入表 |
| `--preset quick` | 小数据连通性测试 |
| `--preset standard` | 日常推荐主分析 |
| `--preset full` | 按需启用耗时专题模块 |
| `--mode downstream` | 从已有上游结果重跑下游 |

使用 `--background` 后从结果目录日志查看状态。失败修复后在原结果目录加 `--resume`。Motif、footprinting 和高分辨率 profile 属于按需分析，不必每次全开。
