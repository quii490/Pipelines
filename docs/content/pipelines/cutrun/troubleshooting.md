# CUT&RUN 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

## Preflight 失败

按报告逐项修复 manifest、FASTQ、control、reference、gzip 或磁盘问题。不要先使用 `--skip-preflight` 绕过。

## 下游模块显示 SKIP

检查 `module_status_latest.tsv` 的 reason。可选依赖缺失允许 SKIP；核心 BAM、tracks、peaks 或 counts 失败则不能继续解释结果。

## TE 外部工具没有结果

确认工具安装状态、输入格式与参考 build。流程应记录 SKIP/FAIL，不会用空图替代真实结果。

## 内存或磁盘不足

切换 `--resource-tier small/standard`，限制 locus 数量或分辨率；不要直接启用全量 5 bp 分析。
