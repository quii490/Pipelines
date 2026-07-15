# 可复现性与结果归档

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Active | Pipeline maintainers | 2026-07-16 |

## 每个正式 run 至少记录

- Git commit 或 release；
- 完整命令和低层 `--extra`；
- samplesheet/manifest、metadata、contrasts；
- environment/profile 与软件版本；
- reference assembly、annotation release 和配置快照；
- pipeline/module status、日志和 QC；
- 关键 raw counts/regions、完整差异表和绘图脚本。

## Immutable 与 latest

`latest` 文件方便查看，但可能被后续 run 覆盖。正式结论应指向带 run ID/时间戳的不可变产物。修改关键方法时创建新结果目录，不覆盖旧 run。

## Resume 与清理

resume 依赖 work/cache 和任务 fingerprint。清理前确认：任务已结束、结果已验证、以后不需要 resume、关键 logs/version 已归档。不要把 work 当作正式结果长期交付，也不要运行中删除它。

## 图的可复现性

每张最终图至少关联：输入 table/matrix、region selection、sample order、normalization、filter/cutoff、脚本 commit 和输出文件。仅保留 PDF 不能重建分析。

## 公开与内部归档分开

内部归档可保存受控的真实 sample mapping 和服务器配置；公开仓库只保存脱敏示例。`git history` 仍会保留曾提交的秘密，因此发现密钥要立即撤销并按安全流程清除历史。
