# 内容状态与完善顺序

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Documentation maintainers | 2026-07-16 | `main` |

| 模块 | 当前阶段 | 内容来源 | 下一次审查重点 |
|---|---|---|---|
| 维护指南 | Active | MkDocs 配置和实际发布流程 | 新成员能否独立新增页面 |
| RNA-seq | Draft review | 当前代码、完整指南、参数手册 | 命令、默认值、输出、QC |
| ATAC-seq | Draft review | 当前代码、pipeline/script 手册 | PE/SE、QC、专题工具 |
| CUT&RUN | Draft / under development | 当前代码、新手 SOP、输出和故障手册 | pipeline 未稳定；接口、默认值和结果继续验证 |
| 工具与脚本 | Draft review | 各工具 README、源码和 `--help` | 用真实小数据验证每个示例 |
| 跨流程专题 | Draft review | 三套流程的共同定义 | 实验室推荐值和解释底线 |
| 教程 | Draft review | 场景化通用示例 | 增加公开 fixtures 和稳定预期结果 |

当前完善优先级为 RNA-seq → ATAC-seq → 入门/工具/专题/教程；CUT&RUN 等代码稳定后再完成最终科学验收。页面只有经过技术检查、真实小数据运行和实验室科学审核后才从 Draft 转为 Active。
