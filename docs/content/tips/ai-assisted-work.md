# 使用 AI 辅助生信工作

AI 适合加速理解、生成初稿和审阅，但不能替代数据授权、软件文档、代码测试、统计设计和科学判断。

## 适合交给 AI 的任务

- 解释一段脱敏报错并列出诊断顺序；
- 根据真实 `--help` 生成命令模板或参数表；
- 把一次成功的交互命令整理成带 `set -euo pipefail` 的脚本；
- 审阅 samplesheet 的列名、重复 sample ID 和明显路径问题；
- 比较两版脚本 diff，指出可能改变结果的参数；
- 编写 smoke test、检查输出存在性或汇总日志；
- 将代码和已有结果说明整理成文档初稿。

不适合直接相信 AI 的事项：凭记忆给出软件最新参数、替你决定 contrast/协变量、虚构 QC 阈值、根据一张图排除样本、解释未检查的真实临床数据。

## 一个有效的提问模板

```text
目标：从已有 clean BAM 重新 call ATAC peaks，不重跑 FASTQ。
环境：Linux + Slurm；仓库 commit <COMMIT>；MACS3 <VERSION>。
输入：PE、hg38、coordinate-sorted BAM，已去 duplicates/blacklist。
事实来源：下面是脚本 --help（已脱敏）。
约束：不要删除或覆盖原结果；先给 dry-run/检查命令。
希望输出：1) 推荐命令；2) 每个非默认参数解释；3) 验证清单。
```

“帮我跑一下 ATAC”信息太少；输入阶段、assembly、layout、已有过滤、目标输出和不能做的操作越清楚，答案越可靠。

## 安全边界

发送给外部 AI 前删除或替换：患者/受试者信息、未发表 sample 名称、真实服务器地址、用户名、内部路径、accession 受限数据、token、SSH key、云密钥和完整 `.env`。

```text
REAL_SERVER_PATH  -> /path/to/project
REAL_SAMPLE_ID    -> SAMPLE_A
REAL_HOSTNAME     -> login.example.org
```

不要让 AI 输出或保存 secret。若 secret 曾进入对话、日志或 Git history，按泄露处理并立即轮换；仅删除当前文件不够。

## 让 AI 先读事实来源

优先提供：

```bash
tool --version
tool --help
git rev-parse --short HEAD
head -n 5 samplesheet.csv
tail -n 80 failing.log
```

对于本仓库，可要求 AI 先检查入口脚本、实际 process/output 和已有 README，再改文档。不要让它仅凭工具名猜参数。

## 验证 AI 生成的命令

按顺序检查：

1. 命令中的工具和参数在当前 `--help` 中存在；
2. 输入路径、assembly、PE/SE、strand、sample IDs 正确；
3. 输出写到新目录，没有 `rm`、覆盖、上传或公开操作；
4. 用一个样本、`--dry-run`、preview 或小 fixture 测试；
5. 查看 exit code、日志、输出非空和关键统计；
6. 用 `git diff` 审查代码改动，再提交。

```bash
bash -n proposed_script.sh
git diff --check
git diff -- proposed_script.sh
```

R/Python/Bash 语法通过不代表科学逻辑正确。还需核对 feature universe、normalization、contrast 方向、重复设计和结果口径。

## 推荐的协作方式

让 AI 每次做一个可验证的小任务：

```text
先只审阅，不修改：列出脚本真实输入、输出和风险参数。
然后修改文档，不能改分析默认值。
最后运行 --help、smoke test、链接检查和严格构建，给出 diff 摘要。
```

要求它区分：已从代码确认的事实、推断、需要实验室维护者确认的科学决策。保留 commit 和测试记录，使 AI 帮助也成为可审计工作流的一部分。
