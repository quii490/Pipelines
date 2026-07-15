# 通用故障排查

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Active | Pipeline maintainers | 2026-07-16 |

## 第一步：找到首个失败点

Nextflow 最后一行通常只是汇总。记录完整命令、Git commit、输入表、失败时间，然后查看 automation log、Nextflow log 和第一个失败 task 的：

```text
.command.sh
.command.err
.command.out
.command.log
```

不要一开始就删除 work、修改多个参数或重复提交。

## 快速检查

```bash
df -h /path/to/results
test -r /path/to/input && echo input-readable
test -w /path/to/results && echo results-writable
java -version
nextflow -version
conda info --envs
```

## 现象到原因

| 现象 | 先检查 | 不建议 |
|---|---|---|
| `Permission denied` | 结果父目录 owner/permissions | 用 sudo 跑整个 pipeline |
| `No space left` | results/work/conda cache 和 inode | 随机删除运行中 work |
| task 被 kill / 137 | OOM、scheduler、并发 | 反复 resume 而不调资源 |
| 找不到 command/package | 激活环境、profile、实际 PATH | 在多个环境混装到能运行 |
| chromosome not found | assembly 与 `chr` style | 只删除 chr 前缀后假设兼容 |
| 输出为空但 exit 0 | 输入是否空、状态与必需文件检查 | 当成阴性生物学结果 |
| resume 全部重跑 | work 路径、session、输入/参数变化 | 使用新 work 后期待旧缓存 |
| 图与旧版不同 | normalization、regions、排序、软件版本 | 只比较 PDF 外观 |

## Resume 决策

```text
只修复环境/资源，输入和关键参数没变？ → 原目录 + 原 work + --resume
改 metadata/contrast，仅下游？          → downstream-only + 新输出目录
改 reference/alignment/filter/TE 策略？ → 新结果目录，重新计算受影响阶段
work 已清理？                           → 无法完整 resume，评估从何阶段重建
```

## 报告问题时提供

- 脱敏后的入口命令与 Git commit；
- assay、species/build、PE/SE、strand/control；
- 输入表 header + 1–2 行匿名示例；
- 第一个失败模块和 50–100 行相关日志；
- 可用磁盘/内存与执行环境；
- 你已经尝试的处理，避免重复诊断。

专题问题见 [RNA-seq](../pipelines/rnaseq/troubleshooting.md)、[ATAC-seq](../pipelines/atacseq/troubleshooting.md)和 [CUT&RUN](../pipelines/cutrun/troubleshooting.md)。
