# 命令行最小知识

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Active | Documentation maintainers | 2026-07-16 |

不需要先学完整 Linux；先掌握下面这些操作即可安全使用流程。

## 路径与文件检查

```bash
pwd                         # 当前目录
ls -lh /path/to/project     # 查看文件和大小
find /path/to/fastq -type f -name '*.fastq.gz'
test -r /path/to/file && echo readable
test -w /path/to/results && echo writable
df -h /path/to/results      # 可用磁盘
```

有空格的路径必须整体加引号。正式分析更建议避免空格和中文标点。

## 变量让命令更清楚

```bash
PROJECT=/path/to/project
FASTQ="$PROJECT/fastq"
RESULTS="$PROJECT/results"
```

先用 `echo "$RESULTS"` 检查变量，再传入流程。不要复制含真实服务器路径的命令到公开文档。

## 查看日志而不停止任务

```bash
tail -n 100 /path/to/run.log
tail -f /path/to/run.log
```

`tail -f` 中按 `Ctrl-C` 只停止查看日志，不会停止后台任务。

## glob 为什么要加引号

```bash
--bw-glob "/path/to/bw/*.bw"
```

引号让脚本自己处理 glob。若不加引号，shell 可能提前展开成多个参数，导致“unknown option”或只读取第一个文件。

## 成功不只看 exit code

命令返回 0 后还要检查：必需文件非空、样本数正确、状态表没有核心失败、QC 合理、参数和版本已记录。`SKIP`、`FAIL`、空表与“没有生物学信号”含义不同。

## 不确定时的安全顺序

```text
--help → 检查输入 → dry-run/preview → 小规模测试 → 正式运行 → QC → 解释
```
