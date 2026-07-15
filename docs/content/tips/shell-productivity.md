# Shell、日志与文件传输技巧

## 先确认自己在哪里

```bash
hostname
pwd
whoami
date
df -h /path/to/project
```

很多“文件不存在”其实是当前主机、目录或挂载点不对。

## 快速寻找文件和文本

```bash
# 仓库内找文件
rg --files | rg 'samplesheet|manifest|QC_REPORT'

# 找参数或错误文字
rg -n -- '--te-k|ERROR|Traceback' /path/to/project

# 查看最近日志
tail -n 100 /path/to/run.log
tail -f /path/to/run.log
```

搜索以 `-` 开头的字符串时用 `--` 分隔选项，例如 `rg -n -- '--resume'`。

## 检查文件而不打开全部内容

```bash
ls -lh /path/to/file
wc -l /path/to/table.tsv
head -n 5 /path/to/table.tsv
samtools quickcheck -v /path/to/sample.bam
samtools flagstat /path/to/sample.bam
```

对 gzip FASTQ 可用 `gzip -t file.fastq.gz` 检查完整性，但大文件会完整读取，先确认 I/O 和作业节点政策。

## 安全传输

```bash
# 单个文件
scp local.txt lab-cluster:/path/to/project/

# 可恢复目录同步；先 dry-run
rsync -avhn --partial /local/data/ lab-cluster:/path/to/data/
rsync -avh --partial /local/data/ lab-cluster:/path/to/data/
```

注意源路径末尾 `/` 会影响 rsync 目录层级。不要在未检查方向时使用 `--delete`。

传输后核对 checksum：

```bash
shasum -a 256 local.file
sha256sum remote.file
```

## 日志过滤

```bash
rg -n -i 'error|fatal|killed|oom|no space|segmentation' logs/
rg -n 'Completed at|Execution complete|exit status' .nextflow.log
```

只截取 error 行可能丢失上下文；向同事或 AI 求助时提供报错前后约 30–50 行、运行命令、工具版本和输入类型，并先脱敏。

## 不要把临时命令变成不可复现历史

探索命令确认有效后，写入 `.sh` 文件：

```bash
#!/usr/bin/env bash
set -euo pipefail

readonly INPUT="${1:?usage: script.sh INPUT OUTDIR}"
readonly OUTDIR="${2:?usage: script.sh INPUT OUTDIR}"
mkdir -p "$OUTDIR"
```

对包含空格或 glob 的路径正确加引号；不要把密码、token 或真实内部地址写进脚本。
