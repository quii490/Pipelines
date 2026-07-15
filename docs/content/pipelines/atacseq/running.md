# 运行 ATAC-seq

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

## 模式

| 模式 | 用途 |
|---|---|
| `auto` | 依次完成可执行的上游和下游步骤 |
| `init` | 发现样本并生成输入文件，不做重计算 |
| `upstream` | FASTQ 到 BAM、track、peak 和 QC |
| `downstream` | 使用已有结果做矩阵、差异和可视化 |

## 推荐顺序

```bash
# 初始化
bash ATAC-seq/run_auto_atacseq.sh --mode init \
  --fastq-dir /path/to/fastq --species hg38 --outdir /path/to/results

# 正式运行
bash ATAC-seq/run_auto_atacseq.sh --mode upstream \
  --fastq-dir /path/to/fastq --species hg38 --outdir /path/to/results \
  --preset standard --resume

# 下游
bash ATAC-seq/run_auto_atacseq.sh --mode downstream \
  --outdir /path/to/results --metadata-csv /path/to/metadata.csv \
  --contrast KO,WT --levels both
```

## 后台与恢复

交互终端可能断开时使用 `--background`，并记录日志位置。失败后先定位第一个失败任务和原因；修复输入或环境后使用 `--resume`。不要删除 work/cache 后再声称是恢复运行。

## 资源控制

`--queue-size`、`--max-cpus`、`--max-memory` 控制全局上限；`--fastp-max-forks` 和 `--fastp-timeout` 用于限制 fastp 并发和时长。资源不足时优先降低并发，不要盲目重复提交。

## 分阶段补跑

已有可靠 BAM 时可只做 downstream。已有上游结果并不代表 metadata 与 contrasts 一定正确；每次补跑仍应保存实际命令、代码 commit、输入文件和日志。
