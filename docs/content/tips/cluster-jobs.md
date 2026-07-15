# Slurm、后台任务与资源

## 最小作业脚本

```bash
#!/usr/bin/env bash
#SBATCH --job-name=atac_test
#SBATCH --partition=compute
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/%x.%j.out
#SBATCH --error=logs/%x.%j.err

set -euo pipefail
cd /path/to/project
source /path/to/conda.sh
conda activate atac_core

bash /path/to/Pipelines/ATAC-seq/run_auto_atacseq.sh \
  --mode auto \
  --fastq-dir /path/to/fastq \
  --species hg38 \
  --outdir /path/to/results \
  --preset standard --resume
```

提交前先建日志目录：

```bash
mkdir -p logs
sbatch run_atac.sbatch
```

## 查看、取消和诊断

```bash
squeue -u "$USER"
scontrol show job JOB_ID
sacct -j JOB_ID --format=JobID,JobName,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
scancel JOB_ID
```

不要因为暂时没有日志就重复提交同一分析。先检查 job state、queue reason、输出路径和运行进程。

## 怎样申请资源

- CPU：从 pipeline 推荐值开始；线程数增加不一定线性加速，I/O-heavy 工具尤其如此。
- Memory：参考历史 `MaxRSS`，下一次留合理余量。OOM 后只增 CPU 通常无效。
- Time：包括排队后真正运行时间；超时退出可用 `-resume` 的前提是 work/cache 未删除。
- 并发：多个 peak/样本任务会同时占内存和 I/O；总资源不是单任务资源乘以一个很小系数就一定安全。

先用一个样本或 quick preset 做 smoke test，再提交全数据。

## `nohup`、`tmux` 与 scheduler

在普通服务器且管理员允许时：

```bash
nohup bash run.sh > run.log 2>&1 &
echo $! > run.pid
```

或使用 `tmux` 保持交互会话。HPC 的正式长任务仍应使用 scheduler；`nohup` 不能替代资源调度，也不能让 login node 上的大任务变得合规。

## Nextflow 恢复

`-resume` 依赖 `.nextflow/`、work 目录和相同输入/参数 hash。不要为了“清理”提前删除 work。参数改变后，Nextflow 会重跑受影响任务；不能假设所有步骤都复用。

每次运行保存：job ID、commit、完整命令、config、manifest、日志和最终完成状态。
