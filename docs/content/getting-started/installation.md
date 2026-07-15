# 系统与安装

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-15 | `main` |

## 支持环境

- Linux 或 HPC/容器环境；macOS 主要用于维护文档和代码。
- Bash、Python 3、Conda/Mamba 与 Nextflow。
- 建议单项目至少 16 CPU、64 GB memory；实际需求随测序深度和并发数变化。

## 获取代码与检查入口

```bash
git clone https://github.com/quii490/lab-pipelines.git
cd lab-pipelines

bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
bash ATAC-seq/run_auto_atacseq.sh --help
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

分析环境按各 pipeline 的 `envs/`、`nextflow.config` 和安装说明建立。不要把个人 Conda 绝对路径写入公共配置；使用环境变量或本地 `*.local.config`。
