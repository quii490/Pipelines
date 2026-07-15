# 系统与安装

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-16 | `main` |

## 支持环境

- Linux、HPC 或持久化容器用于分析；macOS 主要用于维护代码和文档；
- Bash、Python 3、Conda/Mamba、Java 与 Nextflow；
- 建议起步资源 16 CPU、64 GB memory；实际需求取决于深度、样本数、TE 分支和并发；
- 结果与 Nextflow work 需要持久化存储。临时 Pod 被删除时，非持久化 work 无法 resume。

## 获取仓库

```bash
git clone https://github.com/quii490/Pipelines.git
cd Pipelines
```

## 最小入口检查

```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
bash ATAC-seq/run_auto_atacseq.sh --help
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```

这只证明入口可读取，不证明所有依赖已安装。正式运行前还应检查：

```bash
java -version
nextflow -version
conda --version
samtools --version
```

## 环境原则

三个 pipeline 的依赖不同，不建议把所有软件强塞进一个环境。优先使用各目录的 `envs/`、`nextflow.config` 和 profile。独立工具若通过 wrapper 调用另一个环境，应记录 wrapper 和实际 R/Python 环境。

## 本地预览文档

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
mkdocs serve -f docs/mkdocs.yml
```

浏览器打开终端显示的本地地址。退出预览按 `Ctrl-C`。
