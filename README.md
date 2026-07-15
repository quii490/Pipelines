# Lab Bioinformatics Pipelines

可复现的 RNA-seq、ATAC-seq 与 CUT&RUN/CUT&Tag/ChIP-seq 分析流程，以及配套的中文使用文档。

## 项目入口

| 模块 | 推荐入口 | 说明 |
|---|---|---|
| RNA-seq | `RNA-seq/rnaseq/run_auto_rnaseq.sh` | FASTQ/manifest 到 gene、TE 与下游分析 |
| ATAC-seq | `ATAC-seq/run_auto_atacseq.sh` | FASTQ 到 peaks、tracks、QC 与差异分析 |
| CUT&RUN | `CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh` | CUT&RUN、CUT&Tag、ChIP-seq 标准与 TE 分支 |

完整文档：<https://quii490.github.io/Pipelines/>

本地预览：

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r docs/requirements.txt
mkdocs serve -f docs/mkdocs.yml
```

## 公开仓库约束

本仓库不得包含真实测序数据、患者信息、未发表样本标识、密钥、内网地址或实验室私有资源路径。提交前请阅读 [CONTRIBUTING.md](CONTRIBUTING.md) 并运行：

```bash
bash scripts/validate_public_repo.sh
```

## 许可证

代码与文档按 [MIT License](LICENSE) 发布。第三方软件、数据库和参考资源仍受其各自许可证约束。
