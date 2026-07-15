# CUT&RUN 参数

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | CUT&RUN maintainers | 2026-07-15 | `main` |

| 参数 | 默认/选择 | 说明 |
|---|---|---|
| `--fastq-dir` | 必需 | FASTQ 根目录 |
| `--species` | `hg38/mm39` | 参考 build |
| `--assay` | `cutrun` | `chipseq/cutrun/cuttag` |
| `--outdir` | 自动推断 | 结果目录 |
| `--layout` | `auto` | PE/SE |
| `--resource-tier` | `standard` | 下游规模 |
| `--te-k` | `25` | TE 分支 multi-mapping 数 |
| `--te-duplicate-policy` | `mark` | TE duplicate 策略 |
| `--run-downstream` | false | 运行附加分析 |

完整参数：

```bash
bash CUTnRUN/pipelines/chipseq_auto_nf/run_auto_chipseq.sh --help
```
