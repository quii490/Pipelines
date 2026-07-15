# ATAC-seq 参数

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

| 参数 | 默认/选择 | 说明 |
|---|---|---|
| `--fastq-dir` | 必需（上游） | FASTQ 根目录，递归扫描 |
| `--species` | `hg38/mm10/mm39` | 参考基因组 |
| `--outdir` | 自动推断 | 结果目录 |
| `--layout` | `auto` | PE/SE |
| `--preset` | `standard` | `quick/standard/full` |
| `--metadata-csv` | 可选 | condition 与 replicate |
| `--contrast` | 可重复 | `CASE,CONTROL` |
| `--levels` | `both` | downstream 的 peak/bin 层级 |
| `--resume` | 推荐 | 复用 Nextflow 缓存 |

完整参数：

```bash
bash ATAC-seq/run_auto_atacseq.sh --help
```
