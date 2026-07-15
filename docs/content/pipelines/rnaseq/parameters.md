# RNA-seq 参数

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

| 参数 | 必需 | 默认/选择 | 说明 |
|---|---|---|---|
| `--fastq-dir` / `--manifest` | 上游必需其一 | — | 输入来源 |
| `--results-dir` | 推荐 | 自动推断 | 主结果目录 |
| `--species` | 是 | `hg38`、`mm10`、`mm39` | 参考物种 |
| `--strand` | 是 | `unstranded/forward/reverse` | gene 与 TE strandedness |
| `--layout` | 否 | `auto` | PE/SE 检测 |
| `--aligner` | 否 | `star` | gene aligner |
| `--max-cpus` | 否 | 环境决定 | 最大 CPU |
| `--max-memory` | 否 | 环境决定 | 例如 `64 GB` |
| `--dry-run` | 推荐首次使用 | false | 只预检和打印计划 |

完整事实来源：

```bash
bash RNA-seq/rnaseq/run_auto_rnaseq.sh --help
```

TE、多重比对和无重复差异参数属于高级选项，修改前应记录原因并在测试数据比较结果。
