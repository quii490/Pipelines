# RNA-seq 输入准备

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

入口接受本地 `--fastq-dir` 或 `--manifest`。也可提供：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

contrast 文件：

```csv
case,control
KO,WT
```

`KO_vs_WT` 表示 KO 相对 WT。样本名必须唯一；PE 样本必须同时存在 R1/R2；开始分析前确认 strandedness，错误的 strand 会直接影响 gene 和 TE counts。
