# ATAC-seq 输入准备

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

推荐 metadata：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

生成的 samplesheet：

```csv
sample,layout,condition,replicate,r1,r2
WT_1,PE,WT,1,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz
```

contrast：

```csv
case,control
KO,WT
```

PE/SE layout 必须正确；混合 layout 会在计数阶段停止并要求拆分运行。
