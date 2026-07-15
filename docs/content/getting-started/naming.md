# 样本命名与目录规范

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-15 | `main` |

样本名使用字母、数字、下划线和短横线，不使用空格、斜杠或中文标点。推荐：

```text
project/
├── fastq/
│   ├── WT_1_R1.fastq.gz
│   ├── WT_1_R2.fastq.gz
│   ├── KO_1_R1.fastq.gz
│   └── KO_1_R2.fastq.gz
├── metadata/
└── results/
```

同一分析中的 `sample` 必须唯一。PE 数据的 R1/R2 basename 必须一致；PE 和 SE 不在同一个 count matrix 中比较。结果目录不要放进 Git。
