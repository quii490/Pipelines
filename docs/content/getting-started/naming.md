# 样本命名与目录规范

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-16 | `main` |

## 推荐目录

```text
project_name/
├── fastq/                 # 原始数据，只读
├── metadata/
│   ├── samples.csv
│   └── contrasts.csv
├── results/               # 当前正式 run
├── results_test/          # dry-run/小测试
├── deliverables/          # 脱敏后的交付文件
└── notes/                 # 实验设计、版本和决策记录
```

原始 FASTQ 不要在多个项目目录中复制修改。结果目录不进入 Git。

## 样本命名

推荐 `<condition>_<replicate>`，例如 `WT_1`、`KO_2`。只使用 ASCII 字母、数字、下划线和连字符；不用空格、斜杠、括号或中文标点。

```text
WT_1_R1.fastq.gz
WT_1_R2.fastq.gz
KO_1_R1.fastq.gz
KO_1_R2.fastq.gz
```

- 同一分析中 `sample` 唯一；
- R1/R2 basename 必须匹配；
- lane 不是 biological replicate；
- 不要仅从文件名推断 condition/control；以 metadata/manifest 为准；
- PE 和 SE 的计数单位不同，不放入同一 count matrix。

## 版本化结果目录

关键方法改变时使用新目录：

```text
results_v1_hg38_default/
results_v2_hg38_te_k50/
results_v3_mm39/
```

不要在同一目录混合不同 genome build、strandedness、TE multi-mapping 或 peak 策略。
