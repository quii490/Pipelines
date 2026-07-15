# ATAC-seq 输入准备

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | ATAC-seq maintainers | 2026-07-16 | `main` |

## FASTQ

推荐命名：`<sample>_R1.fastq.gz` 与 `<sample>_R2.fastq.gz`。样本名只使用字母、数字、下划线和连字符；同一 PE 样本必须同时有 R1/R2。

```text
WT_1_R1.fastq.gz  WT_1_R2.fastq.gz
WT_2_R1.fastq.gz  WT_2_R2.fastq.gz
KO_1_R1.fastq.gz  KO_1_R2.fastq.gz
```

流程可递归扫描 `--fastq-dir`。正式运行前要排除 lane 拆分、重复命名或把其他 assay FASTQ 混入目录的情况。

## Metadata

建议显式维护：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

`sample` 必须与流程识别的样本 ID 完全一致。`condition` 决定统计设计，`replicate` 是生物学重复编号，不要把 lane 当成生物学重复。

## Contrast

命令行格式为 `CASE,CONTROL`：

```bash
--contrast KO,WT
```

多个比较可重复传入，或使用 `--contrast-file`。运行前把方向写入分析计划；正的 log2 fold change 表示 CASE 相对 CONTROL 增加。

## PE、SE 与参考基因组

- `--layout auto` 会尝试判断布局，但关键项目应人工复核；
- PE 按 fragment 计数，SE 按 read 计数，不能混合构建同一矩阵；
- `mm39` 必须使用匹配版本的 FASTA、index、chrom sizes、blacklist 和 GTF；缺少正确注释时宁可跳过相关模块，也不能回退到 mm10 注释。
