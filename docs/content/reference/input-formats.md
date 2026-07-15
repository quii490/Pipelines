# Samplesheet、Manifest 与 Contrast

## 通用规则

- UTF-8、comma-separated、首行为 header；
- sample ID 唯一，只用字母、数字、下划线和连字符；
- 路径在实际执行节点可读；
- condition 与 contrast 大小写完全一致；
- replicate 是 biological replicate，不是 lane；
- 编辑后用命令行或文本编辑器检查，防止电子表格改变 ID/日期。

## RNA-seq 自动生成 samplesheet

```csv
sample,layout,condition,replicate,r1,r2
WT_1,PE,WT,1,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz
KO_1,PE,KO,1,/path/to/KO_1_R1.fastq.gz,/path/to/KO_1_R2.fastq.gz
```

RNA 下载 manifest 可接受本地 FASTQ、URL、SRR/ERR/DRR 和部分 GSM，字段别名较多；解析后的 samplesheet/metadata 才是运行记录。正式项目还应记录 strandedness 和 batch。

## ATAC-seq samplesheet

```csv
sample,layout,condition,replicate,r1,r2
WT_1,PE,WT,1,/path/to/WT_1_R1.fastq.gz,/path/to/WT_1_R2.fastq.gz
KO_1,PE,KO,1,/path/to/KO_1_R1.fastq.gz,/path/to/KO_1_R2.fastq.gz
```

同一 sample 可以在 `r1/r2` 中记录 lane 文件列表，但必须保证配对顺序和数量一致。PE/SE 不进入同一 count matrix。

外部 metadata 简表：

```csv
sample,condition,replicate
WT_1,WT,1
KO_1,KO,1
```

## CUT&RUN manifest（Draft pipeline）

```csv
sample,species,assay,group,replicate,igg,is_igg,layout,fastq_1,fastq_2
H3K27ac_WT_1,hg38,cutrun,H3K27ac_WT,1,IgG_WT_1,false,PE,/path/to/R1.fastq.gz,/path/to/R2.fastq.gz
IgG_WT_1,hg38,cutrun,IgG_WT,1,,true,PE,/path/to/IgG_R1.fastq.gz,/path/to/IgG_R2.fastq.gz
```

每个 target 的 `igg` 必须指向表中存在且 `is_igg=true` 的 control。不要修改运行后生成的 `resolved_manifest.csv`。

## Contrast

```csv
case,control
KO,WT
```

统一解释为 CASE vs CONTROL；正 log2FC 表示 CASE 高于 CONTROL。RNA 命令行通常使用 `KO:WT`，ATAC 使用 `KO,WT`，以各入口 `--help` 为准。

## 提交前验证

```bash
head -n 5 /path/to/samplesheet.csv
awk -F, 'NR==1{print "columns=" NF} NR>1{print $1}' /path/to/samplesheet.csv
```

更复杂的 CSV（引号内逗号等）应使用 Python/R CSV parser，不依赖简单 `awk`。
