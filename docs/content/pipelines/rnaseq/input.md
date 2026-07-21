# RNA-seq 输入准备
## 输入模式
入口至少需要以下一种：
1. `--fastq-dir`：递归扫描本地 FASTQ。
2. `--manifest`：读取本地文件、FASTQ URL、SRR/ERR/DRR 或部分 GSM accession。
3. `--downstream`：不跑 FASTQ 上游，只读取已有结果和设计文件。
## FASTQ 命名
默认使用 `_R1` 和 `_R2` 判断 mates，可用 `--r1-pattern`、`--r2-pattern` 覆盖。
```text
WT_1_R1.fastq.gz
WT_1_R2.fastq.gz
WT_2_R1.fastq.gz
WT_2_R2.fastq.gz
```
常见错误：重复 sample ID、R2 缺失、R1/R2 basename 不一致、一个样本出现多套无法区分的 FASTQ、路径指向结果目录而非原始数据目录。
## Metadata
可以通过 `--sample-metadata` 提供：
```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```
也可以重复使用：
```bash
--sample-meta WT_1:WT:1 \
--sample-meta WT_2:WT:2
```
或按组映射：
```bash
--condition-map WT:WT_1,WT_2 \
--condition-map KO:KO_1,KO_2
```
## Contrast
推荐 CSV：
```csv
case,control
KO,WT
```
也可以重复使用 `--contrast KO:WT`。`case` 是分子，正 log2 fold change 表示 case 相对 control 增高。
## Strandedness
`--strand` 同时控制 gene 和 TE counting：

| 值 | 含义 |
|---|---|
| `unstranded` | 不区分链 |
| `forward` | 正向链特异 |
| `reverse` | 反向链特异；当前默认 |

strandedness 错误常表现为 featureCounts assignment 显著偏低。不要根据结果高低随意选择，应依据建库方案或用小规模验证工具确认。
## Reference 与物种
- `--species hg38 --human-ref hg38`：常规人类 hg38。
- `--species hg38 --human-ref t2t`：T2T/CHM13 分支；需相应资源。
- `--species mm10`、`--species mm39`：小鼠对应 build。
FASTA、GTF、STAR/HISAT2 index、TE annotation 和 chromosome naming 必须同 build。禁止只替换其中一个文件。
