# 常用生信文件格式

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Active | Documentation maintainers | 2026-07-16 |

## FASTQ

原始 reads 和 base quality。`.fastq.gz` 是 gzip 压缩，不应解压后提交 Git。PE 的 R1/R2 必须一一配对；lane 文件是否合并取决于 pipeline 的 samplesheet 规则。

## SAM / BAM / CRAM

- SAM：文本 alignment，体积大；
- BAM：SAM 的二进制压缩形式；
- CRAM：依赖参考序列的更高压缩格式。

`sorted BAM` 可能指 coordinate-sorted 或 name-sorted，二者用途不同。`.bai` 是 BAM index，要求 coordinate sort。BAM 本身不告诉你全部处理历史；还要记录 aligner、MAPQ、flags、duplicates、blacklist、multi-mapping 和 reference。

```bash
samtools quickcheck sample.bam
samtools view -H sample.bam | head
samtools flagstat sample.bam
```

## BED / narrowPeak / broadPeak

BED 通常是 0-based、half-open 坐标。至少前三列为 chromosome、start、end；BED6 增加 name、score、strand。narrowPeak/broadPeak 是扩展 BED，额外列含 peak score/significance/summit 等；不要按普通 CSV 随意重排列。

BED 与 GTF 坐标体系不同，转换时要使用可靠工具。染色体命名和 genome build 必须匹配。

## GTF / GFF3 / SAF

- GTF/GFF：分层注释，包含 feature、strand 和 attributes；
- SAF：featureCounts 常用的简化区域表，通常含 `GeneID,Chr,Start,End,Strand`。

GTF 的 exon/gene_id、TE attributes 和版本会改变 counts。SAF 常为 1-based inclusive；不要与 BED 坐标直接互换。

## bedGraph / bigWig

bedGraph 是文本区间信号；bigWig 是带 index 的二进制轨迹，适合 IGV/UCSC 和 deepTools。bigWig 不是 counts matrix，数值含义由 CPM/RPGC/RPKM 等 normalization 和 bin size 决定。

## Count matrix

第一列是 feature ID，后面是样本 raw counts。它用于 DESeq2/edgeR 类 count model。不要把 TPM/FPKM/CPM/VST 当作 raw counts。

## TSV / CSV

TSV 用 tab，CSV 用 comma。用电子表格打开时可能自动把 gene 名转换成日期、丢失前导零或改变科学计数法。正式分析优先用脚本读取，并保留原文件校验值。

## PDF / PNG / SVG / HTML

- PDF/SVG：矢量图，适合文章和编辑；
- PNG：位图，适合预览，注意分辨率；
- HTML：交互/综合报告，可能引用外部资源或包含敏感路径。

交付图的同时保留生成它的 TSV/matrix、脚本和参数。
