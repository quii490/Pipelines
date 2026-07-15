# 参考基因组与注释

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Reference maintainers | 2026-07-16 | `main` |

FASTA、GTF、aligner index、chrom sizes、blacklist、gene annotation 和 TE annotation 必须属于同一 genome build。

| 资源 | 用途 | 常见不匹配现象 |
|---|---|---|
| FASTA | index、motif、footprinting | contig/序列不一致 |
| GTF/GFF/SAF | gene/TE counting、annotation | assigned rate 低、无注释 |
| STAR/HISAT2/Bowtie2 index | alignment | mapping rate 异常或直接失败 |
| chrom sizes | bigWig、bins、profiles | chromosome not found |
| blacklist BED | 去除高噪声区域 | overlap/过滤结果异常 |
| TSS/gene-body BED | profile 与 QC | 曲线平、坐标风格不匹配 |
| RepeatMasker annotation | TE family/subfamily/locus | TE 名称和数量不同 |

## 染色体命名风格

`chr1` 与 `1`、`chrM` 与 `MT` 是常见差异。不能简单删除 `chr` 就假设资源兼容；先确认 assembly、contig 集和坐标体系。

## hg38、mm10、mm39 与 T2T

不同 build 不是可互换别名。RNA-seq 的 `--human-ref t2t`、CUT&RUN RepEnTools 的 CHM13/T2T 和主流程 hg38 分支回答的参考体系不同，结果不能直接混合。mm39 缺少正确注释时应跳过依赖注释的模块，不能回退到 mm10。

公共仓库只保存 `.example` 配置；大型 index 和内部路径保存在仓库外。每个正式 run 归档 reference 名称、release、校验值或配置快照。
