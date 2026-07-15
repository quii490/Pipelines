# 参考基因组

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | Lab pipeline maintainers | 2026-07-15 | `main` |

FASTA、GTF、aligner index、chromosome sizes、blacklist、gene annotation 和 TE annotation 必须来自同一 genome build。支持范围以各 pipeline 参数页为准。

| 资源 | 常见用途 |
|---|---|
| FASTA | reference 检查、motif/footprinting |
| GTF/SAF | gene 与 TE counting、peak annotation |
| STAR/Bowtie2 index | RNA 或表观组比对 |
| chromosome sizes | bigWig、window 和 profile |
| blacklist | 去除高噪声区域 |
| RepeatMasker annotation | TE family/subfamily/locus 分析 |

公共仓库只保存配置模板，不保存大型 index。真实路径通过本地配置或环境变量提供。
