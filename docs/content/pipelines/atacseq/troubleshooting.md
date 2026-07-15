# ATAC-seq 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 | `main` |

排错时保存：完整命令、首个失败任务、日志、软件版本和输入摘要。不要只复制最后一行错误。

## R2 缺失或 PE 配对失败

检查 FASTQ 命名、lane 文件和 `--layout`。不要把缺失 R2 的样本伪装成 SE 后与 PE 样本共同分析。

PE blacklist/mitochondrial filtering 后若出现大量 orphan mates，检查过滤是否按 read name 同步处理两端。PE featureCounts/FRiP 应按 fragment，而不是把两端当两个独立 observations。

## Mapping rate 异常低

确认 `--species`、index 版本、read 长度和污染。随机抽查 FASTQ header 与 FastQC；比较同批样本，而不是立刻放宽 MAPQ。

## mm39 注释相关模块为空

检查 FASTA、chrom names、GTF、blacklist 和 chrom sizes 是否同一 assembly。不得用 mm10 GTF 代替 mm39。

## FRiP 或 peak 数量异常

确认 FRiP 的分母是 reads 还是 fragments、peak 文件是否属于该样本、blacklist/duplicate/filtering 是否一致。结合浏览器轨迹判断是实验背景高还是计算口径错误。

## Downstream matrix 为空

检查 BAM/peak 路径、region coordinates、chromosome naming、metadata sample ID 和文件权限。PE/SE 混合也可能导致单位不一致或流程拒绝。

```bash
head /path/to/metadata.csv
head /path/to/consensus_peaks.bed
samtools idxstats /path/to/sample.clean.bam | head
```

对比 chromosome naming 和实际 sample 列。

## fastp 长时间无进展

查看 fastp process、I/O、输入 gzip 完整性和自动化日志。使用 `--fastp-timeout` 与 `--fastp-max-forks` 控制单任务时长和并发；不要因终端安静就重复提交相同样本。

## bigWig/heatmap 为空

确认 glob 加引号且确实匹配文件；检查 bigWig 与 BED assembly/chromosome style。`skipZeros=true` 可能删除全部区域，`missingDataAsZero` 不能修复坐标不匹配。

## Nucleosome phasing 为 NA

先确认数据为 PE、BAM 保留可靠 fragment length 且深度足够。对 SE 或不适用数据报告 NA 是正确行为。

## `--resume` 没有复用任务

检查是否改变 `--work-dir`、输入路径、参数或清理了缓存。Nextflow resume 依赖原工作目录和任务签名。

## TE track 与常规 track 差异很大

核对 MAPQ、normalization 和 relaxed 策略。这通常是方法差异，不应自动解释为真实 TE 生物学变化。

## TE heatmap 过慢或矩阵很大

先提高 bin size、设置确定性 `--max-regions`，并保持 `--write-values`/大型 TSV 关闭。不要通过 `skipZeros` 仅为了加速，因为它会按结果筛行并改变曲线。

## motif/TOBIAS 失败

确认 genome FASTA、motif MEME、peak BED 和 chromosome naming，检查 HOMER/TOBIAS 是否在当前环境。optional module 失败不应覆盖已完成的 BAM/peaks，但要在报告中注明。
