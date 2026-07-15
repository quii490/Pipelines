# ATAC-seq 重建与补跑

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## 已有主流程 results：只跑 downstream

需要：

```text
07_counts/consensus_peak_counts.txt
07_counts/sample_metadata.csv
06_consensus_peaks/consensus_peaks.bed
```

```bash
bash ATAC-seq/scripts/run_atac_downstream_only.sh \
  --result-dir /path/to/results \
  --species hg38 \
  --contrast-file /path/to/contrasts.csv \
  --outdir /path/to/results/08_downstream_manual
```

适合改 contrast、TE family/name filter 或重画图。若 consensus peaks/counts 本身需要改变，应使用 BAM 重建或主流程，不要只跑 downstream。

## 已有 clean BAM：重建完整 downstream

```bash
bash ATAC-seq/scripts/run_atac_from_bam.sh \
  --bam-dir /path/to/02_align \
  --sample-meta /path/to/metadata.csv \
  --contrast-file /path/to/contrasts.csv \
  --species hg38 \
  --outdir /path/to/atac_from_bam \
  --cores 16
```

默认可生成普通 bigWig、peaks、consensus、counts、downstream、fixed-bin 与 bigWig correlation。按需开启：

```text
--run-tss true
--run-gene-body true
--run-motif true
--run-tobias true
--run-nuc-phasing true
--run-diff-peak-heatmap true
```

输入 BAM 建议为 coordinate-sorted `*.clean.bam` + `.bai`。当前 `run_atac_from_bam.sh` **没有** `--run-te-relaxed-tracks` 或 `--te-source-bam-glob` 参数。clean BAM 已去掉的 multi-mappers 无法恢复；需要 relaxed track 时另行运行：

```bash
bash ATAC-seq/scripts/run_atac_te_tracks_from_bam.sh \
  --bam-glob "/path/to/02_align/*.sorted.bam" \
  --outdir /path/to/te_tracks \
  --species hg38 --mapq 0 --cores 16
```

### 分析开关

这些开关均接受 `true`/`false`：

| 参数 | 模块 | 前置条件 |
|---|---|---|
| `--run-bamcoverage` | clean BAM → standard bigWig | `bamCoverage` |
| `--run-tss` | TSS aggregate profile | TSS BED、bigWig |
| `--run-gene-body` | gene-body profile | gene-body BED、bigWig |
| `--run-te-heatmap` | standard track 上的 TE aggregate | TE BED、bigWig；不是 relaxed branch |
| `--run-motif` | HOMER | motif genome 与 contrast BED |
| `--run-tobias` | footprinting | FASTA、MEME motif、BAM、peaks |
| `--run-nuc-phasing` | fragment/NRL | PE BAM |
| `--run-fixedbin` | fixed-bin counts/DE | chrom sizes、metadata、contrasts |
| `--run-bw-correlation` | bigWig correlation/PCA | 至少两个同策略 tracks |
| `--run-diff-peak-heatmap` | differential peak signal | contrast BED、bigWig |
| `--run-consensus-annotation` | consensus peak gene/TE annotation | GTF、可选 TE BED |

一次只增加少量可选模块。若全部打开，先确认 HOMER、TOBIAS、deepTools、R/Bioconductor 和参考文件都齐全，否则很难区分是核心重建失败还是可选工具缺失。

## 只重新 call peaks

```bash
bash ATAC-seq/scripts/run_callpeak_from_bam.sh \
  --bam-glob "/path/to/clean_bam/*.clean.bam" \
  --outdir /path/to/peaks \
  --species hg38 \
  --format auto \
  --qvalue 0.05 \
  --cores 8
```

PE 数据通常让 `--format auto` 检测并使用 `BAMPE`；只有确认 BAM layout 时才强制格式。SE 的 `--shift -100 --extsize 200` 是常见 ATAC 模型，不能用于 PE fragment calling。`--pvalue` 一旦设置会覆盖 `--qvalue`。`--broad` 适合特殊探索，不是普通 ATAC 的默认选择。

调用前确认 BAM 是否已经 Tn5 shift、duplicate/blacklist/MAPQ 如何处理。不同 peak 参数生成的新 peak set 不应覆盖旧目录。peak 数量增加不自动表示质量提高。

## Fixed-bin 分析

```bash
bash ATAC-seq/scripts/run_fixedbin_from_bam.sh \
  --bam-glob "/path/to/clean_bam/*.clean.bam" \
  --sample-meta /path/to/metadata.csv \
  --contrast-file /path/to/contrasts.csv \
  --species hg38 \
  --outdir /path/to/fixedbin \
  --bin-size 100000 \
  --cores 8
```

Fixed bins 不依赖 peak calling，适合发现分散变化；代价是 feature 数量大、multiple testing 更重。bin size 改变分辨率、稀疏度和计算量，比较项目时必须一致。

常用 `--bin-size` 试验可从 50–100 kb 开始；不要在同一结论中挑选多个 bin size 里最显著的一套而不报告多重尝试。
