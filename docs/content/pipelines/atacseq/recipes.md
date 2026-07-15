# ATAC-seq：按现有文件选择入口

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## 我有 FASTQ

使用 `run_auto_atacseq.sh --mode init` 生成输入，再以 `--mode auto` 或 `upstream` 正式运行。metadata 和 contrast 应显式核对。

## 我有完整主流程 results

如果 `07_counts/consensus_peak_counts.txt`、sample metadata 和 consensus peaks 完整，可直接：

```bash
bash ATAC-seq/scripts/run_atac_downstream_only.sh \
  --result-dir /path/to/results \
  --species hg38 \
  --contrast-file /path/to/contrasts.csv \
  --outdir /path/to/results/08_downstream_v2
```

适合修改 contrast、TE filter 或重新绘图，不会重做 alignment/peaks/counts。

## 我有 clean BAM

使用 `run_atac_from_bam.sh` 可重建 bigWig、peaks、consensus、counts、fixed-bin 和 downstream：

```bash
bash ATAC-seq/scripts/run_atac_from_bam.sh \
  --bam-dir /path/to/clean_bam \
  --sample-meta /path/to/metadata.csv \
  --contrast-file /path/to/contrasts.csv \
  --species hg38 \
  --outdir /path/to/atac_from_bam \
  --cores 16
```

目录需含 `*.clean.bam`，建议同时有 `.bai`。先确认 BAM 的 MAPQ、duplicate、mitochondrial、blacklist、proper-pair 和 Tn5 shift 处理。

!!! warning "TE relaxed track 不能从 strict clean BAM 恢复"

    要保留更多重复区域 reads，必须提供 pre-clean `*.sorted.bam` 给 `--te-source-bam-glob`，或直接使用 FASTQ→TE track 工具。已经在 strict clean 阶段删除的信息无法补回来。

## 我有 BAM，只想做一个模块

- 重新 call peak：`run_callpeak_from_bam.sh`；
- fixed-bin：`run_fixedbin_from_bam.sh`；
- TE relaxed track：`run_atac_te_tracks_from_bam.sh`；
- nucleosome phasing：`run_nuc_phasing.sh`（PE）；
- BAM/track 组合见[ATAC 工具手册](../../tools/atacseq/index.md)。

## 我有 bigWig + BED/peaks

可以运行 TSS、gene-body、region/differential peak/TE heatmap 和 bigWig correlation。必须先确认 track normalization、assembly、bin size 与 BED chromosome naming 一致。

## 我有 differential region 表

可以补 cross-contrast logFC heatmap、region annotation 或 motif。不要用已经按显著性筛选的数据重新估计显著性；这些工具主要用于展示和注释。
