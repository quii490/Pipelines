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

输入 BAM 建议为 coordinate-sorted `*.clean.bam` + `.bai`。如果 BAM 已经去掉 multi-mappers，`--run-te-relaxed-tracks true` 仍需要额外提供 pre-clean BAM：

```bash
--te-source-bam-glob "/path/to/02_align/*.sorted.bam"
```

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
