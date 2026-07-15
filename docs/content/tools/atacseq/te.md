# ATAC-seq TE / L1 工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## 从 pre-clean BAM 生成 relaxed tracks

```bash
bash ATAC-seq/scripts/run_atac_te_tracks_from_bam.sh \
  --bam-glob "/path/to/02_align/*.sorted.bam" \
  --outdir /path/to/te_tracks \
  --species hg38 \
  --mapq 0 \
  --exclude-flags 780 \
  --run-markdup false \
  --remove-mito true \
  --remove-blacklist true \
  --normalization CPM \
  --binsize 10 \
  --cores 16
```

输出 `02_align_te/*.te.bam`、flagstat、逐步 clean counts 和 `04_bw_te/*.te.bw`。先看 `*.te.clean_counts.tsv`，确认 input、MAPQ/flag、mitochondrial、blacklist 到 final 每步保留量。

`exclude-flags 780` 不等于移除所有 duplicate-marked reads；`--run-markdup`、proper-pair、MAPQ 的选择会显著改变信号。

## 从 FASTQ 直接生成 TE tracks

当没有 pre-clean BAM 时使用 `run_atac_te_tracks_from_fastq.sh`。它会生成 sorted BAM 和 relaxed BAM/bigWig，适合独立 TE 诊断；不是完整 ATAC peak/downstream 的替代品。先运行：

```bash
bash ATAC-seq/scripts/run_atac_te_tracks_from_fastq.sh \
  --sample SAMPLE_1 \
  --layout PE \
  --fastq1 /path/to/SAMPLE_1_R1.fastq.gz \
  --fastq2 /path/to/SAMPLE_1_R2.fastq.gz \
  --outdir /path/to/te_tracks_from_fastq \
  --species hg38 \
  --mapq 0 \
  --normalization CPM \
  --binsize 10 \
  --cores 8
```

批量模式的 samplesheet 列为 `sample,layout,r1,r2,condition,replicate`：

```bash
bash ATAC-seq/scripts/run_atac_te_tracks_from_fastq.sh \
  --samplesheet /path/to/samplesheet.csv \
  --outdir /path/to/te_tracks_from_fastq \
  --species hg38 \
  --cores 16
```

输出 pre-clean `02_align_sorted/*.sorted.bam`、relaxed BAM/逐步 counts/flagstat 和 `04_bw_te/*.te.bw`。优先使用 batch 模式避免 PE 配对和 sample ID 错误。fastp 长时间无输出时检查进程、I/O、日志和 `--fastp-timeout`，不要直接删除半成品后重复提交多个任务。

## 单个 TE/L1 heatmap

```bash
bash ATAC-seq/scripts/run_te_heatmap.sh \
  --bw-glob "/path/to/te_tracks/04_bw_te/*.te.bw" \
  --te-bed /path/to/L1_5.5kb.bed \
  --outdir /path/to/L1_heatmap \
  --upstream 2000 \
  --downstream 2000 \
  --body-length 5500 \
  --skip-zeros false \
  --sort-regions keep \
  --cores 16
```

完整 L1 可用 5500 bp body；其他 TE 不应机械沿用。`scale-regions` 把不同长度元素缩放，profile 横轴不是统一真实坐标。

## 批量 global/contrast heatmap

```bash
bash ATAC-seq/scripts/run_te_heatmap_batch.sh \
  --bw-glob "/path/to/04_bw_te/*.te.bw" \
  --te-bed /path/to/te.bed \
  --global-peak-bed /path/to/consensus_peaks.bed \
  --contrast-bed-glob "/path/to/contrast_beds/*/*.bed" \
  --outdir /path/to/te_heatmap_batch \
  --max-regions 5000 \
  --skip-zeros false \
  --sort-regions keep \
  --cores 16
```

## 解释底线

- `skipZeros=true` 会按结果删除无信号元素，可能夸大平均曲线；
- `missingDataAsZero` 与真正测得 0 不同；
- 按 anchor track 排序后展示同一 anchor 会强化视觉梯度；
- relaxed track 不做 EM/fractional assignment，不能自动定位唯一 TE locus；
- 必须报告 annotation、MAPQ、flags、duplicates、normalization、bin size 和 region selection。
