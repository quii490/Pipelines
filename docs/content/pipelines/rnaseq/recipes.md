# RNA-seq：按现有文件选择入口

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | RNA-seq maintainers | 2026-07-16 |

## 我有 FASTQ 或下载 manifest

使用 `run_auto_rnaseq.sh`。FASTQ 已在本地用 `--fastq-dir`；需要从 accession/URL 下载用 `--manifest`，并设置可写的 `--download-dir`。先 `--dry-run`，再正式运行。

## 我有完整 results，只想改分组或重画图

1. 备份并编辑 `condition.csv` 与 `contrast.csv`；
2. 使用 `--downstream`；
3. 用新的 `--plot-outdir` 保留旧图；
4. 用 `--only-tools` 限定模块，先确认模块名来自当前输出。

```bash
cd /path/to/Pipelines/RNA-seq/rnaseq
bash run_auto_rnaseq.sh \
  --results-dir /path/to/results \
  --plot-outdir /path/to/results/plots_v2 \
  --species hg38 \
  --downstream \
  --only-tools gene_featureCounts,TE_TEtranscripts
```

## 我只有 BAM

BAM 不能直接生成 volcano/GSEA。先用 `rnaseq-bam-to-counts` 得到 raw count matrix，再用 `rnaseq-counts-to-de`。使用前确认：

- BAM coordinate-sorted，`samtools quickcheck` 通过；
- BAM 与 gene/TE GTF 属于同一 assembly；
- PE/SE 与 strandedness 已知；
- BAM 是否经过 dedup、是否保留 multi-mapping 有记录；
- TE locus 分析不能从只保留 unique reads 的 BAM“恢复”多重比对信息。

## 我有 raw count matrix

使用 `rnaseq-counts-to-de`。第一列为 feature ID，其他列是 raw integer-like counts；列名必须匹配 sample table。不要输入 TPM、FPKM、CPM、VST 或 batch-corrected matrix。

## 我只有 DE matrix

- 补 volcano/MA/注释：`rnaseq-de-visuals`；
- gene pathway：`rnaseq-pathway-de`；
- TE class/family/age：`rnaseq-te-analysis`。

DE matrix 应包含 `log2FoldChange`、`pvalue/padj` 和 abundance（如 `baseMean`）。若原模型和 contrast 不清楚，不建议只凭表重新解释。

## 我有 bigWig 或两样本矩阵

- 多个 bigWig：`rnaseq-bw-cor`；所有 track 必须使用同一 assembly、bin size 逻辑和 normalization；
- 两列样本表达：`rnaseq-two-sample-scatter`；相关性高不等于不存在系统性差异。

## 我想清理或交付结果

- `rnaseq-clean-work` 先预览空间，确认无任务运行后才 `--confirm`；
- `rnaseq-publish` 只整理报告、图、状态和版本，不应发布 FASTQ/BAM/work。
