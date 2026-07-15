# CUT&RUN 信号、peaks 与 annotation

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

## BAM correlation

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_bam_cor.sh \
  --bam-dir /path/to/04_clean_bam \
  --out-prefix /path/to/qc/clean_bam \
  --threads 8 \
  --bin-size 10000
```

比较同一分支和过滤策略。不要把 clean BAM 与 TE BAM 放在同一次 correlation。

## bigWig correlation

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_bw_cor.sh \
  --bw-dir /path/to/05_tracks \
  --out-prefix /path/to/qc/standard_tracks \
  --threads 8 \
  --bin-size 10000
```

RPGC、RPKM、CPM 不混合。PCA/相关性反映全局 signal，不代替 target-specific enrichment QC。

## 通用 deepTools heatmap

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_dt_heatmap.sh \
  --regions /path/to/regions.bed \
  --signals "/path/to/A.bw /path/to/B.bw" \
  --labels "A B" \
  --out-prefix /path/to/heatmap/regions \
  --mode reference-point \
  --reference-point center \
  --before 3000 --after 3000 \
  --bin-size 50 \
  --threads 8
```

`--skip-zeros` 和 `--missing-data-as-zero` 都是 opt-in；它们会改变行集合或数值含义。`--write-matrix-tsv` 可能生成很大文件，日常保留压缩 matrix 即可。

## Consensus peaks

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_consensus_peaks.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --out-dir /path/to/results/09_downstream/consensus_peaks \
  --min-support 2
```

`min-support=2` 需要至少两个 target replicates；重复不足时应 `SKIP`。broad 与 narrow 分开构建，不合并为一个 universe。

## Peak overlap

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_peak_overlap.py \
  --peak A=/path/to/A.narrowPeak \
  --peak B=/path/to/B.narrowPeak \
  --out-dir /path/to/peak_overlap
```

Overlap 受 peak width 和 threshold 强烈影响。至少报告 A 被 B 覆盖比例、B 被 A 覆盖比例和对称指标，不只报告交集条数。

## HOMER

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_homer.sh \
  --input /path/to/peaks.narrowPeak \
  --ref hg38 \
  --out-dir /path/to/homer \
  --threads 8
```

未安装 HOMER 时 wrapper 写 SKIP。motif enrichment 不证明 direct binding，需结合 target、control、signal 与其他实验。

## Peak annotation

主 downstream 默认使用 ChIPseeker。完整参数、backend 切换、TE filtering 与通路命令见[注释、overlap 与通路](annotation-enrichment.md)。`cutrun_annotate_peaks` 需要对应 R 包与正确 GTF/TxDb。最近基因仅是注释规则，不等于功能靶基因；报告 annotation backend、genome build 和 promoter 定义。
