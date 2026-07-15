# CUT&RUN 注释、Peak overlap 与通路工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft / under development | CUT&RUN maintainers | 2026-07-16 |

## Gene 与 TE annotation

默认 backend 是 ChIPseeker：

```bash
CUTnRUN/tools/cutrun_cli/cutrun_annotate_peaks \
  --input /path/to/sample_peaks.narrowPeak \
  --ref hg38 \
  --gtf /path/to/genes.gtf \
  --out-prefix /path/to/annotation/sample \
  --promoter-up 2000 --promoter-down 500 \
  --te-min-overlap 50
```

| 参数组 | 参数 | 用途 |
|---|---|---|
| 必需 | `--input`、`--ref`、`--out-prefix` | peak BED、assembly、输出前缀 |
| gene | `--gtf`、`--promoter-up/down` | TxDb 和 promoter 定义 |
| stringent TE | `--te-min-length`、`--line-min-length`、`--te-min-overlap` | 过滤短注释和最低 overlap |
| relaxed TE | `--relaxed-window` | 邻近/宽松候选，不等同严格 overlap |
| 输出 | `--no-plots` | 只生成表格 |

如果要复现实验性的 legacy ChIPpeakAnno backend：

```bash
PEAK_ANNOTATOR=chippeakanno \
  CUTnRUN/tools/cutrun_cli/cutrun_annotate_peaks \
  --input /path/to/peaks.broadPeak \
  --ref hg38 --out-prefix /path/to/annotation/sample
```

backend、GTF、assembly、promoter window 和 TE filter 都必须写入方法记录。最近基因、promoter overlap 或 relaxed TE proximity 都不是功能因果证据。

## 多个 peak 集合 overlap

`--peak` 必须是可重复的 `LABEL=PATH`，不能只给裸路径：

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_peak_overlap.py \
  --peak WT=/path/to/WT.narrowPeak \
  --peak KO=/path/to/KO.narrowPeak \
  --peak Rescue=/path/to/Rescue.narrowPeak \
  --out-dir /path/to/peak_overlap
```

输出 `pairwise_overlap.tsv` 同时给出 A 被 B 覆盖、B 被 A 覆盖、各自比例和 overlap bases。它使用至少 1 bp overlap，不做距离扩展、reciprocal fraction 或统计显著性检验。比较前统一 caller、peak type 和 threshold；broad 与 narrow 不应直接混合。

## Replicate-supported consensus

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_consensus_peaks.py \
  --results-dir /path/to/results \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --out-dir /path/to/results/09_downstream/consensus_peaks \
  --min-support 2
```

`--min-support 2` 表示至少两个 target replicates 支持。只有一个重复时应 SKIP 或清楚标记 exploratory；不要把参数改成 1 后仍称为 reproducible peaks。

## GO BP / Reactome enrichment

```bash
Rscript CUTnRUN/tools/cutrun_cli/cutrun_pathway_enrichment.R \
  --annotation /path/to/annotation/sample.gene_structure.tsv \
  --out-prefix /path/to/enrichment/sample \
  --organism hg38 \
  --min-genes 10
```

输入表必须有 `gene_id`。工具运行 GO Biological Process ORA；仅在安装 ReactomePA 且 ID 条件满足时产生 Reactome。少于 `--min-genes`、缺包或 ID 列时写 `SKIP` status，而不是空结果冒充“没有富集”。

ORA 的 gene universe 不应默认为全基因组：正式解释时应使用实验中可检测/可注释的候选 universe，并记录多对一 peak-to-gene 映射规则。当前 helper 是便捷汇总，不替代完整的 enrichment design 审核。

## HOMER

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_homer.sh \
  --input /path/to/peaks.narrowPeak \
  --ref hg38 \
  --out-dir /path/to/homer/sample \
  --threads 8
```

wrapper 会保存命令；未安装 HOMER 时写 SKIP。motif enrichment 受 peak 数量、宽度、GC 和背景集合影响，不能单独证明 target 或 TF 直接结合。
