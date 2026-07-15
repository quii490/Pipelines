# ATAC-seq 报告、汇总与跨对比工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

## 重建 QC 报告：`generate_qc_report.py`

当 `03_qc/atac_qc_summary_with_metrics.tsv` 和 `qc_pass_fail.tsv` 已存在，只需重建中文说明时运行：

```bash
python3 ATAC-seq/scripts/generate_qc_report.py \
  --outdir /path/to/results \
  --species hg38
```

输出 `/path/to/results/QC_REPORT.md`。该脚本只汇总已有表格；如果 FRiP、NRF/PBC 或 fragment metrics 缺失，应先补计算而不是期待报告脚本生成数据。

| 参数 | 必需 | 说明 |
|---|---:|---|
| `--outdir` | 是 | pipeline results 根目录 |
| `--species` | 否 | 报告中显示的 assembly；默认 `unknown` |

## 一次生成多类 profile：`run_atac_profile_heatmaps.sh`

适合已有标准 bigWig 后，一次生成 TSS、gene body、TE 和 L1 aggregate 图：

```bash
bash ATAC-seq/scripts/run_atac_profile_heatmaps.sh \
  --result-dir /path/to/results \
  --species hg38 \
  --features tss,gene_body,te,l1 \
  --preset standard \
  --cores 12
```

| 参数 | 默认 | 选择建议 |
|---|---:|---|
| `--preset` | `standard` | `quick` 每类 1,000 regions；`standard` 3,000；`full` 10,000 |
| `--features` | 脚本默认集合 | 逗号分隔 `tss,gene_body,te,l1` |
| `--te-flank` | `15000` | TE/L1 5′ anchor 两侧范围；短 TE 不应机械使用超大范围 |
| `--body-length` | 脚本默认值 | gene/TE scale-regions 的显示长度，不是真实长度 |
| `--max-*-regions` | 随 preset | 设定 deterministic display cap，影响展示集合 |
| `--write-matrix-tab` | `false` | `true` 会输出很大的明文矩阵 |

也可显式提供 `--bw-glob`、`--samplesheet`、`--tss-bed`、`--gene-body-bed`、`--te-bed` 和 `--l1-bed`。所有 BED 与 bigWig 必须属于同一 assembly，并保存实际抽样后的 region 文件。

## 跨 contrast heatmap：`run_cross_contrast_heatmap.sh`

该工具读取已有 differential region 表，不重新运行 edgeR：

```bash
bash ATAC-seq/scripts/run_cross_contrast_heatmap.sh \
  --level-dir /path/to/results/08_downstream/peak_level \
  --contrasts KO_vs_WT,Rescue_vs_KO \
  --top-n 150 \
  --annotation-mode gene_te \
  --outdir /path/to/results/08_downstream/cross_contrast
```

| 参数 | 默认 | 说明 |
|---|---:|---|
| `--level-dir` | 必需 | 含多个 contrast 差异表的 peak/bin level 目录 |
| `--contrasts` | 必需 | 逗号分隔，名称与目录/表一致且方向固定 |
| `--top-n` | `100` | 用于展示的 region 数量，不改变原始 DE 结果 |
| `--annotation-mode` | `gene_te` | `gene_te`、`gene` 或 `none` |
| `--output-prefix` | 自动 | 自定义结果前缀 |

先确认所有 contrast 的 logFC 方向。例如 `KO_vs_WT` 和 `WT_vs_KO` 不能混在同一图中而不翻转符号。

## 输出整理器

`organize_atac_downstream_outputs.py` 将 R 原始输出复制到稳定的 `figures/`、`results/`、`reports/` 结构，并生成文件索引。它通常由 downstream wrapper 自动调用，不建议普通用户直接运行。原始结果保留在 `legacy/raw_r_output`，因此整理不是统计重算，也不应删除原始目录。

`prepare_manifest_fastq.py`、`plot_atac_qc.R` 和 `run_peak_annotation_chipseeker.R` 同样属于内部后端。调试时先运行其上层 wrapper；只有开发接口时才直接调用，并同步更新 wrapper 的 `--help` 与本站文档。
