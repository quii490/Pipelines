# ATAC-seq 统一入口与命令参考

| 状态 | 事实来源 | 最后验证 |
|---|---|---|
| Draft review | `ATAC-seq/scripts/atacseq` 与各 wrapper `--help` | 2026-07-16 |

## 推荐调用方式

在仓库根目录直接调用最不容易混淆版本：

```bash
bash ATAC-seq/scripts/atacseq help
bash ATAC-seq/scripts/atacseq <subcommand> --help
```

也可以把常用入口链接到当前环境：

```bash
conda activate atac_core
bash ATAC-seq/scripts/install_bin_tools.sh "$CONDA_PREFIX/bin"
command -v atacseq
atacseq help
```

符号链接仍指向仓库脚本；移动或删除仓库后链接会失效。多人共用环境时，不要把个人工作副本链接进公共 `bin/`。

## Subcommand 选择

| 命令 | 实际脚本 | 什么时候用 |
|---|---|---|
| `atacseq run` | `run_auto_atacseq.sh` | 正常 FASTQ 完整运行 |
| `atacseq init` | `run_auto_atacseq.sh --init-only` | 只生成 samplesheet/contrast 模板 |
| `atacseq pipeline` | `run_pipeline.sh` | 已有正式 samplesheet，直接运行 Nextflow |
| `atacseq from-bam` | `run_atac_from_bam.sh` | 已有 clean BAM，重建主要下游 |
| `atacseq downstream` | `run_atac_downstream_only.sh` | counts/peaks 不变，只改 contrast 或下游参数 |
| `atacseq callpeak` | `run_callpeak_from_bam.sh` | 只改 MACS3 参数或从 BAM 补 peaks |
| `atacseq annotate` | `run_peak_annotation.sh` | 对任意 peak BED 做 gene/TE annotation |
| `atacseq heatmap` | `run_region_heatmap.sh` | 在任意 BED 上展示 bigWig 信号 |
| `atacseq te` | `run_te_heatmap_batch.sh` | global/contrast TE heatmap |
| `atacseq te-tracks` | `run_atac_te_tracks_from_bam.sh` | 从 pre-clean BAM 生成 relaxed TE tracks |
| `atacseq tss` | `run_tss_qc.sh` | 重画 TSS aggregate profile |
| `atacseq gene-body` | `run_gene_body_profile.sh` | 重画 gene-body aggregate profile |
| `atacseq motif` | `run_motif_homer.sh` | 对 contrast BED 运行 HOMER |
| `atacseq footprint` | `run_tobias.sh` | TOBIAS footprinting |
| `atacseq nuc` | `run_nuc_phasing.sh` | PE fragment periodicity/NRL |
| `atacseq report` | `generate_qc_report.py` | 从已有 QC 表重建 Markdown 报告 |

## 常用命令

```bash
# 1. FASTQ：先初始化并检查自动识别
atacseq init \
  --fastq-dir /path/to/fastq \
  --species hg38 \
  --outdir /path/to/results

# 2. 正式运行
atacseq run \
  --mode auto \
  --fastq-dir /path/to/fastq \
  --metadata-csv /path/to/metadata.csv \
  --contrast KO,WT \
  --species hg38 \
  --outdir /path/to/results \
  --preset standard --resume

# 3. 修改 contrast 后只重跑下游
atacseq downstream \
  --result-dir /path/to/results \
  --species hg38 --levels both \
  --contrast-file /path/to/contrasts.csv \
  --outdir /path/to/results/08_downstream_rerun

# 4. 对候选 enhancer 画信号
atacseq heatmap \
  --regions /path/to/enhancers.bed \
  --bw-glob "/path/to/results/04_bw/*.bw" \
  --outdir /path/to/results/manual_heatmap \
  --name enhancers --before 2000 --after 2000 --threads 8
```

## `init`、`run` 与 `pipeline` 的区别

- `init` 不开始分析，适合先审查 FASTQ 配对、sample ID、condition 和 contrast。
- `run` 是日常首选自动化入口；它负责初始化、调用 pipeline 和组织结果。
- `pipeline` 假定 samplesheet 已正确，适合高级用户和可复现批处理；不会替你修正实验设计。

## 直接脚本还是统一入口

文档中的完整路径最适合可复现记录：

```bash
bash /path/to/Pipelines/ATAC-seq/scripts/run_region_heatmap.sh ...
```

短命令适合交互操作：

```bash
atacseq heatmap ...
```

无论选择哪一种，都在分析记录中保存：仓库 commit、完整命令、输入路径、输出目录和当前 `--help`。不要只记录“运行了 atacseq”。

## 常见误用

- `from-bam` 要求 clean BAM；想恢复已被过滤的 multi-mappers，必须回到 pre-clean BAM/FASTQ。
- `downstream` 不会重做 peaks/counts；改变 peak calling 后不能只运行它。
- `report` 只整理已有指标，不会重新计算缺失 QC。
- `heatmap`、`tss` 和 `gene-body` 使用 bigWig，结果依赖 track normalization 和过滤策略。
- `pipeline` 参数属于低层入口，不要把 `run` 的参数未经核对直接复制过去。
