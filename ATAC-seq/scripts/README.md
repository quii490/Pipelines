# scripts 工具使用说明

普通使用者只需阅读项目根目录的 `README.md`。本文档是逐脚本调试和开发维护参考，不是运行主流程的必读入口。

本 README 只说明 `scripts/` 目录中的独立工具。主流程入口仍然是项目根目录下的 `run_pipeline.sh` 和 `main.nf`；这里的脚本主要用于：

1. 已经有 BAM/bigWig/peak/count 结果时，单独补跑某个分析；
2. 调试某一步，例如 TSS QC、TE heatmap、motif、TOBIAS、nucleosome phasing；
3. 对 ATAC-seq 的 TE/L1 relaxed track 做专门诊断。

> 重要：`*.before_te_atac` 是修改前备份文件，不是正式工具，正常情况下不要运行。

---

## 0. 基本使用方式

在项目根目录下运行：

```bash
conda activate atac_core   # 或你的 ATAC/ChIP 环境
```

多数工具推荐用 `bash scripts/<tool>.sh ...` 运行，例如：

```bash
bash scripts/run_tss_qc.sh \
  --bw-glob "results/04_bw/*.bw" \
  --species hg38 \
  --outdir results/03_qc/tss \
  --cores 8
```

也可以安装为命令：

```bash
bash scripts/install_bin_tools.sh
```

安装后，一部分脚本可以直接作为命令运行，例如：

```bash
run_tss_qc --help
run_te_heatmap --help
run_atac_te_tracks_from_bam --help
```

当前 `install_bin_tools.sh` 会链接这些命令：

```text
atacseq
run_tss_qc
run_gene_body_profile
run_te_heatmap
run_te_heatmap_batch
run_atac_te_tracks_from_bam
run_motif_homer
run_tobias
run_atac_downstream_only
run_atac_from_bam
run_fixedbin_from_bam
species_config_lib
run_nuc_phasing
```

没有被链接的脚本仍然可以用 `bash scripts/<script>.sh` 运行。

推荐日常优先使用统一入口：

```bash
atacseq run --fastq-dir FASTQ_DIR --species hg38 --outdir RESULTS \
  --metadata-csv metadata.csv --contrast KO,WT --preset standard

atacseq downstream --result-dir RESULTS --species hg38 \
  --contrast-file RESULTS/_automation/inputs/contrasts.csv \
  --outdir RESULTS/08_downstream_manual

atacseq heatmap --regions target.bed --bw-glob "RESULTS/04_bw/*.bw" \
  --outdir RESULTS/manual_heatmap --name target --before 3000 --after 3000

bash scripts/run_cross_contrast_heatmap.sh \
  --level-dir RESULTS/08_downstream/peak_level \
  --contrasts KO_vs_WT,Rescue_vs_KO \
  --top-n 100 \
  --annotation-mode gene_te
```

`run_cross_contrast_heatmap.sh` 是一个轻量的独立下游工具。它直接读取已经生成的 differential region 表和 peak annotation 表，不重新计算 edgeR，也不重新运行 FASTQ/BAM 流程。输出包括带 region 行名和 `peak_class`/变化模式行注释的 PDF/PNG，以及对应的 region annotation、logFC matrix 和 selected contrast 文件。

使用 `--annotation-mode gene` 可只保留 GTF/ChIPseeker 基因注释，不使用 TE 名称、family、class 和 overlap；`--annotation-mode none` 则只保留 peak/bin ID。

---

## 1. 物种配置：`species_config_lib.sh`

### 作用

内部辅助脚本，用于从 `conf/species_refs.config` 中读取参考文件路径，例如：

- `blacklist`
- `tss_bed`
- `gene_body_bed`
- `te_bed`
- `gtf_genes`
- `chrom_sizes`
- `genome_fasta`
- `effective_genome_size`
- `mito_chr`
- `genome_size`

### 直接使用

通常不需要单独运行。其他脚本会自动 source：

```bash
source scripts/species_config_lib.sh
```

示例：

```bash
source scripts/species_config_lib.sh
get_species_param hg38 blacklist conf/species_refs.config
get_species_param hg38 te_bed conf/species_refs.config
```

### 输出

打印对应参数值，例如 blacklist BED 路径。

---

## 2. 安装脚本命令：`install_bin_tools.sh`

### 作用

把常用 shell 脚本链接到 `$CONDA_PREFIX/bin` 或指定目录，使其可以直接作为命令调用。

### 用法

```bash
bash scripts/install_bin_tools.sh
```

指定安装目录：

```bash
bash scripts/install_bin_tools.sh /path/to/bin
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| 第 1 个位置参数 | 否 | `$CONDA_PREFIX/bin` | 命令链接安装目录 |

### 输出

在目标目录生成符号链接，例如：

```text
run_tss_qc -> scripts/run_tss_qc.sh
run_te_heatmap -> scripts/run_te_heatmap.sh
run_atac_te_tracks_from_bam -> scripts/run_atac_te_tracks_from_bam.sh
```

---

## 3. 从主流程结果补跑 downstream：`run_atac_downstream_only.sh`

### 作用

当 ATAC 主流程已经生成 consensus peak counts 后，单独补跑 downstream R 分析。

它会读取：

```text
<result-dir>/07_counts/consensus_peak_counts.txt
<result-dir>/07_counts/sample_metadata.csv
<result-dir>/06_consensus_peaks/consensus_peaks.bed
```

然后调用：

```text
atacseq-downstream/run_downstream_atac.R
```

### 用法

```bash
bash scripts/run_atac_downstream_only.sh \
  --result-dir results_atac \
  --species hg38 \
  --contrast-file contrast.csv \
  --outdir results_atac/08_downstream_manual
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--result-dir` | 是 | 无 | 已完成主流程的结果目录 |
| `--species` | 否 | `hg38` | `hg38`、`mm10`、`mm39` |
| `--contrast-file` | 否 | 空 | 对比文件，列通常为 `case,control` |
| `--outdir` | 否 | `<result-dir>/08_downstream_manual` | downstream 输出目录 |
| `--te-bed` | 否 | 从 species config 读取 | TE 注释 BED/GTF |
| `--gtf` | 否 | 从 species config 读取 | 基因注释 GTF |
| `--te-violin-class` | 否 | 空 | 指定 TE class，例如 `LINE`、`LTR` |
| `--te-violin-top-n` | 否 | `12` | violin plot 展示 top N TE |
| `--te-family-filter` | 否 | 空 | 只分析匹配 family 的 TE |
| `--te-name-filter` | 否 | 空 | 只分析匹配 repName/name 的 TE |

### 输出

输出到 `--outdir`，通常包含：

```text
peak_level/
TE相关统计图/
差异 peak 表格/
PCA/correlation/volcano/heatmap 等 downstream 图表
```

具体文件名由 `run_downstream_atac.R` 决定。

---

## 4. 从已有 clean BAM 重建 ATAC downstream：`run_atac_from_bam.sh`

### 作用

从已经生成的 `*.clean.bam` 开始，重新做 ATAC downstream。适合：

- 不想重新从 FASTQ 跑；
- 已经有 clean BAM，只想重新 peak calling/count/downstream；
- 补跑 bigWig、TSS、gene body、TE heatmap、motif、TOBIAS、fixed bin 等分析；
- 新增 TE/L1 relaxed track 诊断。

### 基本用法

```bash
bash scripts/run_atac_from_bam.sh \
  --bam-dir results/02_align \
  --sample-meta sample_metadata.csv \
  --contrast-file contrast.csv \
  --species hg38 \
  --outdir results_from_bam \
  --cores 16
```

`--bam-dir` 需要包含：

```text
*.clean.bam
*.clean.bam.bai
```

### 启用 TE/L1 relaxed track

TE/L1 relaxed track 必须从 pre-clean sorted/raw BAM 生成，不能从已经 strict clean 的 BAM 恢复。

```bash
bash scripts/run_atac_from_bam.sh \
  --bam-dir results/02_align \
  --sample-meta sample_metadata.csv \
  --contrast-file contrast.csv \
  --species hg38 \
  --outdir results_from_bam_te \
  --run-te-heatmap true \
  --run-te-relaxed-tracks true \
  --te-source-bam-glob "results/02_align/*.sorted.bam" \
  --te-bed /path/to/L1_5.5kb.bed \
  --te-track-normalization CPM \
  --te-mapq 0 \
  --te-exclude-flags 780 \
  --cores 16
```

### 参数

#### 必需参数

| 参数 | 说明 |
|---|---|
| `--bam-dir` | 包含 `*.clean.bam` 的目录 |
| `--sample-meta` | 样本信息，列：`sample,condition,replicate` |
| `--contrast-file` | 对比文件，列：`case,control` |

#### 常规参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--species` | `hg38` | `hg38`、`mm10`、`mm39` |
| `--outdir` | `./atac_from_bam` | 输出目录 |
| `--cores` | `8` | 线程数 |
| `--consensus-half-width` | `250` | summit 扩展到 consensus peak 时的半宽 |
| `--background` | `false` | 后台运行 |

#### 参考文件参数

| 参数 | 说明 |
|---|---|
| `--genome` | genome fasta，用于 TOBIAS 等 |
| `--chrom-sizes` | 染色体大小文件，用于 fixed bin |
| `--tss-bed` | TSS BED |
| `--gene-body-bed` | gene body BED |
| `--te-bed` | TE 注释 BED/GTF，L1 分析时可传 L1 BED |
| `--motif-genome` | HOMER genome 名称，例如 `hg38` |
| `--motif-meme` | motif MEME 文件 |
| `--blacklist` | blacklist BED |

未显式提供时，脚本会尝试从 `conf/species_refs.config` 读取。

#### 分析开关

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--run-bamcoverage` | `true` | 从 clean BAM 生成普通 ATAC bigWig |
| `--run-tss` | `false` | 运行 TSS enrichment QC |
| `--run-gene-body` | `false` | 运行 gene body profile |
| `--run-te-heatmap` | `false` | 运行 TE heatmap |
| `--run-te-relaxed-tracks` | `false` | 从 pre-clean BAM 生成 TE/L1 relaxed BAM+bigWig |
| `--run-motif` | `false` | 运行 HOMER motif |
| `--run-tobias` | `false` | 运行 TOBIAS footprint |
| `--run-nuc-phasing` | `false` | 运行 nucleosome phasing |
| `--run-fixedbin` | `true` | 运行 fixed-bin count/downstream |
| `--run-bw-correlation` | `true` | 运行 bigWig correlation/PCA |
| `--run-diff-peak-heatmap` | `false` | 对差异 peak 画 heatmap |
| `--run-consensus-annotation` | `true` | consensus peak 注释 |

#### TE/L1 relaxed 参数

| 参数 | 默认值 | 说明 |
|---|---|---|
| `--te-source-bam-glob` | 空 | pre-clean sorted/raw BAM glob；启用 `--run-te-relaxed-tracks true` 时必需 |
| `--te-bw-glob` | 空 | 已有 TE/L1 bigWig glob；可跳过生成直接画图 |
| `--te-track-normalization` | `CPM` | `CPM`、`RPGC`、`RPKM`、`BPM` |
| `--te-mapq` | `0` | TE/L1 分支 MAPQ 阈值 |
| `--te-exclude-flags` | `780` | SAM flag 过滤；默认不删除 duplicate-marked reads |
| `--te-run-markdup` | `false` | 是否重新 MarkDuplicates 并删除 duplicates |
| `--te-remove-mito` | `true` | 是否去除线粒体 reads |
| `--te-remove-blacklist` | `true` | 是否去除 blacklist |
| `--te-proper-pair-only` | `true` | PE 数据是否要求 proper pair |
| `--te-binsize` | `10` | TE/L1 bigWig bin size |
| `--te-heatmap-skip-zeros` | `false` | 是否给 computeMatrix 加 `--skipZeros` |
| `--te-heatmap-missing-data-as-zero` | `false` | 是否给 computeMatrix 加 `--missingDataAsZero` |
| `--te-heatmap-sort-regions` | `keep` | `keep`、`descend`、`ascend`、`no` |
| `--te-heatmap-sort-using` | `mean` | `mean`、`median`、`max`、`sum`、`region_length` |
| `--te-heatmap-purpose` | `global` | TE heatmap 目的，常用 `global` |

### 主要输出

```text
<outdir>/03_qc/                         QC 图和统计
<outdir>/04_bw/                         普通 ATAC bigWig
<outdir>/04_bw_te/                      TE/L1 relaxed bigWig
<outdir>/05_peaks/                      MACS3 peaks
<outdir>/06_consensus_peaks/            consensus peak BED/SAF
<outdir>/07_counts/                     consensus peak counts
<outdir>/08_downstream/                 downstream R 结果
<outdir>/02_align_te/                   TE/L1 relaxed BAM、flagstat、clean_counts
```

---

## 5. 从 BAM 生成 TE/L1 relaxed ATAC track：`run_atac_te_tracks_from_bam.sh`

### 作用

从 pre-clean sorted/raw BAM 生成 TE/L1 relaxed BAM 和 bigWig。这个工具专门用于 L1/TE heatmap 诊断。

它的设计目的不是替代普通 ATAC peak calling，而是检查完整 L1 body 的信号是否被 strict clean 流程系统性压低。

### 推荐用法

```bash
bash scripts/run_atac_te_tracks_from_bam.sh \
  --bam-glob "results/02_align/*.sorted.bam" \
  --outdir results_te_tracks \
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

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bam-glob` | 是 | 无 | source BAM glob，必须加引号；推荐 pre-clean sorted/raw BAM |
| `--species` | 否 | `hg38` | `hg38`、`mm10`、`mm39` |
| `--outdir` | 否 | `./atac_te_tracks` | 输出目录 |
| `--blacklist` | 否 | species config | blacklist BED |
| `--mito-chr` | 否 | species config 或 `chrM` | 线粒体染色体名 |
| `--effective-genome-size` | 否 | species config | RPGC 时需要 |
| `--cores` | 否 | `8` | 线程数 |
| `--mapq` | 否 | `0` | MAPQ 阈值 |
| `--exclude-flags` | 否 | `780` | SAM flag 过滤；默认不去 duplicate flag |
| `--run-markdup` | 否 | `false` | 是否删除 duplicates |
| `--remove-mito` | 否 | `true` | 是否去线粒体 |
| `--remove-blacklist` | 否 | `true` | 是否去 blacklist |
| `--proper-pair-only` | 否 | `true` | PE BAM 是否要求 proper pair |
| `--normalization` | 否 | `CPM` | `CPM`、`RPGC`、`RPKM`、`BPM` |
| `--binsize` | 否 | `10` | bigWig bin size |
| `--ignore-for-normalization` | 否 | 空 | RPGC 时传给 deepTools |

### 输出

```text
<outdir>/02_align_te/<sample>.te.bam
<outdir>/02_align_te/<sample>.te.bam.bai
<outdir>/02_align_te/<sample>.te.flagstat.txt
<outdir>/02_align_te/<sample>.te.clean_counts.tsv
<outdir>/04_bw_te/<sample>.te.bw
```

`*.te.clean_counts.tsv` 用于检查每一步保留了多少 alignments：

```text
input
post_markdup_step
post_mapq_flag_filter
post_mito_step
final
```

---

## 6. TE/L1 heatmap：`run_te_heatmap.sh`

### 作用

用 bigWig 在 TE/L1 注释区域上画 deepTools heatmap/profile。支持：

- 全部 TE/L1 区域；
- 与 peak 重叠的 TE/L1 区域；
- 按 TE class/family/name 过滤；
- 控制 `skipZeros`、`missingDataAsZero`、排序方式，避免误删零信号 TE 区域。

### 推荐用于 L1_5.5kb

```bash
bash scripts/run_te_heatmap.sh \
  --bw-glob "results_te_tracks/04_bw_te/*.te.bw" \
  --te-bed /path/to/L1_5.5kb.bed \
  --outdir results_te_tracks/L1_5.5kb_heatmap \
  --cores 16 \
  --upstream 2000 \
  --downstream 2000 \
  --body-length 5500 \
  --skip-zeros false \
  --sort-regions keep
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bw-glob` | 是 | 无 | bigWig glob，必须加引号 |
| `--species` | 否 | `hg38` | 用于从 species config 读取 TE 注释 |
| `--te-bed` | 否 | species config | TE BED/GTF；L1 分析建议显式传 L1 BED |
| `--peak-bed` | 否 | 空 | peak BED；提供后可筛选 peak-overlap TE |
| `--purpose` | 否 | `global` | `global`、`peak-overlap`、`contrast` |
| `--outdir` | 否 | `./te_heatmap` | 输出目录 |
| `--cores` | 否 | `4` | 线程数 |
| `--upstream` | 否 | `2000` | TE 区域上游长度 |
| `--downstream` | 否 | `2000` | TE 区域下游长度 |
| `--body-length` | 否 | `3000` | scale-regions 的 body length；完整 L1 推荐 `5500` |
| `--max-regions` | 否 | `0` | 随机抽样 region 数，`0` 表示不抽样 |
| `--mode` | 否 | `scale-regions` | 保留兼容参数，当前实际使用 scale-regions |
| `--label` | 否 | 空 | batch wrapper 兼容标签 |
| `--te-label-level` | 否 | `locus` | 兼容参数 |
| `--te-class-filter` | 否 | 空 | 根据第 4 列 grep 过滤 TE class/name 字段 |
| `--te-family-filter` | 否 | 空 | 根据第 4 列 grep 过滤 family/name 字段 |
| `--te-name-filter` | 否 | 空 | 根据第 4 列 grep 过滤 repName/name 字段 |
| `--skip-zeros` | 否 | `false` | 是否使用 computeMatrix `--skipZeros` |
| `--missing-data-as-zero` | 否 | `false` | 是否使用 computeMatrix `--missingDataAsZero` |
| `--sort-regions` | 否 | `keep` | `keep`、`descend`、`ascend`、`no` |
| `--sort-using` | 否 | `mean` | `mean`、`median`、`max`、`sum`、`region_length` |

### 输出

```text
<outdir>/te.locus.bed
<outdir>/te.filtered.bed
<outdir>/te.sampled.bed
<outdir>/TE_accessibility_matrix.gz
<outdir>/TE_accessibility_profile.pdf
<outdir>/TE_accessibility_profile.png
<outdir>/TE_accessibility_heatmap.pdf
<outdir>/TE_accessibility_heatmap.png
```

---

## 7. 批量 TE/L1 heatmap：`run_te_heatmap_batch.sh`

### 作用

批量调用 `run_te_heatmap.sh`，适合同时生成：

1. 全局 TE heatmap；
2. 每个 contrast-specific peak BED 上的 TE heatmap。

### 用法

```bash
bash scripts/run_te_heatmap_batch.sh \
  --bw-glob "results/04_bw_te/*.te.bw" \
  --te-bed /path/to/L1_5.5kb.bed \
  --global-peak-bed results/06_consensus_peaks/consensus_peaks.bed \
  --contrast-bed-glob "results/08_downstream/peak_level/contrast_beds/*/*.bed" \
  --outdir results/03_qc/te_heatmap_batch \
  --cores 16 \
  --max-regions 5000 \
  --skip-zeros false \
  --sort-regions keep
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bw-glob` | 是 | 无 | bigWig glob |
| `--te-bed` | 否 | species config | TE BED/GTF |
| `--species` | 否 | `hg38` | 物种 |
| `--global-peak-bed` | 否 | 空 | consensus peak BED；用于 global 子目录 |
| `--contrast-bed-glob` | 否 | 空 | contrast BED glob |
| `--outdir` | 否 | `./te_heatmap_batch` | 输出目录 |
| `--cores` | 否 | `4` | 线程数 |
| `--max-regions` | 否 | `5000` | 每次最多抽样 region 数 |
| `--mode` | 否 | `center` | 兼容参数 |
| `--te-label-level` | 否 | `locus` | 兼容参数 |
| `--te-class-filter` | 否 | 空 | TE class 过滤 |
| `--te-family-filter` | 否 | 空 | TE family 过滤 |
| `--te-name-filter` | 否 | 空 | TE name 过滤 |
| `--skip-zeros` | 否 | `false` | 是否跳过全零 regions |
| `--missing-data-as-zero` | 否 | `false` | 是否把缺失当 0 |
| `--sort-regions` | 否 | `keep` | region 排序方式 |
| `--sort-using` | 否 | `mean` | 排序依据 |

### 输出

```text
<outdir>/global/
<outdir>/contrast/<contrast_label>/
```

每个子目录内部输出同 `run_te_heatmap.sh`。

---

## 8. TSS enrichment QC：`run_tss_qc.sh`

### 作用

用 deepTools 在 TSS 附近计算 ATAC signal profile/heatmap，是 ATAC 质量控制常用指标。

### 用法

```bash
bash scripts/run_tss_qc.sh \
  --bw-glob "results/04_bw/*.bw" \
  --species hg38 \
  --outdir results/03_qc/tss \
  --cores 8 \
  --upstream 3000 \
  --downstream 3000
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bw-glob` | 是 | 无 | bigWig glob |
| `--tss-bed` | 否 | species config | TSS BED |
| `--species` | 否 | `hg38` | 物种 |
| `--outdir` | 否 | `./tss_qc` | 输出目录 |
| `--cores` | 否 | `4` | 线程数 |
| `--upstream` | 否 | `3000` | TSS 上游 |
| `--downstream` | 否 | `3000` | TSS 下游 |

### 输出

```text
<outdir>/<tss>.valid.bed
<outdir>/TSS_enrichment_matrix.gz
<outdir>/TSS_enrichment_profile.pdf
<outdir>/TSS_enrichment_profile.png
<outdir>/TSS_enrichment_heatmap.pdf
<outdir>/TSS_enrichment_heatmap.png
```

---

## 9. Gene body profile：`run_gene_body_profile.sh`

### 作用

用 deepTools 在 gene body 上计算 ATAC/bigWig signal profile/heatmap。

### 用法

```bash
bash scripts/run_gene_body_profile.sh \
  --bw-glob "results/04_bw/*.bw" \
  --species hg38 \
  --outdir results/03_qc/gene_body \
  --cores 8 \
  --upstream 3000 \
  --downstream 3000 \
  --body-length 5000
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bw-glob` | 是 | 无 | bigWig glob |
| `--gene-bed` | 否 | species config | gene body BED |
| `--species` | 否 | `hg38` | 物种 |
| `--outdir` | 否 | `./gene_body_profile` | 输出目录 |
| `--cores` | 否 | `4` | 线程数 |
| `--upstream` | 否 | `3000` | gene body 上游 |
| `--downstream` | 否 | `3000` | gene body 下游 |
| `--body-length` | 否 | `5000` | scale-regions body length |

### 输出

```text
<outdir>/<gene_body>.valid.bed
<outdir>/GeneBody_accessibility_matrix.gz
<outdir>/GeneBody_accessibility_profile.pdf
<outdir>/GeneBody_accessibility_profile.png
<outdir>/GeneBody_accessibility_heatmap.pdf
<outdir>/GeneBody_accessibility_heatmap.png
```

---

## 10. 通用 region heatmap：`run_region_heatmap.sh`

### 作用

对任意 BED 区域画 bigWig heatmap/profile。适合：

- differential peaks；
- motif peaks；
- promoter/enhancer；
- 自定义 TE subset；
- 某个蛋白 peak 或 ATAC peak 集合。

### 用法 1：使用 `--bw-glob`

```bash
bash scripts/run_region_heatmap.sh \
  --regions regions.bed \
  --bw-glob "results/04_bw/*.bw" \
  --outdir results/region_heatmap \
  --name my_regions \
  --mode reference-point \
  --reference-point center \
  --before 3000 \
  --after 3000 \
  --threads 8
```

### 用法 2：使用显式 bigWig 列表

```bash
bash scripts/run_region_heatmap.sh \
  --regions regions.bed \
  --signals "A.bw B.bw C.bw" \
  --labels "A B C" \
  --out-prefix results/region_heatmap/my_regions \
  --mode scale-regions \
  --body-length 5000
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--regions` / `-R` | 是 | 无 | BED/narrowPeak/broadPeak |
| `--signals` / `-S` | 二选一 | 空 | 空格分隔 bigWig 列表，需加引号 |
| `--bw-glob` | 二选一 | 空 | bigWig glob |
| `--out-prefix` / `-o` | 二选一 | 空 | 输出前缀 |
| `--outdir` + `--name` | 二选一 | 空 | 另一种输出命名方式 |
| `--labels` | 否 | 空 | 样本标签 |
| `--mode` | 否 | `reference-point` | `reference-point` 或 `scale-regions` |
| `--reference-point` | 否 | `center` | `center`、`TSS`、`TES` |
| `--before` / `-b` | 否 | `3000` | 上游/左侧 bp |
| `--after` / `-a` | 否 | `3000` | 下游/右侧 bp |
| `--body-length` | 否 | `5000` | scale-regions body length |
| `--sort-regions` | 否 | `descend` | `descend`、`ascend`、`no`、`keep` |
| `--sort-using` | 否 | `mean` | `mean`、`median`、`max`、`sum`、`region_length` |
| `--kmeans` | 否 | 空 | heatmap 聚类数 |
| `--z-min` | 否 | 空 | heatmap zMin |
| `--z-max` | 否 | 空 | heatmap zMax |
| `--threads` / `--cores` | 否 | `$SLURM_CPUS_PER_TASK` 或 `8` | 线程数 |

### 输出

```text
<base>/matrix/<name>.matrix.gz
<base>/matrix/<name>.matrix.tab
<base>/matrix/<name>.regions.sorted.bed
<base>/plots/<name>.heatmap.pdf
<base>/plots/<name>.heatmap.png
<base>/plots/<name>.profile.pdf
<base>/plots/<name>.profile.png
<base>/logs/<name>.computeMatrix.log
```

---

## 11. 差异 peak heatmap：`run_diff_peak_heatmap.sh`

### 作用

扫描 differential peak BED 文件，逐个调用 `run_region_heatmap.sh` 画 heatmap/profile。

### 用法

```bash
bash scripts/run_diff_peak_heatmap.sh \
  --contrast-bed-dir results/08_downstream/peak_level/contrast_beds \
  --bw-glob "results/04_bw/*.bw" \
  --outdir results/region_heatmap/differential_peaks \
  --before 3000 \
  --after 3000 \
  --threads 8
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--contrast-bed-dir` / `--bed-dir` | 是 | 无 | 差异 peak BED 目录 |
| `--bw-glob` | 二选一 | 空 | bigWig glob |
| `--signals` | 二选一 | 空 | 显式 bigWig 列表 |
| `--outdir` | 否 | `region_heatmap/differential_peaks` | 输出目录 |
| `--labels` | 否 | 空 | 样本标签 |
| `--before` / `-b` | 否 | `3000` | peak 中心前范围 |
| `--after` / `-a` | 否 | `3000` | peak 中心后范围 |
| `--threads` / `--cores` | 否 | `$SLURM_CPUS_PER_TASK` 或 `8` | 线程数 |

### 输出

```text
<outdir>/<contrast>/<bed_name>/matrix/
<outdir>/<contrast>/<bed_name>/plots/
<outdir>/<contrast>/<bed_name>/logs/
```

---

## 12. bigWig correlation/PCA：`run_bw_correlation.sh`

### 作用

用 deepTools 计算多个 bigWig 的 bin-level signal 相关性和 PCA。

### 用法

```bash
bash scripts/run_bw_correlation.sh \
  --bw-glob "results/04_bw/*.bw" \
  --out-prefix results/03_qc/bw_correlation/atac \
  --bin-size 10000 \
  --cor-method pearson \
  --threads 8
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bw-dir` | 二选一 | 空 | 扫描目录内 `.bw/.bigWig/.bigwig` |
| `--bw-glob` | 二选一 | 空 | bigWig glob |
| `--out-prefix` | 是 | 无 | 输出前缀 |
| `--bin-size` | 否 | `10000` | multiBigwigSummary bin size |
| `--cor-method` | 否 | `pearson` | `pearson` 或 `spearman` |
| `--threads` | 否 | `$SLURM_CPUS_PER_TASK` 或 `8` | 线程数 |

### 输出

```text
<base>/matrix/<name>.bins.npz
<base>/matrix/<name>.raw_counts.tsv
<base>/matrix/<name>.pearson_correlation.tsv
<base>/matrix/<name>.pca.tsv
<base>/plots/<name>.pearson_heatmap.pdf/png
<base>/plots/<name>.pearson_scatter.pdf/png
<base>/plots/<name>.pca.pdf/png
```

---

## 13. 从 BAM call peak：`run_callpeak_from_bam.sh`

### 作用

从 clean BAM 调用 MACS3 peak，自动识别 PE/SE：

- PE：使用 `BAMPE`
- SE：使用 ATAC shift/extsize model

同时可做 blacklist 过滤、summit 标准化、FRiP 统计。

### 用法

```bash
bash scripts/run_callpeak_from_bam.sh \
  --bam-glob "results/02_align/*.clean.bam" \
  --species hg38 \
  --outdir results_callpeak \
  --qvalue 0.05 \
  --keep-dup all \
  --cores 8
```

单个 BAM：

```bash
bash scripts/run_callpeak_from_bam.sh \
  --bam sample.clean.bam \
  --species hg38 \
  --outdir results_callpeak
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bam-dir` | 三选一 | 空 | 扫描目录下 BAM |
| `--bam-glob` | 三选一 | 空 | BAM glob |
| `--bam` | 三选一 | 空 | 单个 BAM |
| `--species` | 否 | `hg38` | 物种 |
| `--outdir` | 是 | 无 | 输出目录 |
| `--genome-size` | 否 | species config | MACS3 `-g` |
| `--blacklist` | 否 | species config | blacklist BED |
| `--format` | 否 | `auto` | `auto`、`BAMPE`、`BAM` |
| `--qvalue` / `-q` | 否 | `0.05` | MACS3 q-value cutoff |
| `--pvalue` / `-p` | 否 | 空 | 若提供，则覆盖 q-value |
| `--shift` | 否 | `-100` | SE ATAC shift |
| `--extsize` | 否 | `200` | SE ATAC extension |
| `--keep-dup` | 否 | `all` | MACS3 duplicate 策略 |
| `--broad` | 否 | `false` | call broad peaks |
| `--no-summits` | 否 | `false` | 不输出/复制 summit |
| `--cores` / `--threads` | 否 | `$SLURM_CPUS_PER_TASK` 或 `8` | 线程数 |

### 输出

```text
<outdir>/05_peaks/raw_macs3/
<outdir>/05_peaks/filtered/
<outdir>/05_peaks/summits/
<outdir>/05_peaks/qc/frip.tsv
```

---

## 14. Fixed-bin ATAC 分析：`run_fixedbin_from_bam.sh`

### 作用

把基因组切成固定大小 bins，对 BAM 计数，然后调用 downstream R 脚本做差异/可视化。

适合做不依赖 peak calling 的全基因组 ATAC 分析。

### 用法

```bash
bash scripts/run_fixedbin_from_bam.sh \
  --bam-glob "results/02_align/*.clean.bam" \
  --sample-meta sample_metadata.csv \
  --contrast-file contrast.csv \
  --species hg38 \
  --outdir results_fixedbin \
  --bin-size 100000 \
  --cores 8
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bam-glob` | 是 | 无 | BAM glob |
| `--sample-meta` | 是 | 无 | 样本信息，列：`sample,condition,replicate` |
| `--contrast-file` | 是 | 无 | 对比文件，列：`case,control` |
| `--species` | 否 | `hg38` | 物种 |
| `--chrom-sizes` | 否 | species config | 染色体大小文件 |
| `--te-bed` | 否 | species config | TE 注释 |
| `--gtf` | 否 | species config | GTF 注释 |
| `--outdir` | 否 | `./fixedbin_from_bam` | 输出目录 |
| `--bin-size` | 否 | `100000` | bin 大小 |
| `--cores` | 否 | `4` | 线程数 |
| `--te-violin-class` | 否 | 空 | TE violin class |
| `--te-violin-top-n` | 否 | `12` | TE violin top N |

### 输出

```text
<outdir>/counts/fixed_bins_<bin_size>.bed
<outdir>/counts/fixed_bins_<bin_size>.saf
<outdir>/counts/fixedbin_counts.txt
<outdir>/downstream 相关图表和表格
```

---

## 15. HOMER motif：`run_motif_homer.sh`

### 作用

对差异 BED 批量运行 HOMER motif enrichment，并可同时运行 `annotatePeaks.pl`。

脚本同时支持：

```text
<bed-dir>/<contrast>.up.bed
<bed-dir>/<contrast>.down.bed
<bed-dir>/<contrast>/up.bed
<bed-dir>/<contrast>/down.bed
```

例如：

```text
contrast_beds/KO_vs_WT/up.bed
contrast_beds/KO_vs_WT/down.bed
```

### 用法

```bash
bash scripts/run_motif_homer.sh \
  --bed-dir results/08_downstream/peak_level/results/differential/beds \
  --genome hg38 \
  --outdir results/motif_homer \
  --size 200 \
  --run-annotation true
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bed-dir` | 是 | 无 | 包含 flat 或 nested up/down BED 的目录 |
| `--genome` | 是 | 无 | HOMER genome 名称，例如 `hg38` |
| `--outdir` | 否 | `./motif` | 输出目录 |
| `--size` | 否 | `200` | 以 peak 中心为基准的 motif 分析窗口 |
| `--run-annotation` | 否 | `true` | 同时输出 HOMER peak annotation |

### 输出

```text
<outdir>/<contrast>/<direction>/knownResults.txt
<outdir>/<contrast>/<direction>/annotated_peaks.tsv
<outdir>/<contrast>/<direction>/homer.log
<outdir>/motif_summary.tsv
```

---

## 16. TOBIAS footprint：`run_tobias.sh`

### 作用

运行 TOBIAS ATACorrect、ScoreBigwig，并根据 contrast 文件尝试运行 BINDetect。

### 用法

```bash
bash scripts/run_tobias.sh \
  --bam-dir results/02_align \
  --peaks results/06_consensus_peaks/consensus_peaks.bed \
  --sample-meta sample_metadata.csv \
  --contrast-file contrast.csv \
  --genome /path/to/hg38.fa \
  --motif-meme conf/JASPAR2026_CORE_non-redundant_pfms_meme.txt \
  --outdir results/tobias \
  --cores 8
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bam-dir` | 是 | 无 | 包含 `*.clean.bam` 的目录 |
| `--peaks` | 是 | 无 | peak BED |
| `--sample-meta` | 是 | 无 | 样本信息 |
| `--contrast-file` | 否 | 空 | 对比文件；提供后运行 BINDetect |
| `--genome` | 是 | 无 | genome fasta |
| `--motif-meme` | 是 | 无 | motif MEME 文件 |
| `--outdir` | 否 | `./tobias` | 输出目录 |
| `--cores` | 否 | `4` | 线程数 |

### 输出

```text
<outdir>/corrected/<sample>_corrected.bw
<outdir>/footprints/<sample>_footprints.bw
<outdir>/bindetect/bindetect_sample_choice.tsv
<outdir>/bindetect/<case>_vs_<control>/
```

---

## 17. Nucleosome phasing：`run_nuc_phasing.sh`

### 作用

从 ATAC BAM 计算 fragment length phasing，并估计 nucleosome repeat length，输出每个样本的 phasing 图和 NRL 表。

### 用法

```bash
bash scripts/run_nuc_phasing.sh \
  --bam-glob "results/02_align/*.clean.bam" \
  --outdir results/03_qc/nuc_phasing \
  --mapq 30 \
  --max-frag 1000 \
  --cores 8 \
  --conda-env chipseq
```

单个 BAM：

```bash
bash scripts/run_nuc_phasing.sh \
  --bam sample.clean.bam \
  --outdir nuc_phasing_sample
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--bam-glob` | 三选一 | 空 | BAM glob |
| `--bam` | 三选一 | 空 | 单个 BAM |
| `--bam-list` | 三选一 | 空 | 每行一个 BAM 路径 |
| `--outdir` | 是 | 无 | 输出目录 |
| `--labels` | 否 | 文件名推断 | 逗号分隔标签；当前 wrapper help 中保留，但实际未传入 R 参数 |
| `--mapq` | 否 | `30` | MAPQ filter |
| `--max-frag` | 否 | `1000` | 最大 fragment size |
| `--lspan` | 否 | `0.35` | 第一轮 loess span |
| `--rspan` | 否 | `0.1` | residual loess span |
| `--cores` | 否 | `1` | 线程数 |
| `--conda-env` | 否 | `chipseq` | Rscript 所在 conda 环境名 |

### 输出

```text
<outdir>/*_nuc_phasing.pdf
<outdir>/NRL_summary.csv
```

---

## 18. Peak annotation：`run_peak_annotation.sh`

### 作用

用 ChIPseeker 对 peak 做基因组注释，并可额外统计 peak 与 TE 的重叠。

### 用法

```bash
bash scripts/run_peak_annotation.sh \
  --input results/06_consensus_peaks/consensus_peaks.bed \
  --species hg38 \
  --out-prefix results/annotation/consensus \
  --gtf /path/to/gencode.gtf \
  --te-bed /path/to/hg38_te.bed
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--input` / `-i` | 是 | 无 | BED/narrowPeak/broadPeak |
| `--species` | 否 | `hg38` | 物种 |
| `--out-prefix` / `-o` | 是 | 无 | 输出前缀 |
| `--gtf` | 否 | species config | GTF 注释 |
| `--te-bed` / `--te-anno` | 否 | species config | TE BED/GTF |
| `--promoter-up` | 否 | `3000` | promoter 上游定义 |
| `--promoter-down` | 否 | `3000` | promoter 下游定义 |
| `--txdb-cache` | 否 | 自动缓存 | TxDb sqlite 缓存路径 |
| `--no-plots` | 否 | `false` | 只输出表，不画图 |

### 输出

由 `run_peak_annotation_chipseeker.R` 生成，通常包括：

```text
<out-prefix>*.tsv / *.csv
<out-prefix>*annotation*.pdf/png
<out-prefix>*TE*.tsv
```

具体文件名取决于 R 脚本内部输出。

---

## 19. Peak annotation R 脚本：`run_peak_annotation_chipseeker.R`

### 作用

`run_peak_annotation.sh` 的底层 R 脚本。一般不需要直接调用，但可以单独运行。

### 用法

```bash
Rscript scripts/run_peak_annotation_chipseeker.R \
  --input peaks.bed \
  --gtf genes.gtf \
  --out-prefix annotation/sample \
  --species hg38 \
  --te-bed te.bed
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--input` | 是 | 无 | peak BED/narrowPeak/broadPeak |
| `--gtf` | 是 | 无 | GTF 文件 |
| `--out-prefix` | 是 | 无 | 输出前缀 |
| `--species` | 否 | `hg38` | 物种 |
| `--te-bed` | 否 | 空 | TE BED/GTF |
| `--promoter-up` | 否 | `3000` | promoter 上游 |
| `--promoter-down` | 否 | `3000` | promoter 下游 |
| `--txdb-cache` | 否 | 自动 | TxDb sqlite cache |
| `--no-plots` | 否 | `false` | 不生成图 |

### 输出

输出注释表、简化注释统计、TE overlap 表和可视化图。

---

## 20. Peak set overlap：`run_peak_overlap.sh`

### 作用

比较多个 peak set 的重叠情况，输出 bedtools multiinter 表、Jaccard heatmap 和 intersection barplot。

### 用法 1：显式指定 peak

```bash
bash scripts/run_peak_overlap.sh \
  --peaks "WT=WT.narrowPeak KO=KO.narrowPeak OE=OE.narrowPeak" \
  --outdir results/peak_overlap \
  --distance 0
```

### 用法 2：扫描目录

```bash
bash scripts/run_peak_overlap.sh \
  --peak-dir results/05_peaks/filtered \
  --outdir results/peak_overlap
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--peaks` | 二选一 | 空 | 空格分隔 `label=file` 列表 |
| `--peak-dir` | 二选一 | 空 | 扫描 `.bed/.narrowPeak/.broadPeak` |
| `--outdir` | 是 | 无 | 输出目录 |
| `--distance` | 否 | `0` | 比较前 merge 距离 |
| `--top-n` | 否 | 空 | 每个文件只取前 N 行 |

### 输出

```text
<outdir>/prepared/*.sorted.bed
<outdir>/tables/peak_multiinter.tsv
<outdir>/tables/peak_union_clusters.bed
<outdir>/tables/intersection_counts.tsv
<outdir>/plots/peak_jaccard_heatmap.pdf/png
<outdir>/plots/peak_intersection_bar.pdf/png
```

---

## 21. ATAC QC plotting：`plot_atac_qc.R`

### 作用

读取 ATAC QC summary，画 FRiP 和 mitochondrial fraction barplot。

### 用法

```bash
Rscript scripts/plot_atac_qc.R \
  --qc-summary results/03_qc/atac_qc_summary.tsv \
  --outdir results/03_qc
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--qc-summary` | 是 | 无 | 至少包含 `sample` 列；若含 `frip` 列会画 FRiP；若含 `total_clean,mt_reads` 会画 mt fraction |
| `--outdir` | 否 | `.` | 输出目录 |

### 输出

```text
<outdir>/FRiP_barplot.pdf/png
<outdir>/Mito_fraction_barplot.pdf/png
<outdir>/atac_qc_summary_with_metrics.tsv
```

---

## 22. FASTQ samplesheet 转 metadata：`prepare_manifest_fastq.py`

### 作用

从输入 samplesheet 中提取 downstream 所需的样本 metadata。

输入 samplesheet 至少建议包含：

```text
sample,condition,replicate,layout
```

### 用法

```bash
python scripts/prepare_manifest_fastq.py samplesheet.csv sample_metadata.csv
```

### 参数

| 参数 | 必需 | 说明 |
|---|---:|---|
| 第 1 个位置参数 | 是 | 输入 samplesheet CSV |
| 第 2 个位置参数 | 是 | 输出 metadata CSV |

### 输出

```text
sample,condition,replicate,layout,note
```

缺失的 `condition/replicate` 会填 `NA`，缺失的 `layout` 会填 `PE`。

---

## 23. 推荐组合流程

### A. 已有 clean BAM，补完整 downstream

```bash
bash scripts/run_atac_from_bam.sh \
  --bam-dir results/02_align \
  --sample-meta results/07_counts/sample_metadata.csv \
  --contrast-file contrast.csv \
  --species hg38 \
  --outdir results_from_bam \
  --run-tss true \
  --run-gene-body true \
  --run-te-heatmap true \
  --run-motif true \
  --run-tobias false \
  --cores 16
```

### B. 专门研究 L1_5.5kb 的 ATAC accessibility

```bash
bash scripts/run_atac_te_tracks_from_bam.sh \
  --bam-glob "results/02_align/*.sorted.bam" \
  --outdir results_L1_te_tracks \
  --species hg38 \
  --mapq 0 \
  --exclude-flags 780 \
  --run-markdup false \
  --normalization CPM \
  --binsize 10 \
  --cores 16

bash scripts/run_te_heatmap.sh \
  --bw-glob "results_L1_te_tracks/04_bw_te/*.te.bw" \
  --te-bed L1_5.5kb.bed \
  --outdir results_L1_te_tracks/L1_5.5kb_heatmap \
  --body-length 5500 \
  --upstream 2000 \
  --downstream 2000 \
  --skip-zeros false \
  --sort-regions keep \
  --cores 16
```

### C. 只重新画某一组 differential peaks 的 heatmap

```bash
bash scripts/run_diff_peak_heatmap.sh \
  --contrast-bed-dir results/08_downstream/peak_level/contrast_beds \
  --bw-glob "results/04_bw/*.bw" \
  --outdir results/replot_diff_peak_heatmap \
  --threads 16
```

---

## 24. L1/TE 分析注意事项

1. **不要用 strict `*.clean.bam` 试图恢复 TE/L1 信号。** 这些 BAM 已经删除了大量低 MAPQ、duplicate、blacklist 或其他 reads。TE/L1 relaxed track 必须尽量从 pre-clean sorted/raw BAM 开始。

2. **完整 L1 body heatmap 建议：**

```text
--mapq 0
--exclude-flags 780
--run-markdup false
--normalization CPM
--binsize 10
--skip-zeros false
--sort-regions keep
--body-length 5500
```

3. **TE/L1 relaxed bigWig 是诊断 track，不等价于普通 ATAC peak track。** 它适合判断 L1 body 中央低谷是否由 mappability/filtering 造成；真正做 peak calling 和常规 ATAC QC 仍然应使用 strict clean BAM。

4. **如果要比较 target/IgG 或不同处理组，优先使用相同参数生成的 bigWig。** 不要把 strict RPGC bigWig 和 relaxed CPM bigWig 直接混合解释。


---

## 新增：从 FASTQ 直接生成 sorted BAM + TE/L1 relaxed BAM/bigWig：`run_atac_te_tracks_from_fastq.sh`

### 作用

这个脚本用于**不跑 Nextflow 主流程**时，从 ATAC-seq FASTQ 直接生成：

```text
02_align_sorted/*.sorted.bam        # pre-clean sorted BAM，可作为后续 TE source BAM
02_align_sorted/*.sorted.bam.bai
02_align_te/*.te.bam                # TE/L1 relaxed BAM
02_align_te/*.te.bam.bai
02_align_te/*.te.clean_counts.tsv   # 每一步过滤后 reads 数
02_align_te/*.te.flagstat.txt
04_bw_te/*.te.bw                    # TE/L1 relaxed bigWig
logs/*.bowtie2.log
logs/*.fastp.log
```

它适合用来研究 L1/TE heatmap 中央低谷是否来自 MAPQ、duplicate、blacklist、mitochondria 等过滤，而不是用于替代普通 ATAC peak calling 主流程。

### 单样本 PE 用法

```bash
bash scripts/run_atac_te_tracks_from_fastq.sh \
  --sample SAMPLE1 \
  --layout PE \
  --fastq1 SAMPLE1_R1.fastq.gz \
  --fastq2 SAMPLE1_R2.fastq.gz \
  --outdir results_te_fastq \
  --species hg38 \
  --cores 16 \
  --mapq 0 \
  --exclude-flags 780 \
  --run-markdup false \
  --normalization CPM \
  --binsize 10
```

### 单样本 SE 用法

```bash
bash scripts/run_atac_te_tracks_from_fastq.sh \
  --sample SAMPLE1 \
  --layout SE \
  --fastq1 SAMPLE1.fastq.gz \
  --outdir results_te_fastq \
  --species hg38 \
  --cores 16
```

### 批量 samplesheet 用法

samplesheet 使用主流程相同的核心列：

```csv
sample,layout,r1,r2,condition,replicate
S1,PE,/path/S1_R1.fastq.gz,/path/S1_R2.fastq.gz,WT,rep1
S2,PE,/path/S2_R1.fastq.gz,/path/S2_R2.fastq.gz,KO,rep1
```

运行：

```bash
bash scripts/run_atac_te_tracks_from_fastq.sh \
  --samplesheet samplesheet.csv \
  --outdir results_te_fastq \
  --species hg38 \
  --cores 16 \
  --normalization CPM \
  --binsize 10
```

### 参数

| 参数 | 必需 | 默认值 | 说明 |
|---|---:|---|---|
| `--samplesheet` | 二选一 | 空 | 批量 CSV；列至少包含 `sample,layout,r1`，PE 还需要 `r2` |
| `--sample` | 二选一 | 空 | 单样本名称 |
| `--layout` | 单样本必需 | 空 | `PE` 或 `SE` |
| `--fastq1` | 单样本必需 | 空 | R1 或 SE FASTQ |
| `--fastq2` | PE 必需 | 空 | R2 FASTQ |
| `--outdir` | 否 | `./atac_te_tracks_from_fastq` | 输出目录 |
| `--species` | 否 | `hg38` | `hg38`、`mm10`、`mm39` |
| `--bowtie2-index` | 否 | species config | Bowtie2 index prefix |
| `--blacklist` | 否 | species config | blacklist BED |
| `--mito-chr` | 否 | species config 或 `chrM` | 线粒体染色体名 |
| `--effective-genome-size` | RPGC 必需 | species config | deepTools RPGC 用 |
| `--cores` | 否 | `8` | 线程数 |
| `--run-fastp` | 否 | `true` | 是否先跑 fastp |
| `--fastp-timeout` | 否 | `120m` | 单个样本 fastp 最长运行时间；设为 `0` 可关闭 timeout |
| `--max-insert` | 否 | `2000` | PE bowtie2 `-X` |
| `--mapq` | 否 | `0` | TE/L1 relaxed BAM 的 MAPQ cutoff |
| `--exclude-flags` | 否 | `780` | 默认不去 duplicate-marked reads |
| `--run-markdup` | 否 | `false` | 是否在 TE 分支去重复 |
| `--remove-mito` | 否 | `true` | 是否去线粒体 reads |
| `--remove-blacklist` | 否 | `true` | 是否去 blacklist |
| `--proper-pair-only` | 否 | `true` | PE 是否仅保留 proper pair |
| `--normalization` | 否 | `CPM` | `CPM`、`RPGC`、`RPKM`、`BPM` |
| `--binsize` | 否 | `10` | bigWig bin size |
| `--ignore-for-normalization` | 否 | 空 | deepTools RPGC 附加参数 |
| `--keep-intermediate` | 否 | `false` | 是否保留 raw/mapq/nomito 临时 BAM |

### 结果解释

- `02_align_sorted/*.sorted.bam`：仅完成 bowtie2 alignment + sort/index，还没有 MAPQ、duplicate、mito、blacklist 过滤；这是最适合后续重新生成 TE/L1 relaxed track 的 source BAM。
- `02_align_te/*.te.bam`：经过 TE/L1 relaxed 参数过滤后的 BAM，默认 `MAPQ=0`、保留 duplicate-marked reads、不去 `XS:i`。
- `04_bw_te/*.te.bw`：由 `*.te.bam` 生成的 bigWig，推荐用于 L1/TE heatmap 诊断。
- `*.te.clean_counts.tsv`：检查每一步保留 reads 数；如果 `post_mapq_flag_filter` 到 `final_te_bam` 掉得很厉害，需要看是哪一步过滤造成。
- `logs/*.fastp.log`：fastp 标准输出和错误输出；如果 fastp 挂起或超时，先看这个日志。

### fastp 长时间无输出时

旧版本直接前台运行 fastp，没有独立日志和 timeout；如果 fastp 卡在读取 FASTQ 或文件系统 I/O，终端会长时间停在某个样本。当前版本默认给每个样本 fastp 加 `120m` timeout，并把日志写到 `logs/<sample>.fastp.log`。如果输入 FASTQ 已经是主流程清洗后的文件，或只是想快速做 TE/L1 诊断，可以显式跳过 fastp：

```bash
bash scripts/run_atac_te_tracks_from_fastq.sh \
  --samplesheet samplesheet.csv \
  --outdir results_te_fastq \
  --species hg38 \
  --run-fastp false
```

---

## 主流程 sorted BAM 输出说明

现在主流程的 `ALIGN_FILTER_PE` 和 `ALIGN_FILTER_SE` 会同时声明输出：

```text
02_align/*.clean.bam
02_align/*.clean.bam.bai
02_align/*.qc.tsv
02_align/*.sorted.bam
02_align/*.sorted.bam.bai
02_align/*.bowtie2.log
```

其中：

- `*.clean.bam`：用于普通 ATAC peak calling、FRiP、TSS enrichment、标准 bigWig。
- `*.sorted.bam`：pre-clean sorted BAM，适合后续作为 TE/L1 relaxed track 的 source BAM；不要把它当作最终 clean BAM 做 peak calling。

如果你想从主流程结果补做 TE/L1 relaxed tracks，推荐使用：

```bash
bash scripts/run_atac_te_tracks_from_bam.sh \
  --bam-glob "results/02_align/*.sorted.bam" \
  --outdir results_te_tracks \
  --species hg38 \
  --mapq 0 \
  --exclude-flags 780 \
  --run-markdup false \
  --normalization CPM \
  --binsize 10
```
