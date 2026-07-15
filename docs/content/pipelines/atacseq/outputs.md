# ATAC-seq 输出与结果解释

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

```text
results/
├── 00_refcheck/
├── 01_fastp/
├── 02_align/                    # sorted/clean BAM、BAI、alignment stats
├── 02_align_te/                 # 可选 relaxed TE BAM
├── 03_qc/                       # QC、TSS/gene body/TE profiles
├── 04_bw/                       # standard bigWig
├── 04_bw_te/                    # 可选 relaxed TE bigWig
├── 05_peaks/                    # sample peaks/summits
├── 06_consensus_peaks/          # consensus BED
├── 07_counts/
│   ├── consensus_peak_counts.txt
│   ├── sample_metadata.csv
│   └── bin_level/
├── 08_downstream/
│   ├── peak_level/
│   └── bin_level/
├── 09_motif/
├── 10_footprinting/
├── 11_nuc_phasing/
├── QC_REPORT.md
└── _automation/
```

具体目录随 preset/开关变化；不存在的可选目录先看运行计划和日志，不自动视为失败。

## 优先查看顺序

1. `QC_REPORT.md` 与 `03_qc/`：样本是否完整、指标定义是否适用；
2. `02_align/`：mapping、proper pairs、duplicate/mitochondrial/blacklist 过滤；
3. `04_bw/`：浏览器中背景、peak shape 与重复一致性；
4. `05_peaks/`、`06_consensus_peaks/`：peak 数量、宽度和 universe；
5. `07_counts/`：样本列、PE fragment/SE read 单位、matrix 是否非空；
6. `08_downstream/`：PCA/correlation 后再看 differential regions；
7. motif/footprinting/TE：作为下游证据，不覆盖基础 QC。

## 关键文件

### BAM

`*.clean.bam` 是常规分析输入；在使用独立脚本前确认 clean 定义。`*.sorted.bam` 可包含后来被 strict filtering 删除的 reads，适合 TE relaxed 分支，但背景也更高。

### bigWig

用于 IGV、correlation 和 profiles。数值依赖 normalization、effective genome size、bin size 与过滤；跨 run 比较前必须核对。

### Peaks 与 consensus

sample peaks 来自各样本；consensus peak universe 用于多样本 counting。改变 peak caller/cutoff/summit half-width 会改变矩阵 feature 定义，因此旧/new counts 不可直接拼接。

### Count matrix

PE 通常按 fragment，SE 按 read。raw region counts进入差异模型；不要把 bigWig signal 当 raw counts。

### Differential tables

至少保留 coordinates、abundance/baseMean、log2FC、p-value、adjusted p-value、contrast 和 region universe。`log2FC > 0` 的方向由 CASE,CONTROL 决定。

## Peak-level 与 bin-level

- peak-level：聚焦富集区域，统计功效通常较高，但依赖 peak universe；
- bin-level：全基因组固定窗口，不依赖 peak calling，但 feature 更多、稀疏和多重检验更重。

二者是互补分析，不是选择“显著更多”的一种。归档表格、metadata、contrast、命令和软件版本，不只保存 PDF。
