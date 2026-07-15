# ATAC-seq 常用工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

完整项目优先跑主 pipeline；独立脚本用于已有可靠中间文件、补跑某一模块或调试。先按目标选择，不需要从头阅读所有参数。每个命令运行前都建议查看当前版本帮助：

```bash
bash ATAC-seq/scripts/<script>.sh --help
# 或安装统一命令后
atacseq help
```

| 你想做什么 | 页面 | 推荐脚本 |
|---|---|---|
| 从已有 results/BAM 重建结果 | [重建与补跑](rebuild.md) | `run_atac_downstream_only.sh`、`run_atac_from_bam.sh`、`run_callpeak_from_bam.sh`、`run_fixedbin_from_bam.sh` |
| TSS/gene body/相关性/nucleosome QC | [信号与 QC](signal-qc.md) | `run_tss_qc.sh`、`run_gene_body_profile.sh`、`run_bw_correlation.sh`、`run_nuc_phasing.sh` |
| 任意 region heatmap、peak annotation/motif | [区域、注释与 motif](regions-motif.md) | `run_region_heatmap.sh`、`run_diff_peak_heatmap.sh`、`run_peak_annotation.sh`、`run_peak_overlap.sh`、`run_motif_homer.sh`、`run_tobias.sh` |
| TE/L1 relaxed track 与 heatmap | [TE 工具](te.md) | `run_atac_te_tracks_from_bam.sh`、`run_atac_te_tracks_from_fastq.sh`、`run_te_heatmap*.sh` |
| 重建 QC 报告、profile、跨 contrast 图 | [报告、汇总与跨对比](reporting-comparison.md) | `generate_qc_report.py`、`run_atac_profile_heatmaps.sh`、`run_cross_contrast_heatmap.sh` |

## 完整脚本清单

“用户入口”可以直接运行；“内部辅助”通常由主流程或 wrapper 调用，不建议手工执行。

| 脚本 | 类型 | 主要输入 | 作用 |
|---|---|---|---|
| `atacseq` | 用户入口 | 取决于 subcommand | 统一分发 `run/init/from-bam/downstream/callpeak/...` |
| `run_atac_downstream_only.sh` | 用户入口 | 已有 counts、metadata、peaks | 只重跑 edgeR、annotation、TE/gene 下游 |
| `run_atac_from_bam.sh` | 用户入口 | clean BAM、metadata、contrasts | 从 BAM 重建 peaks、counts、tracks 和下游 |
| `run_callpeak_from_bam.sh` | 用户入口 | 一个或多个 BAM | 只重新调用 MACS3 peaks |
| `run_fixedbin_from_bam.sh` | 用户入口 | BAM、metadata、contrasts | 固定窗口计数和差异分析 |
| `run_tss_qc.sh` | 用户入口 | bigWig、TSS BED | TSS aggregate profile/heatmap |
| `run_gene_body_profile.sh` | 用户入口 | bigWig、gene-body BED | scale-regions gene-body profile |
| `run_bw_correlation.sh` | 用户入口 | 同策略 bigWig | correlation、scatter、PCA |
| `run_nuc_phasing.sh` | 用户入口 | PE BAM | fragment distribution 与 NRL |
| `run_region_heatmap.sh` | 用户入口 | BED、bigWig | 任意区域 heatmap/profile |
| `run_diff_peak_heatmap.sh` | 用户入口 | contrast BED 目录、bigWig | 批量差异 peak heatmap |
| `run_peak_annotation.sh` | 用户入口 | peak BED、GTF、TE annotation | ChIPseeker gene/TE annotation |
| `run_peak_overlap.sh` | 用户入口 | 多个 peak BED | 多集合 overlap 与图表 |
| `run_motif_homer.sh` | 用户入口 | contrast BED 目录 | HOMER motif/annotation |
| `run_tobias.sh` | 用户入口 | BAM、peaks、genome、motifs | TOBIAS correction/footprint/BINDetect |
| `run_cross_contrast_heatmap.sh` | 用户入口 | 多个 differential tables | 跨 contrast logFC heatmap |
| `run_atac_profile_heatmaps.sh` | 用户入口 | results/bigWig、feature BED | 一次生成 TSS/gene/TE/L1 profiles |
| `run_atac_te_tracks_from_bam.sh` | 用户入口 | pre-clean BAM | relaxed TE BAM/bigWig |
| `run_atac_te_tracks_from_fastq.sh` | 用户入口 | FASTQ/samplesheet | 从 FASTQ 生成 TE diagnostic tracks |
| `run_te_heatmap.sh` | 用户入口 | TE annotation、bigWig | 单次 TE heatmap |
| `run_te_heatmap_batch.sh` | 用户入口 | global/contrast BED、bigWig | 批量调用 TE heatmap |
| `install_bin_tools.sh` | 安装辅助 | 目标 `bin/` | 建立命令符号链接 |
| `generate_qc_report.py` | 报告辅助 | 已有 `03_qc/` | 重建中文 `QC_REPORT.md` |
| `prepare_manifest_fastq.py` | 内部辅助 | pipeline samplesheet | 生成 downstream metadata |
| `plot_atac_qc.R` | 内部辅助 | QC 汇总表 | 生成 QC 图和 pass/fail 表 |
| `organize_atac_downstream_outputs.py` | 内部辅助 | raw R output | 整理面向用户的目录和索引 |
| `run_peak_annotation_chipseeker.R` | 内部后端 | peaks、GTF、TE | `run_peak_annotation.sh` 的 R 后端 |
| `species_config_lib.sh` | 内部函数库 | species config | 为 shell wrapper 读取参考路径 |

## 通用规则

1. 从仓库根目录用 `bash ATAC-seq/scripts/<script>.sh --help`；统一入口见[命令参考](command-reference.md)；
2. 所有 BAM/bigWig/BED 必须同一 assembly；
3. glob 加引号，例如 `--bw-glob "/path/*.bw"`；
4. 输出到新目录，保存完整命令和 `--help` 对应 commit；
5. PE/SE、read/fragment、normalization、strict/relaxed BAM 必须在图注中写清。

可选安装命令链接：

```bash
bash ATAC-seq/scripts/install_bin_tools.sh /path/to/bin
```

没有被安装的脚本仍可通过完整路径运行。
