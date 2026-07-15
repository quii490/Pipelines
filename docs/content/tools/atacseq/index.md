# ATAC-seq 常用工具

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | ATAC-seq maintainers | 2026-07-16 |

完整项目优先跑主 pipeline；独立脚本用于已有可靠中间文件、补跑某一模块或调试。先按目标选择，不需要从头阅读所有参数。

| 你想做什么 | 页面 | 推荐脚本 |
|---|---|---|
| 从已有 results/BAM 重建结果 | [重建与补跑](rebuild.md) | `run_atac_downstream_only.sh`、`run_atac_from_bam.sh`、`run_callpeak_from_bam.sh`、`run_fixedbin_from_bam.sh` |
| TSS/gene body/相关性/nucleosome QC | [信号与 QC](signal-qc.md) | `run_tss_qc.sh`、`run_gene_body_profile.sh`、`run_bw_correlation.sh`、`run_nuc_phasing.sh` |
| 任意 region heatmap、peak annotation/motif | [区域、注释与 motif](regions-motif.md) | `run_region_heatmap.sh`、`run_diff_peak_heatmap.sh`、`run_peak_annotation.sh`、`run_peak_overlap.sh`、`run_motif_homer.sh`、`run_tobias.sh` |
| TE/L1 relaxed track 与 heatmap | [TE 工具](te.md) | `run_atac_te_tracks_from_bam.sh`、`run_atac_te_tracks_from_fastq.sh`、`run_te_heatmap*.sh` |

## 通用规则

1. 从仓库根目录用 `bash ATAC-seq/scripts/<script>.sh --help`；
2. 所有 BAM/bigWig/BED 必须同一 assembly；
3. glob 加引号，例如 `--bw-glob "/path/*.bw"`；
4. 输出到新目录，保存完整命令和 `--help` 对应 commit；
5. PE/SE、read/fragment、normalization、strict/relaxed BAM 必须在图注中写清。

可选安装命令链接：

```bash
bash ATAC-seq/scripts/install_bin_tools.sh /path/to/bin
```

没有被安装的脚本仍可通过完整路径运行。
