# ATAC-seq 工具

| 工具 | 输入 | 用途 |
|---|---|---|
| `run_atac_from_bam.sh` | BAM | 从已有 BAM 补主分析 |
| `run_callpeak_from_bam.sh` | BAM | 重跑 peak calling |
| `run_atac_downstream_only.sh` | 主流程结果 | 重跑 peak/bin downstream |
| `run_region_heatmap.sh` | bigWig + BED | 自定义区域 heatmap |
| `run_tobias.sh` | BAM/bigWig + motifs | Footprinting |
| `run_nuc_phasing.sh` | PE BAM | Nucleosome phasing |

工具位于 `ATAC-seq/scripts/`。先确认输入与主流程 reference build、normalization 和 layout 一致。
