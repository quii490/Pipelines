# ATAC-seq 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | ATAC-seq maintainers | 2026-07-15 | `main` |

## samplesheet 样本或 R2 错误

使用 `--init-only` 重新生成到临时结果目录，比对 FASTQ basename；确认后再用 `--overwrite-inputs` 更新正式输入。

## PE 与 SE 混合报错

按 layout 拆分成独立运行，不能绕过检查后合并 count matrix。

## peaks/FRiP 异常

先检查 clean BAM、深度、blacklist、duplicate 和 MACS 参数。不要只通过降低 peak cutoff 制造更多 peaks。

## mm39 注释缺失

提供同 build 的 GTF。流程会跳过不兼容注释，而不应静默使用 mm10 坐标。
