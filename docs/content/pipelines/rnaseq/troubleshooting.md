# RNA-seq 故障排查

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Active | RNA-seq maintainers | 2026-07-15 | `main` |

## STAR index 无法打开

确认 `--species`、reference build、index 前缀和文件权限一致；不要混用 hg38、T2T 和 mouse index。

## counts 近乎为空

检查 BAM 是否有 reads、GTF chromosome 命名、feature type 和 `--strand`。先解决 gene counts 再解释 TE 或差异结果。

## 下游提示样本不完整

比较 `condition.csv` 中的样本与对应 count matrix 列名。默认跳过不完整模块；不要为得到图片而直接允许部分输入。

## Nextflow 失败

从结果目录自动化日志和具体模块日志定位首个失败任务，修复输入/资源后在同一结果目录使用 `--resume`。
