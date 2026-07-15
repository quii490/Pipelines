# 分析记录与交付习惯

## 每次运行建立一个不可混淆的目录

```text
project/
├── inputs/          # manifests、metadata、contrasts；不复制大 FASTQ
├── configs/         # 本次冻结配置
├── scripts/         # 真正执行过的脚本
├── results/         # pipeline 输出
├── logs/            # scheduler/pipeline 日志
└── README.md        # 目的、版本、命令、结论边界
```

原始数据只读；新分析输出到新目录。不要在同一 `results/` 中手工替换部分文件而不留下记录。

## 运行前五项检查

- sample IDs 唯一，并与 FASTQ/BAM、metadata 一致；
- assembly、annotation 和 chromosome style 一致；
- CASE/CONTROL 方向写清楚；
- biological replicate 与 batch/paired design 正确；
- 预计磁盘、内存、线程和运行时间足够。

## 保存最小 provenance

```bash
git -C /path/to/Pipelines rev-parse HEAD > pipeline.commit.txt
conda env export --from-history > environment.from-history.yml
cp /path/to/manifest.csv inputs/
cp /path/to/run.config configs/
```

完整 conda export 可用于取证，但跨平台复现通常还需 from-history/lock file。reference FASTA、GTF、TE annotation 应记录版本、来源和 checksum。

## 中间结果能不能删

先分类：

- 必须保留：原始输入引用、manifest/config、核心 count/peak/DE 表、报告、版本与日志；
- 可重建但昂贵：sorted/clean BAM、bigWig、Nextflow work；项目完成前通常保留；
- 可重建且便宜：临时矩阵、重复 PNG、下载缓存；验证后可按政策清理。

删除前确认 `-resume`、补报告和复审是否仍需要 work/BAM。清理脚本先 dry-run，并由项目负责人确认。

## 交付前检查

1. 核心表格不是空文件，列名和 contrast 方向可解释；
2. 图能追溯到对应完整表格和 plotting 参数；
3. QC 的分母、过滤阶段和阈值写清；
4. 没有内部绝对路径、secret、患者身份或未授权数据；
5. README 指明优先查看哪些文件、哪些模块 SKIP/FAIL；
6. 在一个干净路径中验证相对链接和解压/读取方式。

“pipeline exit 0”只是开始。可靠交付需要同时满足：文件完整、QC 可解释、统计设计正确、版本可追溯和公开范围合规。
