# RNA-seq Quick Start

| 状态 | 维护人 | 最后审查 | 适用版本 |
|---|---|---|---|
| Draft | RNA-seq maintainers | 2026-07-16 | `main` |

本页只保留第一次跑通所需步骤。开始前确认 FASTQ、物种、library layout 和 strandedness；不确定 strand 时先向建库/测序方确认。

## 1. 检查入口

```bash
cd /path/to/Pipelines/RNA-seq/rnaseq
python3 --version
nextflow -version
bash run_auto_rnaseq.sh --help
```

## 2. 准备 FASTQ

```text
project/
└── fastq/
    ├── WT_1_R1.fastq.gz
    ├── WT_1_R2.fastq.gz
    ├── WT_2_R1.fastq.gz
    ├── WT_2_R2.fastq.gz
    ├── KO_1_R1.fastq.gz
    └── KO_1_R2.fastq.gz
```

样本名不要包含空格。PE 数据必须同时存在 R1/R2，且 mate basename 一致。

## 3. 先 dry-run

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/project/fastq \
  --results-dir /path/to/project/results \
  --species hg38 \
  --strand reverse \
  --layout auto \
  --dry-run
```

dry-run 只验证输入、参数和计划模块，不提交分析任务。若 reference、环境或 FASTQ 检查失败，先解决后再正式运行。

## 4. 正式运行上游

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/project/fastq \
  --results-dir /path/to/project/results \
  --species hg38 \
  --strand reverse \
  --background \
  --resume \
  --max-cpus 16 \
  --max-memory "64 GB"
```

资源值是示例，不是固定要求。运行日志位于：

```text
results/_automation/logs/
```

## 5. 修改实验设计

上游创建：

```text
results/condition.csv
results/contrast.csv
```

确认 `condition.csv`：

```csv
sample,condition,replicate
WT_1,WT,1
WT_2,WT,2
KO_1,KO,1
KO_2,KO,2
```

确认 `contrast.csv`：

```csv
case,control
KO,WT
```

`KO_vs_WT` 表示 KO 相对 WT。不要在没有备份时使用 `--replace-design` 覆盖已经人工修改的设计。

## 6. 运行下游

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/project/results \
  --species hg38 \
  --downstream
```

## 7. 判断是否成功

- [ ] 自动化状态中没有核心模块失败。
- [ ] 所有预期样本均生成 gene alignment/counts。
- [ ] MultiQC 或对应 QC 文件存在且可打开。
- [ ] `condition.csv` 样本与 count matrix 列一致。
- [ ] PCA/correlation 中没有未解释的标签错误或极端离群。
- [ ] 差异结果方向与 `case,control` 定义一致。

完成后按[输出指南](outputs.md)检查关键文件，再按[QC](qc.md)解释样本质量。
