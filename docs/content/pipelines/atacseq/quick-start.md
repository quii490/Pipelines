# ATAC-seq Quick Start
输入fastq文件，完成ATAC-seq分析和可视化。

**建议路径都选择绝对路径，VScode可以右键复制绝对路径**
## 1. 激活环境/检查入口
```bash
conda activate /home/machicheng/.conda/envs/chipseq
#可选 bash /home/machicheng/ATAC-seq/pipeline_final/run_auto_atacseq.sh --help
```
## 2. 准备 FASTQ
推荐目录如下：

```text
project_name/
├── fastq/                 # 原始数据，只读
	├── WT_1_R1.fastq.gz
	├── WT_1_R2.fastq.gz
	├── KO_1_R1.fastq.gz
	└── KO_1_R2.fastq.gz
└── results/              

```
results可以不用手动建文件夹，运行时会自动生成。
样本名、condition 和 replicate 的详细规则见[输入准备](input.md)。
## 3. 初始化并检查输入

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode init \
  --fastq-dir /path/to/fastq \
  --species hg38 \
  --outdir /path/to/results
```

检查自动生成的 samplesheet；分组信息不明确时，显式提供 `--metadata-csv`。

## 4. 正式运行

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode auto \
  --fastq-dir /path/to/fastq \
  --species hg38 \
  --outdir /path/to/results \
  --preset standard \
  --resume
```

若需比较 KO 与 WT：

```bash
bash ATAC-seq/run_auto_atacseq.sh \
  --mode downstream \
  --outdir /path/to/results \
  --metadata-csv /path/to/metadata.csv \
  --contrast KO,WT \
  --levels both
```

## 5. 判断是否完成

- 命令退出码为 0；
- `QC_REPORT.md` 存在且样本未出现未解释的严重异常；
- `02_align/`、`04_bw/`、`05_peaks/` 有对应样本结果；
- downstream 启用时，差异分析结果非空且分组正确。

!!! tip "最先检查"

    先读 `QC_REPORT.md`，再看 mapping、mitochondrial fraction、peak 数量、FRiP、样本相关性与 PCA，最后解释差异 peaks。
