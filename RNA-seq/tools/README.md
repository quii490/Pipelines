# RNA-seq / TE pipeline 使用说明（面向不懂生信的同学）

这份说明尽量用“怎么操作”的方式来写。

你可以把这条流程理解成两段：

1. **上游**：从 `FASTQ` 开始，做质控、比对、计数，得到 `BAM` 和各种 `raw count matrix`。  
2. **下游**：从 `raw count matrix` 开始，做差异分析和出图，例如 PCA、火山图、热图、GSEA。

> **一句话原则**：  
> - `FASTQ` 用来跑完整 pipeline。  
> - `BAM` 主要是中间结果，不是火山图/GSEA 的直接输入。  
> - **火山图、热图、GSEA 的真正输入是 raw count matrix。**

---

## 1. 这条 pipeline 能做什么

### 上游能做的
- `FASTQ -> cleaned FASTQ -> BAM`
- gene 比对：`STAR` 或 `HISAT2`
- gene 定量：
  - `featureCounts`
  - `Salmon`
  - `StringTie`（默认关闭）
- TE 定量：
  - `TElocal`
  - `TEtranscripts`
  - `TEcount`（默认关闭）
  - `REdiscoverTE`
  - `SalmonTE`
- 自动生成：
  - `condition.csv`
  - `contrast.csv`

### 下游能做的
- PCA
- 样本相关性热图
- 火山图
- MA 图
- 热图
- GO 富集（gene）
- GSEA（gene）
  - Hallmark
  - KEGG
  - Reactome
  - GO_BP
  - GO_CC
  - GO_MF

### 独立工具能做的
这次我另外提供了两个独立脚本：

- `run_plot_from_counts.R`  
  从 **raw count matrix** 直接画图，不需要重跑完整 pipeline。
- `run_gsea_standalone.R`  
  从 **gene raw count matrix** 直接跑 GSEA，不需要重跑完整 pipeline。

另外还提供了安装脚本：

- `install_rnaseq_cli_tools.sh`

安装后可以直接在 conda 环境里调用：

- `rnaseq-bam-to-counts`
- `rnaseq-diff-counts`
- `rnaseq-annotate-de`
- `rnaseq-pathway-de`
- `rnaseq-te-analysis`
- `rnaseq-counts-to-de`
- `rnaseq-de-visuals`
- `rnaseq-plot-counts`
- `rnaseq-gsea`
- `rnaseq-bw-cor`
- `rnaseq-two-sample-scatter`
- `rnaseq-qc-metrics-bam`

当前版本的 R 工具会自动寻找 `../rnaseq-downstream/rnaseq-function.R`，所以从 `tools/` 目录直接运行时，通常不需要额外传 `--function-file`。如果通过 `install_rnaseq_cli_tools.sh` 安装，安装脚本也会把 `rnaseq-function.R` 放进目标 `bin/`，方便复制模式继续工作。

如果命令装在 `rnaseq` 环境，但 R 包在 `downstream` 环境，建议安装成 wrapper：

```bash
cd /path/to/RNA-seq/Pipeline/rnaseq/tools
bash install_rnaseq_cli_tools.sh \
  --prefix /path/to/.conda/envs/rnaseq \
  --rscript-cmd "conda run -n downstream Rscript"
```

这样 `rnaseq-pathway-de`、`rnaseq-counts-to-de`、`rnaseq-de-visuals` 等命令可以在 `(rnaseq)` shell 中直接调用，同时实际 R 运行环境是 `downstream`。

## raw counts 直接出图和 pathway

可以。推荐用 `rnaseq-counts-to-de`，它会从 raw count matrix 生成每个 contrast 的 DE matrix，并同时输出火山图、MA 图、top heatmap、normalized counts、VST/log matrix。gene counts 加上 `--run-go true --run-gsea true` 后也会输出 GO/GSEA pathway 结果：

```bash
rnaseq-counts-to-de \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --species hg38 \
  --outdir /path/to/counts_to_visuals \
  --matrix-name gene_featureCounts \
  --make-plots true \
  --run-go true \
  --run-gsea true
```

如果只想先得到 DE matrix，不跑图：

```bash
rnaseq-diff-counts \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --outdir /path/to/diff_only \
  --matrix-name gene_featureCounts
```

然后再对某个 DE matrix 单独补 pathway：

```bash
rnaseq-pathway-de \
  --de-matrix /path/to/gene_featureCounts_KO_vs_WT.DE_matrix.csv \
  --species hg38 \
  --tx2gene-path /path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf \
  --case KO \
  --control WT \
  --outdir /path/to/pathway/KO_vs_WT
```

更完整的小工具说明见同目录：

```text
TOOLS.md
```

如果要把多个下游模块组织成可复跑流程，使用：

```text
../snakemake-visualization/visualization.smk
../snakemake-visualization/config.example.yaml
../snakemake-visualization/README.md
```

---

## 2. 先理解 4 种常见文件

### 2.1 FASTQ
原始测序数据。完整 pipeline 的标准起点。

### 2.2 BAM
比对后的结果文件。  
BAM 本身**不能直接拿来画火山图或 GSEA**，因为差异分析要用的是“每个基因/TE 的计数矩阵”，不是 read alignment 本身。

### 2.3 raw count matrix
这是下游最关键的输入。  
例如：

- `gene_featureCounts_matrix.csv`
- `gene_Salmon_matrix.csv`
- `TE_TElocal_repName_matrix.csv`
- `TE_TEtranscripts_repClass_matrix.csv`

### 2.4 sample table / contrast table
也就是：

- `condition.csv`
- `contrast.csv`

这两个文件告诉脚本：
- 每个样本属于哪个组
- 你要比较哪两组

---

## 3. 推荐的标准使用方式

## 3.1 第一步：从 FASTQ 跑到矩阵

主入口是：

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/fastq_dir \
  --species hg38
```

跑完后会得到：

- BAM
- gene / TE 的各种 count matrix
- `results-dir/condition.csv`
- `results-dir/contrast.csv`

### 支持的 FASTQ 目录形式

#### 形式 A：所有 FASTQ 平铺在一个目录
```text
project_fastq/
  A1KO_1.fq.gz
  A1KO_2.fq.gz
  A2KO_1.fq.gz
  A2KO_2.fq.gz
  NC1_1.fq.gz
  NC1_2.fq.gz
```

#### 形式 B：每个 sample 在自己的子目录
```text
project_fastq/
  A1KO/
    A1KO_1.fq.gz
    A1KO_2.fq.gz
  A2KO/
    A2KO_1.fq.gz
    A2KO_2.fq.gz
  NC1/
    NC1_1.fq.gz
    NC1_2.fq.gz
```

两种都可以直接把 `project_fastq/` 作为 `--fastq-dir`。

### 支持的单双端命名

#### 双端
- `_1/_2`
- `_R1/_R2`
- `_r1/_r2`
- `.1/.2`
- `-R1/-R2`

#### 单端
- `.fastq.gz`
- `.fq.gz`
- `.fastq`
- `.fq`

---

## 3.2 第二步：手动修改 `condition.csv` 和 `contrast.csv`

### `condition.csv` 示例
```csv
sample,condition,replicate
A1KO,KO,1
A2KO,KO,2
NC1,NC,1
NC2,NC,2
```

### `contrast.csv` 示例
```csv
case,control
KO,NC
```

这里的意思是：
- `condition` 列里有 `KO` 和 `NC`
- 我要比较 `KO vs NC`

> 注意：  
> `contrast.csv` 里的 `case/control` 填的是 **condition 名字**，不是 sample 名字。

---

## 3.3 第三步：只跑下游出图

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream
```

如果想指定绘图输出目录：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream \
  --plot-dir /path/to/rnaseq_results/plots_custom
```

---

# 4. 最常用参数说明

## 4.1 基本参数

### `--fastq-dir`
FASTQ 根目录。完整流程最常用输入。

### `--results-dir`
结果目录。  
如果不写，会自动在 FASTQ 同级目录下生成默认结果目录。

### `--species`
支持：
- `hg38`
- `mm10`
- `mm39`

### `--human-ref`
仅 `species=hg38` 时有意义。

支持：
- `hg38`
- `t2t`

例如：

```bash
--species hg38 --human-ref t2t
```

---

## 4.2 比对和定量参数

### `--aligner`
gene 比对方式：
- `star`
- `hisat2`

例如：

```bash
--aligner hisat2
```

> 注意：  
> 这里主要影响 gene 分支。TE 分支仍然按原流程走。

### `--strandedness`
支持：
- `unstranded`
- `forward`
- `reverse`

常见示例：

```bash
--strandedness reverse
```

### `--layout`
支持：
- `auto`
- `PE`
- `SE`

如果你不确定，就用默认：

```bash
--layout auto
```

### `--r1-pattern` / `--r2-pattern`
如果你的 FASTQ 命名不规则，可以手动指定。

例如：

```bash
--r1-pattern _r1 --r2-pattern _r2
```

---

## 4.3 上游模块开关

这些参数决定上游跑哪些工具。

### 默认比较关键的值
- `--run-fastp true`
- `--run-star-fc true`
- `--run-salmon true`
- `--run-stringtie false`
- `--run-tecount false`
- `--run-telocal true`
- `--run-tetranscripts true`
- `--run-dedup true`
- `--run-rediscoverte auto`
- `--run-salmonte true`

### 推荐的稳妥用法
我更建议你平时常用这一版：

```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/fastq_dir \
  --species hg38 \
  --run-dedup false \
  --run-salmonte false
```

原因：
- TE 主分析通常不建议默认 dedup
- SalmonTE 稳定性一般，不建议当主结果依赖

---

## 4.4 下游绘图参数

### `--padj-cutoff`
默认：`0.05`

### `--lfc-cutoff`
默认：`0.585`

### `--baseMean-min`
默认：`5`

### `--label-top-n`
默认：`40`

### `--heatmap-top-n`
默认：`40`

### `--only-tools`
只运行指定下游模块，逗号分隔。例：`TE_TEtranscripts`

### `--skip-tools`
跳过指定下游模块，逗号分隔。例：`TE_TElocal,TE_TEcount`

可用模块名：

```text
gene_featureCounts
gene_Salmon
TE_TElocal
TE_TEcount
TE_TEtranscripts
TE_REdiscoverTE
TE_SalmonTE
panel_plots
```

示例：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream \
  --padj-cutoff 0.05 \
  --lfc-cutoff 0.585 \
  --label-top-n 30 \
  --heatmap-top-n 50
```

只重跑 TEtranscripts：

```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream \
  --only-tools TE_TEtranscripts
```

---

## 4.5 资源参数

### `--max-cpus`
最大 CPU 数。

### `--max-memory`
最大内存，例如：

```bash
--max-memory "300 GB"
```

### `--queue-size`
Nextflow 并行队列大小。

---

# 5. 输出目录怎么看

完整跑完后，通常会看到类似结构：

```text
rnaseq_results/
  condition.csv
  contrast.csv
  _automation/
    logs/
    work/
    inputs/
  01_fastp/
  02_gene_star/ 或 02_gene_hisat2/
  03_gene_featurecounts/
  04_gene_stringtie/
  05_gene_salmon/
  06_te_star/
  07_tecount/
  08_telocal/
  09_TEtranscripts/
  10_rediscoverte/
  10_rediscoverte_rollup/
  11_salmonte/
  plots/
```

### 你最常会用到的原始矩阵
- `gene_featureCounts_matrix.csv`
- `gene_Salmon_matrix.csv`
- `TE_TElocal_repName_matrix.csv`
- `TE_TElocal_repFamily_matrix.csv`
- `TE_TElocal_repClass_matrix.csv`
- `TE_TEtranscripts_repName_matrix.csv`
- `TE_TEtranscripts_repFamily_matrix.csv`
- `TE_TEtranscripts_repClass_matrix.csv`

---

# 6. 只画火山图 / 热图 / PCA：到底输入什么？

## 6.1 结论
### **输入不是 BAM，而是 raw count matrix。**

原因：
- 火山图来自差异分析结果
- 差异分析需要 count matrix
- BAM 只是 read 比对结果，不是差异分析的直接输入

所以：

### 只画火山图，需要：
- raw count matrix
- `condition.csv`
- `contrast.csv`
- （可选）基因注释或 TE 注释

### 只跑 GSEA，需要：
- **gene raw count matrix**
- `condition.csv`
- `contrast.csv`
- 基因注释（GTF 或 tx2gene）

---

# 7. 独立工具：从 raw count 直接画图

- `run_plot_from_counts.R`

它适用于：
- 你已经有 raw count matrix
- 你不想重跑完整 pipeline
- 只想单独画火山图 / MA / 热图 / PCA / 相关性热图

## 7.1 安装成环境命令

```bash
bash install_rnaseq_cli_tools.sh --prefix /path/to/your/conda/env
```

如果你当前已经 `conda activate your_env`，也可以直接：

```bash
bash install_rnaseq_cli_tools.sh
```

安装后可直接调用：

```bash
rnaseq-plot-counts
rnaseq-gsea
```

---

## 7.2 只画 gene 火山图

```bash
rnaseq-plot-counts \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/only_plots \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --tx2gene-path /path/to/gencode.gtf \
  --padj-cutoff 0.05 \
  --lfc-cutoff 0.58 \
  --label-top-n 30 \
  --volcano-orientation classic \
  --gray-nonsig true \
  --plots volcano
```

### 说明
- 这里的 `matrix` 必须是 **raw count matrix**
- 不是 BAM
- `--plots volcano` 表示只画火山图
- `--volcano-orientation classic` 是默认经典竖版；要恢复旧横版可改成 `horizontal`
- `--gray-nonsig true` 表示非显著点统一灰色；如需按 TE/gene 分组保留颜色，可设为 `false`

---

## 7.3 只画 gene 火山图 + MA 图 + 热图

```bash
rnaseq-plot-counts \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/only_plots \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --tx2gene-path /path/to/gencode.gtf \
  --label-top-n 30 \
  --volcano-orientation classic \
  --gray-nonsig true \
  --plots volcano,ma,heatmap
```

---

## 7.4 画 TE 火山图

```bash
rnaseq-plot-counts \
  --matrix /path/to/TE_TElocal_repName_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/only_plots \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --te-annotation-tsv /path/to/te_annotation.tsv \
  --te-label-level repName \
  --te-color-level repFamily \
  --plots volcano
```

### 常用 TE 参数
- `--te-label-level repName`
- `--te-color-level repFamily`

如果是 family/class 矩阵，也可以改成：
- `--te-label-level repFamily`
- `--te-color-level repClass`

---

## 7.5 只画 PCA 和样本相关性热图

```bash
rnaseq-plot-counts \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/only_plots \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --tx2gene-path /path/to/gencode.gtf \
  --plots pca,corr
```

---

# 8. 独立工具：只跑 GSEA

这次同时提供：

- `run_gsea_standalone.R`
- 安装后的命令名：`rnaseq-gsea`

## 8.1 只跑 GSEA（gene）

```bash
rnaseq-gsea \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/gsea_only \
  --species hg38 \
  --tx2gene-path /path/to/gencode.gtf \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv
```

## 8.2 只跑 GSEA（RDS 矩阵）

```bash
rnaseq-gsea \
  --matrix /path/to/gene_Salmon_matrix.rds \
  --matrix-format rds \
  --outdir /path/to/gsea_only \
  --species hg38 \
  --tx2gene-path /path/to/gencode.gtf \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv
```

## 8.3 GSEA 可调参数

### `--threads`
contrast 并行数。

### `--topn-plot`
每类 GSEA 默认展示多少通路。

### `--max-terms-gseaplot2`
代表性峰谷图最多画多少个通路。

### `--min-gs-size`
最小 gene set 大小。

### `--max-gs-size`
最大 gene set 大小。

### `--pvalue-cutoff`
GSEA 输出 cutoff。

### `--disable-gseaplot2 true`
如果你不想画峰谷图，可以关掉。

示例：

```bash
rnaseq-gsea \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --matrix-format csv \
  --outdir /path/to/gsea_only \
  --species hg38 \
  --tx2gene-path /path/to/gencode.gtf \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --threads 8 \
  --topn-plot 20 \
  --max-terms-gseaplot2 10
```

---

# 9. 从 BAM 到出图：当前该怎么做

## 9.1 最推荐方式
**最推荐还是从 FASTQ 直接跑完整 pipeline。**

因为这样最省事：
- 自动对齐
- 自动计数
- 自动生成 matrix
- 自动生成 condition / contrast 模板
- 自动跑下游

## 9.2 如果你手里已经只有 BAM
当前这套工具里，**没有单独做“BAM -> 所有 matrix -> 所有图”的一个总入口脚本**。  
所以如果你只有 BAM，推荐这么做：

### gene
先用 `featureCounts` 生成 gene raw count matrix

### TE
先用：
- `TElocal`
- `TEtranscripts`
- `TEcount`

生成 TE raw count matrix

然后再用：
- `rnaseq-plot-counts`
- `rnaseq-gsea`

去出图。

### 为什么不能直接拿 BAM 画火山图/GSEA
因为：
- 火山图 = 差异分析结果可视化
- GSEA = 排序后的基因统计结果做富集
- 它们都依赖 **count matrix**，不直接依赖 BAM

---

# 10. 常见问题

## 10.1 `contrast.csv` 里填 sample 还是 condition？
填 **condition**。

例如：

```csv
case,control
KO,NC
```

不是：

```csv
case,control
A1KO,NC1
```

---

## 10.2 `condition.csv` 里的 `sample` 要写什么？
必须和矩阵列名一致。

也就是：
- 如果矩阵列名是 `A1KO_rep1`
- 那 `condition.csv` 里的 `sample` 也必须是 `A1KO_rep1`

---

## 10.3 只想看火山图，要不要 GTF？
### gene 火山图
建议提供 `--tx2gene-path`，这样标签会更友好。

### TE 火山图
建议提供 `--te-annotation-tsv`，这样才能按 repName / family / class 上色和标注。

如果都不提供，也能跑，但标签只会是原始 feature ID。

---

## 10.4 为什么有些比较没有 p 值 / padj？
如果某一组重复太少，代码会自动切到 exploratory 模式。  
这种情况下：
- 仍然会给你 log2FC
- 但 p 值 / padj 可能不会正常计算

---

## 10.5 现在图片输出是什么格式？
当前已改成：

- **只输出 PDF**
- 不再输出 PNG

所以你看到的图主要都应该是 `.pdf`。

---

# 11. 最推荐的三种实际工作流

## 工作流 A：完整标准分析（最推荐）

### 第一步：上游
```bash
bash run_auto_rnaseq.sh \
  --fastq-dir /path/to/fastq_dir \
  --species hg38 \
  --run-dedup false \
  --run-salmonte false
```

### 第二步：改表
改：
- `results-dir/condition.csv`
- `results-dir/contrast.csv`

### 第三步：下游
```bash
bash run_auto_rnaseq.sh \
  --results-dir /path/to/rnaseq_results \
  --species hg38 \
  --downstream
```

---

## 工作流 B：已经有 matrix，只想画火山图

```bash
rnaseq-plot-counts \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --outdir /path/to/plot_only \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv \
  --tx2gene-path /path/to/gencode.gtf \
  --plots volcano
```

---

## 工作流 C：已经有 gene matrix，只想跑 GSEA

```bash
rnaseq-gsea \
  --matrix /path/to/gene_featureCounts_matrix.csv \
  --outdir /path/to/gsea_only \
  --species hg38 \
  --tx2gene-path /path/to/gencode.gtf \
  --sample-table /path/to/condition.csv \
  --contrast-file /path/to/contrast.csv
```

---

# 12. 这次新增的文件

这次补充的文件有：

- `README.md`：这份使用说明
- `run_plot_from_counts.R`：从 raw count matrix 直接画图
- `run_gsea_standalone.R`：从 gene raw count matrix 直接跑 GSEA
- `install_rnaseq_cli_tools.sh`：安装成 conda 环境命令

---

# 13. 最后一句最重要

如果你不确定自己现在手里的文件能干什么，就按下面判断：

- **只有 FASTQ** → 跑完整 pipeline
- **只有 BAM** → 先做计数，变成 raw count matrix
- **已经有 raw count matrix** → 可以直接用独立脚本画图 / 跑 GSEA
