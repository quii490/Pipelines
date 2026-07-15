# RNA-seq Core Downstream

这是 RNA-seq pipeline 的核心下游分析模块。主流程入口  会自动调用本目录的 。

## 文件说明

| 文件 | 用途 |
|---|---|
|  | 下游主入口。从 Nextflow 结果目录扫描各类 count 文件，构建矩阵，运行 DESeq2 差异分析并生成全部可视化。 |
|  | 下游函数库（~2900 行）。所有差异分析、可视化、富集分析函数均在此定义，被 、独立 CLI 工具和 Snakemake 规则共用。 |
|  | Conda 环境定义。包含 DESeq2、clusterProfiler、ggplot2 等 80+ R 包。 |
|  | 辅助脚本：从已有结果目录自动生成 sample metadata CSV。 |
|  | 辅助脚本：将 TE GTF 注释文件转换为下游所需的 TSV 格式。 |

## 支持的输入模块

下游自动扫描 Nextflow 结果目录的以下子目录：

| 模块 | Nextflow 目录 | 文件模式 | 描述 |
|---|---|---|---|
|  |  |  | Gene 水平 counts |
|  |  |  | Salmon 转录本定量 |
|  |  |  | TE locus 水平定量 |
|  |  |  | TE class/family/subfamily |
|  |  |  | TE 转录本定量 |
|  |  |  | REdiscoverTE rollup |
|  |  |  | SalmonTE 表达矩阵 |
|  |  |  | Telescope TE 定量 |
|  | — | 聚合已有 PDF | 跨模块 panel 拼图 |

## 输出结构

运行后会在输出目录生成以下结构：



## 关键参数

详细参数见 。最常用：

-  (默认 0.05): 显著性阈值
-  (默认 0.58): log2 fold change 阈值 (~1.5x)
-  (默认 40): 火山图标注基因数
-  (默认 40): 热图展示基因数
- : 只运行指定模块 (如 )
- : 跳过指定模块
- : 火山图方向
- : 非显著点是否统一灰色
- : contrast 并行数

## 可视化产出

每个 contrast 会生成：
1. **PCA 图** (): 全局样本聚类
2. **Pearson 相关性热图** (): 样本间相关性
3. **火山图** (, ): 差异基因/TE 分布（padj 和 pvalue 两种）
4. **MA 图** (, ): Mean-Average plot
5. **热图** (): Top 差异特征表达热图
6. **GO 富集** (gene 模块): BP/CC/MF 三类
7. **GSEA** (gene 模块): Hallmark, KEGG, Reactome, GO_BP/CC/MF

## 依赖

- R ≥ 4.0
- 必需 R 包: DESeq2, edgeR, ggplot2, pheatmap, clusterProfiler, enrichplot, ggrepel, msigdbr 等
- 环境: 

## 与主流程的关系



独立 CLI 工具（,  等）也共用 ，适合已有矩阵的快速补图。

## Snakemake 批量可视化

如需对多个矩阵和 contrast 做可复跑批量可视化，请使用相邻目录：



详见 。
