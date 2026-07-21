# RNA-seq Quality Control
QC 的目标不是让所有指标都“变绿”，而是识别技术失败、样本错误和会改变结论的偏差。阈值依赖物种、组织、read length、建库和测序深度；本页范围是调查触发条件，最终实验室阈值需维护者确认。
## Reads 层

| 指标 | 正常检查 | 异常时排查 |
|---|---|---|
| Per-base quality/Q30 | read 主体保持较高质量 | 测序批次、末端质量、trim |
| Adapter content | trimming 后明显下降 | adapter 类型、read-through |
| GC distribution | 同类样本大体相近 | 污染、建库偏倚、样本混淆 |
| Overrepresented sequences | 能解释为 rRNA/adapter/高表达转录本 | 污染或建库问题 |
## Alignment 层
*主要看STAR的mapping rate，一般在80%以上为优。若mapping rate低，且multi-mapping高则建库质量可能有问题。*
Mapping rate 应结合 uniquely/multi-mapped、unmapped 原因和同批次样本分布解释。单一样本明显低于其余样本时依次检查：物种/build、FASTQ 完整性、污染、read quality、rRNA 和样本身份。
STAR index、GTF 与 FASTA build 不一致可能仍产生 BAM，但 counts 和 annotation 会失真。
## Assignment 和 strandedness

featureCounts assigned rate 异常低时检查：

1. GTF chromosome 名称是否与 BAM 一致。
2. `--strand` 是否符合 library preparation。
3. feature type/gene attribute 是否存在。
4. BAM 是否来自预期 reference。
5. 大量 reads 是否落在 intronic/intergenic/rRNA 区域。
## Duplicate 和复杂度

高 duplicate rate 可能来自 PCR，也可能来自真实高表达基因。默认只把 duplicate 作为 QC，不自动从普通 RNA-seq counts 删除。比较样本时应结合 library size、unique fragments、expression concentration 和建库批次。
## RNA-specific metrics

- Coding/UTR/intronic/intergenic proportions：解释 RNA 类型、建库和 DNA contamination。
- rRNA proportion：评估 depletion/poly(A) selection。
- 5′/3′ bias 与 coverage：排查降解、建库和 transcript length bias。
## 样本关系
PCA、correlation 和 clustering 用于发现标签错误、batch、离群或样本质量差异。
!!! warning "不能只看 PCA 删除样本"
    排除样本前记录原始 QC、样本身份核对、排除对 contrast/replicates 的影响，以及排除前后结果。避免因为样本“不符合预期结论”而删除。

## 分析完成检查单

- [ ] 所有期望样本进入核心 count matrix。
- [ ] Reads 和 alignment QC 有解释，无静默缺失样本。
- [ ] strandedness 与建库一致。
- [ ] gene counts 非空且分布合理。
- [ ] TE 输出存在或有明确 SKIP 原因。
- [ ] 同组重复关系和 batch 已检查。
- [ ] 样本排除、参数修改和重新运行均有记录。
