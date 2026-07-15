# CUT&RUN 输入与 manifest

| 状态 | 维护人 | 最后审查 |
|---|---|---|
| Draft review | CUT&RUN maintainers | 2026-07-16 |

## FASTQ

```text
H3K27ac_WT_1_R1.fastq.gz
H3K27ac_WT_1_R2.fastq.gz
IgG_WT_1_R1.fastq.gz
IgG_WT_1_R2.fastq.gz
```

文件名只辅助初始化，最终设计以 manifest 为准。

## Manifest 必需列

```csv
sample,species,assay,group,replicate,igg,is_igg,layout,fastq_1,fastq_2
```

| 列 | 含义与检查 |
|---|---|
| `sample` | 唯一 sample ID |
| `species` | `hg38` 或 `mm39`，全表一致 |
| `assay` | `cutrun/cuttag/chipseq` |
| `group` | target/mark 与实验组的稳定标签 |
| `replicate` | biological replicate |
| `is_igg` | control 行为 `true`，target 为 `false` |
| `igg` | target 行指向 manifest 中实际存在且 `is_igg=true` 的 control sample |
| `layout` | `PE/SE` |
| `fastq_1/fastq_2` | 可读路径；PE 必须两列完整 |

不要编辑 `manifest/resolved_manifest.csv`；它是运行产物。编辑 `_automation/inputs/manifest.csv` 后使用该文件重新运行。

## 多 control 设计

不同 condition/target 有不同 IgG/Input 时逐行填写 `igg`。可用 `--no-auto-control` 禁止自动分配。`--control-sample` 适合只有一个共享 control 的简单项目。

## Peak 类型与 control

TF 等尖锐结合通常重点解释 narrow peaks；broad histone marks 使用 broad peaks。主流程生成两类结果便于审阅，但下游不得把 narrow/broad 合并成同一集合。无 matched control 虽可能运行 MACS3，但应在报告明确降级证据等级。

## Reference 一致性

FASTA/index、chrom sizes、blacklist、gene/TE annotation 必须同 assembly。RepEnTools 是 CHM13/T2T FASTQ 工作流，不等于 hg38/mm39 主分支。
