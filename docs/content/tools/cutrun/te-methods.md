# CUT&RUN TE 工具与外部方法

| 状态 | 维护人 | 最后验证 |
|---|---|---|
| Draft review | TE method maintainers | 2026-07-16 |

## 检查 TE multi-mapping

```bash
python3 CUTnRUN/tools/cutrun_cli/cutrun_te_multimap.py \
  --bam-dir /path/to/04_te_bam \
  --output /path/to/09_downstream/te_multimap.tsv \
  --max-records 0
```

`0` 扫描全部 alignments，I/O 较大；诊断时可限制行数。检查 secondary alignments、NH tags 和实际 `te-k`，不能只看 BAM 文件存在。

## Strand-aware TE locus heatmap

```bash
bash CUTnRUN/tools/cutrun_cli/cutrun_te_locus_heatmap.sh \
  --regions /path/to/L1.six_column.bed \
  --anchor /path/to/anchor.te_locus_best.bw \
  --signals "/path/to/anchor.bw /path/to/comparison.bw" \
  --labels "Anchor Comparison" \
  --out-prefix /path/to/heatmap/L1 \
  --mode both \
  --before 3000 --after 3000 \
  --body-length 6000 \
  --max-regions 100000 \
  --bin-size 25 \
  --threads 8
```

BED 必须六列且 strand 有效。行按 anchor 在 strand-aware TE 5′ end 周围的信号排序，然后同一顺序复用于所有 tracks。它适合比较形状，不是独立的 locus-specific significance test。

`--max-regions 0` 和很小 bin 会生成巨型 matrix；只有资源充分且确实需要时使用。`--write-values` 默认关闭。

## T3E 与 Allo

```bash
bash CUTnRUN/tools/te_methods/run_te_methods.sh \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --bam-dir /path/to/results/04_te_bam \
  --out-dir /path/to/results/09_downstream/te_methods \
  --methods t3e,allo \
  --execute
```

先运行 plan/status 模式确认依赖。T3E 默认可限制 BED reads 和 iterations；抽样必须确定性并记录。Allo 命令成功后仍要验证 SAM/BAM/BAI 和 flagstat。

常用 plan-only 命令（没有 `--execute`）：

```bash
bash CUTnRUN/tools/te_methods/run_te_methods.sh \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --bam-dir /path/to/results/04_te_bam \
  --out-dir /path/to/results/09_downstream/te_methods \
  --methods t3e,allo \
  --species hg38 \
  --repeat-bed /path/to/repeatmasker.bed \
  --threads 8
```

T3E 关键参数：`--t3e-python`、`--t3e-dir`、`--t3e-iterations`（默认 100）、`--t3e-alpha`（0.05）、`--t3e-enrichment`（1.0）、`--t3e-max-bed-reads`（默认 1,000,000；0=全部）。`--allo-command` 和 `--repentools-command` 是版本相关模板，必须先人工审查 placeholder 与真实工具版本。

## RepEnTools

RepEnTools 是 CHM13/T2T FASTQ 工作流，需要符合其设计的两个 ChIP 和两个 input。不能把 hg38/mm39 BAM 直接当作 CHM13 结果。

第一次使用先建立独立 reference profile：

```bash
bash CUTnRUN/tools/te_methods/prepare_repentools_reference.sh \
  --fasta /path/to/chm13v2.fa.gz \
  --repeat-bed /path/to/chm13v2.repeatmasker.bed \
  --out-dir /path/to/repentools_reference/chm13v2 \
  --threads 16
```

```bash
bash CUTnRUN/tools/te_methods/run_repentools_fastq.sh \
  --manifest /path/to/results/manifest/resolved_manifest.csv \
  --out-dir /path/to/results/09_downstream/te_methods/repentools_fastq \
  --index-dir /path/to/chm13/indexes \
  --gtf /path/to/chm13/repeatmasker.gtf \
  --ret /path/to/RepEnTools/ret \
  --execute
```

非标准设计显式指定 target groups 和 input samples。缺参考、重复或必需输出时状态应为 `SKIP/FAIL`。

```bash
# 仅 stage/validate，不真正运行 ret
bash CUTnRUN/tools/te_methods/run_repentools_fastq.sh \
  --manifest /path/to/resolved_manifest.csv \
  --out-dir /path/to/repentools_run \
  --index-dir /path/to/chm13/indexes \
  --gtf /path/to/chm13/repeatmasker.gtf \
  --ret /path/to/RepEnTools/ret \
  --target-groups "H3K27ac=S1,S2" \
  --input-samples IgG1,IgG2 \
  --threads 8
```

审查 staged FASTQ、group mapping 和 command 后才加 `--execute`。`--force` 会重跑已有完成 sentinel 的 group，不是普通恢复的默认选项。

## 只重画外部方法汇总图

```bash
python3 CUTnRUN/tools/te_methods/render_te_method_visuals.py \
  --methods-dir /path/to/results/09_downstream/te_methods \
  --output-dir /path/to/results/09_downstream/te_methods/visuals
```

只有通过验证的真实表格才应绘图。查看 `method_status.tsv` 和 `visualization_status.tsv`，不要把空图解释成无信号。

## 内部辅助脚本

- `build_repeat_bed.py --te-saf ... --te-anno ... --output ...`：把 pipeline TE resources 转成外部工具需要的 BED；坐标和 ID 映射必须抽查。
- `te_method_plan.py --manifest ... --output-dir ...`：生成 target/control plan，通常由 adapter 调用。
- `verify_installation.sh [output.tsv]`：当前记录维护者预设路径下的 T3E/Allo/RepEnTools 状态，不是跨机器通用安装器；新用户应以实际环境路径和各工具官方安装检查为准。
