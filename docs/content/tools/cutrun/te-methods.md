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

## RepEnTools

RepEnTools 是 CHM13/T2T FASTQ 工作流，需要符合其设计的两个 ChIP 和两个 input。不能把 hg38/mm39 BAM 直接当作 CHM13 结果。

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

## 只重画外部方法汇总图

```bash
python3 CUTnRUN/tools/te_methods/render_te_method_visuals.py \
  --methods-dir /path/to/results/09_downstream/te_methods \
  --output-dir /path/to/results/09_downstream/te_methods/visuals
```

只有通过验证的真实表格才应绘图。查看 `method_status.tsv` 和 `visualization_status.tsv`，不要把空图解释成无信号。
