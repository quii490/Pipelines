# TE method adapters

This directory adds optional adapters for three complementary repeat-analysis
strategies. They are deliberately isolated from the main `chipseq` environment.

## Methods

- **T3E**: input-derived background for family/subfamily enrichment. It expects
  BAMs retaining secondary mappings and a BED4 RepeatMasker view (the T3E
  parser requires `chrom`, `start`, `end`, `repeat_id`).
- **Allo**: probabilistic/learned allocation of multimappers for locus-level
  regulatory analysis. It should output a corrected BAM before peak calling.
- **RepEnTools**: an independent chm13v2/T2T FASTQ workflow. It is a benchmark,
  not a coordinate-compatible replacement for hg38/mm39. Use
  `run_repentools_fastq.sh` to stage a manifest into its required two-ChIP/two-
  input layout and run the upstream `ret` command.

## Plan/status mode

```bash
bash tools/te_methods/run_te_methods.sh \
  --manifest results/manifest/resolved_manifest.csv \
  --bam-dir results/04_te_bam \
  --out-dir results/09_downstream/te_methods \
  --methods t3e,allo,repentools
```

The adapter writes `method_plan.tsv`, `method_plan.json`, one command plan per
method and `method_status.tsv`. Missing software is reported as `SKIP`; it does
not make the core gene/TE analysis fail.

The T3E environment file follows the upstream tested dependency family. Allo is
kept in a separate TensorFlow environment. RepEnTools must be run with its own
chm13v2 reference and a recorded repository revision.

## Installed environment paths

The standard Pod installation uses:

```text
/path/to/.conda/envs/cutrun-t3e
/path/to/.conda/envs/cutrun-allo
/path/to/.local/share/CUTnRUN-tools/RepEnTools
```

Record the installed revisions and runtime versions before a production run:

```bash
bash tools/te_methods/verify_installation.sh \
  results/09_downstream/te_methods/tool_installation.tsv
```

Run T3E after installing its repository with:

```bash
bash tools/te_methods/run_te_methods.sh \
  --manifest results/manifest/resolved_manifest.csv \
  --bam-dir results/04_te_bam \
  --out-dir results/09_downstream/te_methods \
  --methods t3e \
  --t3e-dir /path/to/.local/share/CUTnRUN-tools/T3E \
  --t3e-python /path/to/.conda/envs/cutrun-t3e/bin/python \
  --t3e-max-bed-reads 1000000 \
  --t3e-iterations 20 \
  --execute
```

T3E 的上游 probability 脚本会按 read length 展开每条多重比对；对大型 TE BAM
会不可扩展。适配器支持确定性的 BED 行数上限并保留 `*.sampled.bed`，参数会写入
命令和日志。`--t3e-max-bed-reads 0` 才使用全部记录；生产运行建议保留上限并报告
采样值。

Allo is auto-detected from the isolated environment and defaults to paired-end
allocation (`allo {input_bam} -seq pe -o {output_bam} -p {threads} --ignore`).
The `--ignore` flag handles valid TE BAMs without an optional PG/collate tag;
pass `--allo-command` to override it. Allo writes SAM, so `{output_bam}` is a
`.sam` path; the adapter additionally validates sorted BAM, BAI and flagstat
outputs when samtools is available. RepEnTools uses a separate CHM13 reference
profile and never mixes coordinates with hg38/mm39. Templates may use
`{sample}`, `{control}`, `{input_bam}`, `{output_bam}`, `{repeat_bed}`,
`{manifest}`, `{out_dir}`, `{species}` and `{threads}`.

RepEnTools FASTQ adapter (LZH defaults shown):

```bash
bash tools/te_methods/run_repentools_fastq.sh \
  --manifest results/manifest/resolved_manifest.csv \
  --out-dir results/09_downstream/te_methods/repentools_fastq \
  --index-dir /path/to/CUTnRUN/resources/RepEnTools/chm13v2/indexes \
  --gtf /path/to/CUTnRUN/resources/RepEnTools/chm13v2/annotation/rmsk.gtf \
  --ret /path/to/.local/share/CUTnRUN-tools/RepEnTools/ret \
  --execute
```

如果不提供 `--target-groups` 和 `--input-samples`，适配器会从 manifest 自动推导：
每个 `group` 取前两个 target replicate，`is_igg=true` 的前两个样本作为 input。
推导结果写入 `manifest_mapping.txt`；如果实验设计不是两组 target + 两个 input，
请显式传入映射。

The adapter exports `HISAT2_INDEXES`, `REPENTOOLS_DIR` and `ADAPTERS` before
calling `ret`; this is required by the upstream script. It validates all eight
index shards and every staged FASTQ before execution.

After the adapters finish, the main pipeline runs
`render_te_method_visuals.py`. It writes `te_methods/visualization_status.tsv`
and method-specific SVG/TSV summaries under `te_methods/visuals/`:

```text
visuals/allo/allo_mapped_fraction.svg
visuals/t3e/t3e_top_features.svg
visuals/repentools_fastq/<group>/repentools/repentools_top_features.svg
```

Plots are produced only from validated outputs. A method with no scientific
result is recorded as `SKIP`; a method marked `PASS` without parseable output
is recorded as `FAIL` so an empty upstream directory cannot look successful.

## Reproducibility rule

Review the generated command plan and run external tools only after recording
the tool revision, reference assembly, RepeatMasker release, mappability source,
and the exact target/control mapping.
