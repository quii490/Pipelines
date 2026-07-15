#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

declare -a shell_tools=(
  "RNA-seq/tools/run_bam_to_counts.sh|docs/content/tools/rnaseq/bam-to-counts.md|rnaseq-bam-to-counts"
  "RNA-seq/snakemake-visualization/run_snakemake_visualization.sh|docs/content/tools/rnaseq/snakemake.md|run_snakemake_visualization.sh"
  "ATAC-seq/scripts/run_atac_downstream_only.sh|docs/content/tools/atacseq/rebuild.md"
  "ATAC-seq/scripts/run_atac_from_bam.sh|docs/content/tools/atacseq/rebuild.md"
  "ATAC-seq/scripts/run_callpeak_from_bam.sh|docs/content/tools/atacseq/rebuild.md"
  "ATAC-seq/scripts/run_fixedbin_from_bam.sh|docs/content/tools/atacseq/rebuild.md"
  "ATAC-seq/scripts/atacseq|docs/content/tools/atacseq/command-reference.md|atacseq"
  "ATAC-seq/scripts/run_tss_qc.sh|docs/content/tools/atacseq/signal-qc.md"
  "ATAC-seq/scripts/run_gene_body_profile.sh|docs/content/tools/atacseq/signal-qc.md"
  "ATAC-seq/scripts/run_bw_correlation.sh|docs/content/tools/atacseq/signal-qc.md"
  "ATAC-seq/scripts/run_nuc_phasing.sh|docs/content/tools/atacseq/signal-qc.md"
  "ATAC-seq/scripts/run_region_heatmap.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_diff_peak_heatmap.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_peak_annotation.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_peak_overlap.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_motif_homer.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_tobias.sh|docs/content/tools/atacseq/regions-motif.md"
  "ATAC-seq/scripts/run_atac_profile_heatmaps.sh|docs/content/tools/atacseq/reporting-comparison.md"
  "ATAC-seq/scripts/run_cross_contrast_heatmap.sh|docs/content/tools/atacseq/reporting-comparison.md"
  "ATAC-seq/scripts/run_atac_te_tracks_from_bam.sh|docs/content/tools/atacseq/te.md"
  "ATAC-seq/scripts/run_atac_te_tracks_from_fastq.sh|docs/content/tools/atacseq/te.md"
  "ATAC-seq/scripts/run_te_heatmap.sh|docs/content/tools/atacseq/te.md"
  "ATAC-seq/scripts/run_te_heatmap_batch.sh|docs/content/tools/atacseq/te.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_annotate_peaks|docs/content/tools/cutrun/annotation-enrichment.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_bam_cor.sh|docs/content/tools/cutrun/signal-peaks.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_bw_cor.sh|docs/content/tools/cutrun/signal-peaks.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_dt_heatmap.sh|docs/content/tools/cutrun/signal-peaks.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_homer.sh|docs/content/tools/cutrun/annotation-enrichment.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_te_locus_heatmap.sh|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/te_methods/prepare_repentools_reference.sh|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/te_methods/run_repentools_fastq.sh|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/te_methods/run_te_methods.sh|docs/content/tools/cutrun/te-methods.md"
)

for item in "${shell_tools[@]}"; do
  tool="${item%%|*}"
  rest="${item#*|}"
  doc="${rest%%|*}"
  token="${rest#*|}"
  [[ "$token" != "$rest" ]] || token="$(basename "$tool")"
  test -f "$tool" || { echo "ERROR: missing tool: $tool" >&2; exit 1; }
  test -f "$doc" || { echo "ERROR: missing tool documentation: $doc" >&2; exit 1; }
  bash "$tool" --help >/dev/null
  grep -F -q -- "$token" "$doc" || {
    echo "ERROR: $doc does not mention $token" >&2
    exit 1
  }
done

declare -a python_tools=(
  "ATAC-seq/scripts/generate_qc_report.py|docs/content/tools/atacseq/reporting-comparison.md"
  "ATAC-seq/scripts/organize_atac_downstream_outputs.py|docs/content/tools/atacseq/reporting-comparison.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_consensus_peaks.py|docs/content/tools/cutrun/annotation-enrichment.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_peak_overlap.py|docs/content/tools/cutrun/annotation-enrichment.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_preflight.py|docs/content/tools/cutrun/validation-status.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_report.py|docs/content/tools/cutrun/recovery-report.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_results_summary.py|docs/content/tools/cutrun/validation-status.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_run_manifest.py|docs/content/tools/cutrun/recovery-report.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_status.py|docs/content/tools/cutrun/validation-status.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_te_multimap.py|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/cutrun_cli/cutrun_te_qc.py|docs/content/tools/cutrun/validation-status.md"
  "CUTnRUN/tools/te_methods/build_repeat_bed.py|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/te_methods/render_te_method_visuals.py|docs/content/tools/cutrun/te-methods.md"
  "CUTnRUN/tools/te_methods/te_method_plan.py|docs/content/tools/cutrun/te-methods.md"
)

for item in "${python_tools[@]}"; do
  tool="${item%%|*}"
  doc="${item#*|}"
  test -f "$tool" || { echo "ERROR: missing tool: $tool" >&2; exit 1; }
  python3 "$tool" --help >/dev/null
  grep -F -q -- "$(basename "$tool")" "$doc" || {
    echo "ERROR: $doc does not mention $(basename "$tool")" >&2
    exit 1
  }
done

echo "Common tool help checks passed."
