#!/usr/bin/env bash
set -euo pipefail

reference="${reference:-hg38}"
outputdir="${outputdir:-}"
gene_bam_dir=""
te_bam_dir=""
gene_counts=""
te_counts=""
manifest_file=""
threads="${threads:-8}"
min_frag="${min_frag:-30}"
max_frag="${max_frag:-1200}"
CHIPSEQ_ENV="${CHIPSEQ_ENV:-chipseq}"
DOWNSTREAM_ENV="${DOWNSTREAM_ENV:-downstream}"
TE_classes_of_interest="${TE_classes_of_interest:-LTR,LINE,SINE}"
MIN_TOTAL_NORMALIZED_COUNTS="${MIN_TOTAL_NORMALIZED_COUNTS:-10}"
N_LABELS_GENES="${N_LABELS_GENES:-20}"
N_LABELS_TE_FAMILY="${N_LABELS_TE_FAMILY:-20}"
N_LABELS_TE_REPNAME="${N_LABELS_TE_REPNAME:-40}"
TE_repname_for_boxplot_and_heatmap="${TE_repname_for_boxplot_and_heatmap:-}"
FORCE_R_REPROCESS="${FORCE_R_REPROCESS:-false}"
SKIP_FEATURECOUNTS="${SKIP_FEATURECOUNTS:-false}"

# Run a tool from either an explicit conda-prefix or a named conda environment.
# The Pod has environments in both /opt/conda/envs and ~/.conda/envs; using
# `conda run -n` with a prefix silently resolves against the wrong conda base.
run_env_command() {
  local env_spec="$1"
  local exe="$2"
  shift 2
  if [[ "$env_spec" == /* && -x "$env_spec/bin/$exe" ]]; then
    "$env_spec/bin/$exe" "$@"
  else
    conda run --no-capture-output -n "$env_spec" "$exe" "$@"
  fi
}

usage() {
  cat <<'EOF'
Usage:
  cutrun_count_draw --reference hg38 --output-dir OUT \
    --gene-bam-dir CLEAN_BAMS --te-bam-dir TE_BAMS

Inputs can be BAM directories or existing count matrices:
  --gene-bam-dir DIR       Strict clean gene BAMs
  --te-bam-dir DIR         TE relaxed BAMs with NH tags
  --gene-counts FILE       Existing gene featureCounts matrix
  --te-counts FILE         Existing TE featureCounts matrix
  --bam-dir DIR            Legacy: use one BAM directory for both branches

Other:
  --reference hg38|mm39
  --output-dir DIR
  --threads INT
  --force-r-reprocess
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --reference) reference="$2"; shift 2 ;;
    --output-dir) outputdir="$2"; shift 2 ;;
    --gene-bam-dir) gene_bam_dir="$2"; shift 2 ;;
    --te-bam-dir) te_bam_dir="$2"; shift 2 ;;
    --gene-counts) gene_counts="$2"; shift 2 ;;
    --te-counts) te_counts="$2"; shift 2 ;;
    --manifest) manifest_file="$2"; shift 2 ;;
    --bam-dir)
      gene_bam_dir="$2"
      te_bam_dir="$2"
      echo "WARNING: --bam-dir is legacy; prefer separate --gene-bam-dir and --te-bam-dir" >&2
      shift 2
      ;;
    --te-classes) TE_classes_of_interest="$2"; shift 2 ;;
    --te-repnames) TE_repname_for_boxplot_and_heatmap="$2"; shift 2 ;;
    --threads) threads="$2"; shift 2 ;;
    --chipseq-env) CHIPSEQ_ENV="$2"; shift 2 ;;
    --downstream-env) DOWNSTREAM_ENV="$2"; shift 2 ;;
    --min-frag) min_frag="$2"; shift 2 ;;
    --max-frag) max_frag="$2"; shift 2 ;;
    --min-total-normalized-counts) MIN_TOTAL_NORMALIZED_COUNTS="$2"; shift 2 ;;
    --n-labels-genes) N_LABELS_GENES="$2"; shift 2 ;;
    --n-labels-te-family) N_LABELS_TE_FAMILY="$2"; shift 2 ;;
    --n-labels-te-repname) N_LABELS_TE_REPNAME="$2"; shift 2 ;;
    --skip-featurecounts) SKIP_FEATURECOUNTS=true; shift ;;
    --force-r-reprocess) FORCE_R_REPROCESS=true; shift ;;
    --filter-samples|--gene-boxplot-samples) shift 2 ;;
    --te-only)
      echo "WARNING: --te-only is deprecated; unified output requires gene and TE inputs" >&2
      shift
      ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$outputdir" ]] || { echo "ERROR: --output-dir is required" >&2; exit 1; }

source_path="$(readlink -f "${BASH_SOURCE[0]}")"
script_dir="$(cd "$(dirname "$source_path")" && pwd)"
repo_root="$(cd "${script_dir}/../.." && pwd)"
anno_dir="${CUTRUN_ANNO_DIR:-${repo_root}/resources/CUTRUN_analysis/anno}"

case "$reference" in
  hg38|hs|human)
    REPEAT_ANNO="${anno_dir}/te_anno_hg38.tsv"
    GENE_ANNO="${anno_dir}/gene_anno_ensembl_hg38.tsv"
    DENY="${anno_dir}/hg38-blacklist.v2.bed"
    GENE_SAF="${anno_dir}/gene_hg38.saf"
    TE_SAF="${anno_dir}/te_hg38.saf"
    ;;
  mm39|mm|mouse)
    REPEAT_ANNO="${anno_dir}/te_anno_mm39.tsv"
    GENE_ANNO="${anno_dir}/gene_anno_ensembl_mm39.tsv"
    DENY="${anno_dir}/mm39-blacklist.v2.bed"
    GENE_SAF="${anno_dir}/gene_mm39.saf"
    TE_SAF="${anno_dir}/te_mm39.saf"
    ;;
  *) echo "ERROR: unsupported reference: $reference" >&2; exit 1 ;;
esac

for path in "$REPEAT_ANNO" "$GENE_ANNO" "$DENY" "$GENE_SAF" "$TE_SAF"; do
  [[ -s "$path" ]] || { echo "ERROR: missing annotation: $path" >&2; exit 1; }
done
if awk 'substr($1,1,4)=="ENSG" || substr($1,1,7)=="ENSMUSG"{found=1; exit} END{exit !found}' "$TE_SAF"; then
  echo "ERROR: TE SAF contains gene IDs: $TE_SAF" >&2
  exit 1
fi

mkdir -p "${outputdir}/featurecount" "${outputdir}/results"
gene_counts="${gene_counts:-${outputdir}/featurecount/featurecounts_gene.txt}"
te_counts="${te_counts:-${outputdir}/featurecount/featurecounts_te.txt}"

collect_bams() {
  local mode="$1"
  local directory="$2"
  local -n result="$3"
  [[ -d "$directory" ]] || { echo "ERROR: BAM directory not found: $directory" >&2; exit 1; }
  shopt -s nullglob
  if [[ "$mode" == gene ]]; then
    result=("$directory"/*_clean.bam)
    [[ ${#result[@]} -gt 0 ]] || result=("$directory"/*.sorted.bam)
  else
    result=("$directory"/*_te.bam)
  fi
  [[ ${#result[@]} -gt 0 ]] || result=("$directory"/*.bam)
  [[ ${#result[@]} -gt 0 ]] || { echo "ERROR: no BAMs found in $directory" >&2; exit 1; }
}

if [[ "$SKIP_FEATURECOUNTS" != true && ! -s "$gene_counts" ]]; then
  [[ -n "$gene_bam_dir" ]] || { echo "ERROR: provide --gene-bam-dir or --gene-counts" >&2; exit 1; }
  collect_bams gene "$gene_bam_dir" GENE_BAMS
  run_env_command "$CHIPSEQ_ENV" featureCounts \
    -a "$GENE_SAF" -o "$gene_counts" -F SAF --ignoreDup \
    -p --countReadPairs \
    -T "$threads" "${GENE_BAMS[@]}"
fi

if [[ "$SKIP_FEATURECOUNTS" != true && ! -s "$te_counts" ]]; then
  [[ -n "$te_bam_dir" ]] || { echo "ERROR: provide --te-bam-dir or --te-counts" >&2; exit 1; }
  collect_bams te "$te_bam_dir" TE_BAMS
  run_env_command "$CHIPSEQ_ENV" featureCounts \
    -a "$TE_SAF" -o "$te_counts" -F SAF --ignoreDup \
    -p --countReadPairs \
    -T "$threads" -M --fraction "${TE_BAMS[@]}"
fi

[[ -s "$gene_counts" ]] || { echo "ERROR: missing gene counts: $gene_counts" >&2; exit 1; }
[[ -s "$te_counts" ]] || { echo "ERROR: missing TE counts: $te_counts" >&2; exit 1; }

csv_to_r_vector() {
  awk -v s="$1" 'BEGIN {
    if (s == "") { printf "character(0)"; exit }
    n = split(s, a, ","); printf "c(";
    for (i = 1; i <= n; i++) {
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", a[i]);
      printf "%s\"%s\"", (i == 1 ? "" : ","), a[i]
    }
    printf ")"
  }'
}

mapfile -t SAMPLE_NAMES < <(
  awk 'BEGIN{FS="\t"} /^Geneid/{for(i=7;i<=NF;i++){n=$i; sub(/^.*\//,"",n); sub(/(_clean_te|_clean|_te|\.sorted_te|\.sorted)\.bam$/,"",n); print n}}' "$gene_counts"
)
if [[ -n "$manifest_file" && ! -s "$manifest_file" ]]; then
  echo "ERROR: manifest not found: $manifest_file" >&2
  exit 1
fi
CONTROL_NAMES=()
TARGET_NAMES=()
for sample in "${SAMPLE_NAMES[@]}"; do
  if [[ "$sample" =~ [Ii][Gg][Gg]|[Ii]nput|[Cc]ontrol|[Cc]trl ]]; then
    CONTROL_NAMES+=("$sample")
  else
    TARGET_NAMES+=("$sample")
  fi
done
[[ ${#CONTROL_NAMES[@]} -gt 0 ]] || { echo "ERROR: no IgG/Input/control sample detected" >&2; exit 1; }

join_by_comma() { local IFS=,; echo "$*"; }
named_map_to_r_vector() {
  if [[ "$#" -eq 0 ]]; then
    printf 'c()'
    return
  fi
  printf 'c('
  local i=0 pair key value
  for pair in "$@"; do
    key="${pair%%=*}"
    value="${pair#*=}"
    [[ "$i" -gt 0 ]] && printf ','
    printf '"%s"="%s"' "$key" "$value"
    i=$((i + 1))
  done
  printf ')'
}
TARGET_CONTROL_PAIRS=()
if [[ -n "$manifest_file" ]]; then
  mapfile -t TARGET_CONTROL_PAIRS < <(
    awk -F',' 'NR > 1 && tolower($7) == "false" && $1 != "" && $6 != "" {print $1 "=" $6}' "$manifest_file"
  )
fi
if [[ ${#TARGET_CONTROL_PAIRS[@]} -eq 0 && ${#CONTROL_NAMES[@]} -gt 1 ]]; then
  echo "WARNING: no manifest target-control map; assigning first control (${CONTROL_NAMES[0]}) to every target" >&2
  for target in "${TARGET_NAMES[@]}"; do
    TARGET_CONTROL_PAIRS+=("${target}=${CONTROL_NAMES[0]}")
  done
fi
analysis_fingerprint="${ANALYSIS_FINGERPRINT:-unified_count_draw_v2}"
CONFIG_R="${outputdir}/featurecount/config_count_draw.R"
cat > "$CONFIG_R" <<EOF
repeat_annotations_file <- "$REPEAT_ANNO"
gene_annotations_file <- "$GENE_ANNO"
deny_list_file <- "$DENY"
gene_saf_file <- "$GENE_SAF"
te_saf_file <- "$TE_SAF"
gene_counts_file <- "$gene_counts"
te_counts_file <- "$te_counts"
output_dir <- "${outputdir}/results"
intermediate_dir <- file.path(output_dir, "intermediate_data")
figures_dir <- file.path(output_dir, "figures")
MIN_TOTAL_NORMALIZED_COUNTS <- $MIN_TOTAL_NORMALIZED_COUNTS
N_LABELS_GENES <- $N_LABELS_GENES
N_LABELS_TE_FAMILY <- $N_LABELS_TE_FAMILY
N_LABELS_TE_REPNAME <- $N_LABELS_TE_REPNAME
TE_CLASSES_OF_INTEREST <- $(csv_to_r_vector "$TE_classes_of_interest")
TE_repname_OI <- $(csv_to_r_vector "$TE_repname_for_boxplot_and_heatmap")
CONTROL_SAMPLES <- $(csv_to_r_vector "$(join_by_comma "${CONTROL_NAMES[@]}")")
TARGET_SAMPLES <- $(csv_to_r_vector "$(join_by_comma "${TARGET_NAMES[@]}")")
TARGET_CONTROL_MAP <- $(named_map_to_r_vector "${TARGET_CONTROL_PAIRS[@]}")
ANALYSIS_FINGERPRINT <- "$analysis_fingerprint"
cache_files <- list(
  gene_counts = file.path(intermediate_dir, "01_gene_counts.rds"),
  te_counts = file.path(intermediate_dir, "01_te_counts.rds"),
  gene_plot_data = file.path(intermediate_dir, "02_gene_plot_data.rds"),
  te_plot_data = file.path(intermediate_dir, "02_te_plot_data.rds"),
  plot_data = file.path(intermediate_dir, "02_plot_data.rds"),
  fingerprint = file.path(intermediate_dir, "analysis_fingerprint.txt")
)
EOF

if [[ "$FORCE_R_REPROCESS" == true ]]; then
  rm -rf "${outputdir}/results/intermediate_data"
fi

run_env_command "$DOWNSTREAM_ENV" Rscript \
  "${repo_root}/pipelines/chipseq_auto_nf/analysis/run.R" "$CONFIG_R"

echo "Unified CUT&RUN analysis finished: ${outputdir}/results"
