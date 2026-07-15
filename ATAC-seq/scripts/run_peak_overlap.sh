#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  run_peak_overlap.sh --peaks "A=a.bed B=b.bed C=c.bed" --outdir out/peak_overlap [options]
  run_peak_overlap.sh --peak-dir 05_peaks/filtered --outdir out/peak_overlap [options]

Purpose:
  Compare peak sets using bedtools multiinter and produce overlap tables/plots.

Options:
  --peaks STR                    Space-separated label=file entries.
  --peak-dir DIR                 Alternative: scan .bed/.narrowPeak/.broadPeak files.
  --outdir DIR                   Output directory.
  --distance INT                 Merge nearby intervals before comparison. Default: 0.
  --top-n INT                    Keep top N rows from each input. Default: all.
USAGE
}

peaks_arg=""
peak_dir=""
outdir=""
distance=0
top_n=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --peaks) peaks_arg="$2"; shift 2 ;;
    --peak-dir) peak_dir="$2"; shift 2 ;;
    --outdir) outdir="$2"; shift 2 ;;
    --distance) distance="$2"; shift 2 ;;
    --top-n) top_n="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown option $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ -n "$outdir" ]] || { usage >&2; exit 1; }
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found in PATH" >&2; exit 1; }

labels=()
files=()
sanitize() {
  echo "$1" | sed -E 's/[[:space:]]+/_/g; s/[^A-Za-z0-9_.-]+/_/g; s/^_+|_+$//g'
}

if [[ -n "$peaks_arg" ]]; then
  read -r -a entries <<< "$peaks_arg"
  for entry in "${entries[@]}"; do
    [[ "$entry" == *=* ]] || { echo "ERROR: --peaks entries must be label=file: $entry" >&2; exit 1; }
    label="$(sanitize "${entry%%=*}")"
    file="${entry#*=}"
    labels+=( "$label" )
    files+=( "$file" )
  done
else
  [[ -d "$peak_dir" ]] || { echo "ERROR: peak dir not found: $peak_dir" >&2; exit 1; }
  shopt -s nullglob
  scan=( "$peak_dir"/*.bed "$peak_dir"/*.narrowPeak "$peak_dir"/*.broadPeak )
  for file in "${scan[@]}"; do
    label="$(basename "$file")"
    label="${label%.bed}"
    label="${label%.narrowPeak}"
    label="${label%.broadPeak}"
    labels+=( "$(sanitize "$label")" )
    files+=( "$file" )
  done
fi

[[ ${#files[@]} -ge 2 ]] || { echo "ERROR: at least two peak files are required" >&2; exit 1; }
for file in "${files[@]}"; do
  [[ -s "$file" ]] || { echo "ERROR: peak file not found: $file" >&2; exit 1; }
done

prepared="$outdir/prepared"
tables="$outdir/tables"
plots="$outdir/plots"
mkdir -p "$prepared" "$tables" "$plots"

prepared_files=()
for i in "${!files[@]}"; do
  label="${labels[$i]}"
  file="${files[$i]}"
  tmp="$prepared/${label}.sorted.bed"
  if [[ -n "$top_n" ]]; then
    awk -v n="$top_n" 'BEGIN{OFS="\t"} $1 !~ /^#/ && NF >= 3 {print $1,$2,$3; c++; if(c>=n) exit}' "$file" | bedtools sort -i - > "$tmp"
  else
    awk 'BEGIN{OFS="\t"} $1 !~ /^#/ && NF >= 3 {print $1,$2,$3}' "$file" | bedtools sort -i - > "$tmp"
  fi
  if [[ "$distance" -gt 0 ]]; then
    merged="$prepared/${label}.merged.bed"
    bedtools merge -d "$distance" -i "$tmp" > "$merged"
    tmp="$merged"
  fi
  prepared_files+=( "$tmp" )
done

multi="$tables/peak_multiinter.tsv"
bedtools multiinter -header -names "${labels[@]}" -i "${prepared_files[@]}" > "$multi"
awk 'BEGIN{OFS="\t"} NR>1 {print $1,$2,$3,"cluster_"NR-1,$4,$5}' "$multi" > "$tables/peak_union_clusters.bed"
awk 'BEGIN{OFS="\t"} NR>1 {a[$5]+=$4} END{print "intersection","count"; for(k in a) print k,a[k]}' "$multi" > "$tables/intersection_counts.tsv"

if command -v Rscript >/dev/null; then
  Rscript - "$multi" "$plots" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
multi <- args[[1]]
plot_dir <- args[[2]]
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
save_both <- function(path, p, w = 7, h = 4.8) {
  ggsave(path, p, width = w, height = h, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", path), p, width = w, height = h, dpi = 300, limitsize = FALSE)
}
x <- read.delim(multi, check.names = FALSE)
label_cols <- setdiff(colnames(x), c("chrom", "start", "end", "num", "list"))
mat <- as.matrix(x[, label_cols, drop = FALSE])
mode(mat) <- "numeric"
jacc <- matrix(NA_real_, nrow = ncol(mat), ncol = ncol(mat), dimnames = list(label_cols, label_cols))
for (i in seq_along(label_cols)) {
  for (j in seq_along(label_cols)) {
    both <- sum(mat[, i] == 1 & mat[, j] == 1)
    either <- sum(mat[, i] == 1 | mat[, j] == 1)
    jacc[i, j] <- ifelse(either == 0, NA, both / either)
  }
}
jdf <- as.data.frame(as.table(jacc))
colnames(jdf) <- c("sample1", "sample2", "jaccard")
p1 <- ggplot(jdf, aes(sample1, sample2, fill = jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", jaccard)), size = 3) +
  scale_fill_viridis_c(na.value = "grey90") +
  coord_equal() +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Peak set Jaccard overlap", x = NULL, y = NULL, fill = "Jaccard")
save_both(file.path(plot_dir, "peak_jaccard_heatmap.pdf"), p1, 6.5, 5.8)
counts <- x |>
  count(list, wt = num, name = "count") |>
  arrange(desc(count)) |>
  slice_head(n = 25)
p2 <- ggplot(counts, aes(reorder(list, count), count)) +
  geom_col(fill = "#4C78A8") +
  coord_flip() +
  theme_bw(base_size = 11) +
  labs(title = "Top peak intersections", x = "Intersection", y = "Cluster count")
save_both(file.path(plot_dir, "peak_intersection_bar.pdf"), p2, 7.5, 5)
RSCRIPT
fi

printf "Finished peak overlap: %s\n" "$outdir"
