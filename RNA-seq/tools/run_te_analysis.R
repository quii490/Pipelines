#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq TE analysis from a TE DE matrix\n\n",
    "Required:\n",
    "  --de-matrix PATH         TE DE matrix CSV\n",
    "  --te-annotation-tsv PATH TE annotation TSV with repName/repFamily/repClass/milliDiv if available\n",
    "  --outdir PATH            output directory\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --prefix STR             output prefix; default basename(de-matrix)\n",
    "  --te-label-level STR     locus_id | repName | repFamily | repClass, default repName\n",
    "  --te-color-level STR     repFamily | repClass, default repFamily\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  --baseMean-min NUM       default 5\n",
    "  --top-n NUM              top groups for plots, default 30\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(te_label_level = "repName", te_color_level = "repFamily", padj_cutoff = "0.05", lfc_cutoff = "0.58", baseMean_min = "5", top_n = "30")
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(key, "--")) stop("Unknown argument: ", key)
    if (i == length(args)) stop("Missing value for: ", key)
    res[[gsub("-", "_", sub("^--", "", key))]] <- args[[i + 1]]
    i <- i + 2
  }
  res
}

read_num <- function(x, default) {
  y <- suppressWarnings(as.numeric(x))
  if (is.na(y)) default else y
}
read_int <- function(x, default) {
  y <- suppressWarnings(as.integer(x))
  if (is.na(y) || y < 1) default else y
}
find_function_file <- function(script_dir, opt) {
  cand <- c(
    if (!is.null(opt$function_file) && nzchar(opt$function_file)) normalizePath(opt$function_file, mustWork = FALSE) else character(0),
    normalizePath(file.path(script_dir, "..", "rnaseq-downstream", "rnaseq-function.R"), mustWork = FALSE),
    file.path(script_dir, "rnaseq-function.R")
  )
  cand <- unique(cand[!is.na(cand) & nzchar(cand)])
  hit <- cand[file.exists(cand)][1]
  if (is.na(hit) || !nzchar(hit)) stop("Cannot find rnaseq-function.R: ", paste(cand, collapse = " | "))
  hit
}

add_te_age <- function(df) {
  md <- if ("milliDiv" %in% colnames(df)) suppressWarnings(as.numeric(df$milliDiv)) else rep(NA_real_, nrow(df))
  df$TE_age <- dplyr::case_when(
    is.na(md) ~ "unknown",
    md < 50 ~ "young (<5%)",
    md < 150 ~ "middle (5-15%)",
    TRUE ~ "old (>=15%)"
  )
  df$TE_age <- factor(df$TE_age, levels = c("young (<5%)", "middle (5-15%)", "old (>=15%)", "unknown"))
  df
}

add_direction <- function(df) {
  has_padj <- "padj" %in% colnames(df) && any(!is.na(df$padj))
  df$direction <- if (has_padj) {
    ifelse(!is.na(df$padj) & df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff,
           ifelse(df$log2FoldChange > 0, "up", "down"), "not_sig")
  } else {
    ifelse(df$baseMean >= baseMean_min & abs(df$log2FoldChange) >= lfc_cutoff,
           ifelse(df$log2FoldChange > 0, "up", "down"), "not_sig")
  }
  df
}

plot_group_counts <- function(df, group_col, outfile, top_n = 30) {
  d <- df %>%
    dplyr::filter(direction %in% c("up", "down")) %>%
    dplyr::mutate(group = as.character(.data[[group_col]])) %>%
    dplyr::filter(!is.na(group), group != "", group != "Unknown") %>%
    dplyr::count(direction, group, name = "n")
  if (nrow(d) == 0) return(invisible(NULL))
  keep <- d %>% dplyr::group_by(group) %>% dplyr::summarise(total = sum(n), .groups = "drop") %>% dplyr::arrange(dplyr::desc(total)) %>% dplyr::slice_head(n = top_n) %>% dplyr::pull(group)
  d <- d %>% dplyr::filter(group %in% keep)
  p <- ggplot(d, aes(x = reorder(group, n), y = n, fill = direction)) +
    geom_col(position = "dodge", width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c(up = "#C95F50", down = "#4F789F")) +
    journal_theme(base_size = 12) +
    labs(title = paste0(group_col, " differential TE counts"), x = NULL, y = "TE features", fill = "Direction")
  h <- max(5.5, 0.25 * length(keep) + 2)
  ggsave(outfile, p, width = 8.5, height = h, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", outfile), p, width = 8.5, height = h, dpi = 300, bg = "white", limitsize = FALSE)
}

plot_mean_lfc_heatmap <- function(df, outfile) {
  cols <- intersect(c("repClass", "repFamily", "TE_age"), colnames(df))
  if (!all(c("repClass", "repFamily") %in% cols)) return(invisible(NULL))
  d <- df %>%
    dplyr::filter(!is.na(repClass), repClass != "", !is.na(repFamily), repFamily != "") %>%
    dplyr::group_by(repClass, repFamily) %>%
    dplyr::summarise(mean_log2FC = mean(log2FoldChange, na.rm = TRUE), n = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(abs(mean_log2FC))) %>%
    dplyr::slice_head(n = 60)
  if (nrow(d) < 2) return(invisible(NULL))
  mat_df <- d %>% dplyr::select(repFamily, repClass, mean_log2FC) %>% tidyr::pivot_wider(names_from = repClass, values_from = mean_log2FC, values_fill = 0)
  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- mat_df$repFamily
  pheatmap::pheatmap(mat, filename = outfile, main = "TE family mean log2FC by class", width = 8, height = max(6, 0.18 * nrow(mat) + 2))
}

opt <- parse_args(args)
required <- c("de_matrix", "te_annotation_tsv", "outdir")
missing <- required[!vapply(required, function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_file <- find_function_file(script_dir, opt)
de_file <- normalizePath(opt$de_matrix, mustWork = FALSE)
te_file <- normalizePath(opt$te_annotation_tsv, mustWork = FALSE)
if (!file.exists(de_file)) stop("DE matrix not found: ", de_file)
if (!file.exists(te_file)) stop("TE annotation not found: ", te_file)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- read_num(opt$padj_cutoff, 0.05)
lfc_cutoff <- read_num(opt$lfc_cutoff, 0.58)
baseMean_min <- read_num(opt$baseMean_min, 5)
label_top_n <- 40L
heatmap_top_n <- 40L
assign("outdir", outdir, envir = .GlobalEnv)
assign("padj_cutoff", padj_cutoff, envir = .GlobalEnv)
assign("lfc_cutoff", lfc_cutoff, envir = .GlobalEnv)
assign("baseMean_min", baseMean_min, envir = .GlobalEnv)
assign("label_top_n", label_top_n, envir = .GlobalEnv)
assign("heatmap_top_n", heatmap_top_n, envir = .GlobalEnv)
source(function_file)

prefix <- if (!is.null(opt$prefix) && nzchar(opt$prefix)) opt$prefix else sub("\\.[^.]+$", "", basename(de_file))
prefix <- sanitize_name(prefix)
res <- utils::read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"feature_id" %in% colnames(res)) colnames(res)[1] <- "feature_id"
te_anno <- read_te_annotation(te_file)
res <- annotate_te_res(
  res, te_anno,
  preferred_match_cols = unique(c(opt$te_label_level, opt$te_color_level, "repName", "repFamily", "repClass", "locus_id", "gene_name", "feature_id")),
  label_level = opt$te_label_level,
  color_level = opt$te_color_level
)
res <- add_te_age(add_direction(res))
utils::write.csv(res, file.path(outdir, paste0(prefix, ".TE_annotated_DE_matrix.csv")), row.names = FALSE)

for (col in intersect(c("repClass", "repFamily", "repName", "TE_age"), colnames(res))) {
  tab <- res %>%
    dplyr::mutate(group = as.character(.data[[col]])) %>%
    dplyr::filter(!is.na(group), group != "") %>%
    dplyr::count(direction, group, name = "n") %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(fraction = n / sum(n)) %>%
    dplyr::ungroup()
  utils::write.csv(tab, file.path(outdir, paste0(prefix, ".", col, ".direction_counts.csv")), row.names = FALSE)
  plot_group_counts(res, col, file.path(outdir, paste0(prefix, ".", col, ".up_down_bar.pdf")), top_n = read_int(opt$top_n, 30))
}
plot_mean_lfc_heatmap(res, file.path(outdir, paste0(prefix, ".TE_family_class_mean_log2FC_heatmap.pdf")))

writeLines(c(
  paste0("de_matrix=", de_file),
  paste0("te_annotation_tsv=", te_file),
  paste0("features=", nrow(res))
), file.path(outdir, paste0(prefix, ".te_analysis.summary.txt")))
message("[run_te_analysis] DONE: ", outdir)
