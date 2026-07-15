#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq DE matrix annotation and visualization\n\n",
    "Required:\n",
    "  --de-matrix PATH         DESeq2-style CSV with feature_id, baseMean, log2FoldChange, padj/pvalue\n",
    "  --outdir PATH            output directory\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --prefix STR             output prefix; default basename(de-matrix)\n",
    "  --annotation-mode STR    auto | gene | te | generic, default auto\n",
    "  --tx2gene-path PATH      gene GTF/tx2gene; adds gene_name/gene_type and gene-type plots\n",
    "  --te-annotation-tsv PATH TE annotation TSV; adds repName/family/class/age plots\n",
    "  --te-label-level STR     locus_id | repName | repFamily | repClass, default repName\n",
    "  --te-color-level STR     repFamily | repClass | TE_age, default repFamily\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  --baseMean-min NUM       default 5\n",
    "  --label-top-n NUM        default 40\n",
    "  --volcano-orientation STR classic | horizontal, default classic\n",
    "  --gray-nonsig BOOL       true|false; color non-significant points gray, default true\n",
    "  --top-n NUM              top entries for bar plots, default 25\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(
    annotation_mode = "auto",
    te_label_level = "repName",
    te_color_level = "repFamily",
    padj_cutoff = "0.05",
    lfc_cutoff = "0.58",
    baseMean_min = "5",
    label_top_n = "40",
    volcano_orientation = "classic",
    gray_nonsig = "true",
    top_n = "25"
  )
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
  if (is.na(y) || y < 0) default else y
}

read_bool <- function(x, default = TRUE) {
  y <- tolower(trimws(as.character(x)))
  if (y %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (y %in% c("false", "f", "no", "n", "0")) return(FALSE)
  default
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

sig_direction <- function(df) {
  has_padj <- "padj" %in% colnames(df) && any(!is.na(df$padj))
  if (has_padj) {
    ifelse(!is.na(df$padj) & df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff,
           ifelse(df$log2FoldChange > 0, "up", "down"), "not_sig")
  } else {
    ifelse(df$baseMean >= baseMean_min & abs(df$log2FoldChange) >= lfc_cutoff,
           ifelse(df$log2FoldChange > 0, "up", "down"), "not_sig")
  }
}

add_te_age <- function(df) {
  if (!"milliDiv" %in% colnames(df)) {
    df$TE_age <- "unknown"
  } else {
    md <- suppressWarnings(as.numeric(df$milliDiv))
    df$TE_age <- dplyr::case_when(
      is.na(md) ~ "unknown",
      md < 50 ~ "young (<5%)",
      md < 150 ~ "middle (5-15%)",
      TRUE ~ "old (>=15%)"
    )
  }
  df$TE_age <- factor(df$TE_age, levels = c("young (<5%)", "middle (5-15%)", "old (>=15%)", "unknown"))
  df
}

plot_direction_bar <- function(df, group_col, outfile, title, top_n = 25) {
  if (!group_col %in% colnames(df)) return(invisible(NULL))
  d <- df %>%
    dplyr::filter(.data$direction %in% c("up", "down")) %>%
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
    labs(title = title, x = NULL, y = "Differential features", fill = "Direction")
  ggsave(outfile, p, width = 8.5, height = max(5.5, 0.25 * length(keep) + 2), limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", outfile), p, width = 8.5, height = max(5.5, 0.25 * length(keep) + 2), dpi = 300, bg = "white", limitsize = FALSE)
}

plot_fraction_bar <- function(df, group_col, outfile, title) {
  if (!group_col %in% colnames(df)) return(invisible(NULL))
  d <- df %>%
    dplyr::filter(.data$direction %in% c("up", "down")) %>%
    dplyr::mutate(group = as.character(.data[[group_col]])) %>%
    dplyr::filter(!is.na(group), group != "", group != "Unknown") %>%
    dplyr::count(direction, group, name = "n") %>%
    dplyr::group_by(direction) %>%
    dplyr::mutate(frac = n / sum(n)) %>%
    dplyr::ungroup()
  if (nrow(d) == 0) return(invisible(NULL))
  p <- ggplot(d, aes(x = direction, y = frac, fill = group)) +
    geom_col(width = 0.7) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = journal_palette(unique(d$group))) +
    journal_theme(base_size = 12) +
    labs(title = title, x = NULL, y = "Fraction", fill = group_col)
  ggsave(outfile, p, width = 8, height = 5.5, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", outfile), p, width = 8, height = 5.5, dpi = 300, bg = "white", limitsize = FALSE)
}

opt <- parse_args(args)
if (is.null(opt$de_matrix) || !nzchar(opt$de_matrix) || is.null(opt$outdir) || !nzchar(opt$outdir)) {
  usage()
  stop("--de-matrix and --outdir are required")
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_file <- find_function_file(script_dir, opt)

de_file <- normalizePath(opt$de_matrix, mustWork = FALSE)
if (!file.exists(de_file)) stop("DE matrix not found: ", de_file)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- read_num(opt$padj_cutoff, 0.05)
lfc_cutoff <- read_num(opt$lfc_cutoff, 0.58)
baseMean_min <- read_num(opt$baseMean_min, 5)
label_top_n <- read_int(opt$label_top_n, 40)
volcano_orientation <- tolower(as.character(opt$volcano_orientation))
if (!volcano_orientation %in% c("classic", "horizontal")) stop("--volcano-orientation must be classic or horizontal")
gray_nonsig <- read_bool(opt$gray_nonsig, TRUE)
heatmap_top_n <- 40L
assign("outdir", outdir, envir = .GlobalEnv)
assign("padj_cutoff", padj_cutoff, envir = .GlobalEnv)
assign("lfc_cutoff", lfc_cutoff, envir = .GlobalEnv)
assign("baseMean_min", baseMean_min, envir = .GlobalEnv)
assign("label_top_n", label_top_n, envir = .GlobalEnv)
assign("heatmap_top_n", heatmap_top_n, envir = .GlobalEnv)
assign("volcano_orientation", volcano_orientation, envir = .GlobalEnv)
assign("gray_nonsig", gray_nonsig, envir = .GlobalEnv)

source(function_file)

prefix <- if (!is.null(opt$prefix) && nzchar(opt$prefix)) opt$prefix else sub("\\.[^.]+$", "", basename(de_file))
prefix <- sanitize_name(prefix)

res <- utils::read.csv(de_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"feature_id" %in% colnames(res)) colnames(res)[1] <- "feature_id"
need <- c("feature_id", "baseMean", "log2FoldChange")
miss <- setdiff(need, colnames(res))
if (length(miss) > 0) stop("DE matrix missing columns: ", paste(miss, collapse = ", "))
res$direction <- sig_direction(res)

mode <- opt$annotation_mode
if (mode == "auto") {
  mode <- if (!is.null(opt$te_annotation_tsv) && nzchar(opt$te_annotation_tsv)) "te" else if (!is.null(opt$tx2gene_path) && nzchar(opt$tx2gene_path)) "gene" else "generic"
}
label_col <- "feature_id"
color_col <- NULL

if (mode == "gene") {
  if (is.null(opt$tx2gene_path) || !nzchar(opt$tx2gene_path)) stop("gene mode requires --tx2gene-path")
  gene_anno <- load_gene_anno_from_gtf(normalizePath(opt$tx2gene_path, mustWork = FALSE), species = opt$species)
  res <- annotate_gene_res(res, gene_anno)
  if ("gene_name_plot" %in% colnames(res)) label_col <- "gene_name_plot"
  type_cols <- intersect(c("gene_type", "gene_biotype", "transcript_type", "transcript_biotype"), colnames(res))
  if (length(type_cols) > 0) color_col <- type_cols[1]
}

if (mode == "te") {
  if (is.null(opt$te_annotation_tsv) || !nzchar(opt$te_annotation_tsv)) stop("TE mode requires --te-annotation-tsv")
  te_anno <- read_te_annotation(normalizePath(opt$te_annotation_tsv, mustWork = FALSE))
  res <- annotate_te_res(
    res, te_anno,
    preferred_match_cols = unique(c(opt$te_label_level, opt$te_color_level, "repName", "repFamily", "repClass", "locus_id", "gene_name", "feature_id")),
    label_level = opt$te_label_level,
    color_level = opt$te_color_level
  )
  res <- add_te_age(res)
  if ("te_label_plot" %in% colnames(res)) label_col <- "te_label_plot"
  color_col <- if (opt$te_color_level %in% colnames(res)) opt$te_color_level else if ("te_color_plot" %in% colnames(res)) "te_color_plot" else NULL
}

utils::write.csv(res, file.path(outdir, paste0(prefix, ".annotated_DE_matrix.csv")), row.names = FALSE)
summary_tbl <- res %>%
  dplyr::count(direction, name = "n") %>%
  dplyr::mutate(fraction = n / sum(n))
utils::write.csv(summary_tbl, file.path(outdir, paste0(prefix, ".direction_summary.csv")), row.names = FALSE)

plot_volcano_simple(res, file.path(outdir, paste0(prefix, ".volcano.pdf")), prefix, label_col = label_col, top_n = label_top_n, color_col = color_col, sig_metric = "padj")
plot_ma_simple(res, file.path(outdir, paste0(prefix, ".MA.pdf")), prefix, label_col = label_col, top_n = label_top_n, color_col = color_col, sig_metric = "padj")

if (mode == "te") {
  for (col in intersect(c("repClass", "repFamily", "repName", "TE_age"), colnames(res))) {
    plot_direction_bar(res, col, file.path(outdir, paste0(prefix, ".", col, ".up_down_bar.pdf")), paste0(prefix, " ", col, " up/down"), top_n = read_int(opt$top_n, 25))
    plot_fraction_bar(res, col, file.path(outdir, paste0(prefix, ".", col, ".up_down_fraction.pdf")), paste0(prefix, " ", col, " composition"))
  }
}

if (mode == "gene") {
  for (col in intersect(c("gene_type", "gene_biotype", "transcript_type", "transcript_biotype"), colnames(res))) {
    plot_direction_bar(res, col, file.path(outdir, paste0(prefix, ".", col, ".up_down_bar.pdf")), paste0(prefix, " ", col, " up/down"), top_n = read_int(opt$top_n, 25))
    plot_fraction_bar(res, col, file.path(outdir, paste0(prefix, ".", col, ".up_down_fraction.pdf")), paste0(prefix, " ", col, " composition"))
  }
}

writeLines(c(
  paste0("de_matrix=", de_file),
  paste0("mode=", mode),
  paste0("padj_cutoff=", padj_cutoff),
  paste0("lfc_cutoff=", lfc_cutoff),
  paste0("outputs=", outdir)
), file.path(outdir, paste0(prefix, ".visuals.summary.txt")))

message("[run_de_matrix_visuals] DONE: ", outdir)
