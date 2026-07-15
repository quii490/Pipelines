#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq DE annotation only\n\n",
    "Required:\n",
    "  --de-matrix PATH         DE matrix CSV with feature_id/baseMean/log2FoldChange\n",
    "  --outdir PATH            output directory\n",
    "  --annotation-mode STR    gene | te | generic\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --prefix STR             output prefix; default basename(de-matrix)\n",
    "  --tx2gene-path PATH      gene GTF/tx2gene for gene labels and gene_type\n",
    "  --te-annotation-tsv PATH TE annotation TSV for repName/family/class/milliDiv\n",
    "  --te-label-level STR     locus_id | repName | repFamily | repClass, default repName\n",
    "  --te-color-level STR     repFamily | repClass | TE_age, default repFamily\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  --baseMean-min NUM       default 5\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(annotation_mode = "generic", te_label_level = "repName", te_color_level = "repFamily", padj_cutoff = "0.05", lfc_cutoff = "0.58", baseMean_min = "5")
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
res <- add_direction(res)

mode <- tolower(opt$annotation_mode)
if (mode == "gene") {
  if (is.null(opt$tx2gene_path) || !nzchar(opt$tx2gene_path)) stop("gene annotation requires --tx2gene-path")
  gene_anno <- load_gene_anno_from_gtf(normalizePath(opt$tx2gene_path, mustWork = FALSE), species = opt$species)
  res <- annotate_gene_res(res, gene_anno)
} else if (mode == "te") {
  if (is.null(opt$te_annotation_tsv) || !nzchar(opt$te_annotation_tsv)) stop("TE annotation requires --te-annotation-tsv")
  te_anno <- read_te_annotation(normalizePath(opt$te_annotation_tsv, mustWork = FALSE))
  res <- annotate_te_res(
    res, te_anno,
    preferred_match_cols = unique(c(opt$te_label_level, opt$te_color_level, "repName", "repFamily", "repClass", "locus_id", "gene_name", "feature_id")),
    label_level = opt$te_label_level,
    color_level = opt$te_color_level
  )
  res <- add_te_age(res)
} else if (mode != "generic") {
  stop("--annotation-mode must be gene, te or generic")
}

utils::write.csv(res, file.path(outdir, paste0(prefix, ".annotated_DE_matrix.csv")), row.names = FALSE)
utils::write.csv(res %>% dplyr::count(direction, name = "n") %>% dplyr::mutate(fraction = n / sum(n)),
                 file.path(outdir, paste0(prefix, ".direction_summary.csv")), row.names = FALSE)
writeLines(c(
  paste0("de_matrix=", de_file),
  paste0("mode=", mode),
  paste0("features=", nrow(res))
), file.path(outdir, paste0(prefix, ".annotation.summary.txt")))
message("[run_annotate_de] DONE: ", outdir)
