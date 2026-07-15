#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq differential analysis only: raw counts -> DE matrices\n\n",
    "Required:\n",
    "  --matrix PATH            raw count matrix CSV/RDS; first column is feature ID\n",
    "  --outdir PATH            output directory\n",
    "  --sample-table PATH      CSV with sample,condition,replicate\n",
    "  --contrast-file PATH     CSV with case,control or group_col,case,control\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --matrix-format STR      auto | csv | rds, default auto\n",
    "  --matrix-name STR        output prefix; default basename(matrix)\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  --baseMean-min NUM       default 5\n",
    "  --exploratory-method STR logCPM_diff | edgeR_fixedBCV; no-replicate method, default logCPM_diff\n",
    "  --exploratory-fixed-bcv NUM edgeR_fixedBCV BCV, default 0.4\n",
    "  --threads NUM            contrast-level workers, default 1\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(matrix_format = "auto", padj_cutoff = "0.05", lfc_cutoff = "0.58", baseMean_min = "5", exploratory_method = "logCPM_diff", exploratory_fixed_bcv = "0.4", threads = "1")
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

read_int <- function(x, default = 1L) {
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

read_input_matrix <- function(path, fmt = "auto") {
  if (fmt == "auto") fmt <- if (tolower(tools::file_ext(path)) == "rds") "rds" else "csv"
  if (fmt == "rds") {
    mat <- extract_rds_matrix(path)
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    return(mat)
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 3) stop("matrix must contain feature ID plus at least two samples")
  mat <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[is.na(mat)] <- 0
  rownames(mat) <- as.character(df[[1]])
  mat
}

read_sample_table <- function(path) {
  st <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("sample", "condition", "replicate")
  miss <- setdiff(need, colnames(st))
  if (length(miss) > 0) stop("sample_table missing columns: ", paste(miss, collapse = ", "))
  st <- unique(st[, need, drop = FALSE])
  st$sample <- as.character(st$sample)
  st$condition <- as.character(st$condition)
  st$replicate <- as.character(st$replicate)
  st
}

read_contrasts <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (all(c("case", "control") %in% colnames(df))) {
    lapply(seq_len(nrow(df)), function(i) c("condition", as.character(df$case[i]), as.character(df$control[i])))
  } else if (all(c("group_col", "case", "control") %in% colnames(df))) {
    lapply(seq_len(nrow(df)), function(i) c(as.character(df$group_col[i]), as.character(df$case[i]), as.character(df$control[i])))
  } else {
    stop("contrast_file needs case,control or group_col,case,control")
  }
}

opt <- parse_args(args)
required <- c("matrix", "outdir", "sample_table", "contrast_file")
missing <- required[!vapply(required, function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_file <- find_function_file(script_dir, opt)
matrix_file <- normalizePath(opt$matrix, mustWork = FALSE)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
sample_table_file <- normalizePath(opt$sample_table, mustWork = FALSE)
contrast_file <- normalizePath(opt$contrast_file, mustWork = FALSE)
for (p in c(matrix_file, sample_table_file, contrast_file)) if (!file.exists(p)) stop("File not found: ", p)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- read_num(opt$padj_cutoff, 0.05)
lfc_cutoff <- read_num(opt$lfc_cutoff, 0.58)
baseMean_min <- read_num(opt$baseMean_min, 5)
label_top_n <- 40L
heatmap_top_n <- 40L
plot_n_cores <- read_int(opt$threads, 1L)
exploratory_method <- as.character(opt$exploratory_method)
if (!exploratory_method %in% c("logCPM_diff", "edgeR_fixedBCV")) stop("--exploratory-method must be logCPM_diff or edgeR_fixedBCV")
exploratory_fixed_bcv <- read_num(opt$exploratory_fixed_bcv, 0.4)
if (is.na(exploratory_fixed_bcv) || exploratory_fixed_bcv <= 0) stop("--exploratory-fixed-bcv must be positive")
assign("outdir", outdir, envir = .GlobalEnv)
assign("padj_cutoff", padj_cutoff, envir = .GlobalEnv)
assign("lfc_cutoff", lfc_cutoff, envir = .GlobalEnv)
assign("baseMean_min", baseMean_min, envir = .GlobalEnv)
assign("label_top_n", label_top_n, envir = .GlobalEnv)
assign("heatmap_top_n", heatmap_top_n, envir = .GlobalEnv)
assign("exploratory_method", exploratory_method, envir = .GlobalEnv)
assign("exploratory_fixed_bcv", exploratory_fixed_bcv, envir = .GlobalEnv)
options(plot_n_cores = plot_n_cores)

source(function_file)

matrix_name <- if (!is.null(opt$matrix_name) && nzchar(opt$matrix_name)) opt$matrix_name else sub("\\.[^.]+$", "", basename(matrix_file))
matrix_name <- sanitize_name(matrix_name)
tool_outdir <- file.path(outdir, matrix_name)
dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)

count_mat <- read_input_matrix(matrix_file, opt$matrix_format)
sample_table <- read_sample_table(sample_table_file)
contrast_list <- read_contrasts(contrast_file)
x <- ensure_sample_order(count_mat, sample_table)
count_mat <- x$mat
sample_table <- x$sample_table
write_matrix_csv(count_mat, file.path(tool_outdir, paste0(matrix_name, ".raw_count_matrix.csv")))

run_one <- function(ct) {
  tag <- sanitize_name(paste(matrix_name, ct[2], "vs", ct[3], sep = "_"))
  cdir <- file.path(tool_outdir, tag)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  fit <- run_deseq_simple(count_mat, sample_table, ct)
  res <- clean_res_for_export(fit$res, tag = tag)
  utils::write.csv(res, file.path(cdir, paste0(tag, ".DE_matrix.csv")), row.names = FALSE)
  if (!is.null(fit$dds)) {
    norm_counts <- DESeq2::counts(fit$dds, normalized = TRUE)
    utils::write.csv(data.frame(feature_id = rownames(norm_counts), norm_counts, check.names = FALSE), file.path(cdir, paste0(tag, ".normalized_counts.csv")), row.names = FALSE)
  }
  utils::write.csv(data.frame(feature_id = rownames(fit$vst_mat), fit$vst_mat, check.names = FALSE), file.path(cdir, paste0(tag, ".vst_or_log_matrix.csv")), row.names = FALSE)
  writeLines(c(
    paste0("tag=", tag),
    paste0("mode=", if (!is.null(fit$mode)) fit$mode else "unknown"),
    paste0("features=", nrow(res))
  ), file.path(cdir, paste0(tag, ".diff.summary.txt")))
}

if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && plot_n_cores > 1) {
  parallel::mclapply(contrast_list, run_one, mc.cores = min(plot_n_cores, length(contrast_list)), mc.preschedule = FALSE)
} else {
  lapply(contrast_list, run_one)
}

writeLines(c(
  paste0("matrix=", matrix_file),
  paste0("matrix_name=", matrix_name),
  paste0("contrast_count=", length(contrast_list)),
  paste0("padj_cutoff=", padj_cutoff),
  paste0("lfc_cutoff=", lfc_cutoff)
), file.path(tool_outdir, paste0(matrix_name, ".diff.summary.txt")))

message("[run_diff_from_counts] DONE: ", tool_outdir)
