#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Standalone DE plot runner for RNA-seq raw count matrix\n\n",
    "用途：从 raw count 矩阵直接做差异分析并只画你想要的图，不跑整套 pipeline。\n",
    "注意：不接受 BAM 作为输入。火山图 / MA / 热图 / PCA / GSEA 都需要 count matrix。\n\n",
    "Required:\n",
    "  --matrix PATH            raw count 矩阵，支持 CSV / RDS\n",
    "  --outdir PATH            输出目录\n",
    "  --sample-table PATH      样本表 CSV，至少包含 sample,condition,replicate\n",
    "  --contrast-file PATH     对比表 CSV，列名: case,control（兼容旧格式 group_col,case,control）\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R 路径；默认自动找 rnaseq-downstream/rnaseq-function.R\n",
    "  --matrix-format STR      auto | csv | rds，默认 auto\n",
    "  --matrix-name STR        输出前缀名；默认取 matrix 文件名\n",
    "  --species STR            hg38 | mm10 | mm39；默认 hg38\n",
    "  --threads NUM            对比并行线程数；默认 1\n",
    "  --plots STR              all | volcano,volcano_pvalue,ma,ma_pvalue,heatmap,pca,corr\n",
    "                         默认 all\n",
    "  --tx2gene-path PATH      基因 GTF / tx2gene 文件；提供后可做 gene label 注释\n",
    "  --te-annotation-tsv PATH TE 注释 TSV；提供后可做 TE label / color 注释\n",
    "  --te-label-level STR     locus_id | repName | repFamily | repClass；默认 repName\n",
    "  --te-color-level STR     repFamily | repClass；默认 repFamily\n",
    "  --padj-cutoff NUM        默认 0.05\n",
    "  --lfc-cutoff NUM         默认 1\n",
    "  --baseMean-min NUM       默认 5\n",
    "  --label-top-n NUM        默认 40\n",
    "  --heatmap-top-n NUM      默认 40\n",
    "  --volcano-orientation STR classic | horizontal，默认 classic\n",
    "  --gray-nonsig BOOL       true|false；非显著点是否统一灰色，默认 true\n",
    "  --exploratory-method STR logCPM_diff | edgeR_fixedBCV；无重复差异算法，默认 logCPM_diff\n",
    "  --exploratory-fixed-bcv NUM edgeR_fixedBCV 的 BCV，默认 0.4\n",
    "  -h, --help               显示帮助\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(
    matrix_format = "auto",
    species = "hg38",
    threads = "1",
    plots = "all",
    te_label_level = "repName",
    te_color_level = "repFamily",
    padj_cutoff = "0.05",
    lfc_cutoff = "1",
    baseMean_min = "5",
    label_top_n = "40",
    heatmap_top_n = "40",
    volcano_orientation = "classic",
    gray_nonsig = "true",
    exploratory_method = "logCPM_diff",
    exploratory_fixed_bcv = "0.4"
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("-h", "--help")) {
      usage()
      quit(save = "no", status = 0)
    }
    if (!startsWith(key, "--")) stop("未知参数: ", key)
    if (i == length(args)) stop("参数缺少值: ", key)
    val <- args[[i + 1]]
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    res[[name]] <- val
    i <- i + 2
  }
  res
}

read_int <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(as.integer(default))
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v < 0) return(as.integer(default))
  v
}

read_num <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) return(default)
  v
}

read_bool <- function(x, default = TRUE) {
  if (is.null(x) || !nzchar(x)) return(default)
  y <- tolower(trimws(as.character(x)))
  if (y %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (y %in% c("false", "f", "no", "n", "0")) return(FALSE)
  default
}

log_msg <- function(...) cat("[run_plot_from_counts]", ..., "\n")

opt <- parse_args(args)
required_keys <- c("matrix", "outdir", "sample_table", "contrast_file")
missing_keys <- required_keys[!vapply(required_keys, function(k) !is.null(opt[[k]]) && nzchar(opt[[k]]), logical(1))]
if (length(missing_keys) > 0) {
  usage()
  stop("缺少必需参数: ", paste(missing_keys, collapse = ", "))
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_candidates <- c(
  if (!is.null(opt$function_file) && nzchar(opt$function_file)) normalizePath(opt$function_file, mustWork = FALSE) else character(0),
  normalizePath(file.path(script_dir, "..", "rnaseq-downstream", "rnaseq-function.R"), mustWork = FALSE),
  file.path(script_dir, "rnaseq-function.R"),
  "/mnt/data/rnaseq-function.R"
)
function_candidates <- unique(function_candidates[!is.na(function_candidates) & nzchar(function_candidates)])
function_file <- function_candidates[file.exists(function_candidates)][1]
if (is.na(function_file) || !nzchar(function_file)) stop("未找到函数文件: ", paste(function_candidates, collapse = " | "))

matrix_file <- normalizePath(opt$matrix, mustWork = FALSE)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
sample_table_file <- normalizePath(opt$sample_table, mustWork = FALSE)
contrast_file <- normalizePath(opt$contrast_file, mustWork = FALSE)
matrix_format <- tolower(opt$matrix_format)
species <- opt$species
threads <- read_int(opt$threads, 1)
padj_cutoff <- read_num(opt$padj_cutoff, 0.05)
lfc_cutoff <- read_num(opt$lfc_cutoff, 1)
baseMean_min <- read_num(opt$baseMean_min, 5)
label_top_n <- read_int(opt$label_top_n, 40)
heatmap_top_n <- read_int(opt$heatmap_top_n, 40)
volcano_orientation <- tolower(as.character(opt$volcano_orientation))
if (!volcano_orientation %in% c("classic", "horizontal")) stop("--volcano-orientation 必须是 classic 或 horizontal")
gray_nonsig <- read_bool(opt$gray_nonsig, TRUE)
exploratory_method <- as.character(opt$exploratory_method)
if (!exploratory_method %in% c("logCPM_diff", "edgeR_fixedBCV")) stop("--exploratory-method 必须是 logCPM_diff 或 edgeR_fixedBCV")
exploratory_fixed_bcv <- read_num(opt$exploratory_fixed_bcv, 0.4)
if (is.na(exploratory_fixed_bcv) || exploratory_fixed_bcv <= 0) stop("--exploratory-fixed-bcv 必须是正数")
plots_raw <- trimws(unlist(strsplit(as.character(opt$plots), ",", fixed = TRUE)))
plots_raw <- plots_raw[nzchar(plots_raw)]
if (length(plots_raw) == 0) plots_raw <- "all"
plots_req <- unique(tolower(plots_raw))
valid_plots <- c("all", "volcano", "volcano_pvalue", "ma", "ma_pvalue", "heatmap", "pca", "corr")
invalid_plots <- setdiff(plots_req, valid_plots)
if (length(invalid_plots) > 0) stop("不支持的 --plots: ", paste(invalid_plots, collapse = ", "))
if ("all" %in% plots_req) plots_req <- setdiff(valid_plots, "all")

for (p in c(matrix_file, sample_table_file, contrast_file)) {
  if (!file.exists(p)) stop("文件不存在: ", p)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 先设置全局 cutoffs，再 source 函数文件
options(plot_n_cores = 1L)
options(deseq_n_cores = 1L)
assign("padj_cutoff", padj_cutoff, envir = .GlobalEnv)
assign("lfc_cutoff", lfc_cutoff, envir = .GlobalEnv)
assign("baseMean_min", baseMean_min, envir = .GlobalEnv)
assign("label_top_n", label_top_n, envir = .GlobalEnv)
assign("heatmap_top_n", heatmap_top_n, envir = .GlobalEnv)
assign("volcano_orientation", volcano_orientation, envir = .GlobalEnv)
assign("gray_nonsig", gray_nonsig, envir = .GlobalEnv)
assign("exploratory_method", exploratory_method, envir = .GlobalEnv)
assign("exploratory_fixed_bcv", exploratory_fixed_bcv, envir = .GlobalEnv)

source(function_file)

read_input_matrix <- function(path, fmt = "auto") {
  if (fmt == "auto") {
    ext <- tolower(tools::file_ext(path))
    fmt <- if (ext == "rds") "rds" else "csv"
  }

  if (fmt == "rds") {
    mat <- extract_rds_matrix(path)
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    return(mat)
  }

  dt <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dt) < 3) stop("输入矩阵至少需要 1 列特征 + 2 列样本: ", basename(path))
  feature_id <- as.character(dt[[1]])
  mat <- as.matrix(dt[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  rownames(mat) <- feature_id
  if (anyNA(mat)) {
    log_msg("matrix contains NA after numeric conversion; set NA to 0")
    mat[is.na(mat)] <- 0
  }
  mat
}

sample_table <- utils::read.csv(sample_table_file, stringsAsFactors = FALSE, check.names = FALSE)
need_sample_cols <- c("sample", "condition", "replicate")
miss_sample_cols <- setdiff(need_sample_cols, colnames(sample_table))
if (length(miss_sample_cols) > 0) stop("sample_table 缺少列: ", paste(miss_sample_cols, collapse = ", "))
sample_table <- unique(sample_table[, need_sample_cols, drop = FALSE])
sample_table$sample <- as.character(sample_table$sample)
sample_table$condition <- as.character(sample_table$condition)
sample_table$replicate <- as.character(sample_table$replicate)

contrast_df <- utils::read.csv(contrast_file, stringsAsFactors = FALSE, check.names = FALSE)
if (all(c("case", "control") %in% colnames(contrast_df))) {
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) c("condition", as.character(contrast_df$case[i]), as.character(contrast_df$control[i])))
} else if (all(c("group_col", "case", "control") %in% colnames(contrast_df))) {
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) c(as.character(contrast_df$group_col[i]), as.character(contrast_df$case[i]), as.character(contrast_df$control[i])))
} else {
  stop("contrast_file 缺少列；需要 case,control，或兼容旧格式 group_col,case,control")
}
if (length(contrast_list) == 0) stop("contrast_file 为空")

matrix_name <- if (!is.null(opt$matrix_name) && nzchar(opt$matrix_name)) {
  opt$matrix_name
} else {
  sub("\\.[^.]+$", "", basename(matrix_file))
}
matrix_name <- sanitize_name(matrix_name)
tool_outdir <- file.path(outdir, matrix_name)
dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)
# safe_run 使用全局 outdir
outdir <- tool_outdir
assign("outdir", outdir, envir = .GlobalEnv)

count_mat <- read_input_matrix(matrix_file, fmt = matrix_format)
log_msg("function_file =", function_file)
log_msg("matrix =", matrix_file)
log_msg("outdir =", outdir)
log_msg("matrix dim =", paste(dim(count_mat), collapse = " x "))
log_msg("species =", species)
log_msg("plots =", paste(plots_req, collapse = ","))
log_msg("threads =", threads)

# 注释加载（可选）
gene_anno <- NULL
te_anno <- NULL
annotation_mode <- "generic"
if (!is.null(opt$tx2gene_path) && nzchar(opt$tx2gene_path)) {
  tx2gene_path <- normalizePath(opt$tx2gene_path, mustWork = FALSE)
  if (!file.exists(tx2gene_path)) stop("tx2gene_path 不存在: ", tx2gene_path)
  gene_anno <- load_gene_anno_from_gtf(tx2gene_path, species = species)
  annotation_mode <- "gene"
  log_msg("gene annotation rows =", nrow(gene_anno))
}
if (!is.null(opt$te_annotation_tsv) && nzchar(opt$te_annotation_tsv)) {
  te_path <- normalizePath(opt$te_annotation_tsv, mustWork = FALSE)
  if (!file.exists(te_path)) stop("te_annotation_tsv 不存在: ", te_path)
  te_anno <- read_te_annotation(te_path)
  annotation_mode <- "te"
  log_msg("TE annotation rows =", nrow(te_anno))
}
if (!is.null(gene_anno) && !is.null(te_anno)) stop("tx2gene-path 和 te-annotation-tsv 只能二选一")

make_title <- function(prefix, ct) paste0(prefix, "_", ct[2], "_vs_", ct[3])

render_one_contrast <- function(ct) {
  tag <- sanitize_name(make_title(matrix_name, ct))
  cdir <- file.path(tool_outdir, tag)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  log_msg("contrast start =", paste(ct, collapse = " | "))

  fit <- run_deseq_simple(count_mat, sample_table, ct)
  res <- fit$res
  this_label_col <- "feature_id"
  this_color_col <- NULL
  this_heatmap_row_anno_col <- NULL

  if (!is.null(gene_anno) && all(c("gene_id", "gene_name") %in% colnames(gene_anno))) {
    res <- annotate_gene_res(res, gene_anno)
    if ("gene_name_plot" %in% colnames(res)) this_label_col <- "gene_name_plot"
  } else if (!is.null(te_anno)) {
    pref_cols <- unique(na.omit(c(
      opt$te_label_level,
      opt$te_color_level,
      "repName", "repFamily", "repClass", "locus_id", "gene_name", "feature_id"
    )))
    res <- annotate_te_res(
      res,
      te_anno,
      preferred_match_cols = pref_cols,
      label_level = opt$te_label_level,
      color_level = opt$te_color_level
    )
    if ("te_label_plot" %in% colnames(res)) this_label_col <- "te_label_plot"
    if ("te_color_plot" %in% colnames(res)) this_color_col <- "te_color_plot"
    if ("te_heatmap_group" %in% colnames(res)) this_heatmap_row_anno_col <- "te_heatmap_group"
  }

  res_export <- clean_res_for_export(res, tag = tag)
  utils::write.csv(res_export, file.path(cdir, paste0(tag, "_DE.csv")), row.names = FALSE)

  if ("volcano" %in% plots_req) {
    plot_volcano_simple(
      res_export,
      file.path(cdir, paste0(tag, "_volcano.png")),
      tag,
      label_col = this_label_col,
      top_n = label_top_n,
      color_col = this_color_col,
      sig_metric = "padj"
    )
  }
  if ("volcano_pvalue" %in% plots_req) {
    plot_volcano_simple(
      res_export,
      file.path(cdir, paste0(tag, "_volcano_pvalue.png")),
      tag,
      label_col = this_label_col,
      top_n = label_top_n,
      color_col = this_color_col,
      sig_metric = "pvalue"
    )
  }
  if ("ma" %in% plots_req) {
    plot_ma_simple(
      res_export,
      file.path(cdir, paste0(tag, "_MA.png")),
      tag,
      label_col = this_label_col,
      top_n = label_top_n,
      color_col = this_color_col,
      sig_metric = "padj"
    )
  }
  if ("ma_pvalue" %in% plots_req) {
    plot_ma_simple(
      res_export,
      file.path(cdir, paste0(tag, "_MA_pvalue.png")),
      tag,
      label_col = this_label_col,
      top_n = label_top_n,
      color_col = this_color_col,
      sig_metric = "pvalue"
    )
  }
  if ("heatmap" %in% plots_req) {
    plot_top_heatmap_simple(
      fit$vst_mat,
      res_export,
      fit$sample_table,
      file.path(cdir, paste0(tag, "_heatmap.png")),
      tag,
      top_n = heatmap_top_n,
      label_col = this_label_col,
      row_anno_col = this_heatmap_row_anno_col
    )
  }

  list(tag = tag, fit = fit)
}

# 全局 QC 只需要画一次；如果要求 pca/corr，则用第一个 contrast 的 fit
if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && threads > 1) {
  log_msg("parallel mode = contrast-level; workers =", min(threads, length(contrast_list)))
  fits <- parallel::mclapply(
    contrast_list,
    render_one_contrast,
    mc.cores = min(threads, length(contrast_list)),
    mc.preschedule = FALSE
  )
} else {
  log_msg("parallel mode = serial")
  fits <- lapply(contrast_list, render_one_contrast)
}

first_fit <- NULL
for (x in fits) {
  if (!is.null(x$fit)) {
    first_fit <- x$fit
    break
  }
}

if (!is.null(first_fit)) {
  qc_dir <- file.path(tool_outdir, "global_qc")
  dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
  if ("pca" %in% plots_req) {
    plot_pca_simple(first_fit$vst_mat, first_fit$sample_table, file.path(qc_dir, paste0(matrix_name, "_PCA.png")), matrix_name)
  }
  if ("corr" %in% plots_req) {
    plot_corr_simple(first_fit$vst_mat, file.path(qc_dir, paste0(matrix_name, "_correlation_heatmap.png")), matrix_name)
  }
}

summary_lines <- c(
  paste0("matrix=", matrix_file),
  paste0("matrix_name=", matrix_name),
  paste0("species=", species),
  paste0("annotation_mode=", annotation_mode),
  paste0("plots=", paste(plots_req, collapse = ",")),
  paste0("contrast_count=", length(contrast_list)),
  paste0("threads=", threads)
)
writeLines(summary_lines, file.path(tool_outdir, "run_plot_from_counts.summary.txt"))
log_msg("done")
