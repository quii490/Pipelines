#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "Standalone GSEA runner for RNA-seq gene matrix\n\n",
    "Required:\n",
    "  --matrix PATH            基因计数矩阵，支持 CSV / RDS\n",
    "  --outdir PATH            输出目录\n",
    "  --species STR            hg38 | mm10 | mm39\n",
    "  --tx2gene-path PATH      基因 GTF / tx2gene 文件\n",
    "  --sample-table PATH      样本表 CSV，至少包含 sample,condition,replicate\n",
    "  --contrast-file PATH     对比表 CSV，列名: case,control（兼容旧格式 group_col,case,control）\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R 路径；默认自动找 rnaseq-downstream/rnaseq-function.R\n",
    "  --matrix-format STR      auto | csv | rds，默认 auto\n",
    "  --matrix-name STR        输出前缀名；默认取 matrix 文件名\n",
    "  --threads NUM            并行线程数；默认 1\n",
    "  --topn-plot NUM          GSEA 图默认 topN；默认 20\n",
    "  --max-terms-gseaplot2 N  gseaplot2 最多通路数；默认 12\n",
    "  --min-gs-size NUM        默认 10\n",
    "  --max-gs-size NUM        默认 500\n",
    "  --pvalue-cutoff NUM      默认 1\n",
    "  --disable-gseaplot2 true|false   默认 false\n",
    "  -h, --help               显示帮助\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(
    matrix_format = "auto",
    threads = "1",
    topn_plot = "20",
    max_terms_gseaplot2 = "12",
    min_gs_size = "10",
    max_gs_size = "500",
    pvalue_cutoff = "1",
    disable_gseaplot2 = "false"
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

read_bool <- function(x, default = FALSE) {
  if (is.null(x) || !nzchar(x)) return(default)
  tolower(as.character(x)) %in% c("1", "true", "t", "yes", "y")
}

read_int <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(as.integer(default))
  v <- suppressWarnings(as.integer(x))
  if (is.na(v) || v < 1) return(as.integer(default))
  v
}

log_msg <- function(...) cat("[run_gsea_standalone]", ..., "\n")

opt <- parse_args(args)
required_keys <- c("matrix", "outdir", "species", "tx2gene_path", "sample_table", "contrast_file")
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
species <- opt$species
tx2gene_path <- normalizePath(opt$tx2gene_path, mustWork = FALSE)
sample_table_file <- normalizePath(opt$sample_table, mustWork = FALSE)
contrast_file <- normalizePath(opt$contrast_file, mustWork = FALSE)
matrix_format <- tolower(opt$matrix_format)
threads <- read_int(opt$threads, 1)
topn_plot <- read_int(opt$topn_plot, 20)
max_terms_gseaplot2 <- read_int(opt$max_terms_gseaplot2, 12)
min_gs_size <- read_int(opt$min_gs_size, 10)
max_gs_size <- read_int(opt$max_gs_size, 500)
pvalue_cutoff <- suppressWarnings(as.numeric(opt$pvalue_cutoff))
if (is.na(pvalue_cutoff)) pvalue_cutoff <- 1
make_gseaplot2 <- !read_bool(opt$disable_gseaplot2, default = FALSE)

for (p in c(matrix_file, tx2gene_path, sample_table_file, contrast_file)) {
  if (!file.exists(p)) stop("文件不存在: ", p)
}
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

source(function_file)

matrix_name <- if (!is.null(opt$matrix_name) && nzchar(opt$matrix_name)) {
  opt$matrix_name
} else {
  sub("\\.[^.]+$", "", basename(matrix_file))
}
matrix_name <- sanitize_name(matrix_name)
tool_outdir <- file.path(outdir, matrix_name)
dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)

if (threads > 1) {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
}
options(plot_n_cores = 1L)
options(deseq_n_cores = 1L)

log_msg("function_file =", function_file)
log_msg("matrix =", matrix_file)
log_msg("outdir =", outdir)
log_msg("tool_outdir =", tool_outdir)
log_msg("species =", species)
log_msg("tx2gene_path =", tx2gene_path)
log_msg("sample_table =", sample_table_file)
log_msg("contrast_file =", contrast_file)
log_msg("threads =", threads)
log_msg("matrix_format =", matrix_format)
log_msg("make_gseaplot2 =", make_gseaplot2)

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
  log_msg("contrast schema = case/control")
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) {
    c("condition", as.character(contrast_df$case[i]), as.character(contrast_df$control[i]))
  })
} else if (all(c("group_col", "case", "control") %in% colnames(contrast_df))) {
  log_msg("contrast schema = legacy group_col/case/control")
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) {
    c(as.character(contrast_df$group_col[i]), as.character(contrast_df$case[i]), as.character(contrast_df$control[i]))
  })
} else {
  stop("contrast_file 缺少列；需要 case,control，或兼容旧格式 group_col,case,control")
}
if (length(contrast_list) == 0) stop("contrast_file 为空")

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

count_mat <- read_input_matrix(matrix_file, fmt = matrix_format)
log_msg("matrix dim =", paste(dim(count_mat), collapse = " x "))

safe_run("gene_annotation_load", {
  gene_anno <- load_gene_anno_from_gtf(tx2gene_path, species = species)
  saveRDS(gene_anno, file.path(tool_outdir, "gene_annotation.rds"))
  utils::write.csv(gene_anno, file.path(tool_outdir, "gene_annotation.preview.csv"), row.names = FALSE)
  log_msg("gene annotation rows =", nrow(gene_anno))
})
if (!exists("gene_anno") || is.null(gene_anno)) stop("gene_anno 加载失败")

worker_fun <- function(ct) {
  tag <- sanitize_name(paste(matrix_name, ct[2], "vs", ct[3], sep = "_"))
  safe_run(paste0(tag, "_gsea"), {
    fit <- run_deseq_simple(count_mat, sample_table, ct)
    res <- fit$res
    if (!is.null(gene_anno) && all(c("gene_id", "gene_name") %in% colnames(gene_anno))) {
      res <- annotate_gene_res(res, gene_anno)
    }
    res_export <- clean_res_for_export(res, tag = tag)
    write.csv(res_export, file.path(tool_outdir, paste0(tag, "_DE.csv")), row.names = FALSE)

    gsea_df <- run_gsea_multi_category(
      fit = fit,
      res_df = res_export,
      contrast = ct,
      prefix = file.path(tool_outdir, tag),
      species = species,
      gene_anno = gene_anno,
      gsea_categories = get_default_gsea_categories(),
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      pvalueCutoff = pvalue_cutoff,
      topN_plot = topn_plot,
      max_terms_gseaplot2 = max_terms_gseaplot2,
      make_gseaplot2 = make_gseaplot2
    )

    list(tag = tag, gsea_df = gsea_df)
  })
}

if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && threads > 1) {
  log_msg("parallel mode = contrast-level; workers =", min(threads, length(contrast_list)))
  res_list <- parallel::mclapply(
    contrast_list,
    worker_fun,
    mc.cores = min(threads, length(contrast_list)),
    mc.preschedule = FALSE
  )
} else {
  log_msg("parallel mode = serial")
  if (threads > 1 && length(contrast_list) == 1) {
    options(deseq_n_cores = min(4L, as.integer(threads)))
    log_msg("single contrast; DESeq2 workers =", getOption("deseq_n_cores", 1L))
  }
  res_list <- lapply(contrast_list, worker_fun)
}

gsea_collect <- lapply(res_list, function(x) {
  if (is.null(x) || is.null(x$gsea_df) || !is.data.frame(x$gsea_df) || nrow(x$gsea_df) == 0) return(NULL)
  x$gsea_df
})
gsea_collect <- gsea_collect[!vapply(gsea_collect, is.null, logical(1))]

if (length(gsea_collect) > 0) {
  gsea_all <- dplyr::bind_rows(gsea_collect)
  write.csv(gsea_all, file.path(tool_outdir, paste0(matrix_name, "_GSEA_all_contrasts_all_categories.csv")), row.names = FALSE)

  tryCatch({
    plot_gsea_across_contrast_heatmap_simple(
      gsea_all,
      file.path(tool_outdir, paste0(matrix_name, "_Hallmark_GSEA_across_contrasts_heatmap.png")),
      category_label = "Hallmark",
      top_n = 25
    )
  }, error = function(e) {
    message("[", matrix_name, "_Hallmark_GSEA_across_contrasts_heatmap] ERROR: ", conditionMessage(e))
  })

  log_msg("all contrasts GSEA rows =", nrow(gsea_all))
} else {
  log_msg("no valid GSEA result generated")
}

writeLines(c(
  paste0("matrix=", matrix_file),
  paste0("outdir=", outdir),
  paste0("tool_outdir=", tool_outdir),
  paste0("species=", species),
  paste0("threads=", threads),
  paste0("contrast_count=", length(contrast_list))
), file.path(tool_outdir, "run_gsea_standalone.summary.txt"))

log_msg("DONE")
