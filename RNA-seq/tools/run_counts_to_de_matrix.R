#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq raw counts -> reusable differential matrix\n\n",
    "Required:\n",
    "  --matrix PATH            raw count matrix CSV/RDS; first column is feature ID\n",
    "  --outdir PATH            output directory\n",
    "  --sample-table PATH      CSV with sample,condition,replicate\n",
    "  --contrast-file PATH     CSV with case,control or group_col,case,control\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --matrix-format STR      auto | csv | rds, default auto\n",
    "  --matrix-name STR        output prefix; default basename(matrix)\n",
    "  --species STR            hg38 | mm10 | mm39, default hg38\n",
    "  --tx2gene-path PATH      gene GTF/tx2gene for gene labels and GO/GSEA\n",
    "  --te-annotation-tsv PATH TE annotation TSV for TE labels/classes/age\n",
    "  --te-label-level STR     locus_id | repName | repFamily | repClass, default repName\n",
    "  --te-color-level STR     repFamily | repClass, default repFamily\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  --baseMean-min NUM       default 5\n",
    "  --label-top-n NUM        labels for volcano/MA, default 40; 0 disables labels\n",
    "  --volcano-orientation STR classic | horizontal, default classic\n",
    "  --gray-nonsig BOOL       true|false; color non-significant points gray, default true\n",
    "  --exploratory-method STR logCPM_diff | edgeR_fixedBCV; no-replicate method, default logCPM_diff\n",
    "  --exploratory-fixed-bcv NUM edgeR_fixedBCV BCV, default 0.4\n",
    "  --threads NUM            contrast-level parallel workers, default 1\n",
    "  --make-plots BOOL        true|false, default true\n",
    "  --run-go BOOL            true|false, default false; gene mode only\n",
    "  --run-gsea BOOL          true|false, default false; gene mode only\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(
    matrix_format = "auto",
    species = "hg38",
    te_label_level = "repName",
    te_color_level = "repFamily",
    padj_cutoff = "0.05",
    lfc_cutoff = "0.58",
    baseMean_min = "5",
    label_top_n = "40",
    volcano_orientation = "classic",
    gray_nonsig = "true",
    exploratory_method = "logCPM_diff",
    exploratory_fixed_bcv = "0.4",
    threads = "1",
    make_plots = "true",
    run_go = "false",
    run_gsea = "false"
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
    val <- args[[i + 1]]
    name <- gsub("-", "_", sub("^--", "", key))
    res[[name]] <- val
    i <- i + 2
  }
  res
}

as_bool <- function(x) {
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

read_int <- function(x, default = 1L) {
  y <- suppressWarnings(as.integer(x))
  if (is.na(y) || y < 0) default else y
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

read_input_matrix <- function(path, fmt = "auto") {
  if (fmt == "auto") {
    fmt <- if (tolower(tools::file_ext(path)) == "rds") "rds" else "csv"
  }
  if (fmt == "rds") {
    mat <- extract_rds_matrix(path)
    mat <- as.matrix(mat)
    storage.mode(mat) <- "numeric"
    return(mat)
  }
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 3) stop("matrix must contain feature ID plus at least two samples")
  ids <- as.character(df[[1]])
  mat <- as.matrix(df[, -1, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  mat[is.na(mat)] <- 0
  rownames(mat) <- ids
  mat
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
label_top_n <- read_int(opt$label_top_n, 40L)
heatmap_top_n <- 40L
plot_n_cores <- read_int(opt$threads, 1L)
volcano_orientation <- tolower(as.character(opt$volcano_orientation))
if (!volcano_orientation %in% c("classic", "horizontal")) stop("--volcano-orientation must be classic or horizontal")
gray_nonsig <- as_bool(opt$gray_nonsig)
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
assign("volcano_orientation", volcano_orientation, envir = .GlobalEnv)
assign("gray_nonsig", gray_nonsig, envir = .GlobalEnv)
assign("exploratory_method", exploratory_method, envir = .GlobalEnv)
assign("exploratory_fixed_bcv", exploratory_fixed_bcv, envir = .GlobalEnv)
options(plot_n_cores = plot_n_cores)

source(function_file)

matrix_name <- if (!is.null(opt$matrix_name) && nzchar(opt$matrix_name)) opt$matrix_name else sub("\\.[^.]+$", "", basename(matrix_file))
matrix_name <- sanitize_name(matrix_name)
tool_outdir <- file.path(outdir, matrix_name)
dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)

sample_table <- utils::read.csv(sample_table_file, stringsAsFactors = FALSE, check.names = FALSE)
need_cols <- c("sample", "condition", "replicate")
miss_cols <- setdiff(need_cols, colnames(sample_table))
if (length(miss_cols) > 0) stop("sample_table missing columns: ", paste(miss_cols, collapse = ", "))
sample_table <- unique(sample_table[, need_cols, drop = FALSE])
sample_table$sample <- as.character(sample_table$sample)
sample_table$condition <- as.character(sample_table$condition)
sample_table$replicate <- as.character(sample_table$replicate)

contrast_df <- utils::read.csv(contrast_file, stringsAsFactors = FALSE, check.names = FALSE)
if (all(c("case", "control") %in% colnames(contrast_df))) {
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) c("condition", as.character(contrast_df$case[i]), as.character(contrast_df$control[i])))
} else if (all(c("group_col", "case", "control") %in% colnames(contrast_df))) {
  contrast_list <- lapply(seq_len(nrow(contrast_df)), function(i) c(as.character(contrast_df$group_col[i]), as.character(contrast_df$case[i]), as.character(contrast_df$control[i])))
} else {
  stop("contrast_file needs case,control or group_col,case,control")
}

count_mat <- read_input_matrix(matrix_file, opt$matrix_format)
x_all <- ensure_sample_order(count_mat, sample_table)
count_mat <- x_all$mat
sample_table <- x_all$sample_table
write_matrix_csv(count_mat, file.path(tool_outdir, paste0(matrix_name, ".raw_count_matrix.csv")))

gene_anno <- NULL
te_anno <- NULL
mode <- "generic"
if (!is.null(opt$tx2gene_path) && nzchar(opt$tx2gene_path)) {
  gene_anno <- load_gene_anno_from_gtf(normalizePath(opt$tx2gene_path, mustWork = FALSE), species = opt$species)
  mode <- "gene"
}
if (!is.null(opt$te_annotation_tsv) && nzchar(opt$te_annotation_tsv)) {
  if (!is.null(gene_anno)) stop("--tx2gene-path and --te-annotation-tsv are mutually exclusive")
  te_anno <- read_te_annotation(normalizePath(opt$te_annotation_tsv, mustWork = FALSE))
  mode <- "te"
}

run_one <- function(ct) {
  tag <- sanitize_name(paste(matrix_name, ct[2], "vs", ct[3], sep = "_"))
  cdir <- file.path(tool_outdir, tag)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  fit <- run_deseq_simple(count_mat, sample_table, ct)
  res <- clean_res_for_export(fit$res, tag = tag)
  label_col <- "feature_id"
  color_col <- NULL
  row_anno_col <- NULL
  if (!is.null(gene_anno)) {
    res <- annotate_gene_res(res, gene_anno)
    label_col <- if ("gene_name_plot" %in% colnames(res)) "gene_name_plot" else "feature_id"
  }
  if (!is.null(te_anno)) {
    res <- annotate_te_res(
      res, te_anno,
      preferred_match_cols = unique(c(opt$te_label_level, opt$te_color_level, "repName", "repFamily", "repClass", "locus_id", "gene_name", "feature_id")),
      label_level = opt$te_label_level,
      color_level = opt$te_color_level
    )
    label_col <- if ("te_label_plot" %in% colnames(res)) "te_label_plot" else "feature_id"
    color_col <- if ("te_color_plot" %in% colnames(res)) "te_color_plot" else NULL
    row_anno_col <- if ("te_heatmap_group" %in% colnames(res)) "te_heatmap_group" else NULL
  }
  utils::write.csv(res, file.path(cdir, paste0(tag, ".DE_matrix.csv")), row.names = FALSE)
  if (!is.null(fit$dds)) {
    norm_counts <- DESeq2::counts(fit$dds, normalized = TRUE)
    utils::write.csv(data.frame(feature_id = rownames(norm_counts), norm_counts, check.names = FALSE), file.path(cdir, paste0(tag, ".normalized_counts.csv")), row.names = FALSE)
  }
  utils::write.csv(data.frame(feature_id = rownames(fit$vst_mat), fit$vst_mat, check.names = FALSE), file.path(cdir, paste0(tag, ".vst_or_log_matrix.csv")), row.names = FALSE)

  if (as_bool(opt$make_plots)) {
    plot_volcano_simple(res, file.path(cdir, paste0(tag, ".volcano.pdf")), tag, label_col = label_col, color_col = color_col, sig_metric = "padj")
    plot_ma_simple(res, file.path(cdir, paste0(tag, ".MA.pdf")), tag, label_col = label_col, color_col = color_col, sig_metric = "padj")
    plot_top_heatmap_simple(fit$vst_mat, res, fit$sample_table, file.path(cdir, paste0(tag, ".top_heatmap.pdf")), tag, label_col = label_col, row_anno_col = row_anno_col)
  }
  if (mode == "gene" && as_bool(opt$run_go)) run_go_simple(res, file.path(cdir, tag), species = opt$species)
  if (mode == "gene" && as_bool(opt$run_gsea)) run_gsea_simple(res, file.path(cdir, tag), species = opt$species, fit = fit, contrast = ct, gene_anno = gene_anno)
  invisible(res)
}

if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && plot_n_cores > 1) {
  parallel::mclapply(contrast_list, run_one, mc.cores = min(plot_n_cores, length(contrast_list)), mc.preschedule = FALSE)
} else {
  lapply(contrast_list, run_one)
}

libsize <- colSums(count_mat, na.rm = TRUE)
libsize[libsize == 0] <- 1
qc_mat <- log2(t(t(count_mat) / libsize * 1e6) + 1)
plot_pca_simple(qc_mat, sample_table, file.path(tool_outdir, paste0(matrix_name, ".PCA.pdf")), paste0(matrix_name, " global"))
plot_corr_simple(qc_mat, file.path(tool_outdir, paste0(matrix_name, ".Pearson.pdf")), paste0(matrix_name, " global"))

writeLines(c(
  paste0("matrix=", matrix_file),
  paste0("matrix_name=", matrix_name),
  paste0("mode=", mode),
  paste0("contrast_count=", length(contrast_list)),
  paste0("padj_cutoff=", padj_cutoff),
  paste0("lfc_cutoff=", lfc_cutoff)
), file.path(tool_outdir, paste0(matrix_name, ".counts_to_de.summary.txt")))

message("[run_counts_to_de_matrix] DONE: ", tool_outdir)
