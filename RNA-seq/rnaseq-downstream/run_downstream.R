#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq downstream automation wrapper\n\n",
    "Required:\n",
    "  --results-dir PATH       上游 Nextflow 结果目录\n",
    "  --outdir PATH            下游输出目录\n",
    "  --species STR            hg38 | mm10 | mm39\n",
    "  --te-annotation-tsv PATH TE 注释 TSV\n",
    "  --tx2gene-path PATH      基因 GTF / tx2gene 文件\n",
    "  --sample-table PATH      样本表 CSV，至少包含 sample,condition,replicate\n",
    "  --contrast-file PATH     对比表 CSV，列名: case,control（兼容旧格式 group_col,case,control）\n\n",
    "Optional:\n",
    "  --tecount-strand STR     unstranded | forward | reverse，默认 unstranded\n",
    "  --padj-cutoff NUM        默认 0.05\n",
    "  --lfc-cutoff NUM         默认 1\n",
    "  --baseMean-min NUM       默认 5\n",
    "  --label-top-n NUM        默认 40\n",
    "  --heatmap-top-n NUM      默认 40\n",
    "  --volcano-orientation STR classic | horizontal，默认 classic\n",
    "  --gray-nonsig BOOL       true|false；非显著点是否统一灰色，默认 true\n",
    "  --exploratory-method STR logCPM_diff | edgeR_fixedBCV；无重复差异算法，默认 logCPM_diff\n",
    "  --exploratory-fixed-bcv NUM  edgeR_fixedBCV 的 BCV，human/TE 常用 0.4，默认 0.4\n",
    "  --plot-threads NUM       contrast 级别并行数，默认 1\n",
    "  --only-tools STR         只运行指定模块，逗号分隔；例如 TE_TEtranscripts\n",
    "  --skip-tools STR         跳过指定模块，逗号分隔；例如 TE_TElocal,TE_TEcount\n",
    "  --partial-input-policy STR  skip | error | allow，默认 skip\n",
    "  --allow-partial-inputs true|false  旧兼容参数；true 等价于 policy=allow\n",
    "  -h, --help               显示帮助\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(
    tecount_strand = "unstranded",
    padj_cutoff = 0.05,
    lfc_cutoff = 0.58,
    baseMean_min = 5,
    label_top_n = 40,
    heatmap_top_n = 40,
    volcano_orientation = "classic",
    gray_nonsig = "true",
    exploratory_method = "logCPM_diff",
    exploratory_fixed_bcv = "0.4",
    plot_threads = "1",
    only_tools = "",
    skip_tools = "",
    partial_input_policy = "skip",
    allow_partial_inputs = NULL
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

opt <- parse_args(args)
required_keys <- c("results_dir", "outdir", "species", "te_annotation_tsv", "tx2gene_path", "sample_table", "contrast_file")
missing_keys <- required_keys[!vapply(required_keys, function(k) !is.null(opt[[k]]) && nzchar(opt[[k]]), logical(1))]
if (length(missing_keys) > 0) {
  usage()
  stop("缺少必需参数: ", paste(missing_keys, collapse = ", "))
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_candidates <- c(
  file.path(script_dir, "rnaseq-function.R"),
  file.path(script_dir, "rnaseq-function_v6.R")
)
function_file <- function_candidates[file.exists(function_candidates)][1]
if (is.na(function_file) || !nzchar(function_file)) stop("未找到函数文件: ", paste(function_candidates, collapse = " | "))

results_dir <- normalizePath(opt$results_dir, mustWork = FALSE)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
species <- opt$species
te_annotation_tsv <- normalizePath(opt$te_annotation_tsv, mustWork = FALSE)
tx2gene_path <- normalizePath(opt$tx2gene_path, mustWork = FALSE)
sample_table_file <- normalizePath(opt$sample_table, mustWork = FALSE)
contrast_file <- normalizePath(opt$contrast_file, mustWork = FALSE)
tecount_strand <- opt$tecount_strand
padj_cutoff <- as.numeric(opt$padj_cutoff)
lfc_cutoff <- as.numeric(opt$lfc_cutoff)
baseMean_min <- as.numeric(opt$baseMean_min)
label_top_n <- as.integer(opt$label_top_n)
heatmap_top_n <- as.integer(opt$heatmap_top_n)
volcano_orientation <- tolower(as.character(opt$volcano_orientation))
if (!volcano_orientation %in% c("classic", "horizontal")) stop("--volcano-orientation 必须是 classic 或 horizontal")
exploratory_method <- as.character(opt$exploratory_method)
if (!exploratory_method %in% c("logCPM_diff", "edgeR_fixedBCV")) stop("--exploratory-method 必须是 logCPM_diff 或 edgeR_fixedBCV")
exploratory_fixed_bcv <- suppressWarnings(as.numeric(opt$exploratory_fixed_bcv))
if (is.na(exploratory_fixed_bcv) || exploratory_fixed_bcv <= 0) stop("--exploratory-fixed-bcv 必须是正数")
plot_threads <- suppressWarnings(as.integer(opt$plot_threads))
if (is.na(plot_threads) || plot_threads < 1) plot_threads <- 1L
parse_bool <- function(x, name) {
  y <- tolower(trimws(as.character(x)))
  if (y %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (y %in% c("false", "f", "no", "n", "0")) return(FALSE)
  stop(name, " 必须是 true 或 false")
}
partial_input_policy <- tolower(trimws(as.character(opt$partial_input_policy)))
if (!is.null(opt$allow_partial_inputs)) {
  partial_input_policy <- if (parse_bool(opt$allow_partial_inputs, "--allow-partial-inputs")) "allow" else "skip"
  message("[run_downstream] WARN: --allow-partial-inputs 已弃用，请改用 --partial-input-policy")
}
if (!partial_input_policy %in% c("skip", "error", "allow")) {
  stop("--partial-input-policy 必须是 skip、error 或 allow")
}
allow_partial_inputs <- identical(partial_input_policy, "allow")
gray_nonsig <- parse_bool(opt$gray_nonsig, "--gray-nonsig")
parse_tool_list <- function(x) {
  if (is.null(x) || !nzchar(x)) return(character(0))
  trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
}
only_tools <- parse_tool_list(opt$only_tools)
skip_tools <- parse_tool_list(opt$skip_tools)
blocked_tools <- character(0)
available_tools <- c(
  "gene_featureCounts",
  "gene_Salmon",
  "TE_TElocal",
  "TE_TEcount",
  "TE_TEtranscripts",
  "TE_REdiscoverTE",
  "TE_SalmonTE",
  "TE_Telescope",
  "panel_plots"
)
unknown_tools <- setdiff(c(only_tools, skip_tools), available_tools)
if (length(unknown_tools) > 0) {
  stop(
    "未知 tool 名称: ", paste(unknown_tools, collapse = ", "),
    "\n可用 tool: ", paste(available_tools, collapse = ", ")
  )
}
should_run_tool <- function(tool) {
  (length(only_tools) == 0 || tool %in% only_tools) && !tool %in% skip_tools && !tool %in% blocked_tools
}

log_msg <- function(...) cat("[run_downstream]", ..., "\n")

module_status_rows <- list()
record_module_status <- function(module, status, reason = "", expected = NA_integer_, found = NA_integer_) {
  module_status_rows[[length(module_status_rows) + 1L]] <<- data.frame(
    module = as.character(module), status = as.character(status), reason = as.character(reason),
    expected_samples = expected, found_samples = found, timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}

for (p in c(results_dir, te_annotation_tsv, tx2gene_path, sample_table_file, contrast_file)) {
  if (!file.exists(p)) stop("文件不存在: ", p)
}
if (!dir.exists(results_dir)) stop("results_dir 不存在: ", results_dir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
summary_dir <- file.path(outdir, "_summary")
matrix_dir <- file.path(outdir, "_matrices")
annotation_dir <- file.path(outdir, "_annotations")
downstream_log_dir <- file.path(outdir, "_logs")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(annotation_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(downstream_log_dir, recursive = TRUE, showWarnings = FALSE)

module_dir <- function(name) {
  p <- file.path(matrix_dir, name)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}
matrix_file <- function(name, ext = "rds", must_exist = FALSE) {
  p <- file.path(module_dir(name), paste0(name, "_matrix.", ext))
  old <- c(
    file.path(outdir, paste0(name, "_matrix.", ext)),
    file.path(outdir, name, paste0(name, "_matrix.", ext))
  )
  if (must_exist && !file.exists(p)) {
    old_hit <- old[file.exists(old)][1]
    if (!is.na(old_hit)) return(old_hit)
  }
  p
}
feature_map_file <- function(name, must_exist = FALSE) {
  p <- file.path(module_dir(name), paste0(name, "_feature_map.csv"))
  old <- file.path(outdir, paste0(name, "_feature_map.csv"))
  if (must_exist && !file.exists(p) && file.exists(old)) return(old)
  p
}
source_file <- function(name, must_exist = FALSE) {
  p <- file.path(module_dir(name), paste0(name, "_source.txt"))
  old <- file.path(outdir, paste0(name, "_source.txt"))
  if (must_exist && !file.exists(p) && file.exists(old)) return(old)
  p
}
te_annotation_rds <- function(must_exist = FALSE) {
  p <- file.path(annotation_dir, "te_annotation.rds")
  old <- file.path(outdir, "te_annotation.rds")
  if (must_exist && !file.exists(p) && file.exists(old)) return(old)
  p
}
te_annotation_preview <- function() file.path(annotation_dir, "te_annotation.preview.csv")

log_msg("results_dir =", results_dir)
log_msg("outdir =", outdir)
log_msg("species =", species)
log_msg("te_annotation_tsv =", te_annotation_tsv)
log_msg("tx2gene_path =", tx2gene_path)
log_msg("sample_table =", sample_table_file)
log_msg("contrast_file =", contrast_file)
log_msg("tecount_strand =", tecount_strand)
log_msg("volcano_orientation =", volcano_orientation)
log_msg("gray_nonsig =", gray_nonsig)
log_msg("exploratory_method =", exploratory_method)
log_msg("exploratory_fixed_bcv =", exploratory_fixed_bcv)
log_msg("plot_threads =", plot_threads)
log_msg("only_tools =", if (length(only_tools) == 0) "ALL" else paste(only_tools, collapse = ","))
log_msg("skip_tools =", if (length(skip_tools) == 0) "NONE" else paste(skip_tools, collapse = ","))
log_msg("partial_input_policy =", partial_input_policy)

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

log_msg("sample_table rows =", nrow(sample_table))
log_msg("contrast count =", length(contrast_list))

assign("volcano_orientation", volcano_orientation, envir = .GlobalEnv)
assign("gray_nonsig", gray_nonsig, envir = .GlobalEnv)
assign("exploratory_method", exploratory_method, envir = .GlobalEnv)
assign("exploratory_fixed_bcv", exploratory_fixed_bcv, envir = .GlobalEnv)
options(plot_n_cores = plot_threads)
source(function_file)

scan_files <- function() {
  list(
    fc_files = find_files(file.path(results_dir, "03_gene_featurecounts"), "\\.featureCounts\\.txt$"),
    stringtie_files = find_files(file.path(results_dir, "04_gene_stringtie"), "\\.stringtie\\.gtf$"),
    salmon_files = find_files(file.path(results_dir, "05_gene_salmon"), "quant\\.sf$"),
    telocal_files = find_files(file.path(results_dir, "08_telocal"), "\\.cntTable$"),
    tecount_class_files = find_files(file.path(results_dir, "07_tecount"), paste0("_tecount_class_", tecount_strand, "\\.count\\.txt$")),
    tecount_family_files = find_files(file.path(results_dir, "07_tecount"), paste0("_tecount_family_", tecount_strand, "\\.count\\.txt$")),
    tecount_subfamily_files = find_files(file.path(results_dir, "07_tecount"), paste0("_tecount_subfamily_", tecount_strand, "\\.count\\.txt$")),
    tetrans_files = find_files(file.path(results_dir, "09_TEtranscripts"), "_tetranscripts\\.cntTable$"),
    rediscover_quant_files = find_files(file.path(results_dir, "10_rediscoverte"), "quant\\.sf$"),
    rediscover_rollup_files = find_files(file.path(results_dir, "10_rediscoverte_rollup"), "\\.(RDS|rds)$"),
    salmonte_files = find_files(file.path(results_dir, "11_salmonte"), "^EXPR\\.csv$"),
    salmonte_clade_files = find_files(file.path(results_dir, "11_salmonte"), "^clades\\.csv$"),
    telescope_files = find_files(file.path(results_dir, "12_telescope"), "[-_]telescope_report.tsv$")
  )
}

counts <- scan_files()
summary_lines <- c(
  paste0("featureCounts=", length(counts$fc_files)),
  paste0("StringTie=", length(counts$stringtie_files)),
  paste0("Salmon=", length(counts$salmon_files)),
  paste0("TElocal=", length(counts$telocal_files)),
  paste0("TEcount_class=", length(counts$tecount_class_files)),
  paste0("TEcount_family=", length(counts$tecount_family_files)),
  paste0("TEcount_subfamily=", length(counts$tecount_subfamily_files)),
  paste0("TEtranscripts=", length(counts$tetrans_files)),
  paste0("REdiscoverTE_quant=", length(counts$rediscover_quant_files)),
  paste0("REdiscoverTE_rollup=", length(counts$rediscover_rollup_files)),
  paste0("SalmonTE_expr=", length(counts$salmonte_files)),
  paste0("SalmonTE_clades=", length(counts$salmonte_clade_files)),
  paste0("Telescope=", length(counts$telescope_files))
)
writeLines(summary_lines, file.path(summary_dir, "downstream_inputs_summary.txt"))
for (x in summary_lines) log_msg(x)

safe_skip <- function(tag, reason) {
  msg <- paste0("[", tag, "] SKIP: ", reason)
  message(msg)
  writeLines(msg, file.path(downstream_log_dir, paste0(tag, "_SKIP.txt")))
  record_module_status(tag, "skipped", reason)
}

check_sample_coverage <- function(tool, files, sample_fun, label) {
  if (!should_run_tool(tool) || allow_partial_inputs || length(files) == 0) return(invisible(TRUE))

  expected <- unique(sample_table$sample)
  found <- unique(sample_fun(files))
  found <- found[nzchar(found)]

  missing <- setdiff(expected, found)
  extra <- setdiff(found, expected)

  if (length(missing) == 0) return(invisible(TRUE))

  msg <- c(
    paste0("[", label, "] 输入样本不完整。"),
    paste0("expected_samples=", length(expected)),
    paste0("found_samples=", length(intersect(expected, found))),
    paste0("missing=", paste(missing, collapse = ",")),
    paste0("extra=", if (length(extra) == 0) "NONE" else paste(extra, collapse = ",")),
    paste0("policy=", partial_input_policy)
  )
  writeLines(msg, file.path(downstream_log_dir, paste0(label, "_INCOMPLETE_INPUTS.txt")))
  record_module_status(label, "incomplete", paste0("missing=", paste(missing, collapse = ",")), length(expected), length(intersect(expected, found)))
  if (identical(partial_input_policy, "error")) stop(paste(msg, collapse = "\n"), call. = FALSE)
  blocked_tools <<- unique(c(blocked_tools, tool))
  message("[", label, "] SKIP whole module; missing samples: ", paste(missing, collapse = ","))
  invisible(FALSE)
}

strip_suffix <- function(x, pattern) sub(pattern, "", basename(x))
sample_from_featurecounts <- function(files) strip_suffix(files, "\\.featureCounts\\.txt$")
sample_from_salmon <- function(files) sub("_salmon$", "", basename(dirname(files)))
sample_from_telocal <- function(files) strip_suffix(files, "_telocal\\.cntTable$")
sample_from_tetranscripts <- function(files) strip_suffix(files, "_tetranscripts\\.cntTable$")
sample_from_tecount <- function(files) sub("_tecount_(class|family|subfamily)_.+$", "", basename(files))
sample_from_rediscoverte <- function(files) sub("_rediscoverte$", "", basename(dirname(files)))

check_sample_coverage("gene_featureCounts", counts$fc_files, sample_from_featurecounts, "gene_featureCounts")
check_sample_coverage("gene_Salmon", counts$salmon_files, sample_from_salmon, "gene_Salmon")
check_sample_coverage("TE_TElocal", counts$telocal_files, sample_from_telocal, "TE_TElocal")
check_sample_coverage("TE_TEtranscripts", counts$tetrans_files, sample_from_tetranscripts, "TE_TEtranscripts")
check_sample_coverage("TE_TEcount", counts$tecount_class_files, sample_from_tecount, "TE_TEcount_class")
check_sample_coverage("TE_TEcount", counts$tecount_family_files, sample_from_tecount, "TE_TEcount_family")
check_sample_coverage("TE_TEcount", counts$tecount_subfamily_files, sample_from_tecount, "TE_TEcount_subfamily")
check_sample_coverage("TE_REdiscoverTE", counts$rediscover_quant_files, sample_from_rediscoverte, "TE_REdiscoverTE")

needs_te_annotation <-
  (should_run_tool("TE_TElocal") && length(counts$telocal_files) >= 2) ||
  (should_run_tool("TE_TEcount") && max(length(counts$tecount_class_files), length(counts$tecount_family_files), length(counts$tecount_subfamily_files)) >= 2) ||
  (should_run_tool("TE_TEtranscripts") && length(counts$tetrans_files) >= 2) ||
  (should_run_tool("TE_REdiscoverTE") && length(counts$rediscover_rollup_files) >= 1) ||
  (should_run_tool("TE_SalmonTE") && length(counts$salmonte_files) >= 2) ||
  (should_run_tool("TE_Telescope") && length(counts$telescope_files) >= 2)

if (needs_te_annotation) {
  safe_run("te_annotation_load", {
    te_anno <- read_te_annotation(te_annotation_tsv)
    saveRDS(te_anno, te_annotation_rds())
    utils::write.csv(te_anno, te_annotation_preview(), row.names = FALSE)
    log_msg("TE annotation rows =", nrow(te_anno))
  })
} else {
  safe_skip("te_annotation_load", "没有已选择且输入完整的 TE 模块，跳过大型 TE 注释加载")
}

if (!should_run_tool("gene_featureCounts")) {
  safe_skip("gene_featureCounts", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$fc_files) >= 2) {
  safe_run("gene_featureCounts_prepare", {
    gene_counts_fc <- read_featurecounts_matrix(counts$fc_files)
    saveRDS(gene_counts_fc, matrix_file("gene_featureCounts", "rds"))
    write_matrix_csv(gene_counts_fc, matrix_file("gene_featureCounts", "csv"))
    log_msg("gene_featureCounts matrix dim =", paste(dim(gene_counts_fc), collapse = " x "))
  })
  safe_run("gene_featureCounts_run", {
    gene_counts_fc <- readRDS(matrix_file("gene_featureCounts", "rds", TRUE))
    gene_anno <- load_gene_anno_from_gtf(tx2gene_path, species = species)
    run_and_plot_deseq(
      gene_counts_fc, sample_table, contrast_list,
      out_prefix = "gene_featureCounts", do_go = TRUE, species = species,
      gene_anno = gene_anno
    )
  })
} else {
  safe_skip("gene_featureCounts", "featureCounts 文件少于2个")
}

if (!should_run_tool("gene_Salmon")) {
  safe_skip("gene_Salmon", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$salmon_files) >= 2) {
  safe_run("gene_salmon_prepare", {
    gene_counts_salmon <- read_salmon_quant(counts$salmon_files, tx2gene_path = tx2gene_path)
    saveRDS(gene_counts_salmon, matrix_file("gene_Salmon", "rds"))
    write_matrix_csv(gene_counts_salmon, matrix_file("gene_Salmon", "csv"))
    log_msg("gene_Salmon matrix dim =", paste(dim(gene_counts_salmon), collapse = " x "))
  })
  safe_run("gene_salmon_run", {
    gene_counts_salmon <- readRDS(matrix_file("gene_Salmon", "rds", TRUE))
    gene_anno <- load_gene_anno_from_gtf(tx2gene_path, salmon_files = counts$salmon_files, species = species)
    run_and_plot_deseq(
      gene_counts_salmon, sample_table, contrast_list,
      out_prefix = "gene_Salmon", do_go = TRUE, species = species,
      gene_anno = gene_anno
    )
  })
} else {
  safe_skip("gene_Salmon", "Salmon quant.sf 文件少于2个")
}

if (!should_run_tool("TE_TElocal")) {
  safe_skip("TE_TElocal", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$telocal_files) >= 2 && file.exists(te_annotation_rds(TRUE))) {
  safe_run("TE_TElocal_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    telocal_all <- read_telocal_matrix(counts$telocal_files)
    saveRDS(telocal_all, matrix_file("TE_TElocal_all", "rds"))
    write_matrix_csv(telocal_all, matrix_file("TE_TElocal_all", "csv"))

    telocal_map <- parse_telocal_feature_map(rownames(telocal_all), te_anno)
    utils::write.csv(telocal_map, feature_map_file("TE_TElocal"), row.names = FALSE)

    telocal_locus <- filter_te_rows_by_map(telocal_all, telocal_map, tag = "TE_TElocal_locus")
    saveRDS(telocal_locus, matrix_file("TE_TElocal_locus", "rds"))
    write_matrix_csv(telocal_locus, matrix_file("TE_TElocal_locus", "csv"))

    telocal_map_te <- attr(telocal_locus, "feature_map")
    telocal_repName <- aggregate_by_feature_map(telocal_locus, telocal_map_te, "repName", tag = "TE_TElocal_repName")
    telocal_repFamily <- aggregate_by_feature_map(telocal_locus, telocal_map_te, "repFamily", tag = "TE_TElocal_repFamily")
    telocal_repClass <- aggregate_by_feature_map(telocal_locus, telocal_map_te, "repClass", tag = "TE_TElocal_repClass")

    saveRDS(telocal_repName, matrix_file("TE_TElocal_repName", "rds"))
    saveRDS(telocal_repFamily, matrix_file("TE_TElocal_repFamily", "rds"))
    saveRDS(telocal_repClass, matrix_file("TE_TElocal_repClass", "rds"))
    write_matrix_csv(telocal_repName, matrix_file("TE_TElocal_repName", "csv"))
    write_matrix_csv(telocal_repFamily, matrix_file("TE_TElocal_repFamily", "csv"))
    write_matrix_csv(telocal_repClass, matrix_file("TE_TElocal_repClass", "csv"))
  })

  for (lev in c("locus", "repName", "repFamily", "repClass")) {
    safe_run(paste0("TE_TElocal_", lev, "_run"), {
      mat_file <- matrix_file(paste0("TE_TElocal_", lev), "rds", TRUE)
      mat <- readRDS(mat_file)
      te_anno <- readRDS(te_annotation_rds(TRUE))
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_TElocal_", lev),
        do_go = FALSE,
        gene_anno = te_anno,
        te_label_level = if (lev == "locus") "locus_id" else lev,
        te_color_level = if (lev %in% c("locus", "repName")) "repFamily" else "repClass"
      )
    })
  }
} else {
  safe_skip("TE_TElocal", "TElocal 文件少于2个或 TE 注释未生成")
}

if (!should_run_tool("TE_TEcount")) {
  safe_skip("TE_TEcount", "被 --only-tools/--skip-tools 跳过")
} else {
  safe_run("TE_TEcount_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    groups <- list(
      class = counts$tecount_class_files,
      family = counts$tecount_family_files,
      subfamily = counts$tecount_subfamily_files
    )
    available <- groups[vapply(groups, length, integer(1)) >= 2]
    if (length(available) == 0) stop("未发现可用于分析的 TEcount 文件（至少每层级需要>=2个样本）")

    for (nm in names(available)) {
      mat <- read_count_files_generic(
        available[[nm]],
        sample_fun = function(f) sample_from_parent(f, "_tecount$"),
        tag = paste0("TEcount_", nm)
      )
      mat <- drop_short_repeat_matrix(mat, tag = paste0("TE_TEcount_", nm))
      saveRDS(mat, matrix_file(paste0("TE_TEcount_", nm), "rds"))
      write_matrix_csv(mat, matrix_file(paste0("TE_TEcount_", nm), "csv"))

      label_level <- switch(nm, class = "repClass", family = "repFamily", subfamily = "repName")
      color_level <- if (nm == "subfamily") "repFamily" else "repClass"
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_TEcount_", nm),
        do_go = FALSE,
        gene_anno = te_anno,
        te_label_level = label_level,
        te_color_level = color_level
      )
    }
  })
  if (!file.exists(matrix_file("TE_TEcount_class", "rds", TRUE)) &&
      !file.exists(matrix_file("TE_TEcount_family", "rds", TRUE)) &&
      !file.exists(matrix_file("TE_TEcount_subfamily", "rds", TRUE))) {
    safe_skip("TE_TEcount", "没有足够的 TEcount 输入")
  }
}

if (!should_run_tool("TE_TEtranscripts")) {
  safe_skip("TE_TEtranscripts", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$tetrans_files) >= 2 && file.exists(te_annotation_rds(TRUE))) {
  safe_run("TE_TEtranscripts_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    tetrans_mat <- read_count_files_generic(
      counts$tetrans_files,
      sample_fun = function(f) sample_from_parent(f, "_tetranscripts$"),
      tag = "TE_TEtranscripts"
    )
    tetrans_mat <- drop_short_repeat_matrix(tetrans_mat, tag = "TE_TEtranscripts")

    tetrans_map <- parse_tetranscripts_feature_map(rownames(tetrans_mat), te_anno)
    utils::write.csv(tetrans_map, feature_map_file("TE_TEtranscripts"), row.names = FALSE)

    tetrans_mat <- filter_te_rows_by_map(tetrans_mat, tetrans_map, tag = "TE_TEtranscripts")
    saveRDS(tetrans_mat, matrix_file("TE_TEtranscripts", "rds"))
    write_matrix_csv(tetrans_mat, matrix_file("TE_TEtranscripts", "csv"))

    tetrans_map_te <- attr(tetrans_mat, "feature_map")
    tetrans_repName <- aggregate_by_feature_map(tetrans_mat, tetrans_map_te, "repName", tag = "TE_TEtranscripts_repName")
    tetrans_repFamily <- aggregate_by_feature_map(tetrans_mat, tetrans_map_te, "repFamily", tag = "TE_TEtranscripts_repFamily")
    tetrans_repClass <- aggregate_by_feature_map(tetrans_mat, tetrans_map_te, "repClass", tag = "TE_TEtranscripts_repClass")

    saveRDS(tetrans_repName, matrix_file("TE_TEtranscripts_repName", "rds"))
    saveRDS(tetrans_repFamily, matrix_file("TE_TEtranscripts_repFamily", "rds"))
    saveRDS(tetrans_repClass, matrix_file("TE_TEtranscripts_repClass", "rds"))
    write_matrix_csv(tetrans_repName, matrix_file("TE_TEtranscripts_repName", "csv"))
    write_matrix_csv(tetrans_repFamily, matrix_file("TE_TEtranscripts_repFamily", "csv"))
    write_matrix_csv(tetrans_repClass, matrix_file("TE_TEtranscripts_repClass", "csv"))
  })

  for (lev in c("repName", "repFamily", "repClass")) {
    f <- matrix_file(paste0("TE_TEtranscripts_", lev), "rds", TRUE)
    if (!file.exists(f)) next

    safe_run(paste0("TE_TEtranscripts_", lev, "_run"), {
      mat <- readRDS(f)
      te_anno <- readRDS(te_annotation_rds(TRUE))
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_TEtranscripts_", lev),
        do_go = FALSE,
        gene_anno = te_anno,
        te_label_level = lev,
        te_color_level = if (lev == "repName") "repFamily" else "repClass"
      )
    })
  }
} else {
  safe_skip("TE_TEtranscripts", "TEtranscripts 文件少于2个或 TE 注释未生成")
}

if (!should_run_tool("TE_REdiscoverTE")) {
  safe_skip("TE_REdiscoverTE", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$rediscover_rollup_files) >= 1 && file.exists(te_annotation_rds(TRUE))) {
  safe_run("TE_REdiscoverTE_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    cand <- counts$rediscover_rollup_files[
      grepl("^RE_", basename(counts$rediscover_rollup_files), ignore.case = TRUE) &
        grepl("raw[_]?counts", basename(counts$rediscover_rollup_files), ignore.case = TRUE)
    ]
    if (length(cand) == 0) {
      cand <- counts$rediscover_rollup_files[
        grepl("^RE_", basename(counts$rediscover_rollup_files), ignore.case = TRUE)
      ]
    }
    if (length(cand) == 0) stop("rollup 中未发现 TE 相关 RDS")
    f <- cand[1]
    mat_raw <- extract_rds_matrix(f)
    mat <- tryCatch(
      filter_te_rows(mat_raw, te_anno, preferred_match_cols = c("repName", "locus_id", "gene_name"), tag = "TE_REdiscoverTE"),
      error = function(e) {
        message("[TE_REdiscoverTE] 注释过滤失败，回退到原始矩阵: ", conditionMessage(e))
        mat_raw
      }
    )
    saveRDS(mat, matrix_file("TE_REdiscoverTE", "rds"))
    write_matrix_csv(mat, matrix_file("TE_REdiscoverTE", "csv"))
    writeLines(c(
      paste("using rollup:", f),
      paste("per-sample quant.sf found:", length(counts$rediscover_quant_files))
    ), source_file("TE_REdiscoverTE"))

    for (lev in c("repName", "repFamily", "repClass")) {
      agg <- tryCatch(
        aggregate_by_annotation(mat, te_anno, lev, preferred_match_cols = c("repName", "locus_id", "gene_name")),
        error = function(e) {
          message("[TE_REdiscoverTE_", lev, "] skip: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(agg)) {
        saveRDS(agg, matrix_file(paste0("TE_REdiscoverTE_", lev), "rds"))
        write_matrix_csv(agg, matrix_file(paste0("TE_REdiscoverTE_", lev), "csv"))
      }
    }
  })

  for (lev in c("repName", "repFamily", "repClass")) {
    f <- matrix_file(paste0("TE_REdiscoverTE_", lev), "rds", TRUE)
    if (!file.exists(f)) next
    safe_run(paste0("TE_REdiscoverTE_", lev, "_run"), {
      mat <- readRDS(f)
      te_anno <- readRDS(te_annotation_rds(TRUE))
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_REdiscoverTE_", lev),
        do_go = FALSE,
        gene_anno = te_anno,
        te_label_level = lev,
        te_color_level = if (lev == "repName") "repFamily" else "repClass"
      )
    })
  }
} else {
  safe_skip("TE_REdiscoverTE", "未发现 REdiscoverTE rollup RDS")
}

if (!should_run_tool("TE_SalmonTE")) {
  safe_skip("TE_SalmonTE", "被 --only-tools/--skip-tools 跳过")
} else if (length(counts$salmonte_files) >= 2 && file.exists(te_annotation_rds(TRUE))) {
  safe_run("TE_SalmonTE_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    clade_anno <- read_salmonte_clades(counts$salmonte_clade_files)
    mat <- read_salmonte_expr(counts$salmonte_files)
    saveRDS(mat, matrix_file("TE_SalmonTE", "rds"))
    write_matrix_csv(mat, matrix_file("TE_SalmonTE", "csv"))
    writeLines(c(counts$salmonte_files, counts$salmonte_clade_files), source_file("TE_SalmonTE"))

    anno_use <- choose_repname_annotation(
      rownames(mat), te_anno, clade_anno,
      alt_match_cols = c("feature_id", "repName"),
      tag = "TE_SalmonTE_anno"
    )

    for (lev in c("repName", "repFamily", "repClass")) {
      agg <- tryCatch(
        aggregate_by_annotation(
          mat, anno_use, lev,
          preferred_match_cols = if (identical(anno_use, te_anno)) c("repName", "locus_id", "gene_name") else c("feature_id", "repName")
        ),
        error = function(e) {
          message("[TE_SalmonTE_", lev, "] skip: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(agg)) {
        saveRDS(agg, matrix_file(paste0("TE_SalmonTE_", lev), "rds"))
        write_matrix_csv(agg, matrix_file(paste0("TE_SalmonTE_", lev), "csv"))
      }
    }
  })

  for (lev in c("repName", "repFamily", "repClass")) {
    f <- matrix_file(paste0("TE_SalmonTE_", lev), "rds", TRUE)
    if (!file.exists(f)) next
    safe_run(paste0("TE_SalmonTE_", lev, "_run"), {
      mat <- readRDS(f)
      te_anno <- readRDS(te_annotation_rds(TRUE))
      clade_anno <- read_salmonte_clades(counts$salmonte_clade_files)
      anno_use <- choose_repname_annotation(
        rownames(readRDS(matrix_file("TE_SalmonTE", "rds", TRUE))), te_anno, clade_anno,
        alt_match_cols = c("feature_id", "repName"),
        tag = paste0("TE_SalmonTE_", lev, "_run_anno")
      )
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_SalmonTE_", lev),
        do_go = FALSE,
        gene_anno = anno_use,
        te_label_level = lev,
        te_color_level = if (lev == "repName") "repFamily" else "repClass"
      )
    })
  }
} else {
  safe_skip("TE_SalmonTE", "SalmonTE EXPR.csv 文件少于2个或 TE 注释未生成")
}


if (!should_run_tool("TE_Telescope")) {
  safe_skip("TE_Telescope", "skipped by --only-tools/--skip-tools")
} else if (length(counts$telescope_files) >= 2 && file.exists(te_annotation_rds(TRUE))) {
  safe_run("TE_Telescope_prepare", {
    te_anno <- readRDS(te_annotation_rds(TRUE))
    tscope_mat <- read_telescope_report(counts$telescope_files)
    tscope_mat <- drop_short_repeat_matrix(tscope_mat, tag = "TE_Telescope")
    saveRDS(tscope_mat, matrix_file("TE_Telescope", "rds"))
    write_matrix_csv(tscope_mat, matrix_file("TE_Telescope", "csv"))
    writeLines(counts$telescope_files, source_file("TE_Telescope"))

    tscope_map <- parse_telescope_feature_map(rownames(tscope_mat), te_anno)
    utils::write.csv(tscope_map, feature_map_file("TE_Telescope"), row.names = FALSE)

    tscope_mat <- filter_te_rows_by_map(tscope_mat, tscope_map, tag = "TE_Telescope")
    tscope_map_te <- attr(tscope_mat, "feature_map")

    tscope_repName <- aggregate_by_feature_map(tscope_mat, tscope_map_te, "repName", tag = "TE_Telescope_repName")
    tscope_repFamily <- aggregate_by_feature_map(tscope_mat, tscope_map_te, "repFamily", tag = "TE_Telescope_repFamily")
    tscope_repClass <- aggregate_by_feature_map(tscope_mat, tscope_map_te, "repClass", tag = "TE_Telescope_repClass")

    saveRDS(tscope_repName, matrix_file("TE_Telescope_repName", "rds"))
    saveRDS(tscope_repFamily, matrix_file("TE_Telescope_repFamily", "rds"))
    saveRDS(tscope_repClass, matrix_file("TE_Telescope_repClass", "rds"))
    write_matrix_csv(tscope_repName, matrix_file("TE_Telescope_repName", "csv"))
    write_matrix_csv(tscope_repFamily, matrix_file("TE_Telescope_repFamily", "csv"))
    write_matrix_csv(tscope_repClass, matrix_file("TE_Telescope_repClass", "csv"))
  })

  for (lev in c("repName", "repFamily", "repClass")) {
    f <- matrix_file(paste0("TE_Telescope_", lev), "rds", TRUE)
    if (!file.exists(f)) next
    safe_run(paste0("TE_Telescope_", lev, "_run"), {
      mat <- readRDS(f)
      te_anno <- readRDS(te_annotation_rds(TRUE))
      run_and_plot_deseq(
        mat, sample_table, contrast_list,
        out_prefix = paste0("TE_Telescope_", lev),
        do_go = FALSE,
        gene_anno = te_anno,
        te_label_level = lev,
        te_color_level = if (lev == "repName") "repFamily" else "repClass"
      )
    })
  }
} else {
  safe_skip("TE_Telescope", "Telescope files < 2 or TE annotation missing")
}

make_panel_pdf <- function(pdf_files, outfile, ncol = 2) {
  if (length(pdf_files) == 0) return(NULL)

  if (!requireNamespace("magick", quietly = TRUE)) {
    log_msg("SKIP panel:", outfile, "because R package 'magick' is not installed")
    return(NULL)
  }

  imgs <- lapply(pdf_files, function(f) {
    tryCatch(
      magick::image_read(f),
      error = function(e) {
        log_msg("WARN image_read failed:", f, "msg =", conditionMessage(e))
        NULL
      }
    )
  })
  imgs <- Filter(Negate(is.null), imgs)
  if (length(imgs) == 0) return(NULL)

  rows <- split(imgs, ceiling(seq_along(imgs) / ncol))
  row_imgs <- lapply(rows, function(x) magick::image_append(do.call(c, x)))
  panel <- magick::image_append(do.call(c, row_imgs), stack = TRUE)
  magick::image_write(panel, path = outfile, format = "pdf")
  log_msg("panel written:", outfile)
  invisible(outfile)
}
if (should_run_tool("panel_plots")) {
safe_run("panel_plots", {
  all_pdf <- list.files(outdir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
  if (length(all_pdf) == 0) {
    log_msg("SKIP panel_plots: no pdf found")
    return(NULL)
  }
  panel_dir <- file.path(outdir, "_panels")
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  all_pdf <- all_pdf[!grepl("/_panels/", all_pdf, fixed = TRUE)]
  all_pdf <- unique(all_pdf)
  panel_tasks <- list()

  add_panel_task <- function(files, outfile, ncol = 2) {
    if (length(files) == 0) return(invisible(NULL))
    panel_tasks[[length(panel_tasks) + 1L]] <<- list(files = files[order(files)], outfile = outfile, ncol = ncol)
    invisible(NULL)
  }

  global_panel_specs <- list(
    PCA = c("_PCA_global.pdf", "_PCA.pdf"),
    Pearson = c("_Pearson_global.pdf", "_Pearson.pdf")
  )
  for (panel_kind in names(global_panel_specs)) {
    suffix_hit <- Reduce(
      `|`,
      lapply(global_panel_specs[[panel_kind]], function(suf) endsWith(basename(all_pdf), suf))
    )
    hit <- all_pdf[suffix_hit]
    if (length(hit) == 0) next

    add_panel_task(hit, file.path(panel_dir, paste0("global_", panel_kind, ".pdf")), ncol = 2)
  }

  suffixes <- c("_PCA.pdf", "_Pearson.pdf", "_volcano.pdf", "_MA.pdf", "_heatmap.pdf")

  for (ct in contrast_list) {
    case_name <- as.character(ct[2])
    ctrl_name <- as.character(ct[3])
    contrast_key <- paste0("_", case_name, "_vs_", ctrl_name, "_")

    pdf_hit <- all_pdf[grepl(contrast_key, basename(all_pdf), fixed = TRUE)]
    if (length(pdf_hit) == 0) {
      log_msg("SKIP panel for contrast =", paste(case_name, "vs", ctrl_name), "because no pdf matched")
      next
    }

    for (suf in suffixes) {
      hit <- pdf_hit[endsWith(pdf_hit, suf)]
      if (length(hit) == 0) next

      hit <- hit[order(hit)]
      panel_name <- paste0(
        gsub("[^A-Za-z0-9._-]+", "_", paste(case_name, "vs", ctrl_name)),
        suf
      )
      add_panel_task(hit, file.path(panel_dir, panel_name), ncol = 2)
    }
  }

  run_panel_task <- function(task) make_panel_pdf(task$files, task$outfile, ncol = task$ncol)
  if (length(panel_tasks) == 0) {
    log_msg("SKIP panel_plots: no panel task built")
  } else if (.Platform$OS.type == "unix" && plot_threads > 1 && length(panel_tasks) > 1) {
    parallel::mclapply(panel_tasks, run_panel_task, mc.cores = min(plot_threads, length(panel_tasks), 4L), mc.preschedule = FALSE)
  } else {
    lapply(panel_tasks, run_panel_task)
  }
})
} else {
  safe_skip("panel_plots", "被 --only-tools/--skip-tools 跳过")
}
all_pdf <- list.files(outdir, pattern = "\\.pdf$", full.names = TRUE, recursive = TRUE)
all_csv <- list.files(outdir, pattern = "\\.(csv|txt)$", full.names = TRUE, recursive = TRUE)
log_msg("PDF count =", length(all_pdf))
log_msg("CSV/TXT count =", length(all_csv))
writeLines(c(
  paste0("pdf=", length(all_pdf)),
  paste0("csv_txt=", length(all_csv))
), file.path(summary_dir, "downstream_outputs_summary.txt"))

if (length(module_status_rows) > 0) {
  module_status <- do.call(rbind, module_status_rows)
  utils::write.csv(module_status, file.path(summary_dir, "module_status.csv"), row.names = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(module_status, file.path(summary_dir, "module_status.json"), pretty = TRUE, na = "null")
  }
}

report_script <- normalizePath(file.path(script_dir, "..", "tools", "build_rnaseq_report.py"), mustWork = FALSE)
if (file.exists(report_script)) {
  report_args <- c(
    report_script, "--plot-dir", outdir,
    "--padj-cutoff", as.character(padj_cutoff),
    "--lfc-cutoff", as.character(lfc_cutoff),
    "--base-mean-min", as.character(baseMean_min)
  )
  report_status <- system2("python3", report_args)
  if (!identical(report_status, 0L)) message("[run_downstream] WARN: HTML report generation failed: exit=", report_status)
}
