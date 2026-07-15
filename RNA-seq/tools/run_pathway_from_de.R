#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
    "RNA-seq pathway analysis from a gene DE matrix\n\n",
    "Required:\n",
    "  --de-matrix PATH         gene DE matrix CSV\n",
    "  --outdir PATH            output directory\n",
    "  --species STR            hg38 | mm10 | mm39\n",
    "  --tx2gene-path PATH      gene GTF/tx2gene for gene symbols\n\n",
    "Optional:\n",
    "  --function-file PATH     rnaseq-function.R; default auto-detect\n",
    "  --prefix STR             output prefix; default basename(de-matrix)\n",
    "  --case STR               contrast case label for titles, default case\n",
    "  --control STR            contrast control label for titles, default control\n",
    "  --run-go BOOL            default true\n",
    "  --run-gsea BOOL          default true\n",
    "  --min-gs-size NUM        default 10\n",
    "  --max-gs-size NUM        default 500\n",
    "  --pvalue-cutoff NUM      default 1\n",
    "  --topn-plot NUM          default 20\n",
    "  --disable-gseaplot2 BOOL default false\n",
    "  --padj-cutoff NUM        default 0.05\n",
    "  --lfc-cutoff NUM         default 0.58\n",
    "  -h, --help               show help\n",
    sep = ""
  )
}

parse_args <- function(args) {
  res <- list(case = "case", control = "control", run_go = "true", run_gsea = "true", min_gs_size = "10", max_gs_size = "500", pvalue_cutoff = "1", topn_plot = "20", disable_gseaplot2 = "false", padj_cutoff = "0.05", lfc_cutoff = "0.58")
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

as_bool <- function(x) tolower(as.character(x)) %in% c("1", "true", "yes", "y")
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

opt <- parse_args(args)
required <- c("de_matrix", "outdir", "species", "tx2gene_path")
missing <- required[!vapply(required, function(x) !is.null(opt[[x]]) && nzchar(opt[[x]]), logical(1))]
if (length(missing) > 0) {
  usage()
  stop("Missing required arguments: ", paste(missing, collapse = ", "))
}

script_file <- sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][1])
script_dir <- dirname(normalizePath(script_file))
function_file <- find_function_file(script_dir, opt)
de_file <- normalizePath(opt$de_matrix, mustWork = FALSE)
tx2gene_path <- normalizePath(opt$tx2gene_path, mustWork = FALSE)
if (!file.exists(de_file)) stop("DE matrix not found: ", de_file)
if (!file.exists(tx2gene_path)) stop("tx2gene-path not found: ", tx2gene_path)
outdir <- normalizePath(opt$outdir, mustWork = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- read_num(opt$padj_cutoff, 0.05)
lfc_cutoff <- read_num(opt$lfc_cutoff, 0.58)
baseMean_min <- 5
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
gene_anno <- load_gene_anno_from_gtf(tx2gene_path, species = opt$species)
res <- annotate_gene_res(res, gene_anno)
utils::write.csv(res, file.path(outdir, paste0(prefix, ".gene_annotated_DE_matrix.csv")), row.names = FALSE)

if (as_bool(opt$run_go)) {
  run_go_simple(res, file.path(outdir, prefix), species = opt$species)
}

gsea_df <- data.frame()
if (as_bool(opt$run_gsea)) {
  gsea_df <- run_gsea_multi_category(
    fit = list(mode = "de_matrix"),
    res_df = res,
    contrast = c("condition", opt$case, opt$control),
    prefix = file.path(outdir, prefix),
    species = opt$species,
    gene_anno = gene_anno,
    gsea_categories = get_default_gsea_categories(),
    minGSSize = read_int(opt$min_gs_size, 10),
    maxGSSize = read_int(opt$max_gs_size, 500),
    pvalueCutoff = read_num(opt$pvalue_cutoff, 1),
    topN_plot = read_int(opt$topn_plot, 20),
    max_terms_gseaplot2 = 12,
    make_gseaplot2 = !as_bool(opt$disable_gseaplot2)
  )
}

writeLines(c(
  paste0("de_matrix=", de_file),
  paste0("species=", opt$species),
  paste0("run_go=", as_bool(opt$run_go)),
  paste0("run_gsea=", as_bool(opt$run_gsea)),
  paste0("gsea_rows=", if (is.data.frame(gsea_df)) nrow(gsea_df) else 0)
), file.path(outdir, paste0(prefix, ".pathway.summary.txt")))
message("[run_pathway_from_de] DONE: ", outdir)
