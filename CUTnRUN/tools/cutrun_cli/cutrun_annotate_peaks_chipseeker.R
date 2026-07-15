#!/usr/bin/env Rscript
CUTRUN_ANNOTATE_NO_MAIN <- TRUE
cmd_file <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))[1]])
script_dir <- dirname(normalizePath(cmd_file, mustWork = TRUE))
source(file.path(script_dir, "cutrun_annotate_peaks_chippeakanno.R"))
suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  if (requireNamespace("txdbmaker", quietly = TRUE)) library(txdbmaker)
  library(AnnotationDbi)
})

help_text <- function() {
  paste0(
"Usage:\n",
"  cutrun_annotate_peaks_chipseeker --input peaks.broadPeak --ref hg38 --out-prefix out/sample [options]\n\n",
"This benchmark version uses ChIPseeker::annotatePeak() for gene/promoter annotation,\n",
"and uses the same GRanges relaxed/stringent TE logic as the ChIPpeakAnno version.\n",
"Peak annotation does not need BAM or bigWig.\n\n",
"Important options:\n",
"  --input FILE                  BED/narrowPeak/broadPeak.\n",
"  --ref hg38|mm39               Reference genome.\n",
"  --out-prefix PREFIX           Output prefix.\n",
"  --gtf FILE                    GTF used to build TxDb.\n",
"  --promoter-up INT             ChIPseeker tssRegion upstream bp. Default 2000.\n",
"  --promoter-down INT           ChIPseeker tssRegion downstream bp. Default 500.\n",
"  --te-min-length INT           Stringent TE minimum length. Default 50.\n",
"  --line-min-length INT         Stringent LINE minimum length. Default 800.\n",
"  --te-min-overlap INT          Stringent peak-TE overlap. Default 50.\n",
"  --relaxed-window INT          Relaxed TE window. Default 500.\n",
"  --include-te-class LIST       Keep only classes, comma-separated.\n",
"  --exclude-te-class-regex REG  Drop TE classes matching regex.\n",
"  --exclude-te-family-regex REG Drop TE families matching regex.\n",
"  --exclude-te-repname-regex REG Drop TE repNames matching regex.\n",
"  --no-plots                    Do not generate pdf/png plots.\n"
  )
}


annotation_to_feature <- function(x) {
  y <- rep("intergenic", length(x))
  y[grepl("Promoter", x, ignore.case = TRUE)] <- "promoter"
  y[grepl("Exon", x, ignore.case = TRUE)] <- "exon"
  y[grepl("Intron", x, ignore.case = TRUE)] <- "intron"
  y[grepl("UTR", x, ignore.case = TRUE)] <- "UTR"
  y[grepl("Downstream", x, ignore.case = TRUE)] <- "downstream"
  y[grepl("Distal Intergenic|Intergenic", x, ignore.case = TRUE)] <- "intergenic"
  y
}

default_txdb_cache <- function(gtf) {
  cache_dir <- file.path(Sys.getenv("HOME"), ".cache", "cutrun")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  info <- file.info(gtf)
  tag <- paste0(gsub("[^A-Za-z0-9_.-]", "_", basename(gtf)), "_", format(info$size, scientific = FALSE), "_", as.numeric(info$mtime))
  file.path(cache_dir, paste0("txdb_", tag, ".sqlite"))
}

write_chipseeker_gene_annotation <- function(peaks, args) {
  cache <- args$txdb_cache
  if (is.null(cache) || is.na(cache) || !nzchar(cache)) cache <- default_txdb_cache(args$gtf)
  if (file.exists(cache)) {
    message_run("loading cached TxDb: ", cache)
    txdb <- AnnotationDbi::loadDb(cache)
  } else {
    message_run("building TxDb from GTF for ChIPseeker")
    if (requireNamespace("txdbmaker", quietly = TRUE)) {
      txdb <- txdbmaker::makeTxDbFromGFF(args$gtf, format = "gtf")
    } else {
      txdb <- GenomicFeatures::makeTxDbFromGFF(args$gtf, format = "gtf")
    }
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    AnnotationDbi::saveDb(txdb, file = cache)
    message_run("saved TxDb cache: ", cache)
  }
  message_run("running ChIPseeker::annotatePeak")
  anno <- ChIPseeker::annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-args$promoter_up, args$promoter_down),
    verbose = FALSE
  )
  anno_df <- as.data.frame(anno)
  raw_path <- paste0(args$out_prefix, ".ChIPseeker_gene_annotation_raw.tsv")
  write.table(anno_df, raw_path, sep = "\t", quote = FALSE, row.names = FALSE)
  feature <- annotation_to_feature(anno_df$annotation)
  peak_id <- mcols(peaks)$peak_id
  gene_id <- if ("geneId" %in% names(anno_df)) as.character(anno_df$geneId) else "NA"
  distance <- if ("distanceToTSS" %in% names(anno_df)) as.integer(anno_df$distanceToTSS) else NA_integer_
  out <- data.frame(
    peak_id = peak_id,
    chrom = as.character(seqnames(peaks)),
    start = start(peaks) - 1L,
    end = end(peaks),
    gene_feature = feature,
    gene_id = gene_id,
    gene_name = if ("SYMBOL" %in% names(anno_df)) as.character(anno_df$SYMBOL) else gene_id,
    gene_type = "NA",
    overlap_bp = NA_integer_,
    distance_to_tss = distance,
    chipseeker_annotation = anno_df$annotation,
    stringsAsFactors = FALSE
  )
  out_path <- paste0(args$out_prefix, ".gene_structure.tsv")
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  table(out$gene_feature)
}

main_chipseeker <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  if (is.null(args$gtf)) args$gtf <- default_gtf(args$ref)
  if (is.null(args$te_anno)) args$te_anno <- default_te(args$ref)
  dir.create(dirname(args$out_prefix), showWarnings = FALSE, recursive = TRUE)
  message_run("input peak: ", args$input)
  message_run("GTF: ", args$gtf)
  message_run("TE annotation: ", args$te_anno)
  peak_style <- first_peak_style(args$input)
  ref_style <- te_ref_style(args$te_anno)
  chrom_action <- style_action(peak_style, ref_style)
  message_run("chromosome style action: ", chrom_action)
  peaks <- read_peaks(args$input, chrom_action)
  te_all <- load_te(args$te_anno, chrom_action, args)
  message_run("loaded peaks: ", length(peaks), "; TE records after manual filters: ", length(te_all))
  gene_counts <- write_chipseeker_gene_annotation(peaks, args)
  message_run("running TE relaxed/stringent annotation")
  te_results <- write_te_annotation(peaks, te_all, args$out_prefix, args)
  summary <- write_summary(args$out_prefix, length(peaks), gene_counts, te_results, te_all, args)
  if (!args$no_plots) {
    plot_dir <- plot_outputs(args$out_prefix, summary, te_results, args)
    message_run("plots: ", plot_dir)
  }
  write_parameters(args)
  message_run("done: ", args$out_prefix)
}

tryCatch(main_chipseeker(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})
