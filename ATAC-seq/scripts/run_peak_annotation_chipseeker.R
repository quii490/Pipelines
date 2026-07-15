#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  pkgs <- c("ChIPseeker", "GenomicFeatures", "AnnotationDbi", "GenomicRanges",
            "IRanges", "rtracklayer", "dplyr", "ggplot2", "tidyr", "stringr")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing R packages: ", paste(missing, collapse = ", "))
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
})

args_raw <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat(
"Usage:
  run_peak_annotation_chipseeker.R --input peaks.bed --gtf genes.gtf --out-prefix out/sample [options]

Options:
  --species STR
  --te-bed FILE
  --promoter-up INT
  --promoter-down INT
  --txdb-cache FILE
  --no-plots
",
  sep = "")
}

parse_args <- function(x) {
  out <- list(
    input = NULL, species = "hg38", out_prefix = NULL, gtf = NULL, te_bed = NULL,
    promoter_up = 3000L, promoter_down = 3000L, txdb_cache = NULL, no_plots = FALSE
  )
  i <- 1
  while (i <= length(x)) {
    key <- x[[i]]
    if (key %in% c("-h", "--help")) {
      usage()
      quit(status = 0)
    }
    if (key == "--no-plots") {
      out$no_plots <- TRUE
      i <- i + 1
      next
    }
    if (i == length(x)) stop("Missing value for ", key)
    val <- x[[i + 1]]
    key2 <- gsub("^-+", "", key)
    key2 <- gsub("-", "_", key2)
    if (!key2 %in% names(out)) stop("Unknown option: ", key)
    if (key2 %in% c("promoter_up", "promoter_down")) val <- as.integer(val)
    out[[key2]] <- val
    i <- i + 2
  }
  if (is.null(out$input) || is.null(out$out_prefix) || is.null(out$gtf)) {
    usage()
    stop("--input, --gtf and --out-prefix are required")
  }
  out
}

msg <- function(...) message("[run_peak_annotation] ", paste0(..., collapse = ""))

save_plot_both <- function(prefix_pdf, plot, width = 7, height = 4.8) {
  dir.create(dirname(prefix_pdf), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(prefix_pdf, plot, width = width, height = height, limitsize = FALSE)
  png_path <- sub("\\.pdf$", ".png", prefix_pdf)
  ggplot2::ggsave(png_path, plot, width = width, height = height, dpi = 300, limitsize = FALSE)
}

txdb_cache_default <- function(gtf) {
  info <- file.info(gtf)
  tag <- paste0(gsub("[^A-Za-z0-9_.-]", "_", basename(gtf)), "_",
                format(as.numeric(info$size), scientific = FALSE, trim = TRUE), "_",
                format(as.numeric(info$mtime), scientific = FALSE, trim = TRUE))
  file.path(Sys.getenv("HOME"), ".cache", "atacseq", paste0("txdb_", tag, ".sqlite"))
}

te_cache_default <- function(te_path) {
  info <- file.info(te_path)
  tag <- paste0(gsub("[^A-Za-z0-9_.-]", "_", basename(te_path)), "_",
                format(as.numeric(info$size), scientific = FALSE, trim = TRUE), "_",
                format(as.numeric(info$mtime), scientific = FALSE, trim = TRUE))
  file.path(Sys.getenv("HOME"), ".cache", "atacseq", paste0("te_", tag, ".rds"))
}

load_packaged_txdb <- function(species) {
  pkg <- switch(
    species,
    hg38 = "TxDb.Hsapiens.UCSC.hg38.knownGene",
    mm10 = "TxDb.Mmusculus.UCSC.mm10.knownGene",
    NULL
  )
  if (is.null(pkg) || !requireNamespace(pkg, quietly = TRUE)) return(NULL)
  msg("load packaged TxDb: ", pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  get(pkg, envir = asNamespace(pkg))
}

load_txdb <- function(gtf, cache, species) {
  if (is.null(cache) || !nzchar(cache)) cache <- txdb_cache_default(gtf)
  if (file.exists(cache)) {
    msg("load TxDb cache: ", cache)
    return(AnnotationDbi::loadDb(cache))
  }
  packaged <- load_packaged_txdb(species)
  if (!is.null(packaged)) return(packaged)
  msg("build TxDb from GTF: ", gtf)
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf")
  dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
  AnnotationDbi::saveDb(txdb, cache)
  msg("saved TxDb cache: ", cache)
  txdb
}

read_peak_gr <- function(path) {
  gr <- rtracklayer::import(path)
  if (!inherits(gr, "GRanges")) stop("Cannot import peak file as GRanges: ", path)
  gr <- GenomeInfoDb::keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- gr[width(gr) > 0]
  if (!"name" %in% colnames(mcols(gr))) {
    mcols(gr)$name <- paste0("peak_", seq_along(gr))
  }
  mcols(gr)$peak_id <- ifelse(is.na(mcols(gr)$name) | mcols(gr)$name == "",
                              paste0("peak_", seq_along(gr)), as.character(mcols(gr)$name))
  gr
}

simplify_annotation <- function(x) {
  y <- rep("Intergenic", length(x))
  y[grepl("Promoter", x, ignore.case = TRUE)] <- "Promoter"
  y[grepl("Exon", x, ignore.case = TRUE)] <- "Exon"
  y[grepl("Intron", x, ignore.case = TRUE)] <- "Intron"
  y[grepl("UTR", x, ignore.case = TRUE)] <- "UTR"
  y[grepl("Downstream", x, ignore.case = TRUE)] <- "Downstream"
  y
}

pick_col <- function(df, patterns) {
  nm <- colnames(df)
  for (pat in patterns) {
    hit <- nm[grepl(pat, nm, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[[1]])
  }
  NULL
}

annotate_te <- function(peaks, te_path) {
  out <- data.frame(
    peak_id = mcols(peaks)$peak_id,
    is_TE = FALSE,
    te_class = NA_character_,
    te_family = NA_character_,
    te_name = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.null(te_path) || !nzchar(te_path) || !file.exists(te_path)) return(out)
  cache <- te_cache_default(te_path)
  if (file.exists(cache)) {
    msg("load TE annotation cache: ", cache)
    te <- tryCatch(readRDS(cache), error = function(e) NULL)
  } else {
    msg("load TE annotation: ", te_path)
    te <- tryCatch(rtracklayer::import(te_path), error = function(e) NULL)
    if (!is.null(te) && length(te) > 0) {
      te <- GenomeInfoDb::keepStandardChromosomes(te, pruning.mode = "coarse")
      dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
      saveRDS(te, cache)
      msg("saved TE annotation cache: ", cache)
    }
  }
  if (is.null(te) || length(te) == 0) return(out)
  hits <- GenomicRanges::findOverlaps(peaks, te, ignore.strand = TRUE)
  if (length(hits) == 0) return(out)

  te_df <- as.data.frame(mcols(te))
  cls_col <- pick_col(te_df, c("^repClass$", "class", "te_type", "gene_biotype"))
  fam_col <- pick_col(te_df, c("^repFamily$", "family"))
  name_col <- pick_col(te_df, c("^repName$", "repname", "name", "gene_id", "transcript_id"))
  value_or_na <- function(col, idx) {
    if (is.null(col)) return(rep(NA_character_, length(idx)))
    as.character(te_df[[col]][idx])
  }
  hit_df <- data.frame(
    peak_i = queryHits(hits),
    te_class = value_or_na(cls_col, subjectHits(hits)),
    te_family = value_or_na(fam_col, subjectHits(hits)),
    te_name = value_or_na(name_col, subjectHits(hits)),
    stringsAsFactors = FALSE
  )
  collapse_vals <- function(x) {
    x <- unique(x[!is.na(x) & x != ""])
    if (length(x) == 0) NA_character_ else paste(head(x, 5), collapse = ";")
  }
  collapsed <- hit_df |>
    dplyr::group_by(.data$peak_i) |>
    dplyr::summarise(
      te_class = collapse_vals(.data$te_class),
      te_family = collapse_vals(.data$te_family),
      te_name = collapse_vals(.data$te_name),
      .groups = "drop"
    )
  out$is_TE[collapsed$peak_i] <- TRUE
  out$te_class[collapsed$peak_i] <- collapsed$te_class
  out$te_family[collapsed$peak_i] <- collapsed$te_family
  out$te_name[collapsed$peak_i] <- collapsed$te_name
  out
}

plot_annotation <- function(tab, fig_dir, title) {
  anno_sum <- tab |>
    dplyr::count(.data$annotation_simple, name = "n") |>
    dplyr::mutate(pct = .data$n / sum(.data$n) * 100)
  p_bar <- ggplot(anno_sum, aes(x = reorder(annotation_simple, n), y = n, fill = annotation_simple)) +
    geom_col(width = 0.75) +
    coord_flip() +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = title, x = NULL, y = "Peak count")
  save_plot_both(file.path(fig_dir, "gene_annotation_bar.pdf"), p_bar, 7, 4.8)

  p_pie <- ggplot(anno_sum, aes(x = "", y = n, fill = annotation_simple)) +
    geom_col(width = 1) +
    coord_polar(theta = "y") +
    theme_void(base_size = 12) +
    labs(title = title, fill = "Annotation")
  save_plot_both(file.path(fig_dir, "gene_annotation_pie.pdf"), p_pie, 6, 5.2)

  dist_tab <- tab |> dplyr::filter(!is.na(.data$distanceToTSS))
  if (nrow(dist_tab) > 0) {
    p_dist <- ggplot(dist_tab, aes(x = distanceToTSS)) +
      geom_histogram(bins = 80, fill = "#4C78A8", color = "white", linewidth = 0.1) +
      theme_bw(base_size = 12) +
      labs(title = "Peak distance to nearest TSS", x = "Distance to TSS (bp)", y = "Peak count")
    save_plot_both(file.path(fig_dir, "distance_to_tss_histogram.pdf"), p_dist, 7, 4.5)
  }
}

plot_te <- function(tab, fig_dir) {
  status <- tab |>
    dplyr::mutate(status = ifelse(.data$is_TE, "TE-overlap", "non-TE")) |>
    dplyr::count(.data$status, name = "n")
  p_status <- ggplot(status, aes(x = status, y = n, fill = status)) +
    geom_col(width = 0.68) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none") +
    labs(title = "TE overlap status", x = NULL, y = "Peak count")
  save_plot_both(file.path(fig_dir, "TE_overlap_status_bar.pdf"), p_status, 6, 4.5)

  for (col in c("te_class", "te_family", "te_name")) {
    sub <- tab |>
      dplyr::filter(.data$is_TE) |>
      tidyr::separate_rows(dplyr::all_of(col), sep = ";") |>
      dplyr::mutate(group = dplyr::coalesce(.data[[col]], "Unannotated"),
                    group = ifelse(.data$group == "", "Unannotated", .data$group)) |>
      dplyr::count(.data$group, name = "n") |>
      dplyr::arrange(dplyr::desc(.data$n)) |>
      dplyr::slice_head(n = 20)
    if (nrow(sub) == 0) next
    p <- ggplot(sub, aes(x = reorder(group, n), y = n, fill = group)) +
      geom_col(width = 0.72) +
      coord_flip() +
      theme_bw(base_size = 12) +
      theme(legend.position = "none") +
      labs(title = paste("TE", gsub("^te_", "", col), "distribution"),
           x = NULL, y = "TE-overlapping peak count")
    save_plot_both(file.path(fig_dir, paste0("TE_", gsub("^te_", "", col), "_bar.pdf")), p, 7, 5)
  }
}

main <- function() {
  opt <- parse_args(args_raw)
  dir.create(dirname(opt$out_prefix), showWarnings = FALSE, recursive = TRUE)
  fig_dir <- file.path(dirname(opt$out_prefix), "plots")
  table_dir <- file.path(dirname(opt$out_prefix), "tables")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

  peaks <- read_peak_gr(opt$input)
  msg("peaks: ", length(peaks))
  txdb <- load_txdb(opt$gtf, opt$txdb_cache, opt$species)
  anno <- ChIPseeker::annotatePeak(
    peaks,
    TxDb = txdb,
    tssRegion = c(-opt$promoter_up, opt$promoter_down),
    verbose = FALSE
  )
  anno_df <- as.data.frame(anno)
  raw_path <- file.path(table_dir, paste0(basename(opt$out_prefix), ".ChIPseeker_raw.tsv"))
  write.table(anno_df, raw_path, sep = "\t", quote = FALSE, row.names = FALSE)

  te_df <- annotate_te(peaks, opt$te_bed)
  tab <- data.frame(
    peak_id = mcols(peaks)$peak_id,
    chrom = as.character(seqnames(peaks)),
    start = start(peaks) - 1L,
    end = end(peaks),
    width = width(peaks),
    annotation = anno_df$annotation,
    annotation_simple = simplify_annotation(anno_df$annotation),
    geneId = if ("geneId" %in% colnames(anno_df)) anno_df$geneId else NA_character_,
    SYMBOL = if ("SYMBOL" %in% colnames(anno_df)) anno_df$SYMBOL else NA_character_,
    distanceToTSS = if ("distanceToTSS" %in% colnames(anno_df)) anno_df$distanceToTSS else NA_integer_,
    stringsAsFactors = FALSE
  )
  tab <- dplyr::left_join(tab, te_df, by = "peak_id")
  out_table <- file.path(table_dir, paste0(basename(opt$out_prefix), ".peak_annotation.tsv"))
  write.table(tab, out_table, sep = "\t", quote = FALSE, row.names = FALSE)

  summary <- tab |>
    dplyr::count(.data$annotation_simple, .data$is_TE, name = "n") |>
    dplyr::mutate(frac = .data$n / sum(.data$n))
  write.table(summary, file.path(table_dir, paste0(basename(opt$out_prefix), ".annotation_summary.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  if (!opt$no_plots) {
    plot_annotation(tab, fig_dir, paste0("Peak annotation: ", basename(opt$out_prefix)))
    plot_te(tab, fig_dir)
  }
  msg("done: ", dirname(opt$out_prefix))
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})
