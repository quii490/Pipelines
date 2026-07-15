#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
  library(ggplot2)
})

parse_args <- function(argv) {
  args <- list(
    input = NULL,
    ref = NULL,
    out_prefix = NULL,
    gtf = NULL,
    te_anno = NULL,
    txdb_cache = "",
    promoter_up = 2000L,
    promoter_down = 500L,
    te_min_length = 50L,
    line_min_length = 800L,
    te_min_overlap = 50L,
    relaxed_window = 500L,
    chippeakanno_output = "overlapping",
    chippeakanno_select = "all",
    chippeakanno_maxgap = NA_integer_,
    chippeakanno_binding_region = "none",
    include_te_class = "",
    exclude_te_class_regex = "",
    exclude_te_family_regex = "",
    exclude_te_repname_regex = "",
    top_n_classes = 12L,
    no_plots = FALSE
  )
  i <- 1L
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (key %in% c("-h", "--help")) {
      cat(help_text())
      quit(status = 0)
    }
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- sub("^--", "", key)
    name <- gsub("-", "_", name)
    if (name == "no_plots") {
      args$no_plots <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(argv)) stop("Missing value for ", key)
    val <- argv[[i + 1L]]
    if (!name %in% names(args)) stop("Unknown option: ", key)
    args[[name]] <- val
    i <- i + 2L
  }
  int_fields <- c("promoter_up", "promoter_down", "te_min_length", "line_min_length", "te_min_overlap", "relaxed_window", "top_n_classes")
  for (x in int_fields) args[[x]] <- as.integer(args[[x]])
  if (is.na(args$chippeakanno_maxgap)) args$chippeakanno_maxgap <- args$relaxed_window
  args$chippeakanno_maxgap <- as.integer(args$chippeakanno_maxgap)
  required <- c("input", "ref", "out_prefix")
  missing <- required[vapply(args[required], is.null, logical(1))]
  if (length(missing)) stop("Missing required option(s): --", paste(missing, collapse = ", --"))
  if (!args$ref %in% c("hg38", "mm39")) stop("--ref must be hg38 or mm39")
  args
}

help_text <- function() {
  paste0(
"Usage:\n",
"  cutrun_annotate_peaks --input peaks.broadPeak --ref hg38 --out-prefix out/sample [options]\n\n",
"Required:\n",
"  --input FILE                  Peak file: BED/narrowPeak/broadPeak. BAM/bigWig are not needed.\n",
"  --ref hg38|mm39               Reference genome.\n",
"  --out-prefix PREFIX           Output prefix, e.g. /path/anno/sample.\n\n",
"Annotation resources:\n",
"  --gtf FILE                    Gene GTF. Default uses local hg38/mm39 GENCODE.\n",
"  --te-anno FILE                TE annotation TSV. Default uses CUTnRUN resources.\n\n",
"Gene structure parameters:\n",
"  --promoter-up INT             Promoter upstream bp. Default 2000.\n",
"  --promoter-down INT           Promoter downstream bp. Default 500.\n\n",
"TE stringent parameters:\n",
"  --te-min-length INT           Keep reference TE length >= INT. Default 50.\n",
"  --line-min-length INT         Keep LINE length >= INT. Default 800.\n",
"  --te-min-overlap INT          Peak-TE overlap required for stringent TE. Default 50.\n",
"  --relaxed-window INT          ChIPpeakAnno/relaxed TE window. Default 500.\n\n",
"Manual TE curation knobs:\n",
"  --include-te-class LIST       Keep only classes, comma-separated, e.g. LINE,SINE,LTR.\n",
"  --exclude-te-class-regex REG  Drop TE classes matching regex.\n",
"  --exclude-te-family-regex REG Drop TE families matching regex.\n",
"  --exclude-te-repname-regex REG Drop TE repNames matching regex, e.g. Simple|Low_complexity.\n\n",
"ChIPpeakAnno knobs:\n",
"  --chippeakanno-output STR     annotatePeakInBatch output. Default overlapping.\n",
"  --chippeakanno-select STR     annotatePeakInBatch select. Default all.\n",
"  --chippeakanno-maxgap INT     maxgap passed to ChIPpeakAnno. Default = relaxed-window.\n",
"  --chippeakanno-binding-region A,B  Optional bindingRegion, e.g. -1,1; default none.\n\n",
"Plotting:\n",
"  --top-n-classes INT           Number of TE classes shown in bar plot. Default 12.\n",
"  --no-plots                    Do not generate pdf/png plots.\n"
  )
}

message_run <- function(...) message(format(Sys.time(), "%F %T"), " | ", ...)

open_text <- function(path) {
  if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
}

default_gtf <- function(ref) {
  if (ref == "hg38") return("/path/to/reference/human_hg38/gencode.v49.primary_assembly.annotation.gtf")
  if (ref == "mm39") return("/path/to/reference/mouse_mm39/gencode.vM38.primary_assembly.annotation.gtf")
}

default_te <- function(ref) {
  file.path("/path/to/CUTnRUN/resources/CUTRUN_analysis/anno", paste0("te_anno_", ref, ".tsv"))
}

normalize_chrom <- function(x, action) {
  if (action == "add_chr") {
    x <- ifelse(startsWith(x, "chr"), x, ifelse(x %in% c("M", "MT"), "chrM", paste0("chr", x)))
  } else if (action == "remove_chr") {
    x <- ifelse(x == "chrM", "MT", sub("^chr", "", x))
  }
  x
}

first_peak_style <- function(path) {
  con <- open_text(path); on.exit(close(con))
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (!length(line)) return("chr")
    if (!grepl("^#|^$", line)) {
      chrom <- strsplit(line, "\\t")[[1]][1]
      return(ifelse(startsWith(chrom, "chr"), "chr", "nochr"))
    }
  }
}

te_ref_style <- function(path) {
  tab <- read.delim(path, nrows = 5, check.names = FALSE)
  if (!"GeneID" %in% names(tab)) stop("TE annotation must contain GeneID column: ", path)
  chrom <- sub(":.*$", "", tab$GeneID[which(!is.na(tab$GeneID))[1]])
  ifelse(startsWith(chrom, "chr"), "chr", "nochr")
}

style_action <- function(peak_style, ref_style) {
  if (identical(peak_style, ref_style)) return("keep")
  if (peak_style == "nochr" && ref_style == "chr") return("add_chr")
  if (peak_style == "chr" && ref_style == "nochr") return("remove_chr")
  "keep"
}

read_peaks <- function(path, chrom_action) {
  dat <- read.delim(path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(dat) < 3) stop("Peak file needs at least 3 BED columns")
  start <- as.integer(dat[[2]])
  end <- as.integer(dat[[3]])
  keep <- !is.na(start) & !is.na(end) & end > start
  dat <- dat[keep, , drop = FALSE]
  start <- start[keep]
  end <- end[keep]
  chrom <- normalize_chrom(as.character(dat[[1]]), chrom_action)
  nm <- if (ncol(dat) >= 4) as.character(dat[[4]]) else rep("", nrow(dat))
  empty <- is.na(nm) | nm == ""
  nm[empty] <- paste0("peak_", which(empty))
  nm <- make.unique(nm)
  score <- if (ncol(dat) >= 5) as.character(dat[[5]]) else "0"
  strand <- if (ncol(dat) >= 6) as.character(dat[[6]]) else "*"
  strand[!strand %in% c("+", "-")] <- "*"
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start + 1L, end = end), strand = strand)
  mcols(gr)$peak_id <- nm
  mcols(gr)$score <- score
  names(gr) <- nm
  gr
}

load_te <- function(path, chrom_action, args) {
  cache_dir <- file.path(Sys.getenv("HOME"), ".cache", "cutrun")
  dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
  info <- file.info(path)
  cache_tag <- paste0(gsub("[^A-Za-z0-9_.-]", "_", basename(path)), "_",
                      format(info$size, scientific = FALSE), "_", as.numeric(info$mtime))
  cache_path <- file.path(cache_dir, paste0("te_annotation_", cache_tag, ".rds"))
  if (file.exists(cache_path)) {
    te <- readRDS(cache_path)
  } else {
    te <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    saveRDS(te, cache_path)
  }
  required <- c("GeneID", "repName", "Class", "Family")
  missing <- setdiff(required, names(te))
  if (length(missing)) stop("TE annotation missing column(s): ", paste(missing, collapse = ", "))
  m <- regexec("^([^:]+):(\\d+)-(\\d+)$", te$GeneID)
  parts <- regmatches(te$GeneID, m)
  ok <- lengths(parts) == 4L
  te <- te[ok, , drop = FALSE]
  parts <- parts[ok]
  chrom <- normalize_chrom(vapply(parts, `[`, character(1), 2), chrom_action)
  start <- as.integer(vapply(parts, `[`, character(1), 3))
  end <- as.integer(vapply(parts, `[`, character(1), 4))
  len <- end - start
  te$repName[is.na(te$repName) | te$repName == ""] <- "NA"
  te$Class[is.na(te$Class) | te$Class == ""] <- "NA"
  te$Family[is.na(te$Family) | te$Family == ""] <- "NA"
  keep <- end > start
  if (nzchar(args$include_te_class)) {
    allowed <- trimws(strsplit(args$include_te_class, ",", fixed = TRUE)[[1]])
    keep <- keep & te$Class %in% allowed
  }
  if (nzchar(args$exclude_te_class_regex)) keep <- keep & !grepl(args$exclude_te_class_regex, te$Class)
  if (nzchar(args$exclude_te_family_regex)) keep <- keep & !grepl(args$exclude_te_family_regex, te$Family)
  if (nzchar(args$exclude_te_repname_regex)) keep <- keep & !grepl(args$exclude_te_repname_regex, te$repName)
  te <- te[keep, , drop = FALSE]
  chrom <- chrom[keep]; start <- start[keep]; end <- end[keep]; len <- len[keep]
  gr <- GRanges(seqnames = chrom, ranges = IRanges(start = start + 1L, end = end))
  mcols(gr)$te_id <- te$GeneID
  mcols(gr)$repName <- te$repName
  mcols(gr)$Class <- te$Class
  mcols(gr)$Family <- te$Family
  mcols(gr)$milliDiv <- if ("milliDiv" %in% names(te)) te$milliDiv else NA
  mcols(gr)$te_length <- len
  names(gr) <- make.unique(te$GeneID)
  gr
}

make_gene_features <- function(gtf, chrom_action, promoter_up, promoter_down) {
  tx <- rtracklayer::import(gtf)
  seqlevels(tx) <- normalize_chrom(seqlevels(tx), chrom_action)
  genes <- tx[tx$type == "gene"]
  exons <- tx[tx$type == "exon"]
  getcol <- function(gr, col, fallback = "NA") {
    if (col %in% names(mcols(gr))) as.character(mcols(gr)[[col]]) else rep(fallback, length(gr))
  }
  gene_id <- getcol(genes, "gene_id")
  gene_name <- if ("gene_name" %in% names(mcols(genes))) getcol(genes, "gene_name") else gene_id
  gene_type <- if ("gene_type" %in% names(mcols(genes))) getcol(genes, "gene_type") else getcol(genes, "gene_biotype")
  tss <- ifelse(as.character(strand(genes)) == "-", end(genes), start(genes))
  pstart <- ifelse(as.character(strand(genes)) == "-", pmax(1L, tss - promoter_down), pmax(1L, tss - promoter_up))
  pend <- ifelse(as.character(strand(genes)) == "-", tss + promoter_up, tss + promoter_down)
  promoters <- GRanges(seqnames = seqnames(genes), ranges = IRanges(pstart, pend), strand = strand(genes))
  mcols(promoters)$feature <- "promoter"
  mcols(promoters)$gene_id <- gene_id
  mcols(promoters)$gene_name <- gene_name
  mcols(promoters)$gene_type <- gene_type
  gene_body <- genes
  mcols(gene_body)$feature <- "gene_body"
  mcols(gene_body)$gene_id <- gene_id
  mcols(gene_body)$gene_name <- gene_name
  mcols(gene_body)$gene_type <- gene_type
  exon_gene_id <- getcol(exons, "gene_id")
  exon_gene_name <- if ("gene_name" %in% names(mcols(exons))) getcol(exons, "gene_name") else exon_gene_id
  exon_gene_type <- if ("gene_type" %in% names(mcols(exons))) getcol(exons, "gene_type") else getcol(exons, "gene_biotype")
  mcols(exons)$feature <- "exon"
  mcols(exons)$gene_id <- exon_gene_id
  mcols(exons)$gene_name <- exon_gene_name
  mcols(exons)$gene_type <- exon_gene_type
  c(promoters[, c("feature", "gene_id", "gene_name", "gene_type")], gene_body[, c("feature", "gene_id", "gene_name", "gene_type")], exons[, c("feature", "gene_id", "gene_name", "gene_type")])
}

write_gene_annotation <- function(peaks, features, out_prefix) {
  hits <- findOverlaps(peaks, features, ignore.strand = TRUE)
  raw_path <- paste0(out_prefix, ".gene_intersections.tsv")
  if (length(hits)) {
    q <- queryHits(hits); s <- subjectHits(hits)
    ov <- width(pintersect(peaks[q], features[s], ignore.strand = TRUE))
    raw <- data.frame(
      peak_id = mcols(peaks)$peak_id[q], peak_chrom = as.character(seqnames(peaks[q])), peak_start = start(peaks[q]) - 1L, peak_end = end(peaks[q]),
      feature = mcols(features)$feature[s], gene_id = mcols(features)$gene_id[s], gene_name = mcols(features)$gene_name[s], gene_type = mcols(features)$gene_type[s], overlap_bp = ov,
      stringsAsFactors = FALSE
    )
    write.table(raw, raw_path, sep = "\t", quote = FALSE, row.names = FALSE)
  } else {
    raw <- data.frame()
    write.table(data.frame(), raw_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  priority <- c(promoter = 1L, exon = 2L, gene_body = 3L)
  best <- data.frame(
    peak_id = mcols(peaks)$peak_id, chrom = as.character(seqnames(peaks)), start = start(peaks) - 1L, end = end(peaks),
    gene_feature = "intergenic", gene_id = "NA", gene_name = "NA", gene_type = "NA", overlap_bp = 0L,
    stringsAsFactors = FALSE
  )
  if (nrow(raw)) {
    raw$priority <- priority[raw$feature]
    raw$priority[is.na(raw$priority)] <- 99L
    raw <- raw[order(raw$peak_id, raw$priority, -raw$overlap_bp), ]
    raw_best <- raw[!duplicated(raw$peak_id), ]
    idx <- match(raw_best$peak_id, best$peak_id)
    best$gene_feature[idx] <- raw_best$feature
    best$gene_id[idx] <- raw_best$gene_id
    best$gene_name[idx] <- raw_best$gene_name
    best$gene_type[idx] <- raw_best$gene_type
    best$overlap_bp[idx] <- raw_best$overlap_bp
  }
  out <- paste0(out_prefix, ".gene_structure.tsv")
  write.table(best, out, sep = "\t", quote = FALSE, row.names = FALSE)
  table(best$gene_feature)
}

write_chippeakanno_raw <- function(peaks, te_all, out_prefix, args) {
  if (!requireNamespace("ChIPpeakAnno", quietly = TRUE)) {
    stop("ChIPpeakAnno is not installed in this R library. Install with: conda install -c conda-forge -c bioconda bioconductor-chippeakanno")
  }
  binding <- NULL
  if (!identical(tolower(args$chippeakanno_binding_region), "none")) {
    vals <- as.integer(strsplit(args$chippeakanno_binding_region, ",", fixed = TRUE)[[1]])
    if (length(vals) != 2L || any(is.na(vals))) stop("--chippeakanno-binding-region must look like -1,1 or none")
    binding <- vals
  }
  raw <- ChIPpeakAnno::annotatePeakInBatch(
    myPeakList = peaks,
    AnnotationData = te_all,
    output = args$chippeakanno_output,
    maxgap = args$chippeakanno_maxgap,
    select = args$chippeakanno_select,
    ignore.strand = TRUE,
    bindingRegion = binding
  )
  out <- paste0(out_prefix, ".ChIPpeakAnno_TE_relaxed_raw.tsv")
  write.table(as.data.frame(raw), out, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}

best_te_table <- function(peaks, te, hits, stringent, min_overlap = 0L) {
  base <- data.frame(
    peak_id = mcols(peaks)$peak_id,
    chrom = as.character(seqnames(peaks)),
    start = start(peaks) - 1L,
    end = end(peaks),
    peak_length = width(peaks),
    stringsAsFactors = FALSE
  )
  if (!length(hits)) {
    base$is_TE <- FALSE; base$repName <- "NA"; base$Class <- "NA"; base$Family <- "NA"; base$milliDiv <- "NA"; base$te_length <- "NA"; base$te_id <- "NA"; base$overlap_bp <- 0L; base$distance_bp <- NA_integer_; base$peak_coverage <- 0
    return(base)
  }
  q <- queryHits(hits); s <- subjectHits(hits)
  ov <- width(pintersect(peaks[q], te[s], ignore.strand = TRUE))
  ov[is.na(ov)] <- 0L
  if (stringent) keep <- ov >= min_overlap else keep <- rep(TRUE, length(q))
  q <- q[keep]; s <- s[keep]; ov <- ov[keep]
  if (!length(q)) {
    base$is_TE <- FALSE; base$repName <- "NA"; base$Class <- "NA"; base$Family <- "NA"; base$milliDiv <- "NA"; base$te_length <- "NA"; base$te_id <- "NA"; base$overlap_bp <- 0L; base$distance_bp <- NA_integer_; base$peak_coverage <- 0
    return(base)
  }
  dist <- distance(peaks[q], te[s], ignore.strand = TRUE)
  dist[is.na(dist)] <- 0L
  cand <- data.frame(q = q, s = s, overlap_bp = ov, distance_bp = dist, stringsAsFactors = FALSE)
  cand <- cand[order(cand$q, cand$distance_bp, -cand$overlap_bp), ]
  cand <- cand[!duplicated(cand$q), ]
  base$is_TE <- FALSE; base$repName <- "NA"; base$Class <- "NA"; base$Family <- "NA"; base$milliDiv <- "NA"; base$te_length <- "NA"; base$te_id <- "NA"; base$overlap_bp <- 0L; base$distance_bp <- NA_integer_; base$peak_coverage <- 0
  idx <- cand$q; subj <- cand$s
  base$is_TE[idx] <- TRUE
  base$repName[idx] <- mcols(te)$repName[subj]
  base$Class[idx] <- mcols(te)$Class[subj]
  base$Family[idx] <- mcols(te)$Family[subj]
  base$milliDiv[idx] <- as.character(mcols(te)$milliDiv[subj])
  base$te_length[idx] <- as.character(mcols(te)$te_length[subj])
  base$te_id[idx] <- mcols(te)$te_id[subj]
  base$overlap_bp[idx] <- cand$overlap_bp
  base$distance_bp[idx] <- cand$distance_bp
  base$peak_coverage[idx] <- cand$overlap_bp / base$peak_length[idx]
  base
}

write_te_annotation <- function(peaks, te_all, out_prefix, args) {
  strict_keep <- mcols(te_all)$te_length >= args$te_min_length
  strict_keep[mcols(te_all)$Class == "LINE"] <- mcols(te_all)$te_length[mcols(te_all)$Class == "LINE"] >= args$line_min_length
  te_strict <- te_all[strict_keep]
  relaxed_hits <- findOverlaps(peaks, te_all, maxgap = args$relaxed_window, minoverlap = 0L, ignore.strand = TRUE)
  strict_hits <- findOverlaps(peaks, te_strict, minoverlap = args$te_min_overlap, ignore.strand = TRUE)
  relaxed <- best_te_table(peaks, te_all, relaxed_hits, stringent = FALSE)
  stringent <- best_te_table(peaks, te_strict, strict_hits, stringent = TRUE, min_overlap = args$te_min_overlap)
  names(relaxed)[names(relaxed) == "is_TE"] <- "is_TE_relaxed"
  names(stringent)[names(stringent) == "is_TE"] <- "is_TE_stringent"
  write.table(relaxed, paste0(out_prefix, ".TE_relaxed_window.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(stringent, paste0(out_prefix, ".TE_stringent.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  list(relaxed = relaxed, stringent = stringent, te_strict = te_strict)
}

write_summary <- function(out_prefix, n_peaks, gene_counts, te_results, te_all, args) {
  relaxed <- te_results$relaxed
  stringent <- te_results$stringent
  te_strict <- te_results$te_strict
  strict_classes <- sort(unique(c(as.character(mcols(te_strict)$Class), stringent$Class[stringent$is_TE_stringent])))
  strict_families <- sort(unique(c(as.character(mcols(te_strict)$Family), stringent$Family[stringent$is_TE_stringent])))
  total_strict_bp <- sum(as.numeric(mcols(te_strict)$te_length), na.rm = TRUE)
  total_strict_peak_n <- sum(stringent$is_TE_stringent)
  rows <- list()
  add <- function(section, class, count, fraction, extra = "NA") {
    rows[[length(rows) + 1L]] <<- data.frame(section = section, class = class, count = count, fraction = fraction, extra = extra, stringsAsFactors = FALSE)
  }
  add_te_level <- function(section, label, ref_values, peak_values) {
    for (x in label) {
      peak_n <- sum(stringent$is_TE_stringent & peak_values == x, na.rm = TRUE)
      genome_bp <- sum(as.numeric(mcols(te_strict)$te_length[ref_values == x]), na.rm = TRUE)
      genome_frac <- ifelse(total_strict_bp > 0, genome_bp / total_strict_bp, NA_real_)
      peak_frac <- ifelse(total_strict_peak_n > 0, peak_n / total_strict_peak_n, 0)
      enrich <- ifelse(!is.na(genome_frac) && genome_frac > 0, peak_frac / genome_frac, NA_real_)
      ovbp <- sum(stringent$overlap_bp[stringent$is_TE_stringent & peak_values == x], na.rm = TRUE)
      add(section, x, peak_n, peak_frac, paste0("genome_bp_fraction=", signif(genome_frac, 6), ";peak_vs_genome_enrichment=", signif(enrich, 6), ";overlap_bp=", ovbp, ";genome_bp=", genome_bp))
    }
  }
  for (nm in sort(names(gene_counts))) add("gene_structure", nm, as.integer(gene_counts[[nm]]), as.integer(gene_counts[[nm]]) / n_peaks)
  add(paste0("TE_relaxed_", args$relaxed_window, "bp"), "ALL", sum(relaxed$is_TE_relaxed), sum(relaxed$is_TE_relaxed) / n_peaks, paste0("ChIPpeakAnno_maxgap=", args$chippeakanno_maxgap))
  add("TE_stringent", "ALL", sum(stringent$is_TE_stringent), sum(stringent$is_TE_stringent) / n_peaks, paste0("overlap>=", args$te_min_overlap, "bp;TE>=", args$te_min_length, "bp;LINE>=", args$line_min_length, "bp"))
  add_te_level("TE_stringent_class", strict_classes, as.character(mcols(te_strict)$Class), stringent$Class)
  add_te_level("TE_stringent_family", strict_families, as.character(mcols(te_strict)$Family), stringent$Family)
  summary <- do.call(rbind, rows)
  write.table(summary, paste0(out_prefix, ".summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  summary
}

plot_outputs <- function(out_prefix, summary, te_results, args) {
  plot_dir <- paste0(out_prefix, ".plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  save_plot <- function(p, name, w = 6, h = 4) {
    ggsave(file.path(plot_dir, paste0(name, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
    ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = w, height = h, dpi = 180)
  }
  fmt_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  fmt_pct <- function(x, total) {
    if (length(total) == 1L) {
      if (is.na(total) || total <= 0) return(rep("NA", length(x)))
      return(sprintf("%.1f%%", 100 * x / total))
    }
    out <- rep("NA", length(x))
    ok <- !is.na(total) & total > 0
    out[ok] <- sprintf("%.1f%%", 100 * x[ok] / total[ok])
    out
  }
  count_label <- function(x, total) paste0(fmt_count(x), " / ", fmt_count(total), " (", fmt_pct(x, total), ")")
  extra_num <- function(extra, key) {
    pat <- paste0(".*", key, "=([^;]+).*")
    out <- suppressWarnings(as.numeric(ifelse(grepl(paste0(key, "="), extra), sub(pat, "\\1", extra), NA)))
    out
  }
  plot_te_background <- function(section_name, file_prefix, level_label, top_n = args$top_n_classes) {
    dat <- subset(summary, section == section_name & count > 0)
    if (!nrow(dat)) return(invisible(NULL))
    dat$count <- as.numeric(dat$count)
    dat$peak_fraction <- as.numeric(dat$fraction)
    dat$genome_fraction <- extra_num(dat$extra, "genome_bp_fraction")
    dat$enrichment <- extra_num(dat$extra, "peak_vs_genome_enrichment")
    dat <- dat[!is.na(dat$genome_fraction) & dat$genome_fraction > 0, , drop = FALSE]
    if (!nrow(dat)) return(invisible(NULL))
    dat <- dat[order(-dat$count), ]
    if (nrow(dat) > top_n) {
      top <- dat[seq_len(top_n), ]
      other <- data.frame(
        section = section_name,
        class = "Other",
        count = sum(dat$count[-seq_len(top_n)]),
        fraction = sum(dat$peak_fraction[-seq_len(top_n)]),
        extra = "collapsed",
        peak_fraction = sum(dat$peak_fraction[-seq_len(top_n)]),
        genome_fraction = sum(dat$genome_fraction[-seq_len(top_n)]),
        enrichment = NA_real_,
        stringsAsFactors = FALSE
      )
      other$enrichment <- ifelse(other$genome_fraction > 0, other$peak_fraction / other$genome_fraction, NA_real_)
      dat <- rbind(top, other)
    }
    dat$class <- factor(dat$class, levels = rev(dat$class))
    long <- rbind(
      data.frame(class = dat$class, source = "Stringent TE peaks", fraction = dat$peak_fraction),
      data.frame(class = dat$class, source = "Reference TE genome bp", fraction = dat$genome_fraction)
    )
    p_frac <- ggplot(long, aes(x = class, y = fraction, fill = source)) +
      geom_col(position = position_dodge(width = 0.75), width = 0.68) +
      geom_text(aes(label = sprintf("%.1f%%", 100 * fraction)), position = position_dodge(width = 0.75), hjust = -0.05, size = 2.8) +
      coord_flip(clip = "off") +
      expand_limits(y = max(long$fraction, na.rm = TRUE) * 1.25) +
      theme_bw(base_size = 12) +
      theme(plot.margin = margin(5.5, 38, 5.5, 5.5)) +
      labs(x = NULL, y = "Fraction", fill = NULL, title = paste0("TE ", level_label, ": peaks vs genome background"), subtitle = "Peak fraction is among stringent TE peaks; genome fraction is bp share of filtered reference TE")
    save_plot(p_frac, paste0(file_prefix, "_peak_vs_genome_fraction_bar"), 8.6, 5)

    dat$enrich_label <- ifelse(is.na(dat$enrichment), "NA", sprintf("%.2fx", dat$enrichment))
    p_enrich <- ggplot(dat, aes(x = class, y = enrichment, fill = class)) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "grey35") +
      geom_col(width = 0.72) +
      geom_text(aes(label = enrich_label), hjust = -0.05, size = 3) +
      coord_flip(clip = "off") +
      expand_limits(y = max(dat$enrichment, na.rm = TRUE) * 1.25) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", plot.margin = margin(5.5, 35, 5.5, 5.5)) +
      labs(x = NULL, y = "Peak fraction / genome bp fraction", title = paste0("TE ", level_label, " enrichment over genome background"), subtitle = "Dashed line = no enrichment relative to filtered reference TE bp distribution")
    save_plot(p_enrich, paste0(file_prefix, "_peak_vs_genome_enrichment_bar"), 8, 5)
    invisible(dat)
  }

  gene <- subset(summary, section == "gene_structure")
  if (nrow(gene)) {
    gene$count <- as.numeric(gene$count)
    gene_total <- sum(gene$count)
    gene <- gene[order(gene$count), ]
    gene$class <- factor(gene$class, levels = gene$class)
    gene$label <- count_label(gene$count, gene_total)
    gene_subtitle <- paste0("Best single gene feature per peak; promoter = TSS -", args$promoter_up, " bp to +", args$promoter_down, " bp; total peaks = ", fmt_count(gene_total))
    p <- ggplot(gene, aes(x = class, y = count, fill = class)) +
      geom_col(width = 0.72) +
      geom_text(aes(label = label), hjust = -0.05, size = 3.2) +
      coord_flip(clip = "off") +
      expand_limits(y = max(gene$count) * 1.28) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", plot.margin = margin(5.5, 35, 5.5, 5.5)) +
      labs(x = NULL, y = "Peaks", title = "Gene structure annotation", subtitle = gene_subtitle)
    save_plot(p, "gene_structure_best_bar", 7.8, 4.8)

    gene$pie_label <- paste0(gene$class, "\n", fmt_count(gene$count), " (", fmt_pct(gene$count, gene_total), ")")
    p2 <- ggplot(gene, aes(x = "", y = count, fill = class)) +
      geom_col(width = 1, color = "white") +
      geom_text(aes(label = pie_label), position = position_stack(vjust = 0.5), size = 3) +
      coord_polar(theta = "y") +
      theme_void(base_size = 12) +
      labs(fill = "Gene feature", title = "Gene structure distribution", subtitle = paste0("Total peaks = ", fmt_count(gene_total)))
    save_plot(p2, "gene_structure_pie", 6, 5.5)
  }

  relaxed_n <- sum(te_results$relaxed$is_TE_relaxed)
  strict_n <- sum(te_results$stringent$is_TE_stringent)
  total_n <- nrow(te_results$stringent)
  status <- data.frame(
    class = c("stringent_TE", "relaxed_only_TE", "non_TE"),
    count = c(strict_n, max(relaxed_n - strict_n, 0), max(total_n - relaxed_n, 0)),
    stringsAsFactors = FALSE
  )
  status$label <- paste0(status$class, "\n", count_label(status$count, total_n))
  p_status <- ggplot(status, aes(x = "", y = count, fill = class)) +
    geom_col(width = 1, color = "white") +
    geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
    coord_polar(theta = "y") +
    theme_void(base_size = 12) +
    labs(
      fill = "Peak class",
      title = "TE annotation status",
      subtitle = paste0("Disjoint display: stringent TE, relaxed-only TE, and non-TE; total peaks = ", fmt_count(total_n))
    )
  save_plot(p_status, "TE_status_pie", 6, 5.5)

  te_nested <- data.frame(
    level = c("All peaks", "Relaxed TE", "Stringent TE"),
    count = c(total_n, relaxed_n, strict_n),
    rule = c(
      "input peaks",
      paste0("peak within ", args$relaxed_window, " bp of any TE"),
      paste0("exact overlap >= ", args$te_min_overlap, " bp; TE length >= ", args$te_min_length, " bp; LINE >= ", args$line_min_length, " bp")
    ),
    stringsAsFactors = FALSE
  )
  te_nested$level <- factor(te_nested$level, levels = rev(te_nested$level))
  te_nested$label <- paste0(count_label(te_nested$count, total_n), "\n", te_nested$rule)
  p <- ggplot(te_nested, aes(x = level, y = count, fill = level)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = label), hjust = -0.03, size = 3.1, lineheight = 0.92) +
    coord_flip(clip = "off") +
    expand_limits(y = max(te_nested$count) * 1.52) +
    theme_bw(base_size = 12) +
    theme(legend.position = "none", plot.margin = margin(5.5, 95, 5.5, 5.5)) +
    labs(x = NULL, y = "Peaks", title = "Nested TE annotation", subtitle = "Stringent TE is a high-confidence subset of the broader TE-proximal peaks")
  save_plot(p, "TE_relaxed_stringent_nested_bar", 9.5, 4.7)

  cls <- subset(summary, section == "TE_stringent_class" & count > 0)
  if (nrow(cls)) {
    cls$count <- as.numeric(cls$count)
    cls <- cls[order(-cls$count), ]
    if (nrow(cls) > args$top_n_classes) {
      top <- cls[seq_len(args$top_n_classes), ]
      other <- data.frame(section = "TE_stringent_class", class = "Other", count = sum(cls$count[-seq_len(args$top_n_classes)]), fraction = sum(cls$fraction[-seq_len(args$top_n_classes)]), extra = "collapsed", stringsAsFactors = FALSE)
      cls <- rbind(top, other)
    }
    cls_total <- sum(cls$count)
    cls$label <- count_label(cls$count, cls_total)
    cls$class <- factor(cls$class, levels = rev(cls$class))
    p <- ggplot(cls, aes(x = class, y = count, fill = class)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = label), hjust = -0.05, size = 3.1) +
      coord_flip(clip = "off") +
      expand_limits(y = max(cls$count) * 1.28) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", plot.margin = margin(5.5, 35, 5.5, 5.5)) +
      labs(x = NULL, y = "Stringent TE peaks", title = "Stringent TE class distribution", subtitle = paste0("Denominator = stringent TE peaks shown; n = ", fmt_count(cls_total)))
    save_plot(p, "TE_stringent_class_bar", 7.6, 4.8)
    p2 <- ggplot(cls, aes(x = "", y = count, fill = class)) + geom_col(width = 1, color = "white") + coord_polar(theta = "y") + theme_void(base_size = 12) + labs(fill = "TE class", title = "Stringent TE classes", subtitle = paste0("n = ", fmt_count(cls_total)))
    save_plot(p2, "TE_stringent_class_pie", 5.5, 5.5)
  }

  fam <- subset(summary, section == "TE_stringent_family" & count > 0)
  if (nrow(fam)) {
    fam$count <- as.numeric(fam$count)
    fam <- fam[order(-fam$count), ]
    if (nrow(fam) > args$top_n_classes) {
      top <- fam[seq_len(args$top_n_classes), ]
      other <- data.frame(section = "TE_stringent_family", class = "Other", count = sum(fam$count[-seq_len(args$top_n_classes)]), fraction = sum(fam$fraction[-seq_len(args$top_n_classes)]), extra = "collapsed", stringsAsFactors = FALSE)
      fam <- rbind(top, other)
    }
    fam_total <- sum(fam$count)
    fam$label <- count_label(fam$count, fam_total)
    fam$class <- factor(fam$class, levels = rev(fam$class))
    p <- ggplot(fam, aes(x = class, y = count, fill = class)) +
      geom_col(width = 0.75) +
      geom_text(aes(label = label), hjust = -0.05, size = 3.1) +
      coord_flip(clip = "off") +
      expand_limits(y = max(fam$count) * 1.28) +
      theme_bw(base_size = 12) +
      theme(legend.position = "none", plot.margin = margin(5.5, 35, 5.5, 5.5)) +
      labs(x = NULL, y = "Stringent TE peaks", title = "Stringent TE family distribution", subtitle = paste0("Denominator = stringent TE peaks shown; n = ", fmt_count(fam_total)))
    save_plot(p, "TE_stringent_family_bar", 7.8, 5)
  }

  plot_te_background("TE_stringent_class", "TE_stringent_class", "class")
  plot_te_background("TE_stringent_family", "TE_stringent_family", "family")
  plot_dir
}

write_parameters <- function(args) {
  path <- paste0(args$out_prefix, ".parameters.txt")
  con <- file(path, "wt"); on.exit(close(con))
  for (nm in names(args)) writeLines(paste(nm, args[[nm]], sep = "="), con)
}

main <- function() {
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
  message_run("running ChIPpeakAnno relaxed annotation")
  raw_path <- write_chippeakanno_raw(peaks, te_all, args$out_prefix, args)
  message_run("ChIPpeakAnno raw output: ", raw_path)
  message_run("building gene structure features")
  gene_features <- make_gene_features(args$gtf, chrom_action, args$promoter_up, args$promoter_down)
  gene_counts <- write_gene_annotation(peaks, gene_features, args$out_prefix)
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

if (!exists("CUTRUN_ANNOTATE_NO_MAIN", envir = .GlobalEnv, inherits = FALSE)) {
  tryCatch(main(), error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(status = 1)
  })
}
