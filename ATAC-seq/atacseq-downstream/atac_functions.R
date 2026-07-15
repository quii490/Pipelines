suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(edgeR)
  library(GenomicRanges)
  library(ChIPseeker)
  library(ggrepel)
  library(rtracklayer)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# Resolve bundled resources from the checked-out pipeline rather than from a
# machine-specific historical path. Wrapper scripts set ATAC_PIPELINE_ROOT;
# the fallback keeps direct Rscript/source use functional during development.
atac_pipeline_root <- Sys.getenv('ATAC_PIPELINE_ROOT', unset = '')
if (!nzchar(atac_pipeline_root)) {
  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) '')
  if (nzchar(source_file)) {
    atac_pipeline_root <- dirname(dirname(normalizePath(source_file, mustWork = FALSE)))
  } else if (dir.exists(file.path(getwd(), 'atacseq-downstream'))) {
    atac_pipeline_root <- getwd()
  } else {
    atac_pipeline_root <- dirname(getwd())
  }
}
atac_pipeline_root <- normalizePath(atac_pipeline_root, mustWork = FALSE)
atac_conf_file <- function(name) file.path(atac_pipeline_root, 'conf', name)

log_msg <- function(...) {
  message('[atac_functions] ', paste(..., collapse = ''))
}

as_png_path <- function(path) sub('\\.pdf$', '.png', path, ignore.case = TRUE)

theme_atac <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 1, color = "grey35"),
      axis.title = element_text(color = "grey15"),
      axis.text = element_text(color = "grey20"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank()
    )
}

save_ggplot_both <- function(path_pdf, plot, width = 7, height = 5, dpi = 180) {
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  ggsave(path_pdf, plot, width = width, height = height)
  ggsave(as_png_path(path_pdf), plot, width = width, height = height, dpi = dpi)
}

save_device_both <- function(path_pdf, width = 7, height = 5, plot_fun, res = 180) {
  dir.create(dirname(path_pdf), recursive = TRUE, showWarnings = FALSE)
  pdf(path_pdf, width = width, height = height)
  plot_fun()
  dev.off()
  png(as_png_path(path_pdf), width = width, height = height, units = 'in', res = res)
  plot_fun()
  dev.off()
  invisible(NULL)
}

# ===== Output structure helpers =====
create_contrast_dirs <- function(croot) {
  dirs <- list(
    figures = file.path(croot, "figures"),
    tables  = file.path(croot, "tables"),
    beds    = file.path(croot, "beds")
  )
  for (d in c(croot, unname(dirs))) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  }
  dirs
}
fig_path <- function(dirs, sub, name) file.path(dirs, sub, name)
tbl_path <- function(dirs, sub, name) file.path(dirs, sub, name)
bed_path <- function(dirs, name) file.path(dirs, name)

contrast_label <- function(case, control) paste0(case, '_vs_', control)

contrast_root <- function(outdir, case, control) {
  file.path(outdir, 'contrasts', paste0(case, '_vs_', control))
}
results_dir <- function(croot) { d <- file.path(croot, 'results'); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
figures_dir <- function(croot, sub = NULL) { d <- file.path(croot, 'figures', if (is.null(sub)) '' else sub); dir.create(d, recursive = TRUE, showWarnings = FALSE); d }

make_sample_labels <- function(meta) {
  label <- dplyr::coalesce(meta$condition, rownames(meta))
  label[label == "" | is.na(label)] <- rownames(meta)[label == "" | is.na(label)]
  if ("replicate" %in% colnames(meta) && any(duplicated(label))) {
    rep <- dplyr::coalesce(meta$replicate, rownames(meta))
    label <- ifelse(is.na(rep) | rep == "" | rep == "NA", label, paste0(label, "_rep", rep))
  }
  if (any(duplicated(label))) {
    label <- paste(label, rownames(meta), sep = "_")
  }
  make.unique(label, sep = "_")
}

safe_repel_layer <- function(df, x_col, y_col, label_col, size = 2.5) {
  if (is.null(df) || nrow(df) == 0) return(NULL)
  keep <- is.finite(df[[x_col]]) & is.finite(df[[y_col]]) & !is.na(df[[label_col]]) & df[[label_col]] != ''
  df2 <- df[keep, , drop = FALSE]
  if (nrow(df2) == 0) return(NULL)
  ggrepel::geom_text_repel(
    data = df2,
    aes_string(x = x_col, y = y_col, label = label_col),
    inherit.aes = FALSE,
    size = size,
    color = "grey12",
    max.overlaps = Inf,
    box.padding = 0.25,
    point.padding = 0.08,
    min.segment.length = 0.08,
    segment.color = "grey70",
    segment.size = 0.18,
    segment.alpha = 0.45,
    force = 1.1,
    force_pull = 0.25,
    max.time = 2
  )
}

scaled_radial_score <- function(x, y) {
  scale_abs <- function(v) {
    lim <- suppressWarnings(stats::quantile(abs(v[is.finite(v)]), 0.99, na.rm = TRUE))
    if (!is.finite(lim) || lim <= 0) lim <- 1
    pmin(abs(v) / lim, 1.5)
  }
  sqrt(scale_abs(x)^2 + scale_abs(y)^2)
}

label_by_distance <- function(df, x_col, y_col, label_col = "label", label_source_col = NULL, top_n = 100) {
  if (is.null(df) || nrow(df) == 0) return(df)
  if (!x_col %in% colnames(df) || !y_col %in% colnames(df)) return(df)
  df[[label_col]] <- ""
  x <- dplyr::coalesce(df[[x_col]], 0)
  y <- dplyr::coalesce(df[[y_col]], 0)
  df$rank_score <- scaled_radial_score(x, y)
  ord <- order(-df$rank_score, na.last = TRUE)
  ord <- ord[is.finite(df$rank_score[ord])]
  top_idx <- head(ord, min(top_n, length(ord)))
  if (length(top_idx) == 0) return(df)
  if (is.null(label_source_col)) {
    labels <- make_region_labels(df)
  } else if (label_source_col %in% colnames(df)) {
    labels <- as.character(df[[label_source_col]])
  } else if ("feature_id" %in% colnames(df)) {
    labels <- as.character(df$feature_id)
  } else {
    labels <- rownames(df)
  }
  labels <- clean_te_label(labels)
  missing <- is.na(labels) | labels == ""
  if (any(missing)) {
    labels[missing] <- if ("feature_id" %in% colnames(df)) as.character(df$feature_id[missing]) else rownames(df)[missing]
  }
  df[[label_col]][top_idx] <- labels[top_idx]
  df
}

make_region_labels <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(character())
  te <- if ("te_name" %in% colnames(df)) clean_te_label(df$te_name) else rep("", nrow(df))
  gene <- if ("SYMBOL" %in% colnames(df)) as.character(df$SYMBOL) else rep("", nrow(df))
  fid <- if ("feature_id" %in% colnames(df)) as.character(df$feature_id) else rownames(df)
  te[te %in% c("NA", "NaN", "NULL")] <- ""
  gene[gene %in% c("NA", "NaN", "NULL")] <- ""
  fid[fid %in% c("NA", "NaN", "NULL")] <- ""
  out <- ifelse(!is.na(te) & te != "", te, ifelse(!is.na(gene) & gene != "", gene, fid))
  make.unique(out, sep = "_")
}

read_featurecounts_matrix <- function(count_file) {
  log_msg('reading count matrix: ', count_file)
  x <- read.delim(count_file, comment.char = '#', check.names = FALSE)
  stopifnot(ncol(x) >= 7)
  anno <- x[, 1:6, drop = FALSE]
  mat <- x[, 7:ncol(x), drop = FALSE]
  colnames(mat) <- basename(colnames(mat)) |> sub('\\.clean\\.bam$', '', x = _)
  rownames(mat) <- anno$Geneid
  rownames(anno) <- anno$Geneid
  list(anno = anno, mat = as.matrix(mat))
}

load_sample_meta <- function(sample_meta, colnames_mat) {
  log_msg('reading sample meta: ', sample_meta)
  meta <- read.csv(sample_meta, stringsAsFactors = FALSE)
  meta$sample <- as.character(meta$sample)
  stopifnot(all(c('sample', 'condition', 'replicate') %in% colnames(meta)))
  meta <- meta[match(as.character(colnames_mat), meta$sample), , drop = FALSE]
  rownames(meta) <- meta$sample
  meta$condition <- as.character(meta$condition)
  meta$replicate <- as.character(meta$replicate)
  if (!"sample_label" %in% colnames(meta)) {
    meta$sample_label <- make_sample_labels(meta)
  } else {
    meta$sample_label <- as.character(meta$sample_label)
    bad <- is.na(meta$sample_label) | meta$sample_label == ""
    if (any(bad)) meta$sample_label[bad] <- make_sample_labels(meta[bad, , drop = FALSE])
    meta$sample_label <- make.unique(meta$sample_label, sep = "_")
  }
  meta
}

load_contrasts <- function(contrast_file, meta) {
  if (is.null(contrast_file) || is.na(contrast_file) || contrast_file == '' || !file.exists(contrast_file)) {
    conds <- unique(meta$condition)
    conds <- conds[!is.na(conds) & conds != 'NA' & conds != '']
    if (length(conds) == 2) {
      return(data.frame(case = conds[2], control = conds[1], stringsAsFactors = FALSE))
    }
    return(data.frame(case = character(), control = character(), stringsAsFactors = FALSE))
  }
  x <- read.csv(contrast_file, stringsAsFactors = FALSE)
  if (all(c('case', 'control') %in% colnames(x))) {
    x <- x[, c('case', 'control'), drop = FALSE]
  } else if (all(c('group_col', 'case', 'control') %in% colnames(x))) {
    log_msg('legacy contrast file detected with group_col column; ignore group_col and continue')
    x <- x[, c('case', 'control'), drop = FALSE]
  } else {
    stop('contrast file must contain columns: case,control')
  }
  x <- x[!is.na(x$case) & !is.na(x$control) & x$case != '' & x$control != '', , drop = FALSE]
  x
}

calc_basic_qc <- function(mat, meta, outdir) {
  dir.create(file.path(outdir, 'figures', 'qc'), recursive = TRUE, showWarnings = FALSE)
  libdf <- data.frame(sample = colnames(mat), sample_label = meta$sample_label, library_size = colSums(mat), condition = meta$condition)
  p1 <- ggplot(libdf, aes(x = sample, y = library_size, fill = condition)) +
    geom_col(width = 0.72, color = "grey20", linewidth = 0.15) +
    scale_x_discrete(labels = setNames(libdf$sample_label, libdf$sample)) +
    coord_flip() + theme_atac() + labs(title = 'Library size', x = '', y = 'Counts')
  save_ggplot_both(file.path(outdir, 'figures', 'qc', 'library_size_barplot.pdf'), p1, width = 7, height = 4)

  lcpm <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  df <- as.data.frame(lcpm) |> tibble::rownames_to_column('feature_id') |>
    tidyr::pivot_longer(-feature_id, names_to = 'sample', values_to = 'logCPM')
  label_map <- setNames(meta$sample_label, rownames(meta))
  df$sample_label <- label_map[df$sample]
  p2 <- ggplot(df, aes(x = sample_label, y = logCPM, fill = sample_label)) +
    geom_boxplot(outlier.size = 0.15, linewidth = 0.25, show.legend = FALSE) +
    theme_atac() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = 'logCPM distribution', x = '', y = 'logCPM')
  save_ggplot_both(file.path(outdir, 'figures', 'qc', 'logCPM_boxplot.pdf'), p2, width = 8, height = 4)
}

plot_pca_cor <- function(mat, meta, outdir) {
  dir.create(file.path(outdir, 'figures', 'overview'), recursive = TRUE, showWarnings = FALSE)
  lcpm <- edgeR::cpm(mat, log = TRUE, prior.count = 1)
  row_var <- apply(lcpm, 1, function(z) stats::var(z, na.rm = TRUE))
  keep <- is.finite(row_var) & !is.na(row_var) & row_var > 0
  log_msg('PCA variable features kept: ', sum(keep), ' / ', length(keep))
  if (sum(keep) >= 2) {
    lcpm_use <- lcpm[keep, , drop = FALSE]
    pca <- prcomp(t(lcpm_use), scale. = TRUE)
    ve <- summary(pca)$importance[2, ]
    df <- data.frame(sample = rownames(pca$x), sample_label = meta$sample_label, PC1 = pca$x[, 1], PC2 = pca$x[, 2], condition = meta$condition)
    p <- ggplot(df, aes(PC1, PC2, color = condition, label = sample_label)) +
      geom_hline(yintercept = 0, color = "grey88", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey88", linewidth = 0.3) +
      geom_point(size = 3.2, alpha = 0.95) +
      ggrepel::geom_text_repel(size = 3.4, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
      theme_atac() +
      labs(title = 'PCA of ATAC-seq consensus peaks', x = sprintf('PC1 (%.1f%%)', 100 * ve[1]), y = sprintf('PC2 (%.1f%%)', 100 * ve[2]))
    save_ggplot_both(file.path(outdir, 'figures', 'overview', 'PCA.pdf'), p, width = 6, height = 5)
  } else {
    log_msg('skip PCA: <2 variable features after filtering')
  }

  cor_mat <- cor(lcpm, method = 'pearson')
  colnames(cor_mat) <- meta$sample_label
  rownames(cor_mat) <- meta$sample_label
  save_device_both(file.path(outdir, 'figures', 'overview', 'Correlation_heatmap.pdf'), width = 6, height = 5,
                   plot_fun = function() pheatmap(cor_mat, border_color = "white", fontsize = 9, main = "Sample correlation"))
}

plot_cross_contrast_overview <- function(result_list, outdir,
                                         padj_cutoff = 0.05,
                                         lfc_cutoff = 1,
                                         baseMean_min = 5,
                                         top_n = 100,
                                         annotation_mode = 'gene_te') {
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  if (length(result_list) == 0) return(invisible(NULL))

  summary_rows <- lapply(names(result_list), function(contrast) {
    x <- result_list[[contrast]]
    exploratory <- all(is.na(x$FDR)) && all(is.na(x$PValue))
    if (exploratory) {
      selected <- !is.na(x$baseMean) & x$baseMean >= baseMean_min &
        !is.na(x$logFC) & abs(x$logFC) >= lfc_cutoff
      mode <- "Exploratory (no replicate)"
    } else {
      selected <- !is.na(x$FDR) & x$FDR < padj_cutoff &
        !is.na(x$logFC) & abs(x$logFC) >= lfc_cutoff
      mode <- "Formal (FDR)"
    }
    data.frame(
      contrast = contrast,
      mode = mode,
      Up = sum(selected & x$logFC > 0, na.rm = TRUE),
      Down = sum(selected & x$logFC < 0, na.rm = TRUE),
      tested = nrow(x),
      stringsAsFactors = FALSE
    )
  })
  summary_df <- dplyr::bind_rows(summary_rows)
  result_dir <- file.path(outdir, "results", "overview")
  figure_dir <- file.path(outdir, "figures", "overview")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(summary_df, file.path(result_dir, "contrast_change_summary.csv"),
            row.names = FALSE, quote = FALSE)

  summary_long <- summary_df |>
    tidyr::pivot_longer(c("Up", "Down"), names_to = "direction", values_to = "regions")
  subtitle <- if (any(grepl("^Exploratory", summary_df$mode))) {
    sprintf("Exploratory: |log2FC| >= %s and baseMean >= %s; replicates are required for FDR",
            lfc_cutoff, baseMean_min)
  } else {
    sprintf("FDR < %s and |log2FC| >= %s", padj_cutoff, lfc_cutoff)
  }
  p_summary <- ggplot(summary_long, aes(x = reorder(contrast, regions), y = regions, fill = direction)) +
    geom_col(position = "dodge", width = 0.72) +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#d73027", "Down" = "#2c7fb8")) +
    theme_atac() +
    labs(
      title = "Accessible regions changed across contrasts",
      subtitle = subtitle,
      x = "", y = "Number of regions", fill = ""
    ) +
    theme(legend.position = "bottom")
  save_ggplot_both(file.path(figure_dir, "Contrast_change_summary.pdf"),
                   p_summary, width = 9, height = max(5, 0.34 * nrow(summary_df) + 2))

  common_ids <- Reduce(intersect, lapply(result_list, function(x) x$feature_id))
  if (length(common_ids) < 2) {
    log_msg("skip cross-contrast heatmaps: fewer than two common regions")
    return(invisible(summary_df))
  }
  lfc_mat <- vapply(result_list, function(x) {
    x$logFC[match(common_ids, x$feature_id)]
  }, numeric(length(common_ids)))
  rownames(lfc_mat) <- common_ids
  lfc_mat[!is.finite(lfc_mat)] <- NA_real_
  lfc_mat <- lfc_mat[rowSums(is.finite(lfc_mat)) == ncol(lfc_mat), , drop = FALSE]
  if (nrow(lfc_mat) < 2) return(invisible(summary_df))

  if (ncol(lfc_mat) > 1) {
    cor_mat <- stats::cor(lfc_mat, method = "spearman", use = "pairwise.complete.obs")
    save_device_both(
      file.path(figure_dir, "Contrast_logFC_correlation_heatmap.pdf"),
      width = max(7, 0.48 * ncol(cor_mat) + 3),
      height = max(6, 0.48 * nrow(cor_mat) + 2),
      plot_fun = function() pheatmap(
        cor_mat, border_color = "white", fontsize = 8,
        cluster_rows = nrow(cor_mat) > 1,
        cluster_cols = ncol(cor_mat) > 1,
        main = "Contrast similarity (Spearman correlation of region log2FC)"
      )
    )
  } else {
    log_msg("skip cross-contrast correlation and variable-region heatmaps: only one contrast")
    return(invisible(summary_df))
  }

  plot_cross_contrast_heatmap(
    result_list,
    outdir,
    top_n = top_n,
    output_prefix = "Cross_contrast",
    title = "Top variable regions across contrasts (log2FC, clipped to +/-3)",
    annotation_mode = annotation_mode
  )
  invisible(summary_df)
}

normalize_cross_contrasts <- function(selected_contrasts) {
  if (is.null(selected_contrasts) || length(selected_contrasts) == 0) return(NULL)
  vals <- unlist(strsplit(paste(selected_contrasts, collapse = ","), "[,;]"))
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  unique(vals)
}

plot_cross_contrast_heatmap <- function(result_list, outdir,
                                        selected_contrasts = NULL,
                                        top_n = 100,
                                        output_prefix = "Cross_contrast_selected",
                                        title = NULL,
                                        annotation_mode = 'gene_te') {
  annotation_mode <- match.arg(annotation_mode, c('gene_te', 'gene', 'none'))
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  if (length(result_list) == 0) stop("no contrast result is available")

  requested <- normalize_cross_contrasts(selected_contrasts)
  if (!is.null(requested)) {
    missing <- setdiff(requested, names(result_list))
    if (length(missing) > 0) {
      stop("requested contrast not found: ", paste(missing, collapse = ", "),
           "; available: ", paste(names(result_list), collapse = ", "))
    }
    result_list <- result_list[requested]
  }
  if (length(result_list) < 2) {
    stop("cross-contrast heatmap requires at least two contrasts")
  }

  common_ids <- Reduce(intersect, lapply(result_list, function(x) as.character(x$feature_id)))
  if (length(common_ids) < 2) stop("fewer than two common regions across selected contrasts")

  lfc_mat <- vapply(result_list, function(x) {
    as.numeric(x$logFC[match(common_ids, x$feature_id)])
  }, numeric(length(common_ids)))
  colnames(lfc_mat) <- names(result_list)
  rownames(lfc_mat) <- common_ids
  lfc_mat[!is.finite(lfc_mat)] <- NA_real_
  lfc_mat <- lfc_mat[rowSums(is.finite(lfc_mat)) == ncol(lfc_mat), , drop = FALSE]
  if (nrow(lfc_mat) < 2) stop("fewer than two regions have complete logFC values")

  row_var <- apply(lfc_mat, 1, stats::var, na.rm = TRUE)
  ranking <- data.frame(
    feature_id = rownames(lfc_mat),
    cross_contrast_variance = as.numeric(row_var),
    max_abs_logFC = apply(abs(lfc_mat), 1, max, na.rm = TRUE),
    mean_abs_logFC = rowMeans(abs(lfc_mat), na.rm = TRUE),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(cross_contrast_variance), dplyr::desc(max_abs_logFC), feature_id)
  top_ids <- head(ranking$feature_id, min(top_n, nrow(ranking)))
  top_mat <- lfc_mat[top_ids, , drop = FALSE]
  top_mat_plot <- pmax(pmin(top_mat, 3), -3)

  ann_cols <- c(
    "annotation", "geneId", "distanceToTSS", "SYMBOL", "is_TE",
    "te_name", "te_family", "te_class", "peak_class", "label_name",
    "Chr", "Start", "End"
  )
  first_res <- result_list[[1]]
  ann <- first_res[match(top_ids, first_res$feature_id),
                   intersect(c("feature_id", ann_cols), colnames(first_res)), drop = FALSE]
  if (!"feature_id" %in% colnames(ann)) ann$feature_id <- top_ids
  ann <- ann[match(top_ids, ann$feature_id), , drop = FALSE]
  ann$feature_id <- top_ids
  if (!"peak_class" %in% colnames(ann)) ann$peak_class <- "Region"
  ann$peak_class <- dplyr::coalesce(as.character(ann$peak_class), "Region")
  ann$peak_class[ann$peak_class %in% c("", "NA", "NaN", "NULL")] <- "Region"

  if (annotation_mode == 'gene') {
    ann <- ann |>
      dplyr::select(-dplyr::any_of(c('is_TE', 'te_name', 'te_family', 'te_class', 'label_name')))
    ann$peak_class <- classify_peak_context(
      ann[['annotation']],
      rep(FALSE, nrow(ann))
    )
  } else if (annotation_mode == 'none') {
    ann <- ann |>
      dplyr::select(feature_id)
    ann$peak_class <- 'Region'
  }
  ann$annotation_mode <- annotation_mode

  signs <- lfc_mat[top_ids, , drop = FALSE]
  ann$contrast_pattern <- apply(signs, 1, function(v) {
    if (all(v > 0)) return("Up_all")
    if (all(v < 0)) return("Down_all")
    if (any(v > 0) && any(v < 0)) return("Mixed")
    "Near_zero"
  })
  ann$region_label <- make_region_labels(ann)
  bad_label <- is.na(ann$region_label) | ann$region_label == ""
  ann$region_label[bad_label] <- as.character(ann$feature_id[bad_label])
  bad_label <- is.na(ann$region_label) | ann$region_label == ""
  ann$region_label[bad_label] <- paste0("region_", which(bad_label))
  ann <- dplyr::left_join(ranking, ann, by = "feature_id") |>
    dplyr::distinct(feature_id, .keep_all = TRUE) |>
    dplyr::slice(match(top_ids, feature_id)) |>
    dplyr::mutate(region_label = make.unique(region_label, sep = "_"))
  bad_label <- is.na(ann$region_label) | ann$region_label == ""
  ann$region_label[bad_label] <- paste0("region_", which(bad_label))

  row_anno <- ann |>
    dplyr::select(region_label, contrast_pattern, peak_class) |>
    as.data.frame()
  rownames(row_anno) <- row_anno$region_label
  row_anno$region_label <- NULL
  rownames(top_mat_plot) <- ann$region_label

  result_dir <- file.path(outdir, "results", "overview")
  figure_dir <- file.path(outdir, "figures", "overview")
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  n_top <- nrow(top_mat_plot)
  stem <- sprintf("%s_top%d_logFC_heatmap", output_prefix, n_top)
  write.csv(ann, file.path(result_dir, paste0(stem, "_region_annotations.csv")),
            row.names = FALSE, quote = TRUE)
  matrix_out <- data.frame(feature_id = top_ids, top_mat, check.names = FALSE)
  write.csv(matrix_out, file.path(result_dir, paste0(stem, "_matrix.csv")),
            row.names = FALSE, quote = TRUE)
  writeLines(names(result_list), file.path(result_dir, paste0(stem, "_contrasts.txt")))

  if (is.null(title)) {
    title <- sprintf("Selected contrasts: %s", paste(names(result_list), collapse = " | "))
  }
  heatmap_width <- max(10, 0.55 * ncol(top_mat_plot) + 5)
  heatmap_height <- max(10, 0.17 * nrow(top_mat_plot) + 4)
  save_device_both(
    file.path(figure_dir, paste0(stem, ".pdf")),
    width = heatmap_width,
    height = heatmap_height,
    plot_fun = function() pheatmap(
      top_mat_plot,
      color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(101),
      breaks = seq(-3, 3, length.out = 102),
      annotation_row = row_anno,
      border_color = NA,
      cluster_rows = nrow(top_mat_plot) > 1,
      cluster_cols = ncol(top_mat_plot) > 1,
      show_rownames = TRUE,
      fontsize = 8,
      fontsize_row = ifelse(nrow(top_mat_plot) > 80, 5, 7),
      main = title
    )
  )
  invisible(list(matrix = top_mat, annotation = ann, ranking = ranking,
                 contrasts = names(result_list)))
}

run_edgeR_contrast <- function(mat, anno, meta, case, control, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  keep_samples <- meta$condition %in% c(case, control)
  meta2 <- meta[keep_samples, , drop = FALSE]
  mat2 <- mat[, rownames(meta2), drop = FALSE]
  meta2$condition <- factor(meta2$condition, levels = c(control, case))

  if (ncol(mat2) < 2) {
    log_msg('skip contrast ', case, ' vs ', control, ': <2 samples')
    return(NULL)
  }

  y <- DGEList(counts = mat2, group = meta2$condition)
  keep <- filterByExpr(y, group = meta2$condition)
  if (!any(keep)) {
    log_msg('skip contrast ', case, ' vs ', control, ': no features after filterByExpr')
    return(NULL)
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  lcpm <- cpm(y, log = TRUE, prior.count = 1)
  norm_cpm <- cpm(y, log = FALSE, prior.count = 1)
  rep_count <- table(meta2$condition)

  if (any(rep_count < 2)) {
    log_msg('contrast ', case, ' vs ', control, ': no/low replicate, use exploratory mode')
    ctrl <- rowMeans(lcpm[, meta2$condition == control, drop = FALSE])
    trt <- rowMeans(lcpm[, meta2$condition == case, drop = FALSE])
    res <- data.frame(
      feature_id = rownames(lcpm),
      logFC = trt - ctrl,
      logCPM = rowMeans(lcpm),
      baseMean = rowMeans(norm_cpm),
      meanCount = rowMeans(y$counts),
      PValue = NA_real_,
      FDR = NA_real_,
      stringsAsFactors = FALSE
    )
  } else {
    design <- model.matrix(~ condition, data = meta2)
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = 2)
    res <- topTags(qlf, n = Inf)$table |> as.data.frame()
    res$baseMean <- rowMeans(norm_cpm[rownames(res), , drop = FALSE])
    res$meanCount <- rowMeans(y$counts[rownames(res), , drop = FALSE])
    res$feature_id <- rownames(res)
    res <- res[, c('feature_id', 'logFC', 'logCPM', 'baseMean', 'meanCount', 'PValue', 'FDR')]
  }

  res <- dplyr::left_join(res, tibble::rownames_to_column(anno, 'feature_id'), by = 'feature_id')
  write.csv(res, file.path(outdir, sprintf('differential_peaks_%s_vs_%s.csv', case, control)), row.names = FALSE, quote = TRUE)
  res
}

auto_label_representatives <- function(res, x_col, y_col, label_top_n = 40) {
  if (is.null(res) || nrow(res) == 0) return(res)
  label_by_distance(res, x_col, y_col, label_col = "label_show", top_n = label_top_n)
}

plot_volcano_ma <- function(res, meta_sub, outdir, case, control,
                              padj_cutoff = 0.05, lfc_cutoff = 1,
                              baseMean_min = 5, label_top_n = 40) {
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  exploratory_mode <- all(is.na(res$FDR)) && all(is.na(res$PValue))

  if (exploratory_mode) {
    res$significant <- abs(res$logFC) >= lfc_cutoff &
      !is.na(res$baseMean) & res$baseMean >= baseMean_min
    res$plot_p <- NA_real_
    log_msg("plot volcano/MA in exploratory mode for ", case, " vs ", control)
  } else {
    res$significant <- !is.na(res$FDR) & res$FDR < padj_cutoff &
      abs(res$logFC) >= lfc_cutoff
    res$plot_p <- dplyr::coalesce(res$FDR, res$PValue)
    res$plot_p[is.na(res$plot_p) | res$plot_p <= 0] <- 1e-300
  }

  res$direction <- dplyr::case_when(
    res$significant & res$logFC > 0  ~ "Up",
    res$significant & res$logFC < 0  ~ "Down",
    TRUE                              ~ "NS"
  )
  res$direction <- factor(res$direction, levels = c("Up", "Down", "NS"))
  res$point_alpha <- ifelse(res$direction == "NS", 0.22, 0.82)

  n_up    <- sum(res$direction == "Up")
  n_down  <- sum(res$direction == "Down")
  n_total <- nrow(res)

  color_map <- c("Up" = "#d73027", "Down" = "#2c7fb8", "NS" = "#bdbdbd")

  # ---- MA plot: log2(baseMean) vs log2FC ----
  ma_df <- res
  ma_df$log2_baseMean <- log2(pmax(dplyr::coalesce(ma_df$baseMean, 0), 1, na.rm = TRUE))
  ma_df <- auto_label_representatives(ma_df, "log2_baseMean", "logFC", label_top_n)

  p_ma <- ggplot(ma_df, aes(x = log2_baseMean, y = logFC)) +
    geom_hline(yintercept = 0, color = "grey78", linewidth = 0.28) +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey65", linewidth = 0.28) +
    geom_point(aes(color = direction, alpha = point_alpha), size = 0.75) +
    scale_color_manual(values = color_map, drop = FALSE) +
    scale_alpha_identity() +
    safe_repel_layer(ma_df, "log2_baseMean", "logFC", "label_show") +
    theme_atac() +
    labs(
      title    = sprintf("MA plot: %s vs %s", case, control),
      subtitle = sprintf("Up: %d  |  Down: %d  |  Total: %d", n_up, n_down, n_total),
      x        = "log2(baseMean)",
      y        = "log2FC",
      color    = ""
    ) +
    theme(legend.position = "bottom")

  save_ggplot_both(file.path(outdir, 'differential', sprintf('MA_%s_vs_%s.pdf', case, control)),
                   p_ma, width = 7, height = 6)

  # ---- Volcano plot: log2FC vs -log10(p-value) ----
  if (!exploratory_mode) {
    vol_df <- res
    vol_df$neg_log10_p <- -log10(vol_df$plot_p)
    vol_df <- auto_label_representatives(vol_df, "logFC", "neg_log10_p", label_top_n)

    p_vol <- ggplot(vol_df, aes(x = logFC, y = neg_log10_p)) +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      geom_point(aes(color = direction, alpha = point_alpha), size = 0.75) +
      scale_color_manual(values = color_map, drop = FALSE) +
      scale_alpha_identity() +
      safe_repel_layer(vol_df, "logFC", "neg_log10_p", "label_show") +
      theme_atac() +
      labs(
        title    = sprintf("Volcano: %s vs %s", case, control),
        subtitle = sprintf("Up: %d  |  Down: %d  |  Total: %d  (FDR < %s, |log2FC| > %s)",
                          n_up, n_down, n_total, padj_cutoff, lfc_cutoff),
        x        = "log2FC",
        y        = expression(-log[10](FDR)),
        color    = ""
      ) +
      theme(legend.position = "bottom")

    save_ggplot_both(file.path(outdir, 'differential', sprintf('Volcano_%s_vs_%s.pdf', case, control)),
                     p_vol, width = 7, height = 6)
  }

  invisible(res)
}

plot_heatmap_top <- function(mat, meta, res, case, control, outdir, top_n = 200) {
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))
  exploratory_mode <- all(is.na(res$FDR)) && all(is.na(res$PValue))
  if (exploratory_mode) {
    sig <- res |>
      dplyr::arrange(dplyr::desc(abs(logFC)), dplyr::desc(baseMean)) |>
      dplyr::slice_head(n = top_n)
  } else {
    sig <- res |>
      dplyr::arrange(dplyr::coalesce(FDR, 1), dplyr::desc(abs(logFC))) |>
      dplyr::slice_head(n = top_n)
  }
  ids <- intersect(sig$feature_id, rownames(mat))
  if (length(ids) < 2) return(invisible(NULL))
  lcpm <- edgeR::cpm(mat[ids, rownames(meta), drop = FALSE], log = TRUE, prior.count = 1)
  ann_col <- data.frame(condition = meta$condition)
  rownames(ann_col) <- meta$sample_label
  colnames(lcpm) <- meta$sample_label
  sig2 <- sig[match(ids, sig$feature_id), , drop = FALSE]
  direction <- dplyr::case_when(
    !is.na(sig2$logFC) & sig2$logFC > 0 ~ "Up",
    !is.na(sig2$logFC) & sig2$logFC < 0 ~ "Down",
    TRUE ~ "NS"
  )
  peak_class <- if ("peak_class" %in% colnames(sig2)) as.character(sig2$peak_class) else "Region"
  ann_row <- data.frame(direction = direction, peak_class = peak_class)
  row_labels <- make_region_labels(sig2)
  rownames(lcpm) <- row_labels
  rownames(ann_row) <- row_labels
  show_rows <- nrow(lcpm) <= 80
  save_device_both(file.path(outdir, sprintf('Heatmap_top_%s_vs_%s.pdf', case, control)), width = 7, height = 9,
                   plot_fun = function() pheatmap(lcpm, scale = 'row', annotation_col = ann_col,
                                                  annotation_row = ann_row,
                                                  show_rownames = show_rows,
                                                  border_color = NA, fontsize = 9,
                                                  main = sprintf('Top variable regions: %s vs %s', case, control)))
}

txdb_cache <- new.env(parent = emptyenv())
region_annotation_cache <- new.env(parent = emptyenv())
te_background_cache <- new.env(parent = emptyenv())

get_anno_db_name <- function(species) {
  if (species == 'hg38') return('org.Hs.eg.db')
  if (species %in% c('mm10', 'mm39')) return('org.Mm.eg.db')
  NULL
}

get_txdb <- function(species, gtf = '') {
  if (!is.null(gtf) && !is.na(gtf) && gtf != '' && file.exists(gtf)) {
    key <- paste0('gtf:', normalizePath(gtf, mustWork = FALSE))
    if (exists(key, envir = txdb_cache, inherits = FALSE)) return(get(key, envir = txdb_cache))
    suppressPackageStartupMessages(library(GenomicFeatures))
    log_msg('building TxDb from GTF: ', gtf)
    txdb <- GenomicFeatures::makeTxDbFromGFF(gtf)
    assign(key, txdb, envir = txdb_cache)
    return(txdb)
  }
  if (species == 'hg38') {
    suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
    return(TxDb.Hsapiens.UCSC.hg38.knownGene)
  }
  if (species == 'mm10') {
    suppressPackageStartupMessages(library(TxDb.Mmusculus.UCSC.mm10.knownGene))
    return(TxDb.Mmusculus.UCSC.mm10.knownGene)
  }
  if (species == 'mm39') {
    log_msg('mm39 annotation requires --gtf; refusing to use mm10 TxDb for mm39 coordinates')
    return(NULL)
  }
  NULL
}

run_annotate_peak <- function(gr, txdb, species) {
  anno_db <- get_anno_db_name(species)
  if (is.null(anno_db)) {
    return(annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), verbose = FALSE))
  }
  annotatePeak(gr, TxDb = txdb, tssRegion = c(-3000, 3000), annoDb = anno_db, verbose = FALSE)
}

get_orgdb <- function(species) {
  if (species == 'hg38') {
    suppressPackageStartupMessages(library(org.Hs.eg.db))
    return(org.Hs.eg.db)
  }
  if (species %in% c('mm10', 'mm39')) {
    suppressPackageStartupMessages(library(org.Mm.eg.db))
    return(org.Mm.eg.db)
  }
  NULL
}

pick_te_meta_col <- function(df, candidates) {
  x <- intersect(candidates, colnames(df))
  if (length(x) == 0) return(NULL)
  x[[1]]
}

clean_te_label <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ''
  x <- sub(';.*$', '', x)
  x <- gsub('([:|_]TE)+$', '', x)
  x[x %in% c('NA', 'NaN', 'NULL', 'TRUE', 'FALSE')] <- ''
  x
}

clean_te_value <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ''
  gsub('([:|_]TE)+(?=;|$)', '', x, perl = TRUE)
}

first_token <- function(x) {
  clean_te_label(x)
}

get_te_genome_background <- function(te_bed = '') {
  if (is.null(te_bed) || is.na(te_bed) || te_bed == '' || !file.exists(te_bed)) return(NULL)
  key <- normalizePath(te_bed, mustWork = FALSE)
  if (exists(key, envir = te_background_cache, inherits = FALSE)) {
    return(get(key, envir = te_background_cache, inherits = FALSE))
  }
  te_gr <- tryCatch(rtracklayer::import(te_bed), error = function(e) NULL)
  if (is.null(te_gr) || length(te_gr) == 0) return(NULL)
  meta <- as.data.frame(mcols(te_gr))
  name_col <- pick_te_meta_col(meta, c('repName', 'gene_name', 'name', 'Name', 'transcript_id', 'gene_id'))
  fam_col <- pick_te_meta_col(meta, c('repFamily', 'family', 'family_id'))
  class_col <- pick_te_meta_col(meta, c('repClass', 'class', 'class_id'))
  extract_col <- function(col) {
    if (is.null(col)) return(rep(NA_character_, length(te_gr)))
    value <- clean_te_label(meta[[col]])
    value[value == ''] <- NA_character_
    value
  }
  out <- data.frame(
    class = extract_col(class_col),
    family = extract_col(fam_col),
    subfamily = extract_col(name_col),
    stringsAsFactors = FALSE
  )
  assign(key, out, envir = te_background_cache)
  out
}

get_te_overlap_table <- function(gr, te_bed = '') {
  out <- data.frame(
    feature_id = as.character(mcols(gr)$feature_id),
    is_TE = FALSE,
    te_name = NA_character_,
    te_family = NA_character_,
    te_class = NA_character_,
    stringsAsFactors = FALSE
  )
  if (is.null(te_bed) || is.na(te_bed) || te_bed == '' || !file.exists(te_bed)) {
    log_msg('TE annotation file missing or not set, skip TE overlap')
    return(out)
  }
  te_gr <- tryCatch(rtracklayer::import(te_bed), error = function(e) NULL)
  if (is.null(te_gr) || length(te_gr) == 0) {
    log_msg('TE annotation import failed or empty: ', te_bed)
    return(out)
  }
  hits <- findOverlaps(gr, te_gr, ignore.strand = TRUE)
  if (length(hits) == 0) {
    log_msg('no peak overlaps found against TE annotation: ', te_bed)
    return(out)
  }
  te_meta <- as.data.frame(mcols(te_gr))
  name_col <- pick_te_meta_col(te_meta, c('repName', 'gene_name', 'name', 'Name', 'transcript_id', 'gene_id'))
  fam_col <- pick_te_meta_col(te_meta, c('repFamily', 'family', 'family_id'))
  class_col <- pick_te_meta_col(te_meta, c('repClass', 'class', 'class_id'))
  qh <- queryHits(hits)
  sh <- subjectHits(hits)
  out$is_TE[unique(qh)] <- TRUE

  assign_collapsed_te_meta <- function(source_col, target_col) {
    if (is.null(source_col)) return(invisible(NULL))
    vals <- as.character(te_meta[sh, source_col])
    vals[is.na(vals)] <- ''
    collapsed <- vapply(split(vals, qh), function(v) {
      v <- unique(v[v != ''])
      if (length(v) == 0) return(NA_character_)
      paste(v, collapse = ';')
    }, character(1))
    idx <- as.integer(names(collapsed))
    keep <- !is.na(collapsed) & collapsed != ''
    out[[target_col]][idx[keep]] <<- collapsed[keep]
    invisible(NULL)
  }
  assign_collapsed_te_meta(name_col, "te_name")
  out$te_name <- clean_te_value(out$te_name)
  assign_collapsed_te_meta(fam_col, 'te_family')
  out$te_family <- clean_te_value(out$te_family)
  assign_collapsed_te_meta(class_col, 'te_class')
  out$te_class <- clean_te_value(out$te_class)
  log_msg('TE overlaps annotated for peaks: ', sum(out$is_TE, na.rm = TRUE))
  out
}

classify_peak_context <- function(annotation, is_te) {
  n <- max(length(annotation), length(is_te), 0)
  if (n == 0) return(character(0))
  if (is.null(annotation)) annotation <- rep(NA_character_, n)
  if (is.null(is_te)) is_te <- rep(FALSE, n)
  if (length(annotation) == 1 && n > 1) annotation <- rep(annotation, n)
  if (length(is_te) == 1 && n > 1) is_te <- rep(is_te, n)
  ann <- dplyr::coalesce(as.character(annotation), '')
  te_flag <- rep(FALSE, n)
  te_flag[seq_len(min(length(is_te), n))] <- as.logical(is_te)[seq_len(min(length(is_te), n))] %in% TRUE
  dplyr::case_when(
    te_flag ~ 'TE',
    grepl('Promoter', ann, ignore.case = TRUE) ~ 'Promoter',
    grepl('Exon|Intron|UTR|Downstream', ann, ignore.case = TRUE) ~ 'GeneBody',
    grepl('Intergenic', ann, ignore.case = TRUE) ~ 'Intergenic',
    TRUE ~ 'Other'
  )
}

annotate_peak_table <- function(res, species, te_bed = '', gtf = '', annotation_mode = 'gene_te') {
  annotation_mode <- match.arg(annotation_mode, c('gene_te', 'gene', 'none'))
  if (is.null(res) || nrow(res) == 0) return(NULL)
  peak_df <- res |>
    dplyr::filter(!is.na(Chr), !is.na(Start), !is.na(End)) |>
    dplyr::select(-dplyr::any_of(c(
      'annotation', 'geneId', 'distanceToTSS', 'SYMBOL',
      'is_TE', 'te_name', 'te_family', 'te_class',
      'peak_class', 'label_name', 'feature_key'
    )))
  if (nrow(peak_df) == 0) return(NULL)

  all_peak_df <- peak_df
  cache_key <- paste(
    annotation_mode,
    species,
    normalizePath(gtf %||% '', mustWork = FALSE),
    normalizePath(te_bed %||% '', mustWork = FALSE),
    sep = '|'
  )
  cached <- if (exists(cache_key, envir = region_annotation_cache, inherits = FALSE)) {
    get(cache_key, envir = region_annotation_cache, inherits = FALSE)
  } else {
    data.frame(feature_id = character(), stringsAsFactors = FALSE)
  }
  peak_df <- all_peak_df[!all_peak_df$feature_id %in% cached$feature_id, , drop = FALSE]
  if (nrow(peak_df) == 0) {
    log_msg('reuse static gene/TE annotations from cache: ', nrow(all_peak_df), ' regions')
    return(dplyr::left_join(all_peak_df, cached, by = 'feature_id'))
  }
  log_msg('annotate uncached regions: ', nrow(peak_df), ' / ', nrow(all_peak_df))

  gr <- GRanges(seqnames = peak_df$Chr, ranges = IRanges(peak_df$Start, peak_df$End), feature_id = peak_df$feature_id)

  txdb <- if (annotation_mode == 'none') NULL else get_txdb(species, gtf = gtf)
  gene_ann <- data.frame(
    feature_id = peak_df$feature_id,
    annotation = NA_character_,
    geneId = NA_character_,
    distanceToTSS = NA_real_,
    SYMBOL = NA_character_,
    stringsAsFactors = FALSE
  )
  if (!is.null(txdb)) {
    ap <- run_annotate_peak(gr, txdb, species)
    ap_df <- as.data.frame(ap)
    keep_cols <- intersect(c('annotation', 'geneId', 'distanceToTSS', 'SYMBOL'), colnames(ap_df))
    chr_col <- intersect(c('seqnames', 'Chr', 'chr'), colnames(ap_df))
    start_col <- intersect(c('start', 'Start'), colnames(ap_df))
    end_col <- intersect(c('end', 'End'), colnames(ap_df))
    if (length(keep_cols) > 0 && length(chr_col) > 0 && length(start_col) > 0 && length(end_col) > 0) {
      peak_df$feature_key <- paste0(as.character(peak_df$Chr), ':', peak_df$Start, '-', peak_df$End)
      ap_df$feature_key <- paste0(as.character(ap_df[[chr_col[[1]]]]), ':', ap_df[[start_col[[1]]]], '-', ap_df[[end_col[[1]]]])
      ap_df <- ap_df |>
        dplyr::select(feature_key, dplyr::all_of(keep_cols)) |>
        dplyr::distinct(feature_key, .keep_all = TRUE)
      gene_ann <- peak_df |>
        dplyr::select(feature_id, feature_key) |>
        dplyr::left_join(ap_df, by = 'feature_key') |>
        dplyr::select(-feature_key)
      log_msg('gene annotation matched peaks: ', sum(!is.na(gene_ann$annotation)), ' / ', nrow(gene_ann))
      peak_df$feature_key <- NULL
    } else {
      log_msg('annotatePeak output missing coordinate/annotation columns, skip gene annotation merge')
    }
  } else {
    log_msg('no TxDb configured for species ', species, ', skip gene annotation')
  }

  te_ann <- if (annotation_mode == 'gene_te') {
    get_te_overlap_table(gr, te_bed = te_bed)
  } else {
    data.frame(
      feature_id = peak_df$feature_id,
      is_TE = FALSE,
      te_name = NA_character_,
      te_family = NA_character_,
      te_class = NA_character_,
      stringsAsFactors = FALSE
    )
  }
  out <- peak_df |>
    dplyr::left_join(gene_ann, by = 'feature_id') |>
    dplyr::left_join(te_ann, by = 'feature_id')
  required_annotation_cols <- c(
    'annotation', 'geneId', 'distanceToTSS', 'SYMBOL',
    'is_TE', 'te_name', 'te_family', 'te_class'
  )
  for (col in required_annotation_cols) {
    if (!col %in% colnames(out)) {
      out[[col]] <- if (col == 'is_TE') FALSE else NA
    }
  }
  out$peak_class <- classify_peak_context(out[['annotation']], out[['is_TE']])
  te_label <- clean_te_label(out[['te_name']])
  out$label_name <- ifelse(
    te_label != '', te_label,
    ifelse(!is.na(out[['SYMBOL']]) & out[['SYMBOL']] != '', out[['SYMBOL']], out[['feature_id']])
  )
  static_cols <- c(
    'feature_id', 'annotation', 'geneId', 'distanceToTSS', 'SYMBOL',
    'is_TE', 'te_name', 'te_family', 'te_class', 'peak_class', 'label_name'
  )
  static_new <- out |>
    dplyr::select(dplyr::all_of(static_cols)) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)
  cache_updated <- dplyr::bind_rows(cached, static_new) |>
    dplyr::distinct(feature_id, .keep_all = TRUE)
  assign(cache_key, cache_updated, envir = region_annotation_cache)
  dplyr::left_join(all_peak_df, cache_updated, by = 'feature_id')
}

simplify_annotation <- function(x) {
  x <- dplyr::coalesce(as.character(x), 'Other')
  x <- sub('\n.*$', '', x)
  x <- sub(' \\(.*$', '', x)
  x[x == ''] <- 'Other'
  x
}

plot_simple_annotation_piebar <- function(ap_df, out_prefix, title_prefix = 'Peak annotation') {
  if (is.null(ap_df) || nrow(ap_df) == 0 || !'annotation' %in% colnames(ap_df)) return(invisible(NULL))
  anno_tab <- ap_df |>
    dplyr::mutate(annotation_simple = simplify_annotation(annotation)) |>
    dplyr::count(annotation_simple, name = 'n') |>
    dplyr::arrange(dplyr::desc(n)) |>
    dplyr::mutate(frac = n / sum(n), label = paste0(annotation_simple, ' (', n, ')'))
  if (nrow(anno_tab) == 0) return(invisible(NULL))

  p_pie <- ggplot(anno_tab, aes(x = 1, y = n, fill = annotation_simple)) +
    geom_col(width = 1, color = 'white') +
    coord_polar(theta = 'y') +
    theme_void() +
    labs(title = title_prefix, fill = 'Annotation')
  save_ggplot_both(paste0(out_prefix, '_pie.pdf'), p_pie, width = 6, height = 5)

  p_bar <- ggplot(anno_tab, aes(x = reorder(annotation_simple, n), y = n, fill = annotation_simple)) +
    geom_col() +
    coord_flip() +
    theme_bw() +
    labs(title = title_prefix, x = 'Genomic regions', y = 'Peak count')
  save_ggplot_both(paste0(out_prefix, '_bar.pdf'), p_bar, width = 7.5, height = 5)
}

plot_te_annotation_distribution <- function(ann_df, outdir, title_prefix, top_n = 15,
                                            te_bed = '', table_outdir = NULL) {
  if (is.null(ann_df) || nrow(ann_df) == 0 || !'is_TE' %in% colnames(ann_df)) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  te <- ann_df |>
    dplyr::filter(is_TE %in% TRUE) |>
    dplyr::mutate(
      te_class1 = first_token(te_class),
      te_family1 = first_token(te_family),
      te_name1 = first_token(te_name)
    )
  total_peaks <- nrow(ann_df)
  te_peaks <- nrow(te)
  status_df <- data.frame(
    status = factor(c('TE-overlap', 'non-TE'), levels = c('TE-overlap', 'non-TE')),
    n = c(te_peaks, total_peaks - te_peaks)
  ) |>
    dplyr::mutate(frac = ifelse(total_peaks > 0, n / total_peaks, 0),
                  label = sprintf('%s\n%s (%.1f%%)', status, scales::comma(n), 100 * frac))
  p_status <- ggplot(status_df, aes(x = status, y = n, fill = status)) +
    geom_col(width = 0.72) +
    geom_text(aes(label = label), vjust = -0.2, size = 3.2) +
    theme_bw() +
    theme(legend.position = 'none') +
    labs(title = paste0(title_prefix, ': TE annotation status'),
         subtitle = paste0('Total peaks = ', scales::comma(total_peaks),
                           '; TE-overlapping peaks = ', scales::comma(te_peaks)),
         x = NULL, y = 'Peak count')
  save_ggplot_both(file.path(outdir, 'TE_annotation_status_bar.pdf'), p_status, width = 6.8, height = 4.6)

  if (te_peaks == 0) return(invisible(NULL))
  genome_bg <- get_te_genome_background(te_bed)
  if (!is.null(table_outdir)) dir.create(table_outdir, recursive = TRUE, showWarnings = FALSE)

  for (level in c('class', 'family', 'subfamily')) {
    coln <- c(class = 'te_class1', family = 'te_family1', subfamily = 'te_name1')[[level]]
    tab <- te |>
      dplyr::mutate(group = dplyr::coalesce(.data[[coln]], 'Unannotated'),
                    group = ifelse(group == '', 'Unannotated', group)) |>
      dplyr::count(group, name = 'n') |>
      dplyr::arrange(dplyr::desc(n))
    if (nrow(tab) == 0) next
    if (nrow(tab) > top_n) {
      tab <- dplyr::bind_rows(
        tab |> dplyr::slice_head(n = top_n),
        data.frame(group = 'Other', n = sum(tab$n[-seq_len(top_n)]), stringsAsFactors = FALSE)
      )
    }
    tab <- tab |>
      dplyr::mutate(frac = n / sum(n),
                    label = sprintf('%s (%.1f%%)', scales::comma(n), 100 * frac))
    p <- ggplot(tab, aes(x = reorder(group, n), y = n, fill = group)) +
      geom_col(show.legend = FALSE) +
      geom_text(aes(label = label), hjust = -0.05, size = 3) +
      coord_flip(clip = 'off') +
      theme_bw() +
      theme(plot.margin = margin(5.5, 36, 5.5, 5.5)) +
      labs(title = paste0(title_prefix, ': TE ', level, ' distribution'),
           subtitle = paste0('Denominator = TE-overlapping peaks; n = ', scales::comma(te_peaks)),
           x = NULL, y = 'Peak count')
    save_ggplot_both(file.path(outdir, sprintf('TE_%s_distribution_bar.pdf', level)), p, width = 8.4, height = 5.2)

    if (!is.null(genome_bg) && level %in% colnames(genome_bg)) {
      observed <- te |>
        dplyr::transmute(group = dplyr::coalesce(.data[[coln]], 'Unannotated')) |>
        dplyr::mutate(group = ifelse(group == '', 'Unannotated', group)) |>
        dplyr::count(group, name = 'observed_n')
      background <- genome_bg |>
        dplyr::transmute(group = dplyr::coalesce(.data[[level]], 'Unannotated')) |>
        dplyr::mutate(group = ifelse(group == '', 'Unannotated', group)) |>
        dplyr::count(group, name = 'background_n')
      enrich <- dplyr::full_join(observed, background, by = 'group') |>
        dplyr::mutate(
          observed_n = dplyr::coalesce(observed_n, 0L),
          background_n = dplyr::coalesce(background_n, 0L)
        )
      n_obs <- sum(enrich$observed_n)
      n_bg <- sum(enrich$background_n)
      k <- nrow(enrich)
      enrich <- enrich |>
        dplyr::mutate(
          observed_fraction = (observed_n + 0.5) / (n_obs + 0.5 * k),
          background_fraction = (background_n + 0.5) / (n_bg + 0.5 * k),
          log2_enrichment = log2(observed_fraction / background_fraction),
          expected_n = n_obs * background_fraction,
          binomial_sd = sqrt(pmax(n_obs * background_fraction * (1 - background_fraction), 0)),
          z_score = ifelse(binomial_sd > 0, (observed_n - expected_n) / binomial_sd, NA_real_),
          background_unit = 'TE_instance_count'
        ) |>
        dplyr::arrange(dplyr::desc(abs(log2_enrichment)))
      if (!is.null(table_outdir)) {
        write.csv(
          enrich,
          file.path(table_outdir, sprintf('TE_%s_enrichment_vs_genome.csv', level)),
          row.names = FALSE, quote = FALSE
        )
      }
      plot_enrich <- enrich |>
        dplyr::filter(observed_n > 0) |>
        dplyr::slice_max(order_by = abs(log2_enrichment), n = top_n, with_ties = FALSE) |>
        dplyr::mutate(
          group = reorder(group, log2_enrichment),
          direction = ifelse(log2_enrichment >= 0, 'Enriched', 'Depleted')
        )
      if (nrow(plot_enrich) > 0) {
        p_enrich <- ggplot(plot_enrich, aes(x = group, y = log2_enrichment, fill = direction)) +
          geom_hline(yintercept = 0, color = 'grey45', linewidth = 0.3) +
          geom_col(width = 0.72) +
          coord_flip() +
          scale_fill_manual(values = c(Enriched = '#d73027', Depleted = '#2c7fb8')) +
          theme_atac() +
          labs(
            title = paste0(title_prefix, ': TE ', level, ' enrichment vs genome'),
            subtitle = 'Background = all annotated TE instances; pseudocount = 0.5',
            x = NULL, y = 'log2(observed fraction / genome fraction)', fill = ''
          )
        save_ggplot_both(
          file.path(outdir, sprintf('TE_%s_enrichment_vs_genome.pdf', level)),
          p_enrich, width = 8.8, height = 5.8
        )
      }
    }
  }
}

annotate_diff_peaks <- function(res, species, outdir, case, control, padj_cutoff = 0.05, lfc_cutoff = 1, baseMean_min = 5, te_bed = '', gtf = '', annotation_mode = 'gene_te') {
  croot <- contrast_root(outdir, case, control)
  ann_all <- annotate_peak_table(res, species, te_bed = te_bed, gtf = gtf, annotation_mode = annotation_mode)
  if (is.null(ann_all) || nrow(ann_all) == 0) {
    log_msg('annotated peak table unavailable for ', case, ' vs ', control)
    return(res)
  }

  subdir <- file.path(croot, 'annotation')
  results_dir(croot)
  invisible(figures_dir(croot, "annotation/sig"))
  invisible(figures_dir(croot, "annotation/up"))
  invisible(figures_dir(croot, "annotation/down"))
  invisible(figures_dir(croot, "differential"))
  write.csv(ann_all, file.path(croot, 'results', 'peak_annotation_all.csv'), row.names = FALSE, quote = TRUE)

  exploratory_mode <- all(is.na(ann_all$FDR)) && all(is.na(ann_all$PValue))
  sig <- ann_all |> dplyr::filter(!is.na(logFC), !is.na(Chr), !is.na(Start), !is.na(End))
  if (exploratory_mode) {
    sig <- sig |> dplyr::filter(abs(logFC) >= lfc_cutoff, !is.na(baseMean), baseMean >= baseMean_min)
  } else {
    sig <- sig |> dplyr::filter(!is.na(FDR), FDR < padj_cutoff, abs(logFC) >= lfc_cutoff)
  }
  write.csv(sig, file.path(croot, 'results', 'peak_annotation_sig.csv'), row.names = FALSE, quote = TRUE)

  cls <- sig |>
    dplyr::count(peak_class, name = 'n') |>
    dplyr::mutate(peak_class = factor(peak_class, levels = c('TE', 'Promoter', 'GeneBody', 'Intergenic', 'Other')))
  if (nrow(cls) > 0) {
    p_cls <- ggplot(cls, aes(x = peak_class, y = n, fill = peak_class)) +
      geom_col() + theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = sprintf('Annotated differential peaks: %s vs %s', case, control), x = '', y = 'Peak count')
    save_ggplot_both(file.path(croot, 'figures', 'annotation', 'sig', 'Peak_class_barplot.pdf'), p_cls, width = 7, height = 4.5)
  }
  plot_te_annotation_distribution(
    sig, file.path(croot, 'figures', 'annotation', 'sig', 'TE_distribution'),
    sprintf('Differential peaks %s vs %s', case, control),
    te_bed = te_bed,
    table_outdir = file.path(croot, 'results', 'TE_distribution', 'sig')
  )

  if (nrow(sig) > 0) {
    plot_simple_annotation_piebar(
      sig,
      file.path(croot, 'figures', 'annotation', 'sig', 'Peak_annotation'),
      sprintf('Peak annotation: %s vs %s', case, control)
    )
  }

  # === UP/DOWN-separated annotation ===
  up_dir <- file.path(subdir, 'up')
  down_dir <- file.path(subdir, 'down')
  sig_up <- sig |> dplyr::filter(logFC > 0)
  sig_down <- sig |> dplyr::filter(logFC < 0)

  if (nrow(sig_up) > 0) {
    write.csv(sig_up, file.path(croot, 'results', 'peak_annotation_up.csv'), row.names = FALSE, quote = TRUE)
    cls_u <- sig_up |> dplyr::count(peak_class, name = 'n') |>
      dplyr::mutate(peak_class = factor(peak_class, levels = c('TE','Promoter','GeneBody','Intergenic','Other')))
    if (nrow(cls_u) > 0) {
      p_u <- ggplot(cls_u, aes(x = peak_class, y = n, fill = peak_class)) +
        geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = sprintf('UP peaks: %s vs %s (n=%d)', case, control, nrow(sig_up)), x = '', y = 'Count')
      save_ggplot_both(file.path(croot, 'figures', 'annotation', 'up', 'Peak_class_barplot.pdf'), p_u, width = 6, height = 4.5)
    }
    plot_te_annotation_distribution(
      sig_up, file.path(croot, 'figures', 'annotation', 'up', 'TE_distribution'),
      sprintf('UP peaks %s vs %s', case, control),
      te_bed = te_bed,
      table_outdir = file.path(croot, 'results', 'TE_distribution', 'up')
    )
    plot_simple_annotation_piebar(
      sig_up,
      file.path(croot, 'figures', 'annotation', 'up', 'Peak_annotation'),
      sprintf('UP peaks: %s vs %s', case, control)
    )
  }

  if (nrow(sig_down) > 0) {
    write.csv(sig_down, file.path(croot, 'results', 'peak_annotation_down.csv'), row.names = FALSE, quote = TRUE)
    cls_d <- sig_down |> dplyr::count(peak_class, name = 'n') |>
      dplyr::mutate(peak_class = factor(peak_class, levels = c('TE','Promoter','GeneBody','Intergenic','Other')))
    if (nrow(cls_d) > 0) {
      p_d <- ggplot(cls_d, aes(x = peak_class, y = n, fill = peak_class)) +
        geom_col() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(title = sprintf('DOWN peaks: %s vs %s (n=%d)', case, control, nrow(sig_down)), x = '', y = 'Count')
      save_ggplot_both(file.path(croot, 'figures', 'annotation', 'down', 'Peak_class_barplot.pdf'), p_d, width = 6, height = 4.5)
    }
    plot_te_annotation_distribution(
      sig_down, file.path(croot, 'figures', 'annotation', 'down', 'TE_distribution'),
      sprintf('DOWN peaks %s vs %s', case, control),
      te_bed = te_bed,
      table_outdir = file.path(croot, 'results', 'TE_distribution', 'down')
    )
    plot_simple_annotation_piebar(
      sig_down,
      file.path(croot, 'figures', 'annotation', 'down', 'Peak_annotation'),
      sprintf('DOWN peaks: %s vs %s', case, control)
    )
  }

  # Up vs Down side-by-side
  if (nrow(sig_up) > 0 && nrow(sig_down) > 0) {
    tryCatch({
      anno_dual <- dplyr::bind_rows(
        dplyr::mutate(sig_up, direction = 'Up'),
        dplyr::mutate(sig_down, direction = 'Down')
      ) |>
        dplyr::mutate(annotation_simple = simplify_annotation(annotation)) |>
        dplyr::count(annotation_simple, direction, name = 'n') |>
        dplyr::group_by(direction) |> dplyr::mutate(pct = n / sum(n) * 100) |> dplyr::ungroup()
      p_dual <- ggplot(anno_dual, aes(x = annotation_simple, y = pct, fill = direction)) +
        geom_col(position = 'dodge') + theme_bw() +
        labs(title = sprintf('Up vs Down annotation: %s vs %s', case, control),
             x = '', y = 'Percentage (%)', fill = 'Direction')
      save_ggplot_both(file.path(croot, 'figures', 'annotation', 'Up_vs_Down_annotation_bar.pdf'), p_dual, width = 8, height = 4.5)
    }, error = function(e) log_msg('Up vs Down annotation comparison failed: ', e$message))
  }

  ann_all
}

plot_annotation_volcano <- function(res, outdir, case, control, label_top_n = 40) {
  if (is.null(res) || nrow(res) == 0 || !'peak_class' %in% colnames(res)) return(invisible(NULL))
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  exploratory_mode <- all(is.na(res$FDR)) && all(is.na(res$PValue))
  res$peak_class <- factor(res$peak_class, levels = c('TE', 'Promoter', 'GeneBody', 'Intergenic', 'Other'))
  res$label <- ''
  res$point_alpha <- 0.25

  if (exploratory_mode) {
    res$plot_y <- log10(pmax(dplyr::coalesce(res$baseMean, 0), 0) + 1)
    res$point_alpha[abs(res$logFC) >= 1 & res$plot_y >= log10(6)] <- 0.8
    res <- label_by_distance(res, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'label_name', top_n = label_top_n)
    p <- ggplot(res, aes(x = logFC, y = plot_y, color = peak_class)) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.28) +
      geom_point(aes(alpha = point_alpha), size = 0.9) +
      scale_alpha_identity() +
      safe_repel_layer(res, 'logFC', 'plot_y', 'label') +
      theme_atac() +
      labs(title = sprintf('Peak-class effect-abundance: %s vs %s', case, control), x = 'log2FC', y = 'log10(baseMean + 1)', color = 'Peak class')
    save_ggplot_both(file.path(outdir, sprintf('Peak_class_effect_%s_vs_%s.pdf', case, control)), p, width = 7, height = 5.5)
  } else {
    res$plot_p <- dplyr::coalesce(res$FDR, res$PValue)
    res$plot_p[is.na(res$plot_p) | res$plot_p <= 0] <- 1
    res$point_alpha[res$FDR < 0.05 & abs(res$logFC) >= 1] <- 0.8
    res$plot_y <- -log10(res$plot_p)
    res <- label_by_distance(res, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'label_name', top_n = label_top_n)
    p <- ggplot(res, aes(x = logFC, y = plot_y, color = peak_class)) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      geom_point(aes(alpha = point_alpha), size = 0.9) +
      scale_alpha_identity() +
      safe_repel_layer(res, 'logFC', 'plot_y', 'label') +
      theme_atac() +
      labs(title = sprintf('Peak-class volcano: %s vs %s', case, control), x = 'log2FC', y = '-log10(FDR/PValue)', color = 'Peak class')
    save_ggplot_both(file.path(outdir, sprintf('Peak_class_volcano_%s_vs_%s.pdf', case, control)), p, width = 7, height = 5.5)
  }

  gene_only <- res |> dplyr::filter(!(is_TE %in% TRUE))
  if (nrow(gene_only) >= 2) {
    gene_only$label <- ''
    gene_only$gene_group <- factor(gene_only$peak_class, levels = c('Promoter', 'GeneBody', 'Intergenic', 'Other'))
    if (exploratory_mode) {
      gene_only$plot_y <- log10(pmax(dplyr::coalesce(gene_only$baseMean, 0), 0) + 1)
      gene_only$direction <- dplyr::case_when(abs(gene_only$logFC) >= 1 & gene_only$logFC > 0 ~ "Up",
                                              abs(gene_only$logFC) >= 1 & gene_only$logFC < 0 ~ "Down",
                                              TRUE ~ "NS")
      gene_only$point_alpha <- ifelse(gene_only$direction == "NS", 0.22, 0.82)
      gene_only <- label_by_distance(gene_only, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'SYMBOL', top_n = label_top_n)
      p_gene <- ggplot(gene_only, aes(x = logFC, y = plot_y, color = direction, alpha = point_alpha)) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.28) +
        geom_point(size = 0.85) +
        scale_color_manual(values = c("Up" = "#d73027", "Down" = "#2c7fb8", "NS" = "#bdbdbd"), drop = FALSE) +
        scale_alpha_identity() +
        safe_repel_layer(gene_only, 'logFC', 'plot_y', 'label') +
        theme_atac() +
        labs(title = sprintf('Gene-only effect-abundance: %s vs %s', case, control), x = 'log2FC', y = 'log10(baseMean + 1)', color = '')
      save_ggplot_both(file.path(outdir, sprintf('Gene_only_effect_%s_vs_%s.pdf', case, control)), p_gene, width = 7, height = 5.5)
    } else {
      gene_only$plot_p <- dplyr::coalesce(gene_only$FDR, gene_only$PValue)
      gene_only$plot_p[is.na(gene_only$plot_p) | gene_only$plot_p <= 0] <- 1
      gene_only$plot_y <- -log10(gene_only$plot_p)
      gene_only$direction <- dplyr::case_when(!is.na(gene_only$FDR) & gene_only$FDR < 0.05 & gene_only$logFC > 0 & abs(gene_only$logFC) >= 1 ~ "Up",
                                              !is.na(gene_only$FDR) & gene_only$FDR < 0.05 & gene_only$logFC < 0 & abs(gene_only$logFC) >= 1 ~ "Down",
                                              TRUE ~ "NS")
      gene_only$point_alpha <- ifelse(gene_only$direction == "NS", 0.22, 0.82)
      gene_only <- label_by_distance(gene_only, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'SYMBOL', top_n = label_top_n)
      p_gene <- ggplot(gene_only, aes(x = logFC, y = plot_y, color = direction, alpha = point_alpha)) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.28) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.28) +
        geom_point(size = 0.85) +
        scale_color_manual(values = c("Up" = "#d73027", "Down" = "#2c7fb8", "NS" = "#bdbdbd"), drop = FALSE) +
        scale_alpha_identity() +
        safe_repel_layer(gene_only, 'logFC', 'plot_y', 'label') +
        theme_atac() +
        labs(title = sprintf('Gene-only volcano: %s vs %s', case, control), x = 'log2FC', y = '-log10(FDR/PValue)', color = '')
      save_ggplot_both(file.path(outdir, sprintf('Gene_only_volcano_%s_vs_%s.pdf', case, control)), p_gene, width = 7, height = 5.5)
    }
  }

  te_only <- res |> dplyr::filter(is_TE %in% TRUE)
  if (nrow(te_only) >= 2) {
    te_only$label <- ''
    if (exploratory_mode) {
      te_only$plot_y <- log10(pmax(dplyr::coalesce(te_only$baseMean, 0), 0) + 1)
      te_only$te_label <- first_token(te_only$te_name)
      te_only <- label_by_distance(te_only, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'te_label', top_n = label_top_n)
      te_only$point_alpha <- ifelse(abs(te_only$logFC) >= 1, 0.82, 0.22)
      p_te <- ggplot(te_only, aes(x = logFC, y = plot_y, color = first_token(te_class), alpha = point_alpha)) +
        geom_vline(xintercept = 0, color = "grey80", linewidth = 0.28) +
        geom_point(size = 0.85) +
        scale_alpha_identity() +
        safe_repel_layer(te_only, 'logFC', 'plot_y', 'label') +
        theme_atac() +
        labs(title = sprintf('TE-only effect-abundance: %s vs %s', case, control), x = 'log2FC', y = 'log10(baseMean + 1)', color = 'TE class')
      save_ggplot_both(file.path(outdir, sprintf('TE_only_effect_%s_vs_%s.pdf', case, control)), p_te, width = 7, height = 5.5)
    } else {
      te_only$plot_p <- dplyr::coalesce(te_only$FDR, te_only$PValue)
      te_only$plot_p[is.na(te_only$plot_p) | te_only$plot_p <= 0] <- 1
      te_only$plot_y <- -log10(te_only$plot_p)
      te_only$te_label <- first_token(te_only$te_name)
      te_only <- label_by_distance(te_only, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'te_label', top_n = label_top_n)
      te_only$point_alpha <- ifelse(!is.na(te_only$FDR) & te_only$FDR < 0.05 & abs(te_only$logFC) >= 1, 0.82, 0.22)
      p_te <- ggplot(te_only, aes(x = logFC, y = plot_y, color = first_token(te_class), alpha = point_alpha)) +
        geom_point(size = 0.85) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50", alpha = 0.5) +
        scale_alpha_identity() +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", alpha = 0.5) +
        safe_repel_layer(te_only, 'logFC', 'plot_y', 'label') +
        theme_atac() +
        labs(title = sprintf('TE-only volcano: %s vs %s', case, control), x = 'log2FC', y = '-log10(FDR/PValue)', color = 'TE class')
      save_ggplot_both(file.path(outdir, sprintf('TE_only_volcano_%s_vs_%s.pdf', case, control)), p_te, width = 7, height = 5.5)
    }
  }
}

run_peak_gene_enrichment <- function(res, species, outdir, case, control, padj_cutoff = 0.05, lfc_cutoff = 1, baseMean_min = 5, te_bed = '', gtf = '', label_top_n = 40) {
  croot <- contrast_root(outdir, case, control)
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))
  suppressPackageStartupMessages(library(clusterProfiler))
  orgdb <- get_orgdb(species)
  if (is.null(orgdb)) {
    log_msg('no OrgDb configured for species ', species, ', skip GO/GSEA')
    return(invisible(NULL))
  }

  ann <- annotate_peak_table(res, species, te_bed = te_bed, gtf = gtf)
  if (is.null(ann) || nrow(ann) == 0 || !'geneId' %in% colnames(ann)) {
    log_msg('peak annotation unavailable, skip GO/GSEA')
    return(invisible(NULL))
  }

  subdir <- file.path(croot, "differential", "gene")
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  results_dir(croot)
  invisible(figures_dir(croot, "gene"))
  invisible(figures_dir(croot, "TE"))
  write.csv(ann, file.path(croot, 'results', 'peak_to_gene_annotation_all.csv'), row.names = FALSE, quote = TRUE)

  exploratory_mode <- all(is.na(res$FDR)) && all(is.na(res$PValue))
  sig <- ann |> dplyr::filter(!is.na(logFC), !is.na(geneId), geneId != '')
  if (exploratory_mode) {
    sig <- sig |> dplyr::filter(abs(logFC) >= lfc_cutoff, !is.na(baseMean), baseMean >= baseMean_min)
  } else {
    sig <- sig |> dplyr::filter(!is.na(FDR), FDR < padj_cutoff, abs(logFC) >= lfc_cutoff)
  }
  sig_genes <- unique(as.character(sig$geneId))
  sig_genes <- sig_genes[sig_genes != '' & !is.na(sig_genes)]

  if (length(sig_genes) > 0) {
    ego <- tryCatch(enrichGO(gene = sig_genes, OrgDb = orgdb, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', readable = TRUE), error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego), file.path(croot, 'results', 'GO_BP_enrichment.csv'), row.names = FALSE, quote = FALSE)
      save_device_both(file.path(croot, 'figures', 'gene', 'GO_BP_dotplot.pdf'), width = 9, height = 6,
                       plot_fun = function() print(dotplot(ego, showCategory = 20) + ggplot2::ggtitle(sprintf('GO BP: %s vs %s', case, control))))
    } else {
      log_msg('GO enrichment returned empty for ', case, ' vs ', control)
    }
  }

  rank_df <- ann |>
    dplyr::filter(!is.na(geneId), geneId != '', !is.na(logFC)) |>
    dplyr::mutate(abs_fc = abs(logFC)) |>
    dplyr::arrange(dplyr::desc(abs_fc), dplyr::desc(baseMean)) |>
    dplyr::group_by(geneId) |>
    dplyr::slice_head(n = 1) |>
    dplyr::ungroup()

  gene_list <- rank_df$logFC
  names(gene_list) <- as.character(rank_df$geneId)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!duplicated(names(gene_list))]
  gene_list <- gene_list[!is.na(gene_list)]

  if (length(gene_list) >= 10 && length(unique(gene_list)) > 1) {
    ggo <- tryCatch(gseGO(geneList = gene_list, OrgDb = orgdb, keyType = 'ENTREZID', ont = 'BP', eps = 0, pvalueCutoff = 1), error = function(e) NULL)
    if (!is.null(ggo) && nrow(as.data.frame(ggo)) > 0) {
      write.csv(as.data.frame(ggo), file.path(croot, 'results', 'GSEA_GO_BP.csv'), row.names = FALSE, quote = FALSE)
      save_device_both(file.path(croot, 'figures', 'gene', 'GSEA_GO_BP_dotplot.pdf'), width = 9, height = 6,
                       plot_fun = function() print(dotplot(ggo, showCategory = 20, split = '.sign') + facet_grid(. ~ .sign) + ggplot2::ggtitle(sprintf('GSEA GO BP: %s vs %s', case, control))))
    } else {
      log_msg('GSEA returned empty for ', case, ' vs ', control)
    }
  } else {
    log_msg('skip GSEA for ', case, ' vs ', control, ': insufficient ranked genes')
  }

  # === Hallmark GSEA ===
  hallmark_ok <- FALSE
  hallmark_res <- tryCatch({
    hallmark_gmt <- atac_conf_file('h.all.v2024.1.Hs.symbols.gmt')
    if (!file.exists(hallmark_gmt)) {
      hallmark_gmt <- file.path(dirname(dirname(subdir)), '..', '..', '..', '..', '..', 'conf', 'h.all.v2024.1.Hs.symbols.gmt')
      hallmark_gmt <- normalizePath(hallmark_gmt, mustWork = FALSE)
    }
    if (!file.exists(hallmark_gmt)) {
      stop('Hallmark GMT file not found')
    }
    hallmark_list <- fgsea::gmtPathways(hallmark_gmt)
    # Convert Entrez IDs to gene symbols (Hallmark GMT uses symbols)
    entrez_to_symbol <- tryCatch({
      clusterProfiler::bitr(names(gene_list), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = orgdb)
    }, error = function(e) NULL)
    if (!is.null(entrez_to_symbol) && nrow(entrez_to_symbol) > 0) {
      gene_list_symbol <- gene_list
      names(gene_list_symbol) <- entrez_to_symbol$SYMBOL[match(names(gene_list), entrez_to_symbol$ENTREZID)]
      gene_list_symbol <- gene_list_symbol[!is.na(names(gene_list_symbol))]
      gene_list_symbol <- gene_list_symbol[names(gene_list_symbol) != ""]
      hallmark_res <- fgsea::fgseaMultilevel(pathways = hallmark_list, stats = gene_list_symbol, minSize = 10, maxSize = 500)
    } else {
      hallmark_res <- data.frame()
    }
    if (nrow(hallmark_res) > 0) {
      hallmark_res <- hallmark_res[order(hallmark_res$padj), ]
      hallmark_res$leadingEdge <- vapply(hallmark_res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
      write.csv(as.data.frame(hallmark_res), file.path(croot, 'results', 'GSEA_Hallmark.csv'), row.names = FALSE, quote = FALSE)
      top_h <- head(hallmark_res[!is.na(hallmark_res$padj) & hallmark_res$padj < 0.5, ], 30)
      if (nrow(top_h) == 0) top_h <- head(hallmark_res[order(hallmark_res$pval), ], 15)
      if (nrow(top_h) > 0) {
        top_h$direction <- ifelse(top_h$NES > 0, 'Up', 'Down')
        p_h <- ggplot(top_h, aes(x = NES, y = reorder(pathway, NES), fill = direction)) +
          geom_col() + theme_bw() +
          labs(title = sprintf('Hallmark GSEA: %s vs %s', case, control),
               x = 'Normalized Enrichment Score', y = '') +
          scale_fill_manual(values = c('Up' = '#e41a1c', 'Down' = '#377eb8'))
        save_ggplot_both(file.path(subdir, "GSEA_Hallmark_barplot.pdf"), p_h, width = 11, height = 7.5)
      }
      # Ridgeplot for GO enrichment (robust alternative to cnetplot)
      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        tryCatch({
          p_ridge <- enrichplot::ridgeplot(ego, showCategory = 20) +
            scale_fill_viridis_c(option = "plasma") +
            theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold")) +
            labs(title = sprintf("GO BP Ridgeplot: %s vs %s", case, control))
          save_ggplot_both(file.path(subdir, "GO_BP_ridgeplot.pdf"), p_ridge, width = 10, height = 7)
        }, error = function(e) log_msg("Ridgeplot failed: ", e))
      }

      hallmark_ok <- TRUE
      log_msg('Hallmark GSEA completed for ', case, ' vs ', control, ': ', nrow(hallmark_res), ' pathways')
    }
    hallmark_res
  }, error = function(e) {
    log_msg('Hallmark GSEA failed for ', case, ' vs ', control, ': ', e$message)
    NULL
  })

  # === Gene-level volcano/MA (peak-to-gene aggregation) ===
  gene_df <- tryCatch({
    suppressPackageStartupMessages(library(org.Hs.eg.db)); suppressPackageStartupMessages(library(org.Mm.eg.db))
    gdf <- ann |> dplyr::filter(!is.na(SYMBOL), SYMBOL != '', !is.na(logFC))
    if (nrow(gdf) < 5) { NULL }
    else {
      gdf_agg <- gdf |>
        dplyr::group_by(SYMBOL) |>
        dplyr::summarise(
          mean_logFC = mean(logFC, na.rm = TRUE),
          best_logFC = logFC[which.max(abs(logFC))],
          best_FDR = dplyr::coalesce(dplyr::first(FDR[which.max(abs(logFC))]), NA_real_),
          best_PValue = dplyr::coalesce(dplyr::first(PValue[which.max(abs(logFC))]), NA_real_),
          n_peaks = dplyr::n(),
          mean_baseMean = mean(baseMean, na.rm = TRUE),
          .groups = 'drop'
        ) |>
        dplyr::arrange(dplyr::desc(abs(mean_logFC)))
      write.csv(gdf_agg, file.path(croot, 'results', 'gene_level_differential.csv'), row.names = FALSE, quote = FALSE)

      gdf_agg$plot_y <- dplyr::coalesce(gdf_agg$best_FDR, gdf_agg$best_PValue)
      gdf_agg$plot_y[is.na(gdf_agg$plot_y) | gdf_agg$plot_y <= 0] <- 1
      gdf_agg$neg_log10_p <- -log10(gdf_agg$plot_y)
      gdf_agg$log2_baseMean <- log2(pmax(gdf_agg$mean_baseMean, 1, na.rm = TRUE))
      exploratory_mode_g <- all(is.na(gdf_agg$best_FDR)) && all(is.na(gdf_agg$best_PValue))
      if (exploratory_mode_g) {
        gdf_agg$direction <- dplyr::case_when(
          abs(gdf_agg$mean_logFC) >= lfc_cutoff & gdf_agg$mean_logFC > 0 & gdf_agg$mean_baseMean >= baseMean_min ~ 'Up',
          abs(gdf_agg$mean_logFC) >= lfc_cutoff & gdf_agg$mean_logFC < 0 & gdf_agg$mean_baseMean >= baseMean_min ~ 'Down',
          TRUE ~ 'NS'
        )
      } else {
        gdf_agg$direction <- dplyr::case_when(
          abs(gdf_agg$mean_logFC) >= lfc_cutoff & !is.na(gdf_agg$best_FDR) & gdf_agg$best_FDR < padj_cutoff & gdf_agg$mean_logFC > 0 ~ 'Up',
          abs(gdf_agg$mean_logFC) >= lfc_cutoff & !is.na(gdf_agg$best_FDR) & gdf_agg$best_FDR < padj_cutoff & gdf_agg$mean_logFC < 0 ~ 'Down',
          TRUE ~ 'NS'
        )
      }
      gdf_agg$point_alpha <- ifelse(gdf_agg$direction == 'NS', 0.22, 0.82)
      gdf_agg$label <- ''

      if (exploratory_mode_g) {
        gdf_agg <- label_by_distance(gdf_agg, 'log2_baseMean', 'mean_logFC', label_col = 'label', label_source_col = 'SYMBOL', top_n = label_top_n)
        p_gm <- ggplot(gdf_agg, aes(x = log2_baseMean, y = mean_logFC)) +
          geom_hline(yintercept = 0, color = 'grey78', linewidth = 0.28) +
          geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 'dashed', color = 'grey65', linewidth = 0.28) +
          geom_point(aes(color = direction, alpha = point_alpha), size = 0.9) +
          scale_color_manual(values = c('Up' = '#d73027', 'Down' = '#2c7fb8', 'NS' = '#bdbdbd'), drop = FALSE) +
          scale_alpha_identity() +
          safe_repel_layer(gdf_agg, 'log2_baseMean', 'mean_logFC', 'label') +
          theme_atac() + theme(legend.position = 'bottom') +
          labs(title = sprintf('Gene-level MA: %s vs %s', case, control), x = 'log2(mean baseMean)', y = 'Mean log2FC')
      } else {
        gdf_agg <- label_by_distance(gdf_agg, 'mean_logFC', 'neg_log10_p', label_col = 'label', label_source_col = 'SYMBOL', top_n = label_top_n)
        p_gm <- ggplot(gdf_agg, aes(x = mean_logFC, y = neg_log10_p)) +
          geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 'dashed', color = 'grey65', linewidth = 0.28) +
          geom_hline(yintercept = -log10(padj_cutoff), linetype = 'dashed', color = 'grey65', linewidth = 0.28) +
          geom_point(aes(color = direction, alpha = point_alpha), size = 0.9) +
          scale_color_manual(values = c('Up' = '#d73027', 'Down' = '#2c7fb8', 'NS' = '#bdbdbd'), drop = FALSE) +
          scale_alpha_identity() +
          safe_repel_layer(gdf_agg, 'mean_logFC', 'neg_log10_p', 'label') +
          theme_atac() + theme(legend.position = 'bottom') +
          labs(title = sprintf('Gene-level volcano: %s vs %s', case, control), x = 'Mean log2FC', y = '-log10(P)')
      }
      save_ggplot_both(file.path(croot, 'figures', 'gene', sprintf('Gene_level_%s_vs_%s.pdf', case, control)),
                       p_gm, width = 7, height = 6)
      gdf_agg
    }
  }, error = function(e) {
    log_msg('Gene-level plot failed for ', case, ' vs ', control, ': ', e$message)
    NULL
  })

  # === Split GO enrichment: Up vs Down genes ===
  if (length(sig_genes) > 0) {
    sig_up_genes <- unique(as.character(sig$geneId[sig$logFC > 0]))
    sig_down_genes <- unique(as.character(sig$geneId[sig$logFC < 0]))
    sig_up_genes <- sig_up_genes[sig_up_genes != '' & !is.na(sig_up_genes)]
    sig_down_genes <- sig_down_genes[sig_down_genes != '' & !is.na(sig_down_genes)]

    for (dir_name in c('up', 'down')) {
      gene_set <- if (dir_name == 'up') sig_up_genes else sig_down_genes
      if (length(gene_set) < 5) next
      ego_dir <- tryCatch(enrichGO(gene = gene_set, OrgDb = orgdb, keyType = 'ENTREZID', ont = 'BP', pAdjustMethod = 'BH', readable = TRUE), error = function(e) NULL)
      if (!is.null(ego_dir) && nrow(as.data.frame(ego_dir)) > 0) {
        dir.create(file.path(subdir, dir_name), recursive = TRUE, showWarnings = FALSE)
        write.csv(as.data.frame(ego_dir), file.path(subdir, dir_name, 'GO_BP_enrichment.csv'), row.names = FALSE, quote = FALSE)
        save_device_both(file.path(subdir, dir_name, 'GO_BP_dotplot.pdf'), width = 9, height = 6,
          plot_fun = function() print(dotplot(ego_dir, showCategory = 20) + ggplot2::ggtitle(sprintf('GO BP %s: %s vs %s', dir_name, case, control))))
      }
    }
  }
}

export_diff_peak_beds <- function(res, outdir, case, control, padj_cutoff = 0.05, lfc_cutoff = 1, baseMean_min = 5) {
  if (is.null(res) || nrow(res) == 0) return(invisible(NULL))
  subdir <- file.path(contrast_root(outdir, case, control), 'results', 'beds')
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  exploratory_mode <- all(is.na(res$FDR)) && all(is.na(res$PValue))
  sig <- res |> dplyr::filter(!is.na(logFC), !is.na(Chr), !is.na(Start), !is.na(End))
  if (exploratory_mode) {
    log_msg('export contrast beds in exploratory mode for ', case, ' vs ', control)
    sig <- sig |> dplyr::filter(abs(logFC) >= lfc_cutoff, !is.na(baseMean), baseMean >= baseMean_min)
  } else {
    sig <- sig |> dplyr::filter(!is.na(FDR), FDR < padj_cutoff, abs(logFC) >= lfc_cutoff)
  }
  up <- sig |> dplyr::filter(logFC > 0) |> dplyr::select(Chr, Start, End, feature_id, logFC)
  down <- sig |> dplyr::filter(logFC < 0) |> dplyr::select(Chr, Start, End, feature_id, logFC)
  write.table(up, file.path(subdir, 'up.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down, file.path(subdir, 'down.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  legacy_dir <- file.path(outdir, 'contrast_beds', contrast_label(case, control))
  dir.create(legacy_dir, recursive = TRUE, showWarnings = FALSE)
  write.table(up, file.path(legacy_dir, 'up.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(down, file.path(legacy_dir, 'down.bed'), sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
  log_msg('exported bed files for ', case, ' vs ', control, ': up= ', nrow(up), ' ; down= ', nrow(down))
}

aggregate_counts_by_label <- function(count_mat, labels) {
  labels <- as.character(labels)
  keep <- !is.na(labels) & labels != ''
  if (!any(keep)) return(NULL)
  mat2 <- count_mat[keep, , drop = FALSE]
  lab2 <- labels[keep]
  split_idx <- split(seq_along(lab2), lab2)
  out <- do.call(rbind, lapply(names(split_idx), function(k) colSums(mat2[split_idx[[k]], , drop = FALSE])))
  rownames(out) <- names(split_idx)
  as.matrix(out)
}

run_group_edgeR <- function(group_counts, meta, case, control) {
  keep_samples <- meta$condition %in% c(case, control)
  meta2 <- meta[keep_samples, , drop = FALSE]
  group_counts <- group_counts[, rownames(meta2), drop = FALSE]
  meta2$condition <- factor(meta2$condition, levels = c(control, case))
  if (nrow(group_counts) == 0 || ncol(group_counts) < 2) return(NULL)
  y <- DGEList(counts = group_counts, group = meta2$condition)
  keep <- rowSums(y$counts) > 0
  y <- y[keep, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) return(NULL)
  y <- calcNormFactors(y)
  lcpm <- cpm(y, log = TRUE, prior.count = 1)
  cpm_lin <- cpm(y, log = FALSE, prior.count = 1)
  rep_count <- table(meta2$condition)
  if (any(rep_count < 2)) {
    ctrl <- rowMeans(lcpm[, meta2$condition == control, drop = FALSE])
    trt <- rowMeans(lcpm[, meta2$condition == case, drop = FALSE])
    res <- data.frame(group_id = rownames(lcpm), logFC = trt - ctrl, logCPM = rowMeans(lcpm), baseMean = rowMeans(cpm_lin), meanCount = rowMeans(y$counts), PValue = NA_real_, FDR = NA_real_, stringsAsFactors = FALSE)
  } else {
    design <- model.matrix(~ condition, data = meta2)
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    qlf <- glmQLFTest(fit, coef = 2)
    res <- topTags(qlf, n = Inf)$table |> as.data.frame()
    res$group_id <- rownames(res)
    res$baseMean <- rowMeans(cpm_lin[rownames(res), , drop = FALSE])
    res$meanCount <- rowMeans(y$counts[rownames(res), , drop = FALSE])
    res <- res[, c('group_id', 'logFC', 'logCPM', 'baseMean', 'meanCount', 'PValue', 'FDR')]
  }
  res
}

plot_group_diff <- function(df, out_pdf, title, label_top_n = 40, color_col = NULL) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  exploratory_mode <- all(is.na(df$FDR)) && all(is.na(df$PValue))
  df$label <- ''
  use_color_group <- !is.null(color_col) && color_col %in% colnames(df)
  color_title <- if (use_color_group) gsub("_", " ", color_col) else ""
  hide_large_color_legend <- use_color_group && length(unique(na.omit(df[[color_col]]))) > 15
  if (exploratory_mode) {
    df$plot_y <- log10(pmax(dplyr::coalesce(df$baseMean, 0), 0) + 1)
    df$direction <- dplyr::case_when(abs(df$logFC) >= 1 & df$logFC > 0 ~ "Up",
                                     abs(df$logFC) >= 1 & df$logFC < 0 ~ "Down",
                                     TRUE ~ "NS")
    df$point_alpha <- ifelse(df$direction == "NS", 0.22, 0.82)
    df <- label_by_distance(df, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'group_id', top_n = label_top_n)
    if (use_color_group) {
      p <- ggplot(df, aes(x = logFC, y = plot_y, color = .data[[color_col]], alpha = point_alpha)) +
        geom_point(size = 1.15) +
        labs(color = color_title)
    } else {
      p <- ggplot(df, aes(x = logFC, y = plot_y, color = direction, alpha = point_alpha)) +
        geom_point(size = 1.15) +
        scale_color_manual(values = c("Up" = "#d73027", "Down" = "#2c7fb8", "NS" = "#bdbdbd"), drop = FALSE) +
        labs(color = "")
    }
    p <- p +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.28) +
      scale_alpha_identity() +
      safe_repel_layer(df, 'logFC', 'plot_y', 'label') +
      theme_atac() +
      labs(title = title, x = 'log2FC', y = 'log10(baseMean + 1)')
  } else {
    df$plot_p <- dplyr::coalesce(df$FDR, df$PValue)
    df$plot_p[is.na(df$plot_p) | df$plot_p <= 0] <- 1
    df$plot_y <- -log10(df$plot_p)
    df$direction <- dplyr::case_when(!is.na(df$FDR) & df$FDR < 0.05 & df$logFC > 0 & abs(df$logFC) >= 1 ~ "Up",
                                     !is.na(df$FDR) & df$FDR < 0.05 & df$logFC < 0 & abs(df$logFC) >= 1 ~ "Down",
                                     TRUE ~ "NS")
    df$point_alpha <- ifelse(df$direction == "NS", 0.22, 0.82)
    df <- label_by_distance(df, 'logFC', 'plot_y', label_col = 'label', label_source_col = 'group_id', top_n = label_top_n)
    if (use_color_group) {
      p <- ggplot(df, aes(x = logFC, y = plot_y, color = .data[[color_col]], alpha = point_alpha)) +
        geom_point(size = 1.15) +
        labs(color = color_title)
    } else {
      p <- ggplot(df, aes(x = logFC, y = plot_y, color = direction, alpha = point_alpha)) +
        geom_point(size = 1.15) +
        scale_color_manual(values = c("Up" = "#d73027", "Down" = "#2c7fb8", "NS" = "#bdbdbd"), drop = FALSE) +
        labs(color = "")
    }
    p <- p +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey65", linewidth = 0.28) +
      scale_alpha_identity() +
      safe_repel_layer(df, 'logFC', 'plot_y', 'label') +
      theme_atac() +
      labs(title = title, x = 'log2FC', y = '-log10(FDR/PValue)')
  }
  if (hide_large_color_legend) {
    p <- p + guides(color = "none") +
      labs(subtitle = sprintf("Points are colored by %s; legend hidden because there are many groups.", color_title))
  }
  save_ggplot_both(out_pdf, p, width = 7, height = 5.5)
}

run_te_group_analysis <- function(res, mat, meta, outdir, case, control, te_violin_class = '', te_violin_top_n = 12, te_family_filter = '', te_name_filter = '', label_top_n = 40) {
  croot <- contrast_root(outdir, case, control)
  if (is.null(res) || nrow(res) == 0 || !'is_TE' %in% colnames(res)) return(invisible(NULL))
  te_res <- res |>
    dplyr::filter(is_TE %in% TRUE) |>
    dplyr::mutate(
      te_name1 = first_token(te_name),
      te_family1 = first_token(te_family),
      te_class1 = first_token(te_class)
    )
  if (nrow(te_res) == 0) {
    log_msg('no TE-annotated peaks available for TE group analysis: ', case, ' vs ', control)
    return(invisible(NULL))
  }

  if (!is.null(te_violin_class) && te_violin_class != '') {
    te_res <- te_res |> dplyr::filter(grepl(te_violin_class, te_class1, ignore.case = TRUE))
  }
  if (!is.null(te_family_filter) && te_family_filter != '') {
    te_res <- te_res |> dplyr::filter(grepl(te_family_filter, te_family1, ignore.case = TRUE))
  }
  if (!is.null(te_name_filter) && te_name_filter != '') {
    te_res <- te_res |> dplyr::filter(grepl(te_name_filter, te_name1, ignore.case = TRUE))
  }
  if (nrow(te_res) == 0) {
    log_msg('no TE peaks remain after TE filters: ', case, ' vs ', control)
    return(invisible(NULL))
  }

  subdir <- file.path(croot, 'differential', 'TE')
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  results_dir(croot)
  invisible(figures_dir(croot, 'TE'))
  write.csv(te_res, file.path(croot, 'results', 'TE_peak_annotation_filtered.csv'), row.names = FALSE, quote = FALSE)

  # 保留原有的 class/family/subfamily 聚合，但它们现在是在用户指定的 TE 子集内汇总，而不是全 TE 全局汇总
  level_map <- list(class = 'te_class1', family = 'te_family1', subfamily = 'te_name1')
  for (lev in names(level_map)) {
    coln <- level_map[[lev]]
    agg <- aggregate_counts_by_label(mat[te_res$feature_id, , drop = FALSE], te_res[[coln]])
    if (is.null(agg) || nrow(agg) == 0) next
    diff_df <- run_group_edgeR(agg, meta, case, control)
    if (is.null(diff_df) || nrow(diff_df) == 0) next
    color_col <- NULL
    if (lev == 'family') {
      color_map_df <- te_res |>
        dplyr::filter(!is.na(te_family1), te_family1 != '') |>
        dplyr::count(te_family1, te_class1, sort = TRUE) |>
        dplyr::group_by(te_family1) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::transmute(group_id = te_family1, TE_class = te_class1)
      diff_df <- dplyr::left_join(diff_df, color_map_df, by = 'group_id')
      color_col <- 'TE_class'
    } else if (lev == 'subfamily') {
      color_map_df <- te_res |>
        dplyr::filter(!is.na(te_name1), te_name1 != '') |>
        dplyr::count(te_name1, te_family1, sort = TRUE) |>
        dplyr::group_by(te_name1) |>
        dplyr::slice_head(n = 1) |>
        dplyr::ungroup() |>
        dplyr::transmute(group_id = te_name1, TE_family = te_family1)
      diff_df <- dplyr::left_join(diff_df, color_map_df, by = 'group_id')
      color_col <- 'TE_family'
    }
    write.csv(diff_df, file.path(croot, 'results', sprintf('TE_%s_differential.csv', lev)), row.names = FALSE, quote = FALSE)
    plot_group_diff(diff_df, file.path(croot, 'figures', 'TE', sprintf('TE_%s_%s_vs_%s.pdf', lev, case, control)),
                    sprintf('TE %s differential: %s vs %s', lev, case, control),
                    label_top_n = label_top_n, color_col = color_col)
  }

  # 关键新增：对筛选后的 TE-overlapping peaks 本身画 effect/volcano，体现“指定某类/某家族 TE 上的 peak 可及性变化”
  te_peak_df <- te_res
  te_peak_df$label <- dplyr::coalesce(te_peak_df$te_name1, te_peak_df$feature_id)
  exploratory_mode <- all(is.na(te_peak_df$FDR)) && all(is.na(te_peak_df$PValue))
  if (exploratory_mode) {
    te_peak_df$plot_y <- log10(pmax(dplyr::coalesce(te_peak_df$baseMean, 0), 0) + 1)
    te_peak_df <- label_by_distance(te_peak_df, 'logFC', 'plot_y', label_col = 'label_show', label_source_col = 'label', top_n = label_top_n)
    p_te_peak <- ggplot(te_peak_df, aes(x = logFC, y = plot_y, color = te_class1)) +
      geom_point(alpha = 0.7) +
      safe_repel_layer(te_peak_df, 'logFC', 'plot_y', 'label_show') +
      theme_bw() +
      labs(title = sprintf('TE-overlapping peak effect: %s vs %s', case, control), x = 'log2FC', y = 'log10(baseMean + 1)', color = 'TE class')
    save_ggplot_both(file.path(croot, 'figures', 'TE', sprintf('TE_peak_effect_%s_vs_%s.pdf', case, control)), p_te_peak, width = 8, height = 6)
  } else {
    te_peak_df$plot_p <- dplyr::coalesce(te_peak_df$FDR, te_peak_df$PValue)
    te_peak_df$plot_p[is.na(te_peak_df$plot_p) | te_peak_df$plot_p <= 0] <- 1
    te_peak_df$plot_y <- -log10(te_peak_df$plot_p)
    te_peak_df <- label_by_distance(te_peak_df, 'logFC', 'plot_y', label_col = 'label_show', label_source_col = 'label', top_n = label_top_n)
    p_te_peak <- ggplot(te_peak_df, aes(x = logFC, y = plot_y, color = te_class1)) +
      geom_point(alpha = 0.7) +
      safe_repel_layer(te_peak_df, 'logFC', 'plot_y', 'label_show') +
      theme_bw() +
      labs(title = sprintf('TE-overlapping peak volcano: %s vs %s', case, control), x = 'log2FC', y = '-log10(FDR/PValue)', color = 'TE class')
    save_ggplot_both(file.path(croot, 'figures', 'TE', sprintf('TE_peak_volcano_%s_vs_%s.pdf', case, control)), p_te_peak, width = 8, height = 6)
  }

  top_te <- te_res |>
    dplyr::filter(!is.na(te_name1), te_name1 != '') |>
    dplyr::group_by(te_name1) |>
    dplyr::summarise(score = max(abs(logFC), na.rm = TRUE), .groups = 'drop') |>
    dplyr::arrange(dplyr::desc(score)) |>
    dplyr::slice_head(n = te_violin_top_n)
  if (nrow(top_te) == 0) return(invisible(NULL))
  sel <- te_res |> dplyr::filter(te_name1 %in% top_te$te_name1)
  ids <- intersect(sel$feature_id, rownames(mat))
  if (length(ids) < 2) return(invisible(NULL))
  lcpm <- edgeR::cpm(mat[ids, rownames(meta), drop = FALSE], log = TRUE, prior.count = 1)
  plot_df <- as.data.frame(lcpm) |>
    tibble::rownames_to_column('feature_id') |>
    dplyr::left_join(sel |> dplyr::select(feature_id, te_name1, te_family1, te_class1) |> dplyr::distinct(), by = 'feature_id') |>
    tidyr::pivot_longer(cols = all_of(rownames(meta)), names_to = 'sample', values_to = 'logCPM') |>
    dplyr::left_join(meta, by = 'sample')

  p_violin <- ggplot(plot_df, aes(x = condition, y = logCPM, fill = condition)) +
    geom_violin(scale = 'width', trim = TRUE) +
    facet_wrap(~ te_name1, scales = 'free_y') +
    theme_bw() +
    labs(title = sprintf('Top changed TE subfamilies: %s vs %s', case, control), x = '', y = 'Peak-level logCPM')
  fn <- sprintf('TE_top%d_violin_%s_vs_%s.pdf', te_violin_top_n, case, control)
  save_ggplot_both(file.path(subdir, fn), p_violin, width = 12, height = 8)
}


# Load visualization upgrades (TE Family coloring, msigdbr, better themes)
vis_upgrade_file <- file.path(dirname(normalizePath(system.file(package = "atacseq"), mustWork = FALSE)), "atac_visuals_upgraded.R")
vis_upgrade_file2 <- file.path(atac_pipeline_root, 'atacseq-downstream', 'atac_visuals_upgraded.R')
vis_upgrade_file3 <- file.path(getwd(), 'atac_visuals_upgraded.R')
for (vf in c(vis_upgrade_file2, vis_upgrade_file, vis_upgrade_file3)) {
  if (file.exists(vf)) { tryCatch(source(vf), error = function(e) message("[visuals] load failed: ", e$message)); break }
}
