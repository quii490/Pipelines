message('[02_normalize_and_process] start')

clean_sample_names <- function(names_in) {
  names_out <- basename(names_in) %>%
    str_remove('(_clean_te|_clean|_te|\\.sorted_te|\\.sorted)\\.bam$')
  if (anyDuplicated(names_out)) {
    stop('[02_normalize_and_process] duplicate sample names after BAM suffix removal')
  }
  names_out
}

to_matrix <- function(counts_df) {
  mat <- counts_df %>%
    column_to_rownames(var = 'Geneid') %>%
    as.matrix()
  storage.mode(mat) <- 'numeric'
  colnames(mat) <- clean_sample_names(colnames(mat))
  mat
}

normalize_gene_counts <- function(mat) {
  col_data <- data.frame(row.names = colnames(mat))
  dds <- DESeqDataSetFromMatrix(
    countData = round(mat),
    colData = col_data,
    design = ~1
  )
  dds <- estimateSizeFactors(dds)
  counts(dds, normalized = TRUE)
}

normalize_fractional_counts <- function(mat) {
  # Vectorized geometric means are essential for TE matrices: LZH has 5.6M
  # repeat rows, so an R callback per row can take hours.  The matrix is only
  # eight samples wide; row-wise log means are equivalent and much faster.
  log_mat <- mat
  log_mat[!is.finite(log_mat) | log_mat <= 0] <- NA_real_
  log_mat <- log(log_mat)
  geo_means <- exp(rowMeans(log_mat, na.rm = TRUE))
  geo_means[!is.finite(geo_means)] <- NA_real_
  usable <- is.finite(geo_means) & geo_means > 0
  usable_mat <- mat[usable, , drop = FALSE]
  usable_geo <- geo_means[usable]
  size_factors <- vapply(seq_len(ncol(usable_mat)), function(j) {
    ratios <- usable_mat[, j] / usable_geo
    ratios <- ratios[is.finite(ratios) & ratios > 0]
    if (length(ratios) == 0) return(NA_real_)
    median(ratios)
  }, numeric(1))
  if (any(!is.finite(size_factors) | size_factors <= 0)) {
    library_sizes <- colSums(mat)
    size_factors <- library_sizes / median(library_sizes[library_sizes > 0])
  }
  sweep(mat, 2, size_factors, '/')
}

calculate_enrichment <- function(norm_counts, feature_type) {
  filtered <- norm_counts[rowSums(norm_counts) > MIN_TOTAL_NORMALIZED_COUNTS, , drop = FALSE]
  long_counts <- as.data.frame(filtered) %>%
    rownames_to_column(var = 'geneID') %>%
    pivot_longer(cols = -geneID, names_to = 'sample', values_to = 'norm_count')

  controls <- intersect(CONTROL_SAMPLES, unique(long_counts$sample))
  if (length(controls) == 0) {
    controls <- grep('igg|input|control|ctrl', unique(long_counts$sample), value = TRUE, ignore.case = TRUE)
  }
  if (length(controls) == 0) {
    stop('[02_normalize_and_process] no control sample found for ', feature_type)
  }

  targets <- setdiff(unique(long_counts$sample), controls)
  target_control_map <- if (exists('TARGET_CONTROL_MAP')) TARGET_CONTROL_MAP else character(0)

  bind_rows(lapply(targets, function(target) {
    control <- unname(target_control_map[target])
    if (length(control) == 0 || is.na(control) || control == '') {
      if (length(controls) == 1) {
        control <- controls[[1]]
      } else {
        stop('[02_normalize_and_process] target ', target,
             ' has no unambiguous control for ', feature_type)
      }
    }
    if (!(control %in% controls)) {
      stop('[02_normalize_and_process] mapped control ', control,
           ' for target ', target, ' is absent from count matrix')
    }

    target_counts <- long_counts %>%
      filter(sample == target) %>%
      select(geneID, sample, norm_count)
    control_counts <- long_counts %>%
      filter(sample == control) %>%
      select(geneID, igg_count = norm_count)

    target_counts %>%
      left_join(control_counts, by = 'geneID') %>%
      mutate(
        control_sample = control,
        lfc_over_igg = log2((norm_count + 1) / (igg_count + 1))
      )
  }))
}

if (!file.exists(cache_files$plot_data)) {
  gene_counts <- readRDS(cache_files$gene_counts)
  te_counts <- readRDS(cache_files$te_counts)
  gene_matrix <- to_matrix(gene_counts)
  te_matrix <- to_matrix(te_counts)
  if (!setequal(colnames(gene_matrix), colnames(te_matrix))) {
    stop(
      '[02_normalize_and_process] gene/TE sample mismatch. gene=',
      paste(colnames(gene_matrix), collapse = ','),
      '; TE=',
      paste(colnames(te_matrix), collapse = ',')
    )
  }
  te_matrix <- te_matrix[, colnames(gene_matrix), drop = FALSE]
  gene_norm <- normalize_gene_counts(gene_matrix)
  te_norm <- normalize_fractional_counts(te_matrix)

  gene_enrichment <- calculate_enrichment(gene_norm, 'Gene')
  te_enrichment <- calculate_enrichment(te_norm, 'TE')

  genes <- read_tsv(gene_annotations_file, show_col_types = FALSE,
                    na = c('', 'NA'))
  repeats <- read_tsv(
    repeat_annotations_file,
    skip = 1,
    col_names = c('GeneID', 'repName', 'Class', 'Family', 'milliDiv'),
    show_col_types = FALSE
  )

  gene_plot_data <- gene_enrichment %>%
    left_join(genes, by = c('geneID' = 'GeneID')) %>%
    mutate(
      GeneName = na_if(GeneName, ''),
      GeneBiotype = na_if(GeneBiotype, ''),
      featureType = 'Gene',
      Class = NA_character_,
      Family = GeneBiotype,
      repName = coalesce(GeneName, geneID),
      milliDiv = NA_real_
    ) %>%
    filter(!is.na(Family)) %>%
    select(geneID, sample, control_sample, lfc_over_igg, featureType, Class, Family, repName, milliDiv)

  te_plot_data <- te_enrichment %>%
    left_join(repeats, by = c('geneID' = 'GeneID')) %>%
    mutate(featureType = 'TE') %>%
    filter(!is.na(Family)) %>%
    select(geneID, sample, control_sample, lfc_over_igg, featureType, Class, Family, repName, milliDiv)

  target_samples <- if (exists('TARGET_SAMPLES')) TARGET_SAMPLES else character(0)
  if (length(target_samples) > 0) {
    gene_plot_data <- gene_plot_data %>% filter(sample %in% target_samples)
    te_plot_data <- te_plot_data %>% filter(sample %in% target_samples)
  }

  plot_data <- bind_rows(gene_plot_data, te_plot_data)
  saveRDS(gene_plot_data, cache_files$gene_plot_data)
  saveRDS(te_plot_data, cache_files$te_plot_data)
  saveRDS(plot_data, cache_files$plot_data)
  write_csv(plot_data, file.path(output_dir, 'processed_data_for_plotting.csv'))

  write_csv(gene_enrichment, file.path(output_dir, 'tables/gene/gene_enrichment.csv'))
  write_csv(te_enrichment, file.path(output_dir, 'tables/te/TE_copy_enrichment.csv'))
  te_plot_data %>%
    group_by(sample, Class, Family) %>%
    summarise(mean_lfc_over_igg = mean(lfc_over_igg), n_elements = n(), .groups = 'drop') %>%
    write_csv(file.path(output_dir, 'tables/te/TE_family_enrichment.csv'))
  te_plot_data %>%
    group_by(sample, Class, Family, repName) %>%
    summarise(mean_lfc_over_igg = mean(lfc_over_igg), n_elements = n(), .groups = 'drop') %>%
    write_csv(file.path(output_dir, 'tables/te/TE_subfamily_enrichment.csv'))
}

message('[02_normalize_and_process] unified plot_data ready')
