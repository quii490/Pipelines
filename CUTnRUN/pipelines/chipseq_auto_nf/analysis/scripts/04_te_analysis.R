# ==============================================================================
#  04_te_analysis.R — TE Subfamily/Family-level Analysis (Gold Standard)
#
#  Strategy:
#    1. Read fractional copy-level TE counts from featureCounts (-M --fraction)
#    2. Aggregate to subfamily (repName) level → near-integer counts
#    3. Aggregate to family level
#    4. DESeq2 at family level for differential enrichment
#    5. Require explicit control and target samples for DESeq2
#    6. Visualizations: MA plot, volcano, heatmap, correlation
# ==============================================================================
message('[04_te_analysis] start')

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c('figures/DESeq2', 'figures/Correlation', 'figures/Heatmaps', 'figures/Boxplots')) {
  dir.create(file.path(output_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# ---- Step 1: Read and clean TE counts ----
message('[04_te_analysis] reading TE counts: ', te_counts_file)
counts_raw <- read_tsv(te_counts_file, skip = 1, col_names = TRUE, show_col_types = FALSE) %>%
  filter(!is.na(Start) & !is.na(End))

counts_mat <- counts_raw %>%
  select(Geneid, 7:last_col()) %>%
  column_to_rownames(var = 'Geneid') %>%
  as.matrix()

# Clean sample names
sample_info <- tibble(sample_full = colnames(counts_mat)) %>%
  mutate(sample_id = basename(sample_full) %>% str_remove('_te\\.bam$') %>% str_remove('\\.bam$'))
if (anyDuplicated(sample_info$sample_id)) {
  stop('[04_te_analysis] duplicate sample names after BAM suffix removal')
}
colnames(counts_mat) <- sample_info$sample_id
message('[04_te_analysis] samples: ', paste(colnames(counts_mat), collapse = ', '))

# ---- Step 2: Load TE annotations ----
message('[04_te_analysis] reading TE annotations: ', repeat_annotations_file)
te_anno <- read_table(
  repeat_annotations_file,
  skip = 1,
  col_names = c('GeneID', 'repName', 'Class', 'Family', 'milliDiv'),
  show_col_types = FALSE
)
te_anno <- te_anno %>%
  mutate(
    milliDiv = suppressWarnings(as.numeric(milliDiv)),
    age_bin = cut(
      milliDiv,
      breaks = c(-Inf, 50, 100, 150, 200, 250, 300, Inf),
      labels = c('0-50', '51-100', '101-150', '151-200', '201-250', '251-300', '>300'),
      right = TRUE
    )
  )

# Filter to TE classes of interest
if (exists('TE_CLASSES_OF_INTEREST') && length(TE_CLASSES_OF_INTEREST) > 0) {
  te_anno <- te_anno %>% filter(Class %in% TE_CLASSES_OF_INTEREST)
  message('[04_te_analysis] TE classes: ', paste(TE_CLASSES_OF_INTEREST, collapse = ', '))
}

# Blacklist filtering — remove TE copies overlapping blacklist
if (exists('deny_list_file') && file.exists(deny_list_file)) {
  message('[04_te_analysis] blacklist filtering')
  saf_raw <- read_tsv(te_saf_file, col_names = c('GeneID','Chr','Start','End','Strand'), show_col_types = FALSE)
  saf_gr <- makeGRangesFromDataFrame(saf_raw, keep.extra.columns = TRUE)
  deny <- read_tsv(deny_list_file, col_names = c('Chr','Start','End'), show_col_types = FALSE)
  deny_gr <- makeGRangesFromDataFrame(deny)
  clean_ids <- unique(saf_gr[!overlapsAny(saf_gr, deny_gr)]$GeneID)
  counts_mat <- counts_mat[rownames(counts_mat) %in% clean_ids, ]
  message('[04_te_analysis] TE copies retained after blacklist: ', nrow(counts_mat))
}

# ---- Step 3: Aggregate to subfamily (repName) level ----
message('[04_te_analysis] aggregating to subfamily level')
counts_df <- as.data.frame(counts_mat) %>%
  rownames_to_column(var = 'GeneID') %>%
  left_join(te_anno %>% select(GeneID, repName, Class, Family, milliDiv, age_bin), by = 'GeneID') %>%
  filter(!is.na(repName))

# Copy-level annotation summary.  This is independent of DESeq2 and remains
# available when a project has too few biological replicates for inference.
annotation_summary <- counts_df %>%
  group_by(Class, Family, repName, age_bin) %>%
  summarise(
    loci = n(),
    total_count = sum(as.matrix(across(all_of(colnames(counts_mat)))), na.rm = TRUE),
    mean_milliDiv = mean(milliDiv, na.rm = TRUE),
    .groups = 'drop'
  )
write_csv(annotation_summary, file.path(output_dir, 'TE_annotation_summary.csv'))

age_summary <- counts_df %>%
  group_by(Class, age_bin) %>%
  summarise(
    loci = n(),
    total_count = sum(as.matrix(across(all_of(colnames(counts_mat)))), na.rm = TRUE),
    .groups = 'drop'
  )
write_csv(age_summary, file.path(output_dir, 'TE_age_summary.csv'))

if (nrow(age_summary) > 0) {
  p_age <- ggplot(age_summary, aes(x = age_bin, y = total_count, fill = Class)) +
    geom_col(position = 'stack') +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = 'RepeatMasker divergence (milliDiv)', y = 'Summed TE count',
         title = 'TE signal by estimated repeat age') +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'bottom')
  ggsave(file.path(output_dir, 'figures/Boxplots/TE_age_signal.pdf'), p_age, width = 9, height = 6)
}

# Subfamily-level aggregation
subfam_counts <- counts_df %>%
  group_by(repName, Class, Family) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop')

# Family-level aggregation
message('[04_te_analysis] aggregating to family level')
fam_counts <- subfam_counts %>%
  group_by(Class, Family) %>%
  summarise(across(where(is.numeric), sum), .groups = 'drop')

# ---- Step 4: DESeq2 at family level ----
message('[04_te_analysis] running DESeq2 at family level')

# Identify control and target samples
control_samples <- if (exists('CONTROL_SAMPLES')) CONTROL_SAMPLES else character(0)
if (length(control_samples) == 0) {
  control_samples <- grep('igg|input|control|ctrl', colnames(counts_mat), value = TRUE, ignore.case = TRUE)
}
target_samples <- setdiff(colnames(counts_mat), control_samples)
if (length(control_samples) == 0 || length(target_samples) == 0) {
  stop('[04_te_analysis] DESeq2 requires at least one control and one target sample')
}
message('[04_te_analysis] controls: ', paste(control_samples, collapse = ', '))
message('[04_te_analysis] targets : ', paste(target_samples, collapse = ', '))

# Build DESeq2 dataset
fam_mat <- fam_counts %>%
  unite('id', Class, Family, sep = '|') %>%
  column_to_rownames(var = 'id') %>%
  select(any_of(colnames(counts_mat))) %>%
  as.matrix()
fam_mat <- round(fam_mat)  # family-level counts are large enough that rounding is safe

coldata <- data.frame(
  row.names = colnames(fam_mat),
  condition = ifelse(colnames(fam_mat) %in% control_samples, 'control', 'target')
)

dds <- DESeqDataSetFromMatrix(
  countData = fam_mat,
  colData = coldata,
  design = ~ condition
)

# Filter low-count families
keep <- rowSums(counts(dds)) >= MIN_TOTAL_NORMALIZED_COUNTS
dds <- dds[keep, ]
message('[04_te_analysis] families after filter: ', nrow(dds))

# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c('condition', 'target', 'control'), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  rownames_to_column(var = 'family_id') %>%
  separate(family_id, into = c('Class', 'Family'), sep = '\\|') %>%
  arrange(padj)

message('[04_te_analysis] significant TE families (padj < 0.05): ', sum(res_df$padj < 0.05, na.rm = TRUE))

# ---- Step 5: Normalized counts + lfc_over_igg (for visualization) ----
message('[04_te_analysis] computing normalized enrichment')

norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>%
  rownames_to_column(var = 'family_id') %>%
  separate(family_id, into = c('Class', 'Family'), sep = '\\|')

igg_mean <- norm_counts %>%
  select(Family, Class, any_of(control_samples)) %>%
  pivot_longer(cols = any_of(control_samples), names_to = 'sample', values_to = 'count') %>%
  group_by(Family, Class) %>%
  summarise(igg_mean = mean(count, na.rm = TRUE), .groups = 'drop')

enrichment <- norm_counts %>%
  pivot_longer(cols = any_of(target_samples), names_to = 'sample', values_to = 'norm_count') %>%
  left_join(igg_mean, by = c('Family', 'Class')) %>%
  mutate(lfc_over_igg = log2((norm_count + 1) / (igg_mean + 1))) %>%
  left_join(res_df %>% select(Class, Family, log2FoldChange, padj), by = c('Class', 'Family'))

# ---- Step 6: Visualization ----
message('[04_te_analysis] generating TE visualizations')

# 6a. MA plot
p_ma <- res_df %>%
  mutate(significant = !is.na(padj) & padj < 0.05) %>%
  ggplot(aes(x = log10(baseMean + 1), y = log2FoldChange, color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  scale_color_manual(values = c('grey70', 'firebrick')) +
  geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) +
  labs(x = 'log10(mean normalized count + 1)', y = 'log2 Fold Change',
       title = 'TE Family — DESeq2 MA Plot') +
  theme_minimal(base_size = 12) +
  theme(legend.position = 'bottom')
ggsave(file.path(output_dir, 'figures/DESeq2/MA_plot.pdf'), p_ma, width = 8, height = 7)

# 6b. Volcano plot
p_volcano <- res_df %>%
  mutate(significant = !is.na(padj) & padj < 0.05,
         label = ifelse(significant, Family, '')) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = significant)) +
  geom_point(alpha = 0.6, size = 1) +
  geom_text_repel(aes(label = label), max.overlaps = 30, size = 3) +
  scale_color_manual(values = c('grey70', 'firebrick')) +
  geom_vline(xintercept = 0, linetype = 'dashed', alpha = 0.5) +
  labs(x = 'log2 Fold Change', y = '-log10(p-value)',
       title = 'TE Family — Volcano Plot') +
  theme_minimal(base_size = 12)
ggsave(file.path(output_dir, 'figures/DESeq2/volcano.pdf'), p_volcano, width = 8, height = 7)

# 6c. Top significant families heatmap
sig_fams <- res_df %>% filter(padj < 0.05) %>% arrange(padj) %>% head(40)
if (nrow(sig_fams) > 1) {
  heat_mat <- enrichment %>%
    filter(Family %in% sig_fams$Family) %>%
    select(Family, sample, lfc_over_igg) %>%
    pivot_wider(names_from = sample, values_from = lfc_over_igg) %>%
    column_to_rownames(var = 'Family') %>%
    as.matrix()

  p_heat <- Heatmap(
    heat_mat,
    name = 'log2(FC over IgG)',
    row_title = 'TE Families',
    column_title = 'Samples',
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_dend = TRUE,
    row_names_gp = gpar(fontsize = 7),
    col = colorRamp2(c(-2, 0, 2), c('steelblue', 'white', 'firebrick'))
  )
  pdf(file.path(output_dir, 'figures/Heatmaps/significant_families.pdf'), width = 10, height = 12)
  draw(p_heat)
  dev.off()
}

# 6d. Family-level correlation scatter (top pairs)
top_fams <- res_df %>% arrange(pvalue) %>% head(10) %>% pull(Family)
enrich_wide <- enrichment %>%
  select(Family, Class, sample, lfc_over_igg) %>%
  pivot_wider(names_from = sample, values_from = lfc_over_igg)

if (ncol(enrich_wide) >= 4) {
  sample_cols <- setdiff(colnames(enrich_wide), c('Family', 'Class'))
  pairs_df <- combn(sample_cols, 2, simplify = FALSE)

  for (pair in pairs_df) {
    p_cor <- enrich_wide %>%
      ggplot(aes(x = .data[[pair[1]]], y = .data[[pair[2]]])) +
      geom_point(alpha = 0.5, size = 1) +
      geom_smooth(method = 'lm', se = TRUE, alpha = 0.2) +
      geom_text_repel(
        data = enrich_wide %>% filter(Family %in% top_fams),
        aes(label = Family), size = 3, max.overlaps = Inf
      ) +
      labs(x = pair[1], y = pair[2],
           title = paste('TE Family Correlation:', pair[1], 'vs', pair[2])) +
      theme_minimal(base_size = 11)
    ggsave(
      file.path(output_dir, 'figures/Correlation',
                paste0(pair[1], '_vs_', pair[2], '.pdf')),
      p_cor, width = 8, height = 7
    )
  }
}

# ---- Step 7: Write results ----
message('[04_te_analysis] writing results tables')

# DESeq2 results
write_csv(res_df, file.path(output_dir, 'TE_family_DESeq2_results.csv'))

# Enrichment table
write_csv(enrichment, file.path(output_dir, 'TE_family_enrichment.csv'))

# Summary of significant hits
sig_summary <- res_df %>%
  filter(padj < 0.05) %>%
  mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) %>%
  group_by(Class, direction) %>%
  summarise(n = n(), .groups = 'drop') %>%
  pivot_wider(names_from = direction, values_from = n, values_fill = 0)
write_csv(sig_summary, file.path(output_dir, 'TE_family_significant_summary.csv'))

message('[04_te_analysis] done')
message('[04_te_analysis] output in: ', output_dir)
