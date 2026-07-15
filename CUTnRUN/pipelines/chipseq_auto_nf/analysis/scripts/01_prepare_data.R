message('[01_prepare_data] start')

# Older count_draw configurations did not include these cache controls.  Keep
# them backward-compatible while making new runs fingerprint-aware.
if (!exists('ANALYSIS_FINGERPRINT')) {
  ANALYSIS_FINGERPRINT <- 'legacy_unfingerprinted'
}
if (!exists('cache_files') || !is.list(cache_files)) {
  stop('[01_prepare_data] cache_files must be a named list')
}
if (!'fingerprint' %in% names(cache_files)) {
  cache_files$fingerprint <- file.path(intermediate_dir, 'analysis_fingerprint.txt')
}

dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, 'tables/gene'), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, 'tables/te'), recursive = TRUE, showWarnings = FALSE)

cached_fingerprint <- if (file.exists(cache_files$fingerprint)) {
  trimws(readLines(cache_files$fingerprint, n = 1, warn = FALSE))
} else {
  ''
}
if (!identical(cached_fingerprint, ANALYSIS_FINGERPRINT)) {
  stale_cache <- unname(unlist(cache_files[names(cache_files) != 'fingerprint']))
  stale_cache <- stale_cache[file.exists(stale_cache)]
  if (length(stale_cache) > 0) {
    message('[01_prepare_data] input/code fingerprint changed; invalidating ', length(stale_cache), ' cache files')
    unlink(stale_cache)
  }
}

read_and_filter_counts <- function(counts_path, saf_path, cache_path, feature_type) {
  if (file.exists(cache_path)) {
    message('[01_prepare_data] cache exists: ', cache_path)
    return(invisible(NULL))
  }

  message('[01_prepare_data] reading ', feature_type, ' counts: ', counts_path)
  counts <- read_tsv(counts_path, skip = 1, col_names = TRUE, show_col_types = FALSE)
  counts_df <- counts %>%
    filter(!is.na(Start), !is.na(End)) %>%
    select(Geneid, 7:last_col())

  saf <- read_tsv(
    saf_path,
    col_names = c('GeneID', 'Chr', 'Start', 'End', 'Strand'),
    show_col_types = FALSE
  )
  if (feature_type == 'TE' && any(str_starts(saf$GeneID, 'ENSG') | str_starts(saf$GeneID, 'ENSMUSG'))) {
    stop('[01_prepare_data] TE SAF contains gene IDs: ', saf_path)
  }

  deny <- read_tsv(
    deny_list_file,
    col_names = c('Chr', 'Start', 'End'),
    show_col_types = FALSE
  )
  saf_gr <- makeGRangesFromDataFrame(saf, keep.extra.columns = TRUE)
  deny_gr <- makeGRangesFromDataFrame(deny)
  clean_ids <- unique(saf_gr[!overlapsAny(saf_gr, deny_gr)]$GeneID)
  counts_cleaned <- counts_df %>% filter(Geneid %in% clean_ids)

  message('[01_prepare_data] ', feature_type, ' features retained: ', nrow(counts_cleaned))
  saveRDS(counts_cleaned, cache_path)
}

read_and_filter_counts(gene_counts_file, gene_saf_file, cache_files$gene_counts, 'Gene')
read_and_filter_counts(te_counts_file, te_saf_file, cache_files$te_counts, 'TE')
writeLines(ANALYSIS_FINGERPRINT, cache_files$fingerprint)
