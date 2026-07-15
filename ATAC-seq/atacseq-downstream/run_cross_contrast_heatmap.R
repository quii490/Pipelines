suppressPackageStartupMessages(library(optparse))

script_path <- normalizePath(
  sub('^--file=', '', commandArgs(trailingOnly = FALSE)[grep('^--file=', commandArgs(trailingOnly = FALSE))])
)
Sys.setenv(ATAC_PIPELINE_ROOT = dirname(dirname(script_path)))
source(file.path(dirname(script_path), 'atac_functions.R'))

option_list <- list(
  make_option('--level-dir', dest = 'level_dir', type = 'character',
              help = 'Existing peak_level or bin_level downstream directory.'),
  make_option('--contrasts', dest = 'contrasts', type = 'character',
              help = 'Comma/semicolon-separated contrast names, e.g. KO_vs_WT,Rescue_vs_KO.'),
  make_option('--top-n', dest = 'top_n', type = 'integer', default = 100),
  make_option('--outdir', dest = 'outdir', type = 'character', default = ''),
  make_option('--output-prefix', dest = 'output_prefix', type = 'character',
              default = 'Cross_contrast_selected'),
  make_option('--annotation-mode', dest = 'annotation_mode', type = 'character',
              default = 'gene_te',
              help = 'gene_te | gene | none. Default: gene_te.')
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$level_dir) || !nzchar(opt$level_dir)) stop('missing --level-dir')
if (is.null(opt$contrasts) || !nzchar(opt$contrasts)) stop('missing --contrasts')
if (opt$top_n < 2) stop('--top-n must be >= 2')
opt$annotation_mode <- match.arg(opt$annotation_mode, c('gene_te', 'gene', 'none'))

level_dir <- normalizePath(opt$level_dir, mustWork = FALSE)
outdir <- if (nzchar(opt$outdir)) normalizePath(opt$outdir, mustWork = FALSE) else level_dir
requested <- normalize_cross_contrasts(opt$contrasts)

read_first_csv <- function(candidates, label) {
  candidates <- unique(candidates[file.exists(candidates)])
  if (length(candidates) == 0) {
    stop(label, ' not found. Checked:\n', paste(candidates, collapse = '\n'))
  }
  log_msg('read ', label, ': ', candidates[[1]])
  read.csv(candidates[[1]], check.names = FALSE, stringsAsFactors = FALSE)
}

read_annotation_csv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) == 0) return(data.frame())
  # Old pipeline tables were unquoted; inspect the header rather than any
  # data row, because annotation text may contain literal quote characters.
  if (grepl('^\\s*"', lines[[1]])) {
    return(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  header <- strsplit(lines[[1]], ',', fixed = TRUE)[[1]]
  ann_idx <- match('annotation', header)
  expected <- length(header)
  if (is.na(ann_idx)) {
    return(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  split_unquoted_csv_line <- function(line) {
    fields <- strsplit(line, ',', fixed = TRUE)[[1]]
    if (grepl(',$', line)) fields <- c(fields, '')
    fields
  }
  rows <- lapply(lines[-1], function(line) {
    fields <- split_unquoted_csv_line(line)
    if (length(fields) > expected) {
      after_n <- expected - ann_idx
      ann_end <- length(fields) - after_n
      fields <- c(
        fields[seq_len(ann_idx - 1)],
        paste(fields[ann_idx:ann_end], collapse = ','),
        fields[(ann_end + 1):length(fields)]
      )
    }
    if (length(fields) < expected) fields <- c(fields, rep('', expected - length(fields)))
    fields[seq_len(expected)]
  })
  out <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  colnames(out) <- header
  out
}

result_candidates <- function(contrast) {
  c(
    file.path(level_dir, 'results', 'differential', 'regions',
              paste0(contrast, '.differential_regions.csv')),
    file.path(level_dir, 'legacy', 'raw_r_output', 'contrasts', contrast,
              'differential', paste0('differential_peaks_', contrast, '.csv')),
    file.path(level_dir, 'contrasts', contrast, 'differential',
              paste0('differential_peaks_', contrast, '.csv'))
  )
}

annotation_candidates <- function(contrast) {
  c(
    file.path(level_dir, 'results', 'differential', 'annotation',
              paste0(contrast, '.peak_annotation_all.csv')),
    file.path(level_dir, 'legacy', 'raw_r_output', 'contrasts', contrast,
              'results', 'peak_annotation_all.csv'),
    file.path(level_dir, 'contrasts', contrast, 'results', 'peak_annotation_all.csv')
  )
}

merge_peak_annotation <- function(res, contrast) {
  ann_paths <- unique(annotation_candidates(contrast)[file.exists(annotation_candidates(contrast))])
  if (length(ann_paths) == 0 || !'feature_id' %in% names(res)) return(res)
  ann <- read_annotation_csv(ann_paths[[1]])
  keep <- intersect(c(
    'feature_id', 'annotation', 'geneId', 'distanceToTSS', 'SYMBOL', 'is_TE',
    'te_name', 'te_family', 'te_class', 'peak_class', 'label_name',
    'Chr', 'Start', 'End'
  ), names(ann))
  if (length(keep) <= 1) return(res)
  static <- ann[, keep, drop = FALSE] |> dplyr::distinct(feature_id, .keep_all = TRUE)
  res |>
    dplyr::select(-dplyr::any_of(setdiff(keep, 'feature_id'))) |>
    dplyr::left_join(static, by = 'feature_id')
}

result_list <- list()
for (contrast in requested) {
  res <- read_first_csv(result_candidates(contrast), paste0('result for ', contrast))
  res <- merge_peak_annotation(res, contrast)
  if (!all(c('feature_id', 'logFC') %in% names(res))) {
    stop('result for ', contrast, ' must contain feature_id and logFC')
  }
  result_list[[contrast]] <- res
}

log_msg('selected contrasts: ', paste(names(result_list), collapse = ', '))
plot_cross_contrast_heatmap(
  result_list,
  outdir = outdir,
  selected_contrasts = names(result_list),
  top_n = opt$top_n,
  output_prefix = opt$output_prefix,
  title = sprintf('Selected contrasts: %s (%s annotation)',
                  paste(names(result_list), collapse = ' | '), opt$annotation_mode),
  annotation_mode = opt$annotation_mode
)
log_msg('done; output root: ', outdir)
