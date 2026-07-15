suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(readr)
  library(dplyr)
  library(tidyr)
})

option_list <- list(
  make_option('--qc-summary', type = 'character', dest = 'qc_summary'),
  make_option('--outdir', type = 'character', dest = 'outdir', default = '.')
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$qc_summary) || !file.exists(opt$qc_summary)) stop('missing --qc-summary')
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
save_both <- function(path_pdf, plot, width = 6, height = 4, dpi = 180) {
  ggplot2::ggsave(path_pdf, plot, width = width, height = height)
  ggplot2::ggsave(sub('\\.pdf$', '.png', path_pdf, ignore.case = TRUE), plot, width = width, height = height, dpi = dpi)
}
message('[plot_atac_qc] reading: ', opt$qc_summary)
qc <- read_tsv(opt$qc_summary, show_col_types = FALSE)
if (!'sample' %in% colnames(qc)) stop('qc summary must contain sample column')

to_num <- function(x) suppressWarnings(as.numeric(x))
status_from_cut <- function(x, ideal = NULL, acceptable = NULL, direction = 'high') {
  x <- to_num(x)
  out <- rep('NA', length(x))
  ok <- !is.na(x)
  if (direction == 'high') {
    out[ok & x >= ideal] <- 'PASS'
    out[ok & x < ideal & x >= acceptable] <- 'WARN'
    out[ok & x < acceptable] <- 'FAIL'
  } else {
    out[ok & x <= ideal] <- 'PASS'
    out[ok & x > ideal & x <= acceptable] <- 'WARN'
    out[ok & x > acceptable] <- 'FAIL'
  }
  out
}

if ('frip' %in% colnames(qc)) {
  p <- ggplot(qc, aes(x = sample, y = frip)) + geom_col(fill = '#4C78A8') + coord_flip() + theme_bw() + labs(title = 'FRiP', x = '', y = 'FRiP')
  save_both(file.path(opt$outdir, 'FRiP_barplot.pdf'), p, width = 6, height = 4)
}
if (all(c('mitochondrial_fragments', 'nuclear_fragments_before_blacklist') %in% colnames(qc))) {
  qc <- qc |>
    mutate(mt_fraction = ifelse(
      !is.na(mitochondrial_fragments) & !is.na(nuclear_fragments_before_blacklist) &
        (mitochondrial_fragments + nuclear_fragments_before_blacklist) > 0,
      mitochondrial_fragments / (mitochondrial_fragments + nuclear_fragments_before_blacklist),
      NA_real_
    ))
  p <- ggplot(qc, aes(x = sample, y = mt_fraction)) + geom_col(fill = '#F58518') + coord_flip() + theme_bw() + labs(title = 'Mito fraction', x = '', y = 'mitochondrial fragments / pre-blacklist fragments')
  save_both(file.path(opt$outdir, 'Mito_fraction_barplot.pdf'), p, width = 6, height = 4)
} else if (all(c('total_clean', 'mt_reads') %in% colnames(qc))) {
  qc <- qc |> mutate(mt_fraction = ifelse(total_clean > 0, mt_reads / total_clean, NA_real_))
  p <- ggplot(qc, aes(x = sample, y = mt_fraction)) + geom_col(fill = '#F58518') + coord_flip() + theme_bw() + labs(title = 'Mito fraction (legacy)', x = '', y = 'mt_reads / total_clean')
  save_both(file.path(opt$outdir, 'Mito_fraction_barplot.pdf'), p, width = 6, height = 4)
}

if (all(c('NRF', 'PBC1', 'PBC2') %in% colnames(qc))) {
  complexity_df <- qc |>
    select(sample, NRF, PBC1, PBC2) |>
    mutate(across(c(NRF, PBC1, PBC2), to_num)) |>
    pivot_longer(c(NRF, PBC1, PBC2), names_to = 'metric', values_to = 'value')
  p <- ggplot(complexity_df, aes(x = sample, y = value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.68) +
    coord_flip() +
    theme_bw() +
    labs(title = 'Library complexity', x = '', y = 'Metric value')
  save_both(file.path(opt$outdir, 'Library_complexity_barplot.pdf'), p, width = 7, height = 4.5)
}

if (all(c('nfr_fragments', 'mono_fragments', 'di_fragments', 'tri_plus_fragments') %in% colnames(qc))) {
  frag_df <- qc |>
    select(sample, nfr_fragments, mono_fragments, di_fragments, tri_plus_fragments) |>
    mutate(across(-sample, to_num)) |>
    pivot_longer(-sample, names_to = 'class', values_to = 'fragments') |>
    filter(!is.na(fragments)) |>
    group_by(sample) |>
    mutate(frac = ifelse(sum(fragments) > 0, fragments / sum(fragments), NA_real_)) |>
    ungroup()
  if (nrow(frag_df) > 0) {
    frag_df$class <- factor(frag_df$class, levels = c('nfr_fragments', 'mono_fragments', 'di_fragments', 'tri_plus_fragments'))
    p <- ggplot(frag_df, aes(x = sample, y = frac, fill = class)) +
      geom_col(width = 0.72) +
      coord_flip() +
      theme_bw() +
      labs(title = 'ATAC fragment classes', x = '', y = 'Fraction of PE fragments', fill = 'Class')
    save_both(file.path(opt$outdir, 'Fragment_class_fraction_barplot.pdf'), p, width = 7.5, height = 4.8)
  }
}

qc_pass <- qc |> select(sample, any_of(c('layout', 'frip', 'mt_fraction', 'NRF', 'PBC1', 'PBC2', 'final_nuclear_fragments')))
if ('frip' %in% colnames(qc_pass)) qc_pass$frip_status <- status_from_cut(qc_pass$frip, ideal = 0.30, acceptable = 0.20, direction = 'high')
if ('mt_fraction' %in% colnames(qc_pass)) qc_pass$mt_status <- status_from_cut(qc_pass$mt_fraction, ideal = 0.20, acceptable = 0.50, direction = 'low')
if ('NRF' %in% colnames(qc_pass)) qc_pass$NRF_status <- status_from_cut(qc_pass$NRF, ideal = 0.90, acceptable = 0.80, direction = 'high')
if ('PBC1' %in% colnames(qc_pass)) qc_pass$PBC1_status <- status_from_cut(qc_pass$PBC1, ideal = 0.90, acceptable = 0.50, direction = 'high')
if ('PBC2' %in% colnames(qc_pass)) qc_pass$PBC2_status <- status_from_cut(qc_pass$PBC2, ideal = 10, acceptable = 1, direction = 'high')
status_cols <- grep('_status$', colnames(qc_pass), value = TRUE)
qc_pass$overall_qc <- apply(qc_pass[, status_cols, drop = FALSE], 1, function(x) {
  if (length(x) == 0 || all(x == 'NA')) return('NA')
  if (any(x == 'FAIL')) return('FAIL')
  if (any(x == 'WARN')) return('WARN')
  'PASS'
})
write.table(qc_pass, file.path(opt$outdir, 'qc_pass_fail.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
write.table(qc, file.path(opt$outdir, 'atac_qc_summary_with_metrics.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
