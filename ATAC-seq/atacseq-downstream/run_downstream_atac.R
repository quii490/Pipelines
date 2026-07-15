suppressPackageStartupMessages(library(optparse))
script_path <- normalizePath(sub('^--file=', '', commandArgs(trailingOnly = FALSE)[grep('^--file=', commandArgs(trailingOnly = FALSE))]))
Sys.setenv(ATAC_PIPELINE_ROOT = dirname(dirname(script_path)))
source(file.path(dirname(script_path), 'atac_functions.R'))

option_list <- list(
  make_option('--count-file', dest = 'count_file', type = 'character'),
  make_option('--sample-meta', dest = 'sample_meta', type = 'character'),
  make_option('--contrast-file', dest = 'contrast_file', type = 'character', default = ''),
  make_option('--peak-bed', dest = 'peak_bed', type = 'character'),
  make_option('--outdir', dest = 'outdir', type = 'character', default = '.'),
  make_option('--analysis-level', dest = 'analysis_level', type = 'character', default = 'peak'),
  make_option('--species', dest = 'species', type = 'character', default = 'hg38'),
  make_option('--padj-cutoff', dest = 'padj_cutoff', type = 'double', default = 0.05),
  make_option('--lfc-cutoff', dest = 'lfc_cutoff', type = 'double', default = 1),
  make_option('--heatmap-top-n', dest = 'heatmap_top_n', type = 'integer', default = 60),
  make_option('--label-top-n', dest = 'label_top_n', type = 'integer', default = 40),
  make_option('--baseMean-min', dest = 'baseMean_min', type = 'double', default = 5),
  make_option('--gtf', dest = 'gtf', type = 'character', default = ''),
  make_option('--te-bed', dest = 'te_bed', type = 'character', default = ''),
  make_option('--annotation-mode', dest = 'annotation_mode', type = 'character', default = 'gene_te'),
  make_option('--te-violin-class', dest = 'te_violin_class', type = 'character', default = ''),
  make_option('--te-violin-top-n', dest = 'te_violin_top_n', type = 'integer', default = 12),
  make_option('--te-family-filter', dest = 'te_family_filter', type = 'character', default = ''),
  make_option('--te-name-filter', dest = 'te_name_filter', type = 'character', default = '')
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$count_file) || is.na(opt$count_file) || opt$count_file == '') stop('missing required --count-file')
if (is.null(opt$sample_meta) || is.na(opt$sample_meta) || opt$sample_meta == '') stop('missing required --sample-meta')
if (!opt$analysis_level %in% c('peak', 'bin')) stop('--analysis-level must be peak or bin')
opt$annotation_mode <- match.arg(opt$annotation_mode, c('gene_te', 'gene', 'none'))
if (is.null(opt$peak_bed) || is.na(opt$peak_bed) || opt$peak_bed == '') log_msg('peak bed not used in current downstream body: empty --peak-bed')

log_msg('downstream args: count_file=', opt$count_file,
        '; sample_meta=', opt$sample_meta,
        '; contrast_file=', opt$contrast_file %||% '',
        '; peak_bed=', opt$peak_bed %||% '',
        '; outdir=', opt$outdir,
        '; analysis_level=', opt$analysis_level,
        '; species=', opt$species,
        '; annotation_mode=', opt$annotation_mode,
        '; gtf=', opt$gtf %||% '',
        '; te_bed=', opt$te_bed %||% '',
        '; te_violin_class=', opt$te_violin_class %||% '',
        '; te_family_filter=', opt$te_family_filter %||% '',
        '; te_name_filter=', opt$te_name_filter %||% '')

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
log_msg('start downstream')
x <- read_featurecounts_matrix(opt$count_file)
meta <- load_sample_meta(opt$sample_meta, colnames(x$mat))
contrasts <- load_contrasts(opt$contrast_file, meta)

calc_basic_qc(x$mat, meta, opt$outdir)
plot_pca_cor(x$mat, meta, opt$outdir)

if (nrow(contrasts) == 0) {
  log_msg('no valid contrasts found, skip differential analysis')
} else {
  contrast_results <- list()
  for (i in seq_len(nrow(contrasts))) {
    case <- contrasts$case[i]
    control <- contrasts$control[i]
    log_msg('running contrast: ', case, ' vs ', control)
    croot <- contrast_root(opt$outdir, case, control)
    contrast_dir <- file.path(croot, 'differential')
    fig_dir <- file.path(croot, 'figures')
    meta_sub <- meta[meta$condition %in% c(case, control), , drop = FALSE]
    res <- run_edgeR_contrast(x$mat, x$anno, meta, case, control, contrast_dir)
    if (opt$analysis_level == 'peak') {
      res <- annotate_diff_peaks(res, opt$species, opt$outdir, case, control,
                                 padj_cutoff = opt$padj_cutoff,
                                 lfc_cutoff = opt$lfc_cutoff,
                                 baseMean_min = opt$baseMean_min,
                                 te_bed = if (opt$annotation_mode == 'gene_te') opt$te_bed else '',
                                 gtf = if (opt$annotation_mode == 'none') '' else opt$gtf,
                                 annotation_mode = opt$annotation_mode)
    }
    contrast_results[[paste0(case, '_vs_', control)]] <- res
    export_diff_peak_beds(res, opt$outdir, case, control,
                          padj_cutoff = opt$padj_cutoff,
                          lfc_cutoff = opt$lfc_cutoff,
                          baseMean_min = opt$baseMean_min)
    plot_volcano_ma(res, meta_sub, fig_dir, case, control, label_top_n = opt$label_top_n)
    plot_heatmap_top(x$mat, meta_sub, res, case, control, fig_dir, top_n = opt$heatmap_top_n)
    if (opt$analysis_level == 'peak') {
      plot_annotation_volcano(res, fig_dir, case, control, label_top_n = opt$label_top_n)
      run_peak_gene_enrichment(res, opt$species, opt$outdir, case, control,
                               padj_cutoff = opt$padj_cutoff,
                               lfc_cutoff = opt$lfc_cutoff,
                               baseMean_min = opt$baseMean_min,
                               te_bed = if (opt$annotation_mode == 'gene_te') opt$te_bed else '',
                               gtf = if (opt$annotation_mode == 'none') '' else opt$gtf,
                               label_top_n = opt$label_top_n)
      if (opt$annotation_mode == 'gene_te') {
        run_te_group_analysis(res, x$mat, meta, opt$outdir, case, control,
                              te_violin_class = opt$te_violin_class,
                              te_violin_top_n = opt$te_violin_top_n,
                              te_family_filter = opt$te_family_filter,
                              te_name_filter = opt$te_name_filter,
                              label_top_n = opt$label_top_n)
      }
    }
  }
  plot_cross_contrast_overview(
    contrast_results, opt$outdir,
    padj_cutoff = opt$padj_cutoff,
    lfc_cutoff = opt$lfc_cutoff,
    baseMean_min = opt$baseMean_min,
    annotation_mode = opt$annotation_mode
  )
}

log_msg('downstream done')
