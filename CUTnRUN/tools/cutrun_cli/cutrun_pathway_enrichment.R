#!/usr/bin/env Rscript

raw <- commandArgs(trailingOnly=TRUE)
get_arg <- function(flag, default=NULL) {
  i <- match(flag, raw)
  if (is.na(i) || i == length(raw)) default else raw[[i + 1L]]
}
annotation <- get_arg('--annotation')
out_prefix <- get_arg('--out-prefix')
organism <- get_arg('--organism', 'hg38')
min_genes <- as.integer(get_arg('--min-genes', '10'))
if (is.null(annotation) || is.null(out_prefix)) stop('--annotation and --out-prefix are required')
args <- list(annotation=annotation, out_prefix=out_prefix, organism=organism, min_genes=min_genes)
dir.create(dirname(args$out_prefix), recursive=TRUE, showWarnings=FALSE)

write_status <- function(status, reason, n=0L) {
  writeLines(sprintf('status\t%s\nreason\t%s\ngenes\t%s\n', status, reason, n), paste0(args$out_prefix, '.STATUS.tsv'))
}
if (!file.exists(args$annotation) || file.info(args$annotation)$size == 0) {
  write_status('SKIP', 'annotation missing')
  quit(status=42)
}
if (!requireNamespace('clusterProfiler', quietly=TRUE)) {
  write_status('SKIP', 'clusterProfiler is not installed')
  quit(status=42)
}
org_pkg <- if (args$organism == 'hg38') 'org.Hs.eg.db' else if (args$organism == 'mm39') 'org.Mm.eg.db' else ''
if (!nzchar(org_pkg) || !requireNamespace(org_pkg, quietly=TRUE)) {
  write_status('SKIP', paste0(org_pkg, ' is not installed'))
  quit(status=42)
}

anno <- tryCatch(read.delim(args$annotation, stringsAsFactors=FALSE, check.names=FALSE), error=function(e) NULL)
if (is.null(anno) || !'gene_id' %in% colnames(anno)) {
  write_status('SKIP', 'gene_id column missing')
  quit(status=42)
}
genes <- unique(as.character(anno$gene_id))
genes <- genes[!is.na(genes) & nzchar(genes) & genes != 'NA']
if (length(genes) < args$min_genes) {
  write_status('SKIP', paste0('fewer than ', args$min_genes, ' genes'), length(genes))
  quit(status=42)
}

suppressPackageStartupMessages(library(org_pkg, character.only=TRUE))
org <- get(org_pkg)
key_type <- if (all(grepl('^[0-9]+$', genes))) 'ENTREZID' else 'ENSEMBL'
if (key_type == 'ENSEMBL') genes <- sub('\\..*$', '', genes)
ego <- tryCatch(clusterProfiler::enrichGO(gene=genes, OrgDb=org, keyType=key_type, ont='BP', readable=TRUE), error=function(e) NULL)
if (is.null(ego)) {
  write_status('FAIL', 'enrichGO failed', length(genes))
  quit(status=1)
}
write.csv(as.data.frame(ego), paste0(args$out_prefix, '.GO_BP.csv'), row.names=FALSE)
if (requireNamespace('ReactomePA', quietly=TRUE) && key_type == 'ENTREZID') {
  er <- tryCatch(ReactomePA::enrichPathway(gene=genes, organism=ifelse(args$organism == 'hg38', 'human', 'mouse'), readable=TRUE), error=function(e) NULL)
  if (!is.null(er)) write.csv(as.data.frame(er), paste0(args$out_prefix, '.Reactome.csv'), row.names=FALSE)
}
write_status('PASS', 'GO biological process enrichment completed', length(genes))
