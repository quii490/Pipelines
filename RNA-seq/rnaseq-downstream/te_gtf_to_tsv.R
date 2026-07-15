suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript te_gtf_to_tsv.R input.gtf output.tsv")
}

infile <- args[1]
outfile <- args[2]

# 更兼容的读取方式
gtf <- read.delim(
  infile,
  sep = "\t",
  header = FALSE,
  quote = "",
  comment.char = "#",
  stringsAsFactors = FALSE
)

if (ncol(gtf) < 9) {
  stop("Input file does not look like a valid GTF with 9 columns.")
}

colnames(gtf) <- c(
  "chr", "source", "feature", "start", "end",
  "score", "strand", "phase", "attribute"
)

extract_attr <- function(x, key) {
  pat <- paste0(key, ' "([^"]*)"')
  m <- regexec(pat, x)
  sapply(regmatches(x, m), function(z) {
    if (length(z) >= 2) z[2] else NA_character_
  })
}

gtf$gene_id       <- extract_attr(gtf$attribute, "gene_id")
gtf$transcript_id <- extract_attr(gtf$attribute, "transcript_id")
gtf$family_id     <- extract_attr(gtf$attribute, "family_id")
gtf$class_id      <- extract_attr(gtf$attribute, "class_id")
gtf$gene_name     <- extract_attr(gtf$attribute, "gene_name")

out <- data.frame(
  locus_id   = gtf$transcript_id,
  chr        = gtf$chr,
  start      = gtf$start,
  end        = gtf$end,
  strand     = gtf$strand,
  width      = gtf$end - gtf$start + 1,
  feature    = gtf$feature,
  source     = gtf$source,
  score      = gtf$score,
  phase      = gtf$phase,
  repName    = gtf$gene_id,
  repFamily  = gtf$family_id,
  repClass   = gtf$class_id,
  gene_name  = gtf$gene_name,
  stringsAsFactors = FALSE
)

na_locus <- is.na(out$locus_id) | out$locus_id == ""
out$locus_id[na_locus] <- paste0(
  out$chr[na_locus], ":",
  out$start[na_locus], "-",
  out$end[na_locus], "(",
  out$strand[na_locus], ")"
)

write.table(
  out,
  file = outfile,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done. Wrote:", outfile, "\n")
cat("Rows:", nrow(out), "\n")