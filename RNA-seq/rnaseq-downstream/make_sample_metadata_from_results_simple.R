args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript make_sample_metadata_from_results_simple.R <results_dir> <output_csv> [sample_split_regex]")
}

results_dir <- args[1]
out_csv <- args[2]
sample_split_regex <- ifelse(length(args) >= 3, args[3], "^(.*)_([^_]+)$")

fc_files <- Sys.glob(file.path(results_dir, "03_gene_featurecounts", "*.featureCounts.txt"))
telocal_files <- Sys.glob(file.path(results_dir, "08_telocal", "*_telocal.cntTable"))

get_sample_from_fc <- function(x) sub("\\.featureCounts\\.txt$", "", basename(x))
get_sample_from_telocal <- function(x) sub("_telocal\\.cntTable$", "", basename(x))

samples <- unique(c(get_sample_from_fc(fc_files), get_sample_from_telocal(telocal_files)))
if (length(samples) == 0) stop("No sample files found under results_dir")

condition <- sub(sample_split_regex, "\\1", samples)
replicate <- sub(sample_split_regex, "\\2", samples)
bad <- condition == samples & replicate == samples
replicate[bad] <- as.character(seq_len(sum(bad)))

meta <- data.frame(
  sample = samples,
  condition = condition,
  replicate = replicate,
  stringsAsFactors = FALSE
)

write.csv(meta, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n")
print(meta)
