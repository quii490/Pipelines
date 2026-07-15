# ==============================================================================
#                      主运行脚本 (run.R)
# ==============================================================================
config_file <- commandArgs(trailingOnly = TRUE)
if (length(config_file) == 0) {
  stop('[run.R] missing config file path')
}

message('[run.R] loading config: ', config_file[1])
source(config_file[1])

message('[run.R] loading packages')
suppressPackageStartupMessages({
  library(tidyverse)
  library(GenomicRanges)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(viridis)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(DESeq2)
  library(patchwork)
  library(jsonlite)
})
message('[run.R] packages loaded')

cmd_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub('--file=', '', grep('--file=', cmd_args, value = TRUE))
if (length(script_path) != 1) {
  stop('[run.R] unable to determine script path')
}
this_dir <- dirname(normalizePath(script_path))

message('[run.R] starting analysis pipeline')
source(file.path(this_dir, 'scripts', '01_prepare_data.R'))
source(file.path(this_dir, 'scripts', '02_normalize_and_process.R'))
source(file.path(this_dir, '..', '..', 'count_draw', 'scripts', '03_generate_visuals.R'))

# Differential TE-family analysis is deliberately gated by replicate count.
# DESeq2 with a single library in either condition is not a valid inferential
# test, so the core count/visual pipeline remains successful while the report
# records an explicit SKIP for that optional module.
te_deseq2_status <- file.path(output_dir, 'TE_family_DESeq2_status.json')
if (exists('RUN_TE_DESEQ2') && isTRUE(RUN_TE_DESEQ2)) {
  controls_n <- if (exists('CONTROL_SAMPLES')) length(CONTROL_SAMPLES) else 0
  targets_n <- if (exists('TARGET_SAMPLES')) length(TARGET_SAMPLES) else 0
  if (controls_n >= 2 && targets_n >= 2) {
    tryCatch({
      source(file.path(this_dir, 'scripts', '04_te_analysis.R'))
      writeLines(jsonlite::toJSON(list(status = 'PASS', controls = controls_n, targets = targets_n, generated_utc = format(Sys.time(), tz = 'UTC')), auto_unbox = TRUE, pretty = TRUE), te_deseq2_status)
    }, error = function(e) {
      writeLines(jsonlite::toJSON(list(status = 'FAIL', reason = conditionMessage(e), generated_utc = format(Sys.time(), tz = 'UTC')), auto_unbox = TRUE, pretty = TRUE), te_deseq2_status)
      message('[run.R] optional TE DESeq2 failed: ', conditionMessage(e))
    })
  } else {
    writeLines(jsonlite::toJSON(list(status = 'SKIP', reason = 'at least two controls and two targets are required', controls = controls_n, targets = targets_n, generated_utc = format(Sys.time(), tz = 'UTC')), auto_unbox = TRUE, pretty = TRUE), te_deseq2_status)
    message('[run.R] optional TE DESeq2 skipped: insufficient replicates')
  }
} else {
  writeLines(jsonlite::toJSON(list(status = 'SKIP', reason = 'RUN_TE_DESEQ2=false', generated_utc = format(Sys.time(), tz = 'UTC')), auto_unbox = TRUE, pretty = TRUE), te_deseq2_status)
}
message('[run.R] analysis finished')
