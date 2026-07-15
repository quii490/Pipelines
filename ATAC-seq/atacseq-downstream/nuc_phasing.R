suppressPackageStartupMessages({
  library(Rsamtools)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
  library(optparse)
})

log_msg <- function(...) message("[nuc_phasing] ", paste(..., collapse = ""))
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

extract_frag_sizes <- function(bam_path, mapq = 30, max_frag = 1000) {
  log_msg("extracting fragment sizes from ", bam_path)
  sbp <- ScanBamParam(
    flag = scanBamFlag(isPaired = TRUE, isProperPair = TRUE,
                       isFirstMateRead = TRUE),
    mapqFilter = mapq,
    what = "isize"
  )
  isize <- scanBam(bam_path, param = sbp)[[1]][["isize"]]
  isize <- abs(isize[!is.na(isize)])
  isize <- isize[isize <= max_frag & isize >= 1]
  log_msg("  extracted ", length(isize), " proper-pair fragments")
  isize
}

compute_nrl <- function(frag_sizes, lspan = 0.35, rspan = 0.1, min_frag = 51) {
  log_msg("computing NRL with lspan=", lspan, ", rspan=", rspan)
  if (length(frag_sizes) == 0) {
    xs <- seq_len(1)
    empty <- rep(NA_real_, length(xs))
    return(list(nrl = NA_real_, fit1 = NULL, fit2 = NULL, xs = xs, freq = empty,
                cnt = integer(length(xs)), pfreq = empty, rfreq = empty,
                pred2 = empty, peaks = c(mono = NA_real_, di = NA_real_, tri = NA_real_),
                keep = rep(FALSE, length(xs)), note = "no_proper_pair_fragments"))
  }
  cnt <- tabulate(frag_sizes)
  xs <- seq_along(cnt)
  freq <- cnt / sum(cnt)

  keep <- xs >= min_frag & cnt > 0
  if (sum(keep) < 20) {
    warning("too few distinct fragment sizes for NRL computation")
    empty <- rep(NA_real_, length(freq))
    return(list(nrl = NA_real_, fit1 = NULL, fit2 = NULL, xs = xs, freq = freq,
                cnt = cnt, pfreq = empty, rfreq = empty, pred2 = empty,
                peaks = c(mono = NA_real_, di = NA_real_, tri = NA_real_),
                keep = keep, note = "too_few_distinct_fragment_sizes"))
  }

  fit1 <- loess(log(freq[keep]) ~ log(xs[keep]), span = lspan, surface = "direct")
  pfreq <- rep(NA_real_, length(freq))
  pfreq[keep] <- exp(predict(fit1, log(xs[keep])))

  rfreq <- freq - pfreq
  fit2 <- loess(rfreq[keep] ~ xs[keep], span = rspan)

  pred2 <- rep(NA_real_, length(freq))
  pred2[keep] <- predict(fit2, xs[keep])

  nucl_ranges <- list(
    mono = c(157, 227),
    di   = c(314, 454),
    tri  = c(471, 681)
  )
  peaks <- sapply(names(nucl_ranges), function(nm) {
    rng <- nucl_ranges[[nm]]
    idx <- which(xs >= rng[1] & xs <= rng[2] & keep)
    if (length(idx) == 0) return(NA_real_)
    idx[which.max(pred2[idx])]
  })
  names(peaks) <- names(nucl_ranges)
  nrl <- if (all(is.na(peaks))) NA_real_ else {
    mean(xs[peaks] / c(1, 2, 3), na.rm = TRUE)
  }

  list(nrl = nrl, fit1 = fit1, fit2 = fit2, xs = xs, freq = freq,
       cnt = cnt, pfreq = pfreq, rfreq = rfreq, pred2 = pred2,
       peaks = peaks, keep = keep)
}

plot_nuc_phasing <- function(nrl_res, out_pdf, sample_name = "",
                              xlims = c(0, 800)) {
  if (is.null(nrl_res) || length(nrl_res$freq) == 0 || all(is.na(nrl_res$freq))) {
    log_msg("no data to plot for ", sample_name)
    return(invisible(NULL))
  }

  xs    <- nrl_res$xs
  freq  <- nrl_res$freq
  keep  <- nrl_res$keep
  pred2 <- nrl_res$pred2
  peaks <- nrl_res$peaks
  nrl   <- nrl_res$nrl

  in_range <- xs >= xlims[1] & xs <= xlims[2]
  xs_plot    <- xs[in_range]
  rfreq_plot <- (freq - nrl_res$pfreq)[in_range] / 1e-4
  pred2_plot <- pred2[in_range] / 1e-4
  if (!any(is.finite(rfreq_plot)) && !any(is.finite(pred2_plot))) {
    log_msg("insufficient finite phasing values to plot for ", sample_name)
    return(invisible(NULL))
  }

  ylims <- range(c(rfreq_plot[is.finite(rfreq_plot)], pred2_plot[is.finite(pred2_plot)]),
                 na.rm = TRUE)
  if (diff(ylims) == 0) ylims <- c(-1, 1)

  subtitle <- if (!is.na(nrl)) sprintf("NRL = %.1f bp", nrl) else "NRL = NA"

  draw_plot <- function() {
    par(mar = c(4.5, 4.5, 3, 1))
    plot(xs_plot, rfreq_plot, type = "n",
         xlab = "Fragment size (bp)",
         ylab = expression("Residual frequency" %*% 10^{-4}),
         main = paste0("Nucleosome phasing: ", sample_name),
         sub  = subtitle,
         xlim = xlims, ylim = ylims,
         xaxt = "n")
    axis(1, at = seq(0, 800, by = 200))
    lines(xs_plot, rfreq_plot, col = "grey60", lwd = 0.8)
    keep_plot <- in_range & keep
    if (any(keep_plot)) {
      xs_fit <- xs[keep_plot]
      pred_fit <- pred2[keep_plot] / 1e-4
      ord <- order(xs_fit)
      lines(xs_fit[ord], pred_fit[ord], col = "#e41a1c", lwd = 1.5)
    }
    if (!is.null(peaks) && any(!is.na(peaks))) {
      peak_pos <- xs[peaks[!is.na(peaks)]]
      peak_val <- pred2[peaks[!is.na(peaks)]] / 1e-4
      points(peak_pos, peak_val, pch = 21, col = "grey40", bg = "#e41a1c", cex = 1.2)
    }
    legend("topright",
           legend = c("Observed", "Loess fit (span=0.1)"),
           col = c("grey60", "#e41a1c"), lwd = c(0.8, 1.5),
           bty = "n", cex = 0.85)
  }

  pdf(out_pdf, width = 7, height = 5)
  draw_plot()
  dev.off()
  png(sub("\\.pdf$", ".png", out_pdf, ignore.case = TRUE), width = 7, height = 5, units = "in", res = 180)
  draw_plot()
  dev.off()
  log_msg("saved ", out_pdf)
  invisible(out_pdf)
}

run_nuc_phasing <- function(bam_paths, outdir, labels = NULL,
                             mapq = 30, max_frag = 1000,
                             lspan = 0.35, rspan = 0.1, cores = 1) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(labels)) labels <- basename(bam_paths)

  nrl_table <- data.frame(
    sample = character(0), nrl = numeric(0),
    mono_peak = numeric(0), di_peak = numeric(0), tri_peak = numeric(0),
    n_fragments = numeric(0), note = character(0),
    stringsAsFactors = FALSE
  )

  for (i in seq_along(bam_paths)) {
    bam  <- bam_paths[i]
    lab  <- labels[i]
    log_msg("processing ", lab)

    frag <- extract_frag_sizes(bam, mapq = mapq, max_frag = max_frag)
    nrl_res <- compute_nrl(frag, lspan = lspan, rspan = rspan)

    out_pdf <- file.path(outdir, paste0(lab, "_nuc_phasing.pdf"))
    plot_nuc_phasing(nrl_res, out_pdf, sample_name = lab)

    nrl_table[i, "sample"] <- lab
    nrl_table[i, "nrl"]    <- if (is.na(nrl_res$nrl)) NA_real_ else round(nrl_res$nrl, 1)
    nrl_table[i, "n_fragments"] <- length(frag)
    nrl_table[i, "note"] <- nrl_res$note %||% ""
    if (!is.null(nrl_res$peaks)) {
      nrl_table[i, "mono_peak"] <- nrl_res$xs[nrl_res$peaks["mono"]]
      nrl_table[i, "di_peak"]   <- nrl_res$xs[nrl_res$peaks["di"]]
      nrl_table[i, "tri_peak"]  <- nrl_res$xs[nrl_res$peaks["tri"]]
    }
  }

  write.csv(nrl_table, file.path(outdir, "NRL_summary.csv"),
            row.names = FALSE, quote = FALSE)
  log_msg("NRL summary written to ", file.path(outdir, "NRL_summary.csv"))
  invisible(nrl_table)
}

# ---- CLI entry ----
if (!interactive()) {
  option_list <- list(
    make_option("--bam", dest = "bam", type = "character", default = "",
                help = "Single BAM path (alternative to --bam-list)"),
    make_option("--bam-list", dest = "bam_list", type = "character", default = "",
                help = "File with one BAM path per line"),
    make_option("--bam-glob", dest = "bam_glob", type = "character", default = "",
                help = "Glob pattern for BAM files, e.g. /path/*.clean.bam"),
    make_option("--outdir", dest = "outdir", type = "character", default = "."),
    make_option("--mapq", dest = "mapq", type = "integer", default = 30),
    make_option("--max-frag", dest = "max_frag", type = "integer", default = 1000),
    make_option("--lspan", dest = "lspan", type = "double", default = 0.35),
    make_option("--rspan", dest = "rspan", type = "double", default = 0.1),
    make_option("--cores", dest = "cores", type = "integer", default = 1)
  )
  opt <- parse_args(OptionParser(option_list = option_list))

  bam_paths <- c()
  if (opt$bam != "") bam_paths <- opt$bam
  if (opt$bam_list != "" && file.exists(opt$bam_list)) {
    bam_paths <- c(bam_paths, readLines(opt$bam_list))
  }
  if (opt$bam_glob != "") {
    bam_paths <- c(bam_paths, Sys.glob(opt$bam_glob))
  }
  bam_paths <- bam_paths[bam_paths != "" & file.exists(bam_paths)]

  if (length(bam_paths) == 0) stop("No BAM files found. Use --bam, --bam-list, or --bam-glob.")

  run_nuc_phasing(bam_paths, opt$outdir,
                   mapq = opt$mapq, max_frag = opt$max_frag,
                   lspan = opt$lspan, rspan = opt$rspan,
                   cores = opt$cores)
}
