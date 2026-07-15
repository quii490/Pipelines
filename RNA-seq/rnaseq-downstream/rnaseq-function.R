# ==== 1. 载入包与通用函数 ====
required_pkgs <- c(
  "DESeq2","edgeR","ggplot2","pheatmap","RColorBrewer","data.table","dplyr","tibble",
  "stringr","readr","tidyr","tibble","ggrepel","matrixStats","clusterProfiler",
  "enrichplot","GenomicRanges","rtracklayer","SummarizedExperiment",
  "org.Hs.eg.db","org.Mm.eg.db","AnnotationDbi","msigdbr","ggridges"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("缺少 R 包: ", paste(missing_pkgs, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readr)
  library(tidyr)
  library(ggrepel)
  library(matrixStats)
  library(clusterProfiler)
  library(enrichplot)
  library(GenomicRanges)
  library(rtracklayer)
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  library(ggridges)
})

plot_n_cores <- getOption("plot_n_cores", 1L)
if (!is.numeric(plot_n_cores) || is.na(plot_n_cores) || plot_n_cores < 1) plot_n_cores <- 1L
plot_n_cores <- as.integer(plot_n_cores)
message("[init] plot_n_cores=", plot_n_cores)

safe_run <- function(tag, expr) {
  log_dir <- get0("downstream_log_dir", envir = .GlobalEnv, ifnotfound = outdir)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  log_file <- file.path(log_dir, paste0(tag, "_ERROR.txt"))
  if (file.exists(log_file)) unlink(log_file)
  message("[", tag, "] START")
  ok <- FALSE
  res <- tryCatch(
    {
      val <- force(expr)
      ok <- TRUE
      val
    },
    error = function(e) {
      msg <- paste0("[", tag, "] ERROR: ", conditionMessage(e))
      calls <- tryCatch(vapply(sys.calls(), function(x) paste(deparse(x), collapse = " "), character(1)), error = function(...) character(0))
      payload <- c(msg, if (length(calls) > 0) c("--- sys.calls ---", calls) else NULL)
      message(msg)
      message("[", tag, "] error log: ", log_file)
      writeLines(payload, log_file)
      if (exists("record_module_status", envir = .GlobalEnv)) {
        get("record_module_status", envir = .GlobalEnv)(tag, "failed", conditionMessage(e))
      }
      return(NULL)
    }
  )
  if (ok) {
    message("[", tag, "] DONE")
    if (exists("record_module_status", envir = .GlobalEnv)) {
      get("record_module_status", envir = .GlobalEnv)(tag, "success", "")
    }
  }
  res
}

find_files <- function(path, pattern, recursive = TRUE) {
  if (!dir.exists(path)) return(character(0))
  files <- list.files(path, pattern = pattern, full.names = TRUE, recursive = recursive)
  files[!dir.exists(files)]
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  x
}

as_pdf_path <- function(path) {
  sub("\\.png$", ".pdf", path, ignore.case = TRUE)
}

sample_from_parent <- function(path, suffix_pattern = NULL) {
  x <- basename(dirname(path))
  if (!is.null(suffix_pattern)) x <- sub(suffix_pattern, "", x)
  x
}

normalize_id_vec <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', "", x)
  x <- sub("^\\s+|\\s+$", "", x)

  ens_like <- grepl("^(ENS[A-Z0-9]+|ENSG[0-9]+|ENST[0-9]+|ENSMUSG[0-9]+|ENSMUST[0-9]+)\\.[0-9]+$", x)
  x[ens_like] <- sub("\\.[0-9]+$", "", x[ens_like])

  tolower(x)
}


is_short_repeat_name <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', "", x)
  x <- sub("^\\s+|\\s+$", "", x)
  grepl("^\\([^()]+\\)n$", x)
}

drop_short_repeat_annotation <- function(df, tag = "te_annotation") {
  if (!"repName" %in% colnames(df)) return(df)
  bad <- !is.na(df$repName) & is_short_repeat_name(df$repName)
  if (any(bad)) {
    message("[", tag, "] drop short-repeat annotation rows: ", sum(bad))
    df <- df[!bad, , drop = FALSE]
  }
  df
}

drop_short_repeat_matrix <- function(mat, tag = "matrix") {
  rn <- rownames(mat)
  if (is.null(rn)) return(mat)
  bad <- grepl("^\\([^()]+\\)n($|:|_)", rn)
  if (any(bad)) {
    message("[", tag, "] drop short-repeat rows: ", sum(bad))
    mat <- mat[!bad, , drop = FALSE]
  }
  mat
}

finalize_feature_matrix <- function(df, tag = "matrix") {
  if (!"feature_id" %in% colnames(df)) stop("[", tag, "] 缺少 feature_id 列")
  df$feature_id <- as.character(df$feature_id)
  df$feature_id <- sub("^\\s+|\\s+$", "", df$feature_id)

  bad <- is.na(df$feature_id) | df$feature_id == ""
  if (any(bad)) {
    message("[", tag, "] drop empty feature_id rows: ", sum(bad))
    df <- df[!bad, , drop = FALSE]
  }
  if (nrow(df) == 0) stop("[", tag, "] 过滤后 feature 行数为0")

  count_cols <- setdiff(colnames(df), "feature_id")
  if (length(count_cols) == 0) stop("[", tag, "] 没有计数列")

  for (nm in count_cols) {
    df[[nm]] <- suppressWarnings(as.numeric(df[[nm]]))
    df[[nm]][is.na(df[[nm]])] <- 0
  }

  dup_n <- sum(duplicated(df$feature_id))
  if (dup_n > 0) {
    message("[", tag, "] collapse duplicated feature_id rows: ", dup_n)
    df <- df %>%
      group_by(feature_id) %>%
      summarise(across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
  }

  df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(df) <- df$feature_id
  mat <- as.matrix(df[, count_cols, drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

matrix_to_export_df <- function(mat, id_col = "feature_id") {
  df <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)
  if (is.null(rownames(df)) || all(rownames(df) %in% as.character(seq_len(nrow(df))))) {
    message("[matrix_to_export_df] warning: rownames 缺失或为默认数字，导出时将保留当前行号")
  }
  df <- tibble::rownames_to_column(df, var = id_col)
  rownames(df) <- NULL
  df
}

write_matrix_csv <- function(mat, file, id_col = "feature_id") {
  df <- matrix_to_export_df(mat, id_col = id_col)
  write.csv(df, file, row.names = FALSE)
}

read_count_files_generic <- function(files, sample_fun, id_col = 1, count_col = NULL, tag = "generic") {
  mats <- lapply(files, function(f) {
    dt <- fread(f, data.table = FALSE)
    if (ncol(dt) < 2) stop("文件列数不足: ", basename(f))
    num_idx <- which(vapply(dt, is.numeric, logical(1)))
    if (length(num_idx) == 0) stop("文件无数值列: ", basename(f))
    if (is.null(count_col)) {
      cnt_col <- colnames(dt)[num_idx[length(num_idx)]]
    } else if (is.numeric(count_col)) {
      cnt_col <- colnames(dt)[count_col]
    } else {
      if (!count_col %in% colnames(dt)) stop("缺少计数列 ", count_col, ": ", basename(f))
      cnt_col <- count_col
    }
    feature_col <- colnames(dt)[id_col]
    sample_name <- sample_fun(f)
    out <- dt[, c(feature_col, cnt_col), drop = FALSE]
    colnames(out) <- c("feature_id", sample_name)
    out
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  finalize_feature_matrix(merged, tag = tag)
}

read_salmonte_expr <- function(files) {
  mats <- lapply(files, function(f) {
    dt <- fread(f, data.table = FALSE)
    if (ncol(dt) < 2) stop("SalmonTE EXPR.csv 列数不足: ", basename(f))
    feature_col <- colnames(dt)[1]
    num_idx <- which(vapply(dt, is.numeric, logical(1)))
    if (length(num_idx) == 0) stop("SalmonTE EXPR.csv 无数值列: ", basename(f))
    cnt_col <- colnames(dt)[num_idx[1]]
    sample_name <- sample_from_parent(f, "_SalmonTE_output$")
    out <- dt[, c(feature_col, cnt_col), drop = FALSE]
    colnames(out) <- c("feature_id", sample_name)
    out
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  finalize_feature_matrix(merged, tag = "SalmonTE")
}

ensure_sample_order <- function(mat, sample_table) {
  keep <- intersect(sample_table$sample, colnames(mat))
  if (length(keep) < 2) stop("矩阵中可匹配样本少于2个")
  st <- sample_table[match(keep, sample_table$sample), , drop = FALSE]
  mat <- mat[, keep, drop = FALSE]
  list(mat = mat, sample_table = st)
}

run_exploratory_simple <- function(count_mat, sample_table, contrast) {
  exploratory_method <- get0("exploratory_method", envir = .GlobalEnv, ifnotfound = "logCPM_diff")
  exploratory_method <- as.character(exploratory_method)[1]
  if (!exploratory_method %in% c("logCPM_diff", "edgeR_fixedBCV")) {
    stop("[run_exploratory_simple] unsupported exploratory_method: ", exploratory_method)
  }
  fixed_bcv <- suppressWarnings(as.numeric(get0("exploratory_fixed_bcv", envir = .GlobalEnv, ifnotfound = 0.4)))
  if (is.na(fixed_bcv) || fixed_bcv <= 0) fixed_bcv <- 0.4

  # 只保留当前 contrast 的样本
  sample_table <- sample_table[sample_table[[contrast[1]]] %in% contrast[2:3], , drop = FALSE]

  x <- ensure_sample_order(count_mat, sample_table)
  count_mat <- as.matrix(x$mat)
  storage.mode(count_mat) <- "numeric"
  count_mat[is.na(count_mat)] <- 0
  sample_table <- x$sample_table

  sample_table[[contrast[1]]] <- factor(sample_table[[contrast[1]]], levels = c(contrast[3], contrast[2]))
  sample_table[[contrast[1]]] <- droplevels(sample_table[[contrast[1]]])

  # 这里 sample_table 的 rownames 可能还没设，保险起见：
  rownames(sample_table) <- sample_table$sample
  g1 <- sample_table$sample[sample_table[[contrast[1]]] == contrast[2]]
  g0 <- sample_table$sample[sample_table[[contrast[1]]] == contrast[3]]

  if (length(g1) == 0 || length(g0) == 0) {
    stop(
      "[run_exploratory_simple] contrast 在表达矩阵中匹配失败: ",
      contrast[2], " matched=", length(g1),
      "; ", contrast[3], " matched=", length(g0)
    )
  }
  message(
    "[run_exploratory_simple] matched samples: ",
    contrast[2], "=", length(g1), "; ",
    contrast[3], "=", length(g0)
  )
  message(
    "[run_exploratory_simple] exploratory_method=", exploratory_method,
    if (identical(exploratory_method, "edgeR_fixedBCV")) paste0("; fixed_BCV=", fixed_bcv) else ""
  )

  if (identical(exploratory_method, "edgeR_fixedBCV")) {
    dge <- edgeR::DGEList(counts = round(count_mat))
    dge <- edgeR::calcNormFactors(dge, method = "TMM")
    eff_lib <- dge$samples$lib.size * dge$samples$norm.factors
    eff_lib[!is.finite(eff_lib) | eff_lib <= 0] <- 1
    norm_mat <- sweep(dge$counts, 2, eff_lib / mean(eff_lib), "/")
    log_mat <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

    group_samples <- c(g0, g1)
    dge_sub <- dge[, group_samples, keep.lib.sizes = FALSE]
    dge_sub$samples$group <- factor(
      c(rep(contrast[3], length(g0)), rep(contrast[2], length(g1))),
      levels = c(contrast[3], contrast[2])
    )
    dge_sub <- edgeR::calcNormFactors(dge_sub, method = "TMM")
    et <- edgeR::exactTest(dge_sub, pair = c(contrast[3], contrast[2]), dispersion = fixed_bcv^2)
    tt <- edgeR::topTags(et, n = Inf, sort.by = "none")$table
    tt <- tt[rownames(count_mat), , drop = FALSE]

    baseMean <- rowMeans(norm_mat[, group_samples, drop = FALSE], na.rm = TRUE)
    log2FoldChange <- as.numeric(tt$logFC)
    pvalue <- as.numeric(tt$PValue)
    padj <- p.adjust(pvalue, method = "fdr")
  } else {
    # 简单 library size 标准化
    libsize <- colSums(count_mat, na.rm = TRUE)
    libsize[libsize == 0] <- 1
    norm_mat <- t(t(count_mat) / libsize * 1e6)   # CPM
    log_mat <- log2(norm_mat + 1)

    mean1 <- if (length(g1) == 1) log_mat[, g1] else rowMeans(log_mat[, g1, drop = FALSE], na.rm = TRUE)
    mean0 <- if (length(g0) == 1) log_mat[, g0] else rowMeans(log_mat[, g0, drop = FALSE], na.rm = TRUE)

    baseMean <- rowMeans(norm_mat[, c(g0, g1), drop = FALSE], na.rm = TRUE)
    log2FoldChange <- mean1 - mean0
    pvalue <- NA_real_
    padj <- NA_real_
  }

  res <- data.frame(
    baseMean = baseMean,
    log2FoldChange = log2FoldChange,
    lfcSE = NA_real_,
    stat = NA_real_,
    pvalue = pvalue,
    padj = padj,
    feature_id = rownames(count_mat),
    stringsAsFactors = FALSE
  )
  res$exploratory_method <- exploratory_method
  if (identical(exploratory_method, "edgeR_fixedBCV")) res$fixed_BCV <- fixed_bcv

  # exploratory 模式下定义候选差异
  res$exploratory_sig <- ifelse(
    res$baseMean >= 5 & abs(res$log2FoldChange) >= lfc_cutoff,
    "sig", "ns"
  )

  res <- res %>% arrange(desc(abs(log2FoldChange)), desc(baseMean))

  list(
    dds = NULL,
    res = res,
    vst_mat = log_mat,
    sample_table = sample_table,
    mode = "exploratory"
  )
}

run_deseq_simple <- function(count_mat, sample_table, contrast) {
  # 只保留当前 contrast 的样本
  sample_table2 <- sample_table[sample_table[[contrast[1]]] %in% contrast[2:3], , drop = FALSE]

  # 统计每组样本数
  tab <- table(sample_table2[[contrast[1]]])
  n_case <- ifelse(contrast[2] %in% names(tab), unname(tab[contrast[2]]), 0)
  n_ctrl <- ifelse(contrast[3] %in% names(tab), unname(tab[contrast[3]]), 0)

  # 只要任一组 < 2，就自动切换 exploratory 模式
  if (n_case < 2 || n_ctrl < 2) {
    msg_method <- get0("exploratory_method", envir = .GlobalEnv, ifnotfound = "logCPM_diff")
    msg_tail <- if (identical(as.character(msg_method)[1], "edgeR_fixedBCV")) {
      "；使用 edgeR fixedBCV 输出探索性 pvalue/padj"
    } else {
      "；不计算p值/padj"
    }
    message("[", paste(contrast, collapse = " | "), "] replicate不足，自动切换为 exploratory 模式", msg_tail)
    return(run_exploratory_simple(count_mat, sample_table, contrast))
  }

  x <- ensure_sample_order(count_mat, sample_table2)
  count_mat <- round(as.matrix(x$mat))
  sample_table2 <- x$sample_table
  message(
    "[run_deseq_simple] contrast=", paste(contrast, collapse = " | "),
    "; matched matrix samples=", paste(colnames(count_mat), collapse = ",")
  )

  sample_table2[[contrast[1]]] <- factor(sample_table2[[contrast[1]]], levels = c(contrast[3], contrast[2]))
  sample_table2[[contrast[1]]] <- droplevels(sample_table2[[contrast[1]]])

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = sample_table2,
    design = as.formula(paste0("~ ", contrast[1]))
  )

  keep <- rowSums(counts(dds) >= 5) >= 2
  dds <- dds[keep, ]
  if (nrow(dds) < 2) stop("过滤后特征数太少，无法做 DESeq2")

  dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = contrast)
    res <- as.data.frame(res)
    res$feature_id <- rownames(res)
    res <- res %>% arrange(padj, desc(abs(log2FoldChange)))

    # --- 修复代码开始 ---
    vst_mat <- tryCatch({
      if (nrow(dds) < 1000) {
        assay(varianceStabilizingTransformation(dds, blind = TRUE))
      } else {
        assay(vst(dds, blind = TRUE))
      }
    }, error = function(e) {
      message("[run_deseq_simple] vst failed, fallback to log2: ", conditionMessage(e))
      log2(counts(dds, normalized = TRUE) + 1)
    })
    # --- 修复代码结束 ---

    list(dds = dds, res = res, vst_mat = vst_mat, sample_table = sample_table2, mode = "deseq2")
}

run_deseq_native_full <- function(count_mat, sample_table, contrast) {
  if (!(contrast[1] %in% colnames(sample_table))) {
    stop("sample_table 缺少 contrast 列: ", contrast[1])
  }

  x <- ensure_sample_order(count_mat, sample_table)
  count_mat <- round(as.matrix(x$mat))
  sample_table2 <- x$sample_table

  col_data <- data.frame(
    sample = sample_table2$sample,
    condition = as.character(sample_table2[[contrast[1]]]),
    stringsAsFactors = FALSE
  )
  rownames(col_data) <- col_data$sample
  col_data <- col_data[, "condition", drop = FALSE]
  col_data$condition <- factor(col_data$condition)

  if (!all(c(contrast[2], contrast[3]) %in% levels(col_data$condition))) {
    stop("contrast 水平不存在于 sample_table 中: ", paste(contrast[2:3], collapse = ", "))
  }

  dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = col_data,
    design = ~ condition
  )

  keep <- rowSums(counts(dds)) > 10
  message("[run_deseq_native_full] keep rows by rowSums>10: ", sum(keep), "/", length(keep))
  dds <- dds[keep, ]
  if (nrow(dds) < 2) stop("过滤后特征数太少，无法做 native REdiscoverTE DESeq2")

  dds <- tryCatch(
    DESeq(dds, quiet = TRUE),
    error = function(e) {
      msg <- conditionMessage(e)
      message("[run_deseq_native_full] DESeq fallback: ", msg)

      if (grepl("invalid 'x'|newsplit: out of vertex space", msg)) {
        return(DESeq(dds, fitType = "mean", quiet = TRUE))
      }
      if (grepl("all gene-wise dispersion estimates", msg, ignore.case = TRUE)) {
        dds <- estimateSizeFactors(dds)
        dds <- estimateDispersionsGeneEst(dds)
        dispersions(dds) <- mcols(dds)$dispGeneEst
        dds <- nbinomWaldTest(dds)
        return(dds)
      }
      stop(e)
    }
  )

  res <- as.data.frame(results(dds, contrast = c("condition", contrast[2], contrast[3])))
  res$feature_id <- rownames(res)
  res <- res %>% arrange(padj, desc(abs(log2FoldChange)))

  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_df <- data.frame(feature_id = rownames(norm_counts), norm_counts, check.names = FALSE, stringsAsFactors = FALSE)

  vst_mat <- tryCatch(
    assay(vst(dds, blind = TRUE)),
    error = function(e) {
      message("[run_deseq_native_full] vst failed, fallback to log2 normalized counts: ", conditionMessage(e))
      log2(as.matrix(norm_counts) + 1)
    }
  )

  list(
    dds = dds,
    res = res,
    norm_counts = norm_counts_df,
    vst_mat = vst_mat,
    sample_table = sample_table2,
    mode = "deseq2_native"
  )
}

clean_res_for_export <- function(res_df, tag = "DE_export") {
  bad <- is.na(res_df$baseMean) |
    (is.na(res_df$log2FoldChange) &
     is.na(res_df$lfcSE) &
     is.na(res_df$stat) &
     is.na(res_df$pvalue) &
     is.na(res_df$padj))

  if (any(bad, na.rm = TRUE)) {
    message("[", tag, "] drop all-NA DE rows: ", sum(bad, na.rm = TRUE))
  }

  res_df[!bad, , drop = FALSE]
}

run_and_plot_deseq_native <- function(count_mat, sample_table, contrast_list, out_prefix, gene_anno = NULL, te_label_level = NULL, te_color_level = NULL) {
  tool_outdir <- file.path(outdir, out_prefix)
  dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)
  figure_dir <- file.path(tool_outdir, "figures")
  de_dir <- file.path(tool_outdir, "de_matrices")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

  # 先画一次全局 PCA / Pearson（每个工具只画一次）
  safe_run(paste0(out_prefix, "_global_qc"), {
    x_all <- ensure_sample_order(count_mat, sample_table)
    count_mat_all <- round(as.matrix(x_all$mat))
    sample_table_all <- x_all$sample_table

    if (ncol(count_mat_all) >= 2) {
      vst_mat_all <- NULL

      # native 这里尽量复用 DESeq2 风格的 vst；失败时回退到 log2(count+1)
      col_data_all <- data.frame(
        sample = sample_table_all$sample,
        condition = as.character(sample_table_all$condition),
        stringsAsFactors = FALSE
      )
      rownames(col_data_all) <- col_data_all$sample
      col_data_all <- col_data_all[, "condition", drop = FALSE]
      col_data_all$condition <- factor(col_data_all$condition)

      dds_all <- DESeqDataSetFromMatrix(
        countData = count_mat_all,
        colData = col_data_all,
        design = ~ condition
      )

      keep_all <- rowSums(counts(dds_all)) > 10
      dds_all <- dds_all[keep_all, ]

      if (nrow(dds_all) >= 2) {
        vst_mat_all <- tryCatch(
          assay(vst(dds_all, blind = TRUE)),
          error = function(e) {
            message("[", out_prefix, "_global_qc] vst fallback: ", conditionMessage(e))
            log2(counts(dds_all, normalized = FALSE) + 1)
          }
        )

        safe_plot <- function(plot_tag, expr_fun, outfile) {
          tryCatch(
            expr_fun(),
            error = function(e) {
              message("[", plot_tag, "] ERROR: ", conditionMessage(e))
              writeLines(paste0("[", plot_tag, "] ERROR: ", conditionMessage(e)), paste0(outfile, ".txt"))
              save_placeholder_plot(outfile, plot_tag, paste("plot failed:", conditionMessage(e)))
              NULL
            }
          )
        }

        safe_plot(
          paste0(out_prefix, "_PCA_global"),
          function() plot_pca_simple(
            vst_mat_all,
            sample_table_all,
            file.path(figure_dir, paste0(out_prefix, "_PCA_global.pdf")),
            paste0(out_prefix, "_global")
          ),
          file.path(figure_dir, paste0(out_prefix, "_PCA_global.pdf"))
        )

        safe_plot(
          paste0(out_prefix, "_Pearson_global"),
          function() plot_corr_simple(
            vst_mat_all,
            file.path(figure_dir, paste0(out_prefix, "_Pearson_global.pdf")),
            paste0(out_prefix, "_global")
          ),
          file.path(figure_dir, paste0(out_prefix, "_Pearson_global.pdf"))
        )
      } else {
        message("[", out_prefix, "_global_qc] too few rows after filtering, skip global PCA/Pearson")
      }
    } else {
      message("[", out_prefix, "_global_qc] too few samples, skip global PCA/Pearson")
    }
  })

  worker_fun <- function(ct) {
    tag <- sanitize_name(paste(out_prefix, ct[2], "vs", ct[3], sep = "_"))
    safe_run(tag, {
      fit <- run_deseq_native_full(count_mat, sample_table, ct)
      res <- fit$res

      this_label_col <- "feature_id"
      this_color_col <- NULL
      this_heatmap_row_anno_col <- NULL

      if (!is.null(gene_anno) && any(c("repName","repFamily","repClass","locus_id") %in% colnames(gene_anno))) {
        pref_cols <- unique(na.omit(c(
          te_label_level,
          te_color_level,
          "feature_id", "repName", "repFamily", "repClass", "locus_id", "gene_name"
        )))
        res <- annotate_te_res(
          res, gene_anno,
          preferred_match_cols = pref_cols,
          label_level = te_label_level,
          color_level = te_color_level
        )
        if ("te_label_plot" %in% colnames(res)) this_label_col <- "te_label_plot"
        if ("te_color_plot" %in% colnames(res)) this_color_col <- "te_color_plot"
        if ("te_heatmap_group" %in% colnames(res)) this_heatmap_row_anno_col <- "te_heatmap_group"
      }

      res_export <- clean_res_for_export(res, tag = tag)
      write.csv(res_export, file.path(de_dir, paste0(tag, "_DE.csv")), row.names = FALSE)
      write.csv(fit$norm_counts, file.path(de_dir, paste0(tag, "_normalized_counts.csv")), row.names = FALSE)

      safe_plot <- function(plot_tag, expr_fun, outfile) {
        tryCatch(
          expr_fun(),
          error = function(e) {
            message("[", plot_tag, "] ERROR: ", conditionMessage(e))
            writeLines(paste0("[", plot_tag, "] ERROR: ", conditionMessage(e)), paste0(outfile, ".txt"))
            save_placeholder_plot(outfile, plot_tag, paste("plot failed:", conditionMessage(e)))
            NULL
          }
        )
      }

      # 去掉每个 contrast 重复画 PCA / Pearson
      volcano_metric <- plot_metric_for_res(res_export, "padj")
      volcano_pvalue_metric <- plot_metric_for_res(res_export, "pvalue")
      plot_tasks <- list(
        list(
          tag = paste0(tag, "_volcano"),
          outfile = file.path(figure_dir, paste0(tag, "_volcano.pdf")),
          fun = function(outfile) plot_volcano_simple(res_export, outfile, tag, label_col = this_label_col, top_n = label_top_n, color_col = this_color_col, sig_metric = volcano_metric)
        ),
        list(
          tag = paste0(tag, "_volcano_pvalue"),
          outfile = file.path(figure_dir, paste0(tag, "_volcano_pvalue.pdf")),
          fun = function(outfile) plot_volcano_simple(res_export, outfile, tag, label_col = this_label_col, top_n = label_top_n, color_col = this_color_col, sig_metric = volcano_pvalue_metric)
        ),
        list(
          tag = paste0(tag, "_MA"),
          outfile = file.path(figure_dir, paste0(tag, "_MA.pdf")),
          fun = function(outfile) plot_ma_simple(res_export, outfile, tag, label_col = this_label_col, top_n = label_top_n, color_col = this_color_col, sig_metric = volcano_metric)
        ),
        list(
          tag = paste0(tag, "_MA_pvalue"),
          outfile = file.path(figure_dir, paste0(tag, "_MA_pvalue.pdf")),
          fun = function(outfile) plot_ma_simple(res_export, outfile, tag, label_col = this_label_col, top_n = label_top_n, color_col = this_color_col, sig_metric = volcano_pvalue_metric)
        ),
        list(
          tag = paste0(tag, "_heatmap"),
          outfile = file.path(figure_dir, paste0(tag, "_heatmap.pdf")),
          fun = function(outfile) plot_top_heatmap_simple(fit$vst_mat, res_export, fit$sample_table, outfile, tag, top_n = heatmap_top_n, label_col = this_label_col, row_anno_col = this_heatmap_row_anno_col)
        )
      )
      run_one_plot <- function(task) safe_plot(task$tag, function() task$fun(task$outfile), task$outfile)
      if (length(contrast_list) == 1 && .Platform$OS.type == "unix" && plot_n_cores > 1) {
        parallel::mclapply(plot_tasks, run_one_plot, mc.cores = min(plot_n_cores, length(plot_tasks)), mc.preschedule = FALSE)
      } else {
        lapply(plot_tasks, run_one_plot)
      }
    })
    NULL
  }

  if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && plot_n_cores > 1) {
    parallel::mclapply(
      contrast_list,
      worker_fun,
      mc.cores = min(plot_n_cores, length(contrast_list)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(contrast_list, worker_fun)
  }
}


plot_pca_simple <- function(vst_mat, sample_table, outfile, title) {
  outfile <- as_pdf_path(outfile)
  pca <- prcomp(t(vst_mat), scale. = FALSE)
  df <- data.frame(
    sample = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    stringsAsFactors = FALSE
  ) %>% left_join(sample_table, by = "sample")

  p <- ggplot(df, aes(PC1, PC2, color = .data[[names(sample_table)[2]]], label = sample)) +
    geom_point(size = 3, alpha = 0.88) +
    geom_text_repel(size = 3.6, max.overlaps = 20, segment.color = "#A6A6A6", segment.size = 0.25) +
    scale_color_manual(values = journal_palette(unique(df[[names(sample_table)[2]]]))) +
    journal_theme(base_size = 13) +
    labs(title = title)
  ggsave(outfile, p, width = 6.5, height = 5.5)
  invisible(p)
}

save_placeholder_plot <- function(outfile, title, msg = "not enough data") {
  outfile <- as_pdf_path(outfile)
  p <- ggplot() +
    annotate("text", x = 1, y = 1, label = msg, size = 6) +
    theme_void(base_size = 13) +
    labs(title = title)
  ggsave(outfile, p, width = 7, height = 5.5)
  invisible(NULL)
}

plot_corr_simple <- function(vst_mat, outfile, title) {
  outfile <- as_pdf_path(outfile)
  if (is.null(dim(vst_mat)) || ncol(vst_mat) < 2) {
    save_placeholder_plot(outfile, title, "too few samples for correlation")
    return(invisible(NULL))
  }

  cor_mat <- suppressWarnings(cor(vst_mat, method = "pearson", use = "pairwise.complete.obs"))
  if (is.null(dim(cor_mat))) {
    cor_mat <- matrix(1, nrow = ncol(vst_mat), ncol = ncol(vst_mat))
    colnames(cor_mat) <- colnames(vst_mat)
    rownames(cor_mat) <- colnames(vst_mat)
  }
  cor_mat[is.na(cor_mat)] <- 1

  pheatmap(
    cor_mat,
    main = title,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    filename = outfile,
    width = 7,
    height = 6
  )
}
collapse_plot_groups <- function(x, max_groups = 12) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"

  tb <- sort(table(x), decreasing = TRUE)
  if (length(tb) > max_groups) {
    keep <- names(tb)[seq_len(max_groups - 1)]
    x[!x %in% keep] <- "Other"
  }
  x
}

plot_bool <- function(x, default = TRUE) {
  if (missing(x) || is.null(x) || length(x) == 0 || is.na(x)) return(default)
  y <- tolower(trimws(as.character(x)[1]))
  if (y %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (y %in% c("false", "f", "no", "n", "0")) return(FALSE)
  default
}

get_plot_option <- function(name, default) {
  val <- get0(name, envir = .GlobalEnv, ifnotfound = default)
  if (is.null(val) || length(val) == 0 || is.na(val[1])) default else val[1]
}

journal_theme <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "plain", size = base_size + 1, hjust = 0),
      axis.title = element_text(color = "#1F2933"),
      axis.text = element_text(color = "#1F2933"),
      axis.line = element_line(color = "#202020", linewidth = 0.55),
      axis.ticks = element_line(color = "#202020", linewidth = 0.45),
      legend.position = "right",
      legend.title = element_text(size = base_size - 3, color = "#222222"),
      legend.text = element_text(size = base_size - 3, color = "#222222"),
      legend.key.height = grid::unit(0.34, "cm"),
      legend.key.width  = grid::unit(0.34, "cm"),
      plot.margin = margin(8, 16, 8, 8)
    )
}

journal_palette <- function(groups) {
  groups <- as.character(groups)
  fallback <- c(
    "#4E79A7", "#E15759", "#59A14F", "#F28E2B", "#B07AA1", "#76B7B2",
    "#EDC948", "#9C755F", "#FF9DA7", "#8CD17D", "#499894", "#D37295",
    "#B6992D", "#86BCB6", "#FABFD2", "#A0CBE8", "#FFBE7D", "#C85200",
    "#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08",
    "#E17C05", "#CC503E", "#94346E", "#6F4070", "#994E95", "#666666",
    "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"
  )
  pal <- setNames(rep(fallback, length.out = length(groups)), groups)
  fixed <- c(
    NotSig = "#C9CCD3", ns = "#C9CCD3", `No change` = "#C9CCD3",
    Unknown = "#BFC3CB", Other = "#C9CCD3",
    Up = "#C95F50", Down = "#4F789F",
    LINE = "#D55E00", SINE = "#009E73", LTR = "#0072B2", DNA = "#CC79A7",
    Satellite = "#8CD17D", Simple_repeat = "#9C755F", Low_complexity = "#B8A6A0",
    `ERV class I` = "#4E79A7", `ERV class II` = "#F28E2B", `ERV class III` = "#59A14F",
    ERV1 = "#4E79A7", ERVK = "#5F4690", ERVL = "#1D6996", `ERVL-MaLR` = "#38A6A5",
    L1 = "#D55E00", L2 = "#E17C05", CR1 = "#CC503E", RTE = "#94346E",
    Alu = "#009E73", MIR = "#73AF48", B2 = "#0F8554", B4 = "#8CD17D",
    hAT = "#B07AA1", TcMar = "#994E95", PiggyBac = "#6F4070", tRNA = "#B6992D"
  )
  for (nm in intersect(names(fixed), groups)) pal[nm] <- fixed[nm]
  pal[grepl("^LINE($|/|:)", groups, ignore.case = TRUE)] <- "#D55E00"
  pal[grepl("^SINE($|/|:)", groups, ignore.case = TRUE)] <- "#009E73"
  pal[grepl("^LTR($|/|:)", groups, ignore.case = TRUE)] <- "#0072B2"
  pal[grepl("^DNA($|/|:)", groups, ignore.case = TRUE)] <- "#CC79A7"
  pal[grepl("satellite", groups, ignore.case = TRUE)] <- "#A3A86B"
  pal[grepl("simple|low", groups, ignore.case = TRUE)] <- "#B8A6A0"
  pal
}

add_plot_groups <- function(df, color_col = NULL, gray_nonsig = TRUE, max_groups = 24) {
  df$fc_direction <- dplyr::case_when(
    is.na(df$log2FoldChange) ~ "No change",
    df$log2FoldChange >= lfc_cutoff ~ "Up",
    df$log2FoldChange <= -lfc_cutoff ~ "Down",
    TRUE ~ "No change"
  )

  if (!is.null(color_col) && color_col %in% colnames(df)) {
    df$plot_group <- collapse_plot_groups(df[[color_col]], max_groups = max_groups)
    n_other <- sum(df$plot_group == "Other", na.rm = TRUE)
    if (n_other > 0) message("[plot] collapse small groups to Other: ", n_other)
    if (gray_nonsig) df$plot_group[df$sig == "ns"] <- "NotSig"
  } else {
    df$plot_group <- ifelse(df$sig == "sig", df$fc_direction, if (gray_nonsig) "NotSig" else df$fc_direction)
  }

  df$plot_alpha <- ifelse(df$sig == "sig", 0.82, 0.22)
  df$plot_size <- ifelse(df$sig == "sig", 1.75, 1.25)
  df <- df %>% dplyr::arrange(.data$sig == "sig")
  df$plot_group <- factor(df$plot_group, levels = unique(df$plot_group))
  df
}

select_extreme_labels <- function(df, top_n = label_top_n, id_col = "feature_id", y_col = "neglog10y") {
  top_n <- max(0L, as.integer(top_n))
  if (top_n == 0 || nrow(df) == 0) return(df[0, , drop = FALSE])

  df$.label_id <- if (id_col %in% colnames(df)) as.character(df[[id_col]]) else as.character(seq_len(nrow(df)))
  df$.label_y_score <- if (y_col %in% colnames(df)) {
    suppressWarnings(as.numeric(df[[y_col]]))
  } else if ("neglog10y" %in% colnames(df)) {
    suppressWarnings(as.numeric(df$neglog10y))
  } else if ("baseMean" %in% colnames(df)) {
    log10(pmax(suppressWarnings(as.numeric(df$baseMean)), 0) + 1)
  } else {
    rep(0, nrow(df))
  }
  df$.label_base_mean <- if ("baseMean" %in% colnames(df)) {
    suppressWarnings(as.numeric(df$baseMean))
  } else {
    rep(0, nrow(df))
  }

  candidates <- df %>%
    dplyr::filter(.data$sig == "sig") %>%
    dplyr::filter(!is.na(.data$log2FoldChange), !is.na(.data[[".label_y_score"]]))
  if (nrow(candidates) == 0) {
    candidates <- df %>%
      dplyr::filter(!is.na(.data$log2FoldChange), !is.na(.data[[".label_y_score"]]))
  }
  if (nrow(candidates) == 0) return(df[0, , drop = FALSE])

  n_y <- ceiling(top_n * 0.45)
  n_abs <- top_n - n_y
  lab_y <- candidates %>%
    dplyr::arrange(dplyr::desc(.data[[".label_y_score"]]), dplyr::desc(abs(.data$log2FoldChange)), dplyr::desc(.data[[".label_base_mean"]])) %>%
    utils::head(n_y)
  lab_abs <- candidates %>%
    dplyr::arrange(dplyr::desc(abs(.data$log2FoldChange)), dplyr::desc(.data[[".label_y_score"]]), dplyr::desc(.data[[".label_base_mean"]])) %>%
    utils::head(n_abs)

  dplyr::bind_rows(lab_y, lab_abs) %>%
    dplyr::distinct(.data[[".label_id"]], .keep_all = TRUE) %>%
    dplyr::arrange(dplyr::desc(abs(.data$log2FoldChange)), dplyr::desc(.data[[".label_y_score"]])) %>%
    utils::head(top_n)
}

is_exploratory_res <- function(res_df) {
  "exploratory_method" %in% colnames(res_df) && any(!is.na(res_df$exploratory_method) & res_df$exploratory_method != "")
}

plot_metric_for_res <- function(res_df, preferred = "padj") {
  if (is_exploratory_res(res_df)) "baseMean" else preferred
}

get_sig_logic <- function(df, metric = c("padj", "pvalue", "baseMean", "auto"), plot_type = c("volcano", "ma")) {
  metric <- match.arg(metric)
  plot_type <- match.arg(plot_type)

  has_padj <- "padj" %in% colnames(df) && any(!is.na(df$padj))
  has_pvalue <- "pvalue" %in% colnames(df) && any(!is.na(df$pvalue))
  df$baseMean <- suppressWarnings(as.numeric(df$baseMean))
  df$log2FoldChange <- suppressWarnings(as.numeric(df$log2FoldChange))

  if (metric == "auto") {
    metric <- if (is_exploratory_res(df)) "baseMean" else if (has_padj) "padj" else if (has_pvalue) "pvalue" else "baseMean"
  }

  if (metric == "padj" && has_padj) {
    df$sig <- ifelse(!is.na(df$padj) & df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff, "sig", "ns")
    cutoff_line <- -log10(padj_cutoff)
    ylab <- "-log10(padj)"
    if (plot_type == "volcano") df$neglog10y <- -log10(pmax(df$padj, 1e-300))
    title_suffix <- ""
  } else if (metric == "pvalue" && has_pvalue) {
    df$sig <- ifelse(!is.na(df$pvalue) & df$pvalue < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff, "sig", "ns")
    cutoff_line <- -log10(padj_cutoff)
    ylab <- "-log10(pvalue)"
    if (plot_type == "volcano") df$neglog10y <- -log10(pmax(df$pvalue, 1e-300))
    title_suffix <- " (pvalue)"
  } else {
    df$sig <- ifelse(df$baseMean >= baseMean_min & abs(df$log2FoldChange) >= lfc_cutoff, "sig", "ns")
    cutoff_line <- if (plot_type == "volcano") log10(baseMean_min + 1) else NA_real_
    ylab <- if (plot_type == "volcano") "log10(baseMean + 1)" else "log10(baseMean)"
    if (plot_type == "volcano") df$neglog10y <- log10(pmax(df$baseMean + 1, 1))
    title_suffix <- " (exploratory)"
  }

  list(df = df, cutoff_line = cutoff_line, ylab = ylab, title_suffix = title_suffix)
}

plot_volcano_simple <- function(res_df, outfile, title, label_col = "feature_id", top_n = label_top_n, color_col = NULL, sig_metric = "auto", orientation = get_plot_option("volcano_orientation", "classic"), gray_nonsig = plot_bool(get_plot_option("gray_nonsig", TRUE), TRUE)) {
  outfile <- as_pdf_path(outfile)
  orientation <- match.arg(tolower(as.character(orientation)), c("classic", "horizontal"))
  res_df <- as.data.frame(res_df)
  sig_metric <- plot_metric_for_res(res_df, sig_metric)
  sig_info <- get_sig_logic(res_df, metric = sig_metric, plot_type = "volcano")
  df <- sig_info$df
  ylab_sig <- sig_info$ylab
  title <- paste0(title, sig_info$title_suffix)
  df <- add_plot_groups(df, color_col = color_col, gray_nonsig = gray_nonsig, max_groups = 24)

  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, title, "no features for volcano")
    return(invisible(NULL))
  }

  lab_df <- select_extreme_labels(df, top_n = top_n, id_col = label_col, y_col = "neglog10y")

  plot_groups <- unique(df$plot_group)
  pal <- journal_palette(plot_groups)

  if (orientation == "horizontal") {
    df$x_plot <- df$neglog10y
    df$y_plot <- df$log2FoldChange
    lab_df$x_plot <- lab_df$neglog10y
    lab_df$y_plot <- lab_df$log2FoldChange
    x_lab <- ylab_sig
    y_lab <- expression(log[2](fold~change))
    width <- 9.8
    height <- 6.4
  } else {
    df$x_plot <- df$log2FoldChange
    df$y_plot <- df$neglog10y
    lab_df$x_plot <- lab_df$log2FoldChange
    lab_df$y_plot <- lab_df$neglog10y
    x_lab <- expression(log[2](fold~change))
    y_lab <- ylab_sig
    width <- 7.8
    height <- 6.6
  }

  p <- ggplot(df, aes(x = x_plot, y = y_plot, color = plot_group, alpha = sig, size = sig)) +
    geom_point() +
    scale_color_manual(values = pal) +
    scale_alpha_manual(values = c(ns = 0.18, sig = 0.78), guide = "none") +
    scale_size_manual(values = c(ns = 0.85, sig = 1.35), guide = "none") +
    journal_theme(base_size = 13) +
    guides(color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(size = 2.6, alpha = 1)
    )) +
    labs(title = title, x = x_lab, y = y_lab, color = ifelse(is.null(color_col), "Direction", color_col))

  if (orientation == "horizontal") {
    p <- p +
      geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, linewidth = 0.45, color = "#8A8A8A") +
      geom_vline(xintercept = 0, linewidth = 0.35, color = "#202020")
    if (!is.na(sig_info$cutoff_line)) {
      p <- p + geom_vline(xintercept = sig_info$cutoff_line, linetype = 2, linewidth = 0.45, color = "#8A8A8A")
    }
  } else {
    p <- p +
      geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, linewidth = 0.45, color = "#8A8A8A") +
      geom_hline(yintercept = 0, linewidth = 0.35, color = "#202020")
    if (!is.na(sig_info$cutoff_line)) {
      p <- p + geom_hline(yintercept = sig_info$cutoff_line, linetype = 2, linewidth = 0.45, color = "#8A8A8A")
    }
  }

  if (nrow(lab_df) > 0 && label_col %in% colnames(df)) {
    p <- p + geom_text_repel(
      data = lab_df,
      aes(x = x_plot, y = y_plot, label = .data[[label_col]], color = plot_group),
      size = 2.7,
      max.overlaps = 50,
      show.legend = FALSE,
      segment.size = 0.22,
      segment.color = "#A6A6A6",
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15
    )
  }

  ggsave(outfile, p, width = width, height = height, limitsize = FALSE)
  invisible(p)
}


plot_ma_simple <- function(res_df, outfile, title, label_col = "feature_id", top_n = label_top_n, color_col = NULL, sig_metric = "auto", gray_nonsig = plot_bool(get_plot_option("gray_nonsig", TRUE), TRUE)) {
  outfile <- as_pdf_path(outfile)
  df <- as.data.frame(res_df)
  sig_metric <- plot_metric_for_res(df, sig_metric)
  df$baseMean <- pmax(df$baseMean, 1e-8)

  sig_info <- get_sig_logic(df, metric = sig_metric, plot_type = "ma")
  df <- sig_info$df
  title <- paste0(title, sig_info$title_suffix)

  df <- add_plot_groups(df, color_col = color_col, gray_nonsig = gray_nonsig, max_groups = 24)

  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, title, "no features for MA")
    return(invisible(NULL))
  }

  lab_df <- select_extreme_labels(df, top_n = top_n, id_col = label_col, y_col = "baseMean")

  plot_groups <- unique(df$plot_group)
  pal <- journal_palette(plot_groups)

  p <- ggplot(df, aes(x = log10(baseMean + 1), y = log2FoldChange, color = plot_group, alpha = sig, size = sig)) +
    geom_point() +
    scale_color_manual(values = pal) +
    scale_alpha_manual(values = c(ns = 0.18, sig = 0.78), guide = "none") +
    scale_size_manual(values = c(ns = 0.85, sig = 1.35), guide = "none") +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "#202020") +
    geom_hline(yintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, linewidth = 0.45, color = "#8A8A8A") +
    journal_theme(base_size = 13) +
    guides(color = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(size = 2.6, alpha = 1)
    )) +
    labs(title = title, x = "log10(baseMean + 1)", y = expression(log[2](fold~change)), color = ifelse(is.null(color_col), "Direction", color_col))

  if (nrow(lab_df) > 0 && label_col %in% colnames(df)) {
    p <- p + geom_text_repel(
      data = lab_df,
      aes(label = .data[[label_col]], color = plot_group),
      size = 2.7,
      max.overlaps = 50,
      show.legend = FALSE,
      segment.size = 0.22,
      segment.color = "#A6A6A6",
      min.segment.length = 0,
      box.padding = 0.25,
      point.padding = 0.15
    )
  }

  ggsave(outfile, p, width = 7.8, height = 6.6, limitsize = FALSE)
  invisible(p)
}

plot_top_heatmap_simple <- function(vst_mat, res_df, sample_table, outfile, title, top_n = heatmap_top_n, label_col = "feature_id", row_anno_col = NULL) {
  outfile <- as_pdf_path(outfile)
  has_padj <- "padj" %in% colnames(res_df) && any(!is.na(res_df$padj))

  if (has_padj) {
    sel_df <- res_df %>%
      filter(!is.na(padj)) %>%
      arrange(padj, desc(abs(log2FoldChange))) %>%
      head(top_n)
  } else {
    sel_df <- res_df %>%
      filter(baseMean >= baseMean_min) %>%
      arrange(desc(abs(log2FoldChange)), desc(baseMean)) %>%
      head(top_n)
    title <- paste0(title, " (exploratory)")
  }

  sel <- intersect(sel_df$feature_id, rownames(vst_mat))
  if (length(sel) < 2) {
    writeLines("not enough features for heatmap", paste0(outfile, ".txt"))
    save_placeholder_plot(outfile, title, "not enough features for heatmap")
    return(invisible(NULL))
  }

  mat <- vst_mat[sel, , drop = FALSE]
  row_ann <- NULL

  if (label_col %in% colnames(sel_df)) {
    lab_map <- sel_df[[label_col]]
    names(lab_map) <- sel_df$feature_id
    rownames(mat) <- ifelse(!is.na(lab_map[rownames(mat)]) & lab_map[rownames(mat)] != "", lab_map[rownames(mat)], rownames(mat))
  }

  if (!is.null(row_anno_col) && row_anno_col %in% colnames(sel_df)) {
    row_ann <- data.frame(group = sel_df[[row_anno_col]][match(sel, sel_df$feature_id)], stringsAsFactors = FALSE)
    rownames(row_ann) <- rownames(mat)
    colnames(row_ann) <- row_anno_col
  }

  ann <- data.frame(condition = sample_table$condition, row.names = sample_table$sample, stringsAsFactors = FALSE)

  pheatmap(
    mat,
    scale = "row",
    annotation_col = ann,
    annotation_row = row_ann,
    show_rownames = TRUE,
    main = title,
    filename = outfile,
    width = 8,
    height = 8
  )
}

run_go_simple <- function(res_df, prefix, species = "hg38") {
  sig <- res_df %>%
    dplyr::filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff)

  if (grepl("^hg|human", species, ignore.case = TRUE)) {
    OrgDb <- org.Hs.eg.db
  } else {
    OrgDb <- org.Mm.eg.db
  }

  run_one_go <- function(df_sub, suffix) {
    if (nrow(df_sub) < 5) {
      msg <- paste0("too few ", suffix, " sig genes for GO")
      message("[run_go_simple] ", msg)
      writeLines(msg, paste0(prefix, "_GO_", suffix, ".txt"))
      return(invisible(NULL))
    }

    gene_ids <- unique(sub("\\..*$", "", df_sub$feature_id))
    gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
    if (length(gene_ids) < 5) {
      msg <- paste0("too few valid ", suffix, " gene IDs for GO")
      message("[run_go_simple] ", msg)
      writeLines(msg, paste0(prefix, "_GO_", suffix, ".txt"))
      return(invisible(NULL))
    }

    ek <- suppressMessages(
      tryCatch(
        enrichGO(
          gene = gene_ids,
          OrgDb = OrgDb,
          keyType = "ENSEMBL",
          ont = "BP",
          pAdjustMethod = "BH",
          readable = TRUE
        ),
        error = function(e) {
          message("[run_go_simple] GO ", suffix, " failed: ", conditionMessage(e))
          NULL
        }
      )
    )

    if (is.null(ek) || nrow(as.data.frame(ek)) == 0) {
      msg <- paste0("no ", suffix, " GO enrichment")
      message("[run_go_simple] ", msg)
      writeLines(msg, paste0(prefix, "_GO_", suffix, ".txt"))
      return(invisible(NULL))
    }

    utils::write.csv(as.data.frame(ek), paste0(prefix, "_GO_", suffix, ".csv"), row.names = FALSE)
    p <- enrichplot::dotplot(ek, showCategory = 15) + ggtitle(paste0(basename(prefix), " GO ", suffix))
    ggsave(paste0(prefix, "_GO_", suffix, "_dotplot.pdf"), p, width = 8, height = 6)
    invisible(ek)
  }

  sig_up <- sig %>% dplyr::filter(log2FoldChange > 0)
  sig_down <- sig %>% dplyr::filter(log2FoldChange < 0)

  run_one_go(sig_up, "up")
  run_one_go(sig_down, "down")
  invisible(NULL)
}


get_orgdb_by_species <- function(species = "hg38") {
  if (grepl("^hg|human", species, ignore.case = TRUE)) org.Hs.eg.db else org.Mm.eg.db
}

get_msig_species_name <- function(species = "hg38") {
  if (grepl("^hg|human", species, ignore.case = TRUE)) "Homo sapiens" else "Mus musculus"
}

get_default_gsea_categories <- function() {
  list(
    list(category = "H",  subcategory = NULL,            label = "Hallmark"),
    list(category = "C2", subcategory = "CP:KEGG_LEGACY", label = "KEGG"),
    list(category = "C2", subcategory = "CP:REACTOME",   label = "Reactome"),
    list(category = "C5", subcategory = "GO:BP",         label = "GO_BP"),
    list(category = "C5", subcategory = "GO:CC",         label = "GO_CC"),
    list(category = "C5", subcategory = "GO:MF",         label = "GO_MF")
  )
}

wrap_term_label <- function(x, width = 45) {
  vapply(as.character(x), function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

get_msig_term2gene <- function(species = "hg38", category = "H", subcategory = NULL) {
  msig_species <- get_msig_species_name(species)

  load_msig <- function(cat, subcat) {
    fml <- names(formals(msigdbr::msigdbr))
    args <- list(species = msig_species)
    if ("collection" %in% fml) {
      args$collection <- cat
      if (!is.null(subcat) && nzchar(subcat)) args$subcollection <- subcat
    } else {
      args$category <- cat
      if (!is.null(subcat) && nzchar(subcat)) args$subcategory <- subcat
    }
    suppressMessages(do.call(msigdbr::msigdbr, args))
  }

  sub_try <- if (is.null(subcategory) || !nzchar(subcategory)) {
    list(NULL)
  } else if (subcategory == "CP:KEGG") {
    list("CP:KEGG_LEGACY", "CP:KEGG_MEDICUS", "CP:KEGG")
  } else {
    list(subcategory)
  }
  msig <- NULL
  last_error <- NULL
  for (subcat in sub_try) {
    msig <- tryCatch(
      load_msig(category, subcat),
      error = function(e) {
        last_error <<- e
        NULL
      }
    )
    if (!is.null(msig) && nrow(msig) > 0) break
  }
  if (is.null(msig) && !is.null(last_error)) stop(conditionMessage(last_error))

  if (is.null(msig) || nrow(msig) == 0) {
    return(data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE))
  }

  term_col <- if ("gs_name" %in% colnames(msig)) "gs_name" else colnames(msig)[1]
  gene_col <- if ("gene_symbol" %in% colnames(msig)) {
    "gene_symbol"
  } else if ("human_gene_symbol" %in% colnames(msig)) {
    "human_gene_symbol"
  } else {
    stop("msigdbr 返回结果中未找到 gene_symbol 列")
  }

  out <- msig[, c(term_col, gene_col), drop = FALSE]
  colnames(out) <- c("term", "gene")
  out$term <- as.character(out$term)
  out$gene <- as.character(out$gene)
  out <- out[!is.na(out$term) & out$term != "" & !is.na(out$gene) & out$gene != "", , drop = FALSE]
  unique(out)
}

map_feature_to_gene_symbol <- function(feature_ids, gene_anno = NULL, species = "hg38") {
  feature_ids <- as.character(feature_ids)
  gene_id <- sub("\\..*$", "", feature_ids)
  gene_symbol <- gene_id

  if (!is.null(gene_anno) && all(c("gene_id", "gene_name") %in% colnames(gene_anno))) {
    anno_df <- unique(gene_anno[, c("gene_id", "gene_name"), drop = FALSE])
    anno_df$gene_id <- sub("\\..*$", "", as.character(anno_df$gene_id))
    anno_df$gene_name <- as.character(anno_df$gene_name)
    hit <- anno_df$gene_name[match(gene_id, anno_df$gene_id)]
    ok <- !is.na(hit) & hit != "" & !grepl("^ENS", hit)
    gene_symbol[ok] <- hit[ok]
  }

  need_map <- grepl("^(ENSG|ENSMUSG|ENS[A-Z0-9]+)", gene_id) &
    (is.na(gene_symbol) | gene_symbol == "" | grepl("^ENS", gene_symbol))

  if (any(need_map)) {
    db_obj <- get_orgdb_by_species(species)
    sym_map <- AnnotationDbi::mapIds(
      db_obj,
      keys = unique(gene_id[need_map]),
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    sym_map <- as.character(sym_map)
    hit2 <- sym_map[gene_id[need_map]]
    ok2 <- !is.na(hit2) & hit2 != ""
    idx2 <- which(need_map)
    if (any(ok2)) gene_symbol[idx2[ok2]] <- hit2[ok2]
  }

  data.frame(
    feature_id = feature_ids,
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    stringsAsFactors = FALSE
  )
}

collapse_rank_to_gene_symbol <- function(rank_vec, gene_anno = NULL, species = "hg38") {
  rank_vec <- rank_vec[is.finite(rank_vec)]
  if (length(rank_vec) == 0) return(rank_vec)

  score_df <- data.frame(
    feature_id = as.character(names(rank_vec)),
    score = as.numeric(rank_vec),
    stringsAsFactors = FALSE
  )

  map_df <- map_feature_to_gene_symbol(score_df$feature_id, gene_anno = gene_anno, species = species)
  score_df <- score_df %>% left_join(map_df, by = "feature_id")

  miss <- is.na(score_df$gene_symbol) | score_df$gene_symbol == ""
  score_df$gene_symbol[miss] <- score_df$gene_id[miss]
  score_df <- score_df[!is.na(score_df$gene_symbol) & score_df$gene_symbol != "", , drop = FALSE]
  if (nrow(score_df) == 0) return(setNames(numeric(0), character(0)))

  score_df <- score_df %>%
    group_by(gene_symbol) %>%
    summarise(score = score[which.max(abs(score))], .groups = "drop")

  out <- score_df$score
  names(out) <- score_df$gene_symbol
  out <- out[is.finite(out)]
  out <- out + stats::runif(length(out), min = -1e-10, max = 1e-10)
  sort(out, decreasing = TRUE)
}

make_rank_from_fit <- function(fit, res_df, contrast) {
  if (!is.null(fit$mode) && identical(fit$mode, "exploratory") && !is.null(fit$vst_mat)) {
    st <- fit$sample_table
    rownames(st) <- st$sample
    g1 <- st$sample[st[[contrast[1]]] == contrast[2]]
    g0 <- st$sample[st[[contrast[1]]] == contrast[3]]

    mean1 <- if (length(g1) == 1) fit$vst_mat[, g1] else rowMeans(fit$vst_mat[, g1, drop = FALSE], na.rm = TRUE)
    mean0 <- if (length(g0) == 1) fit$vst_mat[, g0] else rowMeans(fit$vst_mat[, g0, drop = FALSE], na.rm = TRUE)

    v <- mean1 - mean0
    names(v) <- rownames(fit$vst_mat)
    message("[GSEA] use logCPM difference ranking (exploratory / no replicate)")
    return(sort(v[is.finite(v)], decreasing = TRUE))
  }

  if ("stat" %in% colnames(res_df) && any(is.finite(res_df$stat))) {
    v <- res_df$stat
    names(v) <- res_df$feature_id
    message("[GSEA] use DESeq2 stat ranking")
    return(sort(v[is.finite(v)], decreasing = TRUE))
  }

  if ("log2FoldChange" %in% colnames(res_df) && any(is.finite(res_df$log2FoldChange))) {
    v <- res_df$log2FoldChange
    names(v) <- res_df$feature_id
    message("[GSEA] fallback to log2FoldChange ranking")
    return(sort(v[is.finite(v)], decreasing = TRUE))
  }

  stop("GSEA 无法构建 ranking：缺少可用 stat / log2FC / logCPM 差值")
}

plot_gsea_bar_simple <- function(df, outfile, title, top_n = 20) {
  outfile <- as_pdf_path(outfile)
  df2 <- df %>%
    filter(!is.na(NES)) %>%
    arrange(p.adjust, desc(abs(NES))) %>%
    head(top_n) %>%
    arrange(NES)

  if (nrow(df2) == 0) {
    save_placeholder_plot(outfile, title, "no GSEA bar data")
    return(invisible(NULL))
  }

  df2$Description_wrap <- wrap_term_label(df2$Description)
  df2$Description_wrap <- factor(df2$Description_wrap, levels = df2$Description_wrap)

  p <- ggplot(df2, aes(x = Description_wrap, y = NES, fill = NES)) +
    geom_col() +
    coord_flip() +
    scale_fill_gradient2(low = "#4F789F", mid = "#F4F4F1", high = "#C95F50") +
    journal_theme(base_size = 12) +
    labs(title = title, x = NULL, y = "NES")

  ggsave(outfile, p, width = 10, height = 7)
  invisible(p)
}

plot_gsea_nes_bubble_simple <- function(gsea_tbl, outfile, title, top_n_each_sign = 10) {
  outfile <- as_pdf_path(outfile)
  df <- gsea_tbl %>%
    dplyr::filter(!is.na(p.adjust), is.finite(NES)) %>%
    dplyr::mutate(Sign = ifelse(NES > 0, 'Positive NES', 'Negative NES')) %>%
    dplyr::group_by(Sign) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_each_sign) %>%
    dplyr::ungroup()
  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, title, 'no significant GSEA terms for NES bubble')
    return(invisible(NULL))
  }
  df$Description_wrap <- wrap_term_label(df$Description, width = 36)
  p <- ggplot(df, aes(x = NES, y = reorder(Description_wrap, NES), color = p.adjust, size = setSize)) +
    geom_point(alpha = 0.85) +
    scale_color_gradientn(
    colors = c("#d50745","#cf4935","#FDDBC7","#D1E5F0","#2686c1", "#0050ac"),
    name = "FDR"
    ) +
    facet_grid(Sign ~ ., scales = 'free_y', space = 'free_y') +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5), axis.text.y = element_text(size = 8)) +
    labs(title = title, x = 'Normalized Enrichment Score (NES)', y = NULL, size = 'Gene set size')
  ggsave(outfile, p, width = 7.5, height = 8.5, limitsize = FALSE)
  invisible(p)
}

plot_gsea_nes_fdr_scatter_simple <- function(gsea_tbl, outfile, title, top_n_each_sign = 8) {
  outfile <- as_pdf_path(outfile)
  df <- gsea_tbl %>%
    dplyr::filter(!is.na(p.adjust), is.finite(NES)) %>%
    dplyr::mutate(Sign = ifelse(NES >= 0, 'Positive NES', 'Negative NES'))

  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, title, 'no significant GSEA terms for NES-FDR scatter')
    return(invisible(NULL))
  }

  pick_df <- df %>%
    dplyr::group_by(Sign) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_each_sign) %>%
    dplyr::ungroup()

  pick_df$label_color <- ifelse(pick_df$NES >= 0, '#D62728', '#2C7FB8')
  pick_df$Description_wrap <- wrap_term_label(pick_df$Description, width = 34)

  xmax <- max(df$p.adjust[df$p.adjust <= 0.2], na.rm = TRUE)
  if (!is.finite(xmax)) xmax <- max(df$p.adjust, na.rm = TRUE)
  xmax <- max(0.05, min(xmax * 1.1, 1))

  p <- ggplot(df, aes(x = p.adjust, y = NES)) +
    geom_point(color = '#2B3036', alpha = 0.65, size = 1.8) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.45, color = "#8A8A8A") +
    coord_cartesian(xlim = c(0, xmax)) +
    journal_theme(base_size = 12) +
    labs(title = title, x = 'FDR q value', y = 'NES')

  if (nrow(pick_df) > 0) {
    p <- p +
      geom_point(data = pick_df, aes(x = p.adjust, y = NES), color = pick_df$label_color, size = 2.2, inherit.aes = FALSE) +
      geom_text_repel(
        data = pick_df,
        aes(x = p.adjust, y = NES, label = Description_wrap),
        color = pick_df$label_color,
        size = 3.2,
        segment.size = 0.3,
        box.padding = 0.25,
        point.padding = 0.15,
        min.segment.length = 0,
        show.legend = FALSE,
        inherit.aes = FALSE,
        max.overlaps = Inf
      )
  }

  ggsave(outfile, p, width = 10.5, height = 6.2, limitsize = FALSE)
  invisible(p)
}

plot_gsea_bubble_simple <- function(gsea_all_tbl, outfile, title, top_n_per_cat = 6) {
  outfile <- as_pdf_path(outfile)
  df <- gsea_all_tbl %>%
    dplyr::filter(!is.na(p.adjust), is.finite(NES)) %>%
    dplyr::group_by(category_label) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES)), .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n_per_cat) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      logFDR = -log10(p.adjust + 1e-300),
      y_group = paste0(as.character(category_label), ' | ', Description)
    )
  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, title, 'no significant GSEA terms for collection overview')
    return(invisible(NULL))
  }
  y_levels <- df %>%
    dplyr::arrange(category_label, p.adjust, dplyr::desc(abs(NES))) %>%
    dplyr::pull(y_group) %>% unique()
  df$y_group <- factor(df$y_group, levels = rev(y_levels))
  p <- ggplot(df, aes(x = category_label, y = y_group, size = logFDR, color = NES)) +
    geom_point(alpha = 0.9) +
    scale_y_discrete(labels = function(x) sub('^[^|]+ \\| ', '', x)) +
    scale_color_gradient2(low = '#4F789F', mid = '#F4F4F1', high = '#C95F50', midpoint = 0) +
    journal_theme(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 8)) +
    labs(title = title, x = NULL, y = NULL, size = '-log10(FDR)', color = 'NES')
  ggsave(outfile, p, width = 11.5, height = max(8, 0.25 * length(y_levels) + 3), limitsize = FALSE)
  invisible(p)
}

select_gsea_representative_ids <- function(df, max_terms = 8) {
  ids_up <- df %>%
    dplyr::filter(!is.na(NES), NES > 0) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = ceiling(max_terms / 2)) %>%
    dplyr::pull(ID)
  ids_down <- df %>%
    dplyr::filter(!is.na(NES), NES < 0) %>%
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) %>%
    dplyr::slice_head(n = floor(max_terms / 2)) %>%
    dplyr::pull(ID)
  unique(c(ids_up, ids_down))
}

plot_gsea_across_contrast_heatmap_simple <- function(gsea_df, outfile, category_label = "Hallmark", top_n = 25) {
  outfile <- as_pdf_path(outfile)
  df <- gsea_df %>%
    filter(.data$category_label == category_label, !is.na(NES))

  if (nrow(df) == 0) {
    save_placeholder_plot(outfile, paste0(category_label, " across contrasts"), "no GSEA heatmap data")
    return(invisible(NULL))
  }

  term_keep <- df %>%
    group_by(Description) %>%
    summarise(max_abs_nes = max(abs(NES), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_nes)) %>%
    slice_head(n = top_n) %>%
    pull(Description)

  mat_df <- df %>%
    filter(Description %in% term_keep) %>%
    select(contrast, Description, NES) %>%
    distinct() %>%
    tidyr::pivot_wider(names_from = contrast, values_from = NES, values_fill = 0)

  if (nrow(mat_df) < 2 || ncol(mat_df) < 2) {
    save_placeholder_plot(outfile, paste0(category_label, " across contrasts"), "not enough contrasts/terms")
    return(invisible(NULL))
  }

  mat <- as.matrix(mat_df[, -1, drop = FALSE])
  rownames(mat) <- wrap_term_label(mat_df$Description, width = 40)

  pheatmap(
    mat,
    main = paste0(category_label, " NES across contrasts"),
    filename = outfile,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    width = 9,
    height = max(6, 0.24 * nrow(mat) + 2)
  )
  invisible(mat)
}

run_gsea_multi_category <- function(fit, res_df, contrast, prefix, species = "hg38", gene_anno = NULL,
                                    gsea_categories = NULL, minGSSize = 10, maxGSSize = 500,
                                    pvalueCutoff = 1, topN_plot = 20, max_terms_gseaplot2 = 12,
                                    make_gseaplot2 = TRUE) {
  if (is.null(gsea_categories)) gsea_categories <- get_default_gsea_categories()

  contrast_name <- paste0(contrast[2], "_vs_", contrast[3])
  out_dir <- paste0(prefix, "_GSEA")
  overview_dir <- file.path(out_dir, "overview")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)

  rank_vec <- make_rank_from_fit(fit, res_df, contrast)
  rank_vec <- collapse_rank_to_gene_symbol(rank_vec, gene_anno = gene_anno, species = species)

  if (length(rank_vec) < 20) {
    msg <- paste0("too few ranked genes after symbol collapse for GSEA: ", length(rank_vec))
    message("[run_gsea_multi_category] ", msg)
    writeLines(msg, file.path(overview_dir, "GSEA.txt"))
    return(data.frame())
  }

  utils::write.csv(data.frame(gene_symbol = names(rank_vec), score = as.numeric(rank_vec), stringsAsFactors = FALSE),
                   file.path(overview_dir, "ranked_gene_list.csv"), row.names = FALSE)

  all_tables <- list()

  for (cfg in gsea_categories) {
    cat_label <- cfg$label

    term2gene <- tryCatch(
      get_msig_term2gene(species = species, category = cfg$category, subcategory = cfg$subcategory),
      error = function(e) {
        message("[run_gsea_multi_category] term2gene load failed for ", cat_label, ": ", conditionMessage(e))
        return(data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE))
      }
    )

    if (nrow(term2gene) == 0) {
      message("[run_gsea_multi_category] skip empty gene set category: ", cat_label)
      next
    }

    overlap_n <- sum(unique(term2gene$gene) %in% names(rank_vec))
    message("[run_gsea_multi_category] ", contrast_name, " | ", cat_label, " overlap genes = ", overlap_n)
    if (overlap_n < 10) {
      message("[run_gsea_multi_category] skip ", cat_label, " because overlap < 10")
      next
    }

    gse <- tryCatch(
      clusterProfiler::GSEA(
        geneList = rank_vec,
        TERM2GENE = term2gene,
        minGSSize = minGSSize,
        maxGSSize = maxGSSize,
        pvalueCutoff = pvalueCutoff,
        verbose = FALSE,
        eps = 0
      ),
      error = function(e) {
        message("[run_gsea_multi_category] ", contrast_name, " | ", cat_label, " failed: ", conditionMessage(e))
        return(NULL)
      }
    )

    if (is.null(gse)) next
    df <- as.data.frame(gse)
    if (nrow(df) == 0) {
      message("[run_gsea_multi_category] ", contrast_name, " | ", cat_label, " returned 0 pathways")
      next
    }

    df$contrast <- contrast_name
    df$category <- cfg$category
    df$subcategory <- if (is.null(cfg$subcategory)) NA_character_ else cfg$subcategory
    df$category_label <- cat_label
    all_tables[[cat_label]] <- df

    cdir <- file.path(out_dir, sanitize_name(cat_label))
    dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(df, file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_all.csv')), row.names = FALSE)
    utils::write.csv(df %>% dplyr::filter(!is.na(p.adjust), p.adjust < 0.05),
                     file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_sig.csv')), row.names = FALSE)

    tryCatch({
      plot_gsea_nes_bubble_simple(
        df,
        file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_bubble.pdf')),
        paste0(contrast_name, ' ', cat_label, ' GSEA'),
        top_n_each_sign = max(6, min(topN_plot, 10))
      )
    }, error = function(e) {
      message('[run_gsea_multi_category] NES bubble failed: ', conditionMessage(e))
    })

    tryCatch({
      plot_gsea_nes_fdr_scatter_simple(
        df,
        file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_NES_vs_FDR.pdf')),
        paste0(contrast_name, ' ', cat_label, ' NES vs FDR'),
        top_n_each_sign = max(6, min(topN_plot, 10))
      )
    }, error = function(e) {
      message('[run_gsea_multi_category] NES vs FDR scatter failed: ', conditionMessage(e))
    })

    if (make_gseaplot2) {
      ids <- select_gsea_representative_ids(df, max_terms = max_terms_gseaplot2)
      if (length(ids) > 0) {
        gdir <- file.path(cdir, 'running_score')
        dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
        for (id in ids) {
          row1 <- df[df$ID == id, , drop = FALSE]
          if (nrow(row1) == 0) next
          ttl <- paste0(
            contrast_name, ' | ', cat_label, ' | ', row1$Description[1],
            ' | NES=', round(row1$NES[1], 3),
            ' | FDR=', signif(row1$p.adjust[1], 3)
          )
          tryCatch({
            p <- enrichplot::gseaplot2(gse, geneSetID = id, title = ttl, pvalue_table = TRUE)
            ggsave(as_pdf_path(file.path(gdir, paste0('gseaplot2_', sanitize_name(id), '.pdf'))), p, width = 10, height = 6)
          }, error = function(e) {
            message('[run_gsea_multi_category] gseaplot2 failed: ', conditionMessage(e))
          })
        }
      }
    }
  }

  if (length(all_tables) == 0) {
    writeLines('no multi-category GSEA enrichment', file.path(overview_dir, 'GSEA.txt'))
    return(data.frame())
  }

  out_df <- bind_rows(all_tables)
  utils::write.csv(out_df, file.path(overview_dir, 'GSEA_all_categories.csv'), row.names = FALSE)

  tryCatch({
    plot_gsea_bubble_simple(
      out_df,
      file.path(overview_dir, 'GSEA_collection_overview_bubble.pdf'),
      paste0(contrast_name, ' | multi-category GSEA'),
      top_n_per_cat = 6
    )
  }, error = function(e) {
    message('[run_gsea_multi_category] overview bubble plot failed: ', conditionMessage(e))
  })

  out_df
}

run_gsea_simple <- function(res_df, prefix, species = "hg38", fit = NULL, contrast = NULL, gene_anno = NULL) {
  if (is.null(fit) || is.null(contrast)) {
    message("[run_gsea_simple] missing fit/contrast, fallback disabled")
    writeLines("missing fit/contrast for GSEA", paste0(prefix, "_GSEA.txt"))
    return(invisible(data.frame()))
  }

  run_gsea_multi_category(
    fit = fit,
    res_df = res_df,
    contrast = contrast,
    prefix = prefix,
    species = species,
    gene_anno = gene_anno,
    gsea_categories = get_default_gsea_categories(),
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    topN_plot = 20,
    max_terms_gseaplot2 = 12,
    make_gseaplot2 = TRUE
  )
}

run_gsea_from_counts_pipeline <- function(counts_csv, out_dir, control_level, species = "hg38",
                                          ranking_method = c("logCPM_diff", "edgeR_fixedBCV"),
                                          fixed_BCV = 0.2, gsea_categories = NULL,
                                          minGSSize = 10, maxGSSize = 500,
                                          pvalueCutoff = 1, topN_plot = 20,
                                          make_gseaplot2 = TRUE, max_terms_gseaplot2 = 30,
                                          gene_anno = NULL) {
  ranking_method <- match.arg(ranking_method)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("[run_gsea_from_counts_pipeline] counts_csv=", counts_csv)
  message("[run_gsea_from_counts_pipeline] out_dir=", out_dir)
  message("[run_gsea_from_counts_pipeline] control_level=", control_level)
  message("[run_gsea_from_counts_pipeline] ranking_method=", ranking_method)

  expr <- read.csv(counts_csv, check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(expr) < 3) stop("输入 counts 至少需要 1 列基因 + 2 列样本")

  genes <- as.character(expr[[1]])
  mat <- as.matrix(expr[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- genes

  if (!control_level %in% colnames(mat)) stop("control_level 不在样本列中: ", control_level)

  if (any(duplicated(rownames(mat)))) {
    message("[run_gsea_from_counts_pipeline] collapse duplicated rows by sum")
    df_tmp <- as.data.frame(mat, check.names = FALSE) %>% tibble::rownames_to_column("gene")
    df_tmp <- df_tmp %>%
      group_by(gene) %>%
      summarise(across(everything(), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")
    mat <- as.matrix(df_tmp[, -1, drop = FALSE])
    storage.mode(mat) <- "numeric"
    rownames(mat) <- df_tmp$gene
  }

  keep <- rowSums(mat, na.rm = TRUE) >= 10
  mat <- mat[keep, , drop = FALSE]
  message("[run_gsea_from_counts_pipeline] matrix after filter = ", nrow(mat), " x ", ncol(mat))

  dge <- edgeR::DGEList(counts = round(mat))
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

  make_rank <- function(ko) {
    if (identical(ranking_method, "logCPM_diff")) {
      v <- logCPM[, ko] - logCPM[, control_level]
      names(v) <- rownames(logCPM)
      return(sort(v[is.finite(v)], decreasing = TRUE))
    }

    phi <- fixed_BCV^2
    x <- round(mat[, c(control_level, ko), drop = FALSE])
    grp <- factor(c("CTRL", "KO"))
    y <- edgeR::DGEList(counts = x, group = grp)
    y <- edgeR::calcNormFactors(y, method = "TMM")
    et <- edgeR::exactTest(y, dispersion = phi)
    tt <- edgeR::topTags(et, n = Inf)$table
    v <- tt$logFC
    names(v) <- rownames(tt)
    sort(v[is.finite(v)], decreasing = TRUE)
  }

  if (is.null(gsea_categories)) gsea_categories <- get_default_gsea_categories()

  ko_levels <- setdiff(colnames(mat), control_level)
  all_gsea <- list()

  for (ko in ko_levels) {
    contrast_name <- paste0(ko, "_vs_", control_level)
    message("[run_gsea_from_counts_pipeline] running ", contrast_name)

    rank_vec <- make_rank(ko)
    rank_vec <- collapse_rank_to_gene_symbol(rank_vec, gene_anno = gene_anno, species = species)
    if (length(rank_vec) < 20) {
      message("[run_gsea_from_counts_pipeline] skip ", contrast_name, " because ranked genes < 20")
      next
    }

    prefix <- file.path(out_dir, sanitize_name(contrast_name))
    overview_dir <- file.path(prefix, 'overview')
    dir.create(prefix, recursive = TRUE, showWarnings = FALSE)
    dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(data.frame(gene_symbol = names(rank_vec), score = as.numeric(rank_vec), stringsAsFactors = FALSE),
                     file.path(overview_dir, 'ranked_gene_list.csv'), row.names = FALSE)

    tables_one <- list()
    for (cfg in gsea_categories) {
      cat_label <- cfg$label
      term2gene <- tryCatch(
        get_msig_term2gene(species = species, category = cfg$category, subcategory = cfg$subcategory),
        error = function(e) {
          message("[run_gsea_from_counts_pipeline] term2gene failed for ", cat_label, ": ", conditionMessage(e))
          return(data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE))
        }
      )
      if (nrow(term2gene) == 0) next

      overlap_n <- sum(unique(term2gene$gene) %in% names(rank_vec))
      message("[run_gsea_from_counts_pipeline] ", contrast_name, " | ", cat_label, " overlap genes = ", overlap_n)
      if (overlap_n < 10) next

      gse <- tryCatch(
        clusterProfiler::GSEA(
          geneList = rank_vec,
          TERM2GENE = term2gene,
          minGSSize = minGSSize,
          maxGSSize = maxGSSize,
          pvalueCutoff = pvalueCutoff,
          verbose = FALSE,
          eps = 0
        ),
        error = function(e) {
          message("[run_gsea_from_counts_pipeline] GSEA failed for ", contrast_name, " | ", cat_label, ": ", conditionMessage(e))
          return(NULL)
        }
      )
      if (is.null(gse)) next

      df <- as.data.frame(gse)
      if (nrow(df) == 0) next

      df$contrast <- contrast_name
      df$category <- cfg$category
      df$subcategory <- if (is.null(cfg$subcategory)) NA_character_ else cfg$subcategory
      df$category_label <- cat_label
      tables_one[[cat_label]] <- df

      cdir <- file.path(prefix, sanitize_name(cat_label))
      dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
      utils::write.csv(df, file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_all.csv')), row.names = FALSE)
      utils::write.csv(df %>% dplyr::filter(!is.na(p.adjust), p.adjust < 0.05),
                       file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_sig.csv')), row.names = FALSE)

      tryCatch({
        plot_gsea_nes_bubble_simple(
          df,
          file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_bubble.pdf')),
          paste0(contrast_name, ' ', cat_label, ' GSEA'),
          top_n_each_sign = max(6, min(topN_plot, 10))
        )
      }, error = function(e) {
        message('[run_gsea_from_counts_pipeline] NES bubble failed: ', conditionMessage(e))
      })

      tryCatch({
        plot_gsea_nes_fdr_scatter_simple(
          df,
          file.path(cdir, paste0('GSEA_', sanitize_name(cat_label), '_NES_vs_FDR.pdf')),
          paste0(contrast_name, ' ', cat_label, ' NES vs FDR'),
          top_n_each_sign = max(6, min(topN_plot, 10))
        )
      }, error = function(e) {
        message('[run_gsea_from_counts_pipeline] NES vs FDR scatter failed: ', conditionMessage(e))
      })

      if (make_gseaplot2) {
        ids <- select_gsea_representative_ids(df, max_terms = max_terms_gseaplot2)
        if (length(ids) > 0) {
          gdir <- file.path(cdir, 'running_score')
          dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
          for (id in ids) {
            row1 <- df[df$ID == id, , drop = FALSE]
            if (nrow(row1) == 0) next
            ttl <- paste0(
              contrast_name, ' | ', cat_label, ' | ', row1$Description[1],
              ' | NES=', round(row1$NES[1], 3),
              ' | FDR=', signif(row1$p.adjust[1], 3)
            )
            tryCatch({
              p <- enrichplot::gseaplot2(gse, geneSetID = id, title = ttl, pvalue_table = TRUE)
              ggsave(as_pdf_path(file.path(gdir, paste0('gseaplot2_', sanitize_name(id), '.pdf'))), p, width = 10, height = 6)
            }, error = function(e) {
              message('[run_gsea_from_counts_pipeline] gseaplot2 failed: ', conditionMessage(e))
            })
          }
        }
      }
    }

    if (length(tables_one) > 0) {
      gsea_df <- bind_rows(tables_one)
      utils::write.csv(gsea_df, file.path(overview_dir, 'GSEA_all_categories.csv'), row.names = FALSE)

      tryCatch({
        plot_gsea_bubble_simple(
          gsea_df,
          file.path(overview_dir, 'GSEA_collection_overview_bubble.pdf'),
          paste0(contrast_name, ' | multi-category GSEA'),
          top_n_per_cat = 6
        )
      }, error = function(e) {
        message('[run_gsea_from_counts_pipeline] overview bubble plot failed: ', conditionMessage(e))
      })

      all_gsea[[contrast_name]] <- gsea_df
    }
  }

  if (length(all_gsea) > 0) {
    gsea_all <- bind_rows(all_gsea)
    utils::write.csv(gsea_all, file.path(out_dir, 'GSEA_all_contrasts_all_categories.csv'), row.names = FALSE)
    return(gsea_all)
  }

  invisible(data.frame())
}

read_featurecounts_matrix <- function(files) {
  mats <- lapply(files, function(f) {
    dt <- fread(f, data.table = FALSE)
    if (ncol(dt) < 7) stop("featureCounts 文件列数异常: ", basename(f))
    id_col <- colnames(dt)[1]
    cnt_col <- tail(colnames(dt), 1)
    sample_name <- sub("\\.(featureCounts|featurecounts)\\.txt$", "", basename(f), ignore.case = TRUE)
    out <- dt[, c(id_col, cnt_col)]
    colnames(out) <- c("feature_id", sample_name)
    out
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  finalize_feature_matrix(merged, tag = "featureCounts")
}

read_telocal_matrix <- function(files) {
  mats <- lapply(files, function(f) {
    dt <- fread(f, data.table = FALSE)
    if (ncol(dt) < 2) stop("TElocal 文件列数异常: ", basename(f))
    id_col <- colnames(dt)[1]
    num_idx <- which(vapply(dt, is.numeric, logical(1)))
    if (length(num_idx) == 0) stop("TElocal 文件无数值列: ", basename(f))
    cnt_col <- colnames(dt)[num_idx[length(num_idx)]]
    sample_name <- sub("_telocal.*$", "", basename(f))
    out <- dt[, c(id_col, cnt_col), drop = FALSE]
    colnames(out) <- c("feature_id", sample_name)
    out
  })
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  finalize_feature_matrix(merged, tag = "TElocal")
}

parse_telocal_feature_map <- function(feature_ids, anno_df = NULL) {
  x <- as.character(feature_ids)
  x <- gsub('^"|"$', "", x)
  x <- sub("^\\s+|\\s+$", "", x)
  spl <- strsplit(x, ":", fixed = TRUE)

  out <- data.frame(
    feature_id = x,
    locus_id = NA_character_,
    repName = NA_character_,
    repFamily = NA_character_,
    repClass = NA_character_,
    gene_id = NA_character_,
    stringsAsFactors = FALSE
  )

  is_te <- lengths(spl) >= 4
  if (any(is_te)) {
    out$locus_id[is_te] <- vapply(spl[is_te], function(z) z[1], character(1))
    out$repName[is_te] <- vapply(spl[is_te], function(z) z[2], character(1))
    out$repFamily[is_te] <- vapply(spl[is_te], function(z) z[3], character(1))
    out$repClass[is_te] <- vapply(spl[is_te], function(z) paste(z[4:length(z)], collapse = ":"), character(1))
  }

  gene_like <- grepl("^ENS(G|MUSG|T|MUST)", x)
  out$gene_id[!is_te & gene_like] <- sub("\\.[0-9]+$", "", x[!is_te & gene_like])

  short_bad <- !is.na(out$repName) & is_short_repeat_name(out$repName)
  if (any(short_bad)) {
    message("[parse_telocal_feature_map] drop short-repeat TE rows: ", sum(short_bad))
    out[short_bad, c("locus_id", "repName", "repFamily", "repClass")] <- NA_character_
  }

  if (!is.null(anno_df) && "locus_id" %in% colnames(anno_df)) {
    locus_norm <- unique(normalize_id_vec(anno_df$locus_id))
    te_ok <- is_te & normalize_id_vec(out$locus_id) %in% locus_norm
    out[is_te & !te_ok, c("locus_id", "repName", "repFamily", "repClass")] <- NA_character_
    message("[parse_telocal_feature_map] compound TE rows=", sum(is_te), "/", length(x), "; locus matched=", sum(te_ok), "/", sum(is_te))
  } else {
    message("[parse_telocal_feature_map] compound TE rows=", sum(is_te), "/", length(x))
  }

  out
}

parse_tetranscripts_feature_map <- function(feature_ids, anno_df = NULL) {
  x <- as.character(feature_ids)
  x <- gsub('^"|"$', "", x)
  x <- sub("^\\s+|\\s+$", "", x)
  spl <- strsplit(x, ":", fixed = TRUE)

  out <- data.frame(
    feature_id = x,
    repName = NA_character_,
    repFamily = NA_character_,
    repClass = NA_character_,
    gene_id = NA_character_,
    stringsAsFactors = FALSE
  )

  is_te <- lengths(spl) >= 3
  if (any(is_te)) {
    out$repName[is_te] <- vapply(spl[is_te], function(z) z[1], character(1))
    out$repFamily[is_te] <- vapply(spl[is_te], function(z) z[2], character(1))
    out$repClass[is_te] <- vapply(spl[is_te], function(z) paste(z[3:length(z)], collapse = ":"), character(1))
  }

  gene_like <- grepl("^ENS(G|MUSG|T|MUST)", x)
  out$gene_id[!is_te & gene_like] <- sub("\\.[0-9]+$", "", x[!is_te & gene_like])

  short_bad <- !is.na(out$repName) & is_short_repeat_name(out$repName)
  if (any(short_bad)) {
    message("[parse_tetranscripts_feature_map] drop short-repeat TE rows: ", sum(short_bad))
    out[short_bad, c("repName", "repFamily", "repClass")] <- NA_character_
  }

  if (!is.null(anno_df) && "repName" %in% colnames(anno_df)) {
    rep_norm <- unique(normalize_id_vec(anno_df$repName))
    te_ok <- is_te & normalize_id_vec(out$repName) %in% rep_norm
    out[is_te & !te_ok, c("repName", "repFamily", "repClass")] <- NA_character_
    message("[parse_tetranscripts_feature_map] compound TE rows=", sum(is_te), "/", length(x), "; repName matched=", sum(te_ok), "/", sum(is_te))
  } else {
    message("[parse_tetranscripts_feature_map] compound TE rows=", sum(is_te), "/", length(x))
  }

  out
}

filter_te_rows_by_map <- function(count_mat, feature_map, tag = "TE_filter_by_map") {
  keep <- rep(FALSE, nrow(feature_map))
  if ("locus_id" %in% colnames(feature_map)) {
    keep <- keep | (!is.na(feature_map$locus_id) & feature_map$locus_id != "")
  }
  if ("repName" %in% colnames(feature_map)) {
    keep <- keep | (!is.na(feature_map$repName) & feature_map$repName != "")
  }
  if ("repName" %in% colnames(feature_map)) {
    bad_sr <- !is.na(feature_map$repName) & is_short_repeat_name(feature_map$repName)
    if (any(bad_sr)) message("[", tag, "] drop short-repeat parsed TE rows: ", sum(bad_sr))
    keep <- keep & !bad_sr
  }
  message("[", tag, "] keep parsed TE rows: ", sum(keep), "/", nrow(feature_map))
  if (sum(keep) == 0) stop("[", tag, "] 0 个 TE 行被识别到")
  out <- count_mat[feature_map$feature_id[keep], , drop = FALSE]
  attr(out, "feature_map") <- feature_map[keep, , drop = FALSE]
  out
}

aggregate_by_feature_map <- function(count_mat, feature_map, key_col, tag = "aggregate_by_feature_map") {
  if (!key_col %in% colnames(feature_map)) stop("[", tag, "] feature_map 中缺少列: ", key_col)

  cnt_df <- data.frame(feature_id = rownames(count_mat), count_mat, check.names = FALSE, stringsAsFactors = FALSE)
  df <- cnt_df %>%
    dplyr::left_join(feature_map[, unique(c("feature_id", key_col)), drop = FALSE], by = "feature_id") %>%
    dplyr::filter(!is.na(.data[[key_col]]), .data[[key_col]] != "")

  if (nrow(df) == 0) stop("[", tag, "] 无法按 ", key_col, " 聚合：匹配后 0 行")

  count_cols <- setdiff(colnames(df), c("feature_id", key_col))
  agg <- df %>%
    group_by(.data[[key_col]]) %>%
    summarise(across(all_of(count_cols), ~ sum(as.numeric(.x), na.rm = TRUE)), .groups = "drop")

  agg <- as.data.frame(agg, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(agg) <- agg[[key_col]]
  mat <- as.matrix(agg[, setdiff(colnames(agg), key_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

choose_repname_annotation <- function(feature_ids, te_anno, alt_anno = NULL, alt_match_cols = c("feature_id", "repName"), tag = "choose_repname_annotation") {
  count_hits <- function(map_df) {
    if (is.null(map_df)) return(-1L)
    annot_cols <- intersect(c("locus_id","repName","repFamily","repClass","gene_name"), colnames(map_df))
    if (length(annot_cols) == 0) return(0L)
    sum(rowSums(!is.na(map_df[, annot_cols, drop = FALSE]) & map_df[, annot_cols, drop = FALSE] != "") > 0)
  }

  te_map <- tryCatch(match_annotation_table(feature_ids, te_anno, preferred_match_cols = c("repName", "locus_id", "gene_name")), error = function(e) NULL)
  te_n <- count_hits(te_map)

  alt_n <- -1L
  if (!is.null(alt_anno)) {
    alt_map <- tryCatch(match_annotation_table(feature_ids, alt_anno, preferred_match_cols = alt_match_cols), error = function(e) NULL)
    alt_n <- count_hits(alt_map)
  }

  message("[", tag, "] te_anno matches=", te_n, "/", length(feature_ids),
          if (!is.null(alt_anno)) paste0("; alt matches=", alt_n, "/", length(feature_ids)) else "")

  if (!is.null(alt_anno) && alt_n > te_n) alt_anno else te_anno
}

match_annotation_table <- function(feature_ids, anno_df, preferred_match_cols = NULL) {
  keep_cols <- intersect(c("locus_id","repName","repFamily","repClass","gene_name","transcript_id","feature_id","gene_id","milliDiv"), colnames(anno_df))
  if (length(keep_cols) == 0) stop("annotation 中没有可用列")

  feat_df <- data.frame(
    feature_id = as.character(feature_ids),
    feature_norm = normalize_id_vec(feature_ids),
    stringsAsFactors = FALSE
  )

  candidate_cols <- if (is.null(preferred_match_cols)) {
    keep_cols
  } else {
    intersect(preferred_match_cols, keep_cols)
  }
  if (length(candidate_cols) == 0) stop("annotation 中没有可用匹配列")

  best_df <- feat_df
  best_n <- -1
  best_col <- NA_character_

  for (mcol in candidate_cols) {
    anno2 <- anno_df %>% dplyr::select(all_of(unique(c(mcol, keep_cols))))
    anno2$match_norm <- normalize_id_vec(anno2[[mcol]])
    anno2 <- anno2 %>%
      dplyr::filter(!is.na(match_norm), match_norm != "") %>%
      dplyr::distinct(match_norm, .keep_all = TRUE)

    if ("feature_id" %in% colnames(anno2)) {
      colnames(anno2)[colnames(anno2) == "feature_id"] <- "anno_feature_id"
    }

    df_try <- feat_df %>%
      dplyr::left_join(anno2[, c("match_norm", setdiff(colnames(anno2), "match_norm")), drop = FALSE], by = c("feature_norm" = "match_norm"))

    df_try$feature_id <- feat_df$feature_id
    annot_cols <- intersect(c("locus_id","repName","repFamily","repClass","gene_name"), colnames(df_try))
    n_match <- if (length(annot_cols) == 0) 0 else sum(rowSums(!is.na(df_try[, annot_cols, drop = FALSE]) & df_try[, annot_cols, drop = FALSE] != "") > 0)

    if (n_match > best_n) {
      best_df <- df_try
      best_n <- n_match
      best_col <- mcol
    }
  }

  best_df$feature_id <- feat_df$feature_id
  best_df$match_col <- best_col
  message("[match_annotation_table] matched via ", best_col, "; n=", best_n, "/", nrow(best_df))
  best_df
}

filter_te_rows <- function(count_mat, anno_df, preferred_match_cols = NULL, tag = "TE_filter") {
  map_df <- match_annotation_table(rownames(count_mat), anno_df, preferred_match_cols = preferred_match_cols)
  annot_cols <- intersect(c("locus_id","repName","repFamily","repClass","gene_name"), colnames(map_df))
  keep <- if (length(annot_cols) == 0) rep(FALSE, nrow(map_df)) else rowSums(!is.na(map_df[, annot_cols, drop = FALSE]) & map_df[, annot_cols, drop = FALSE] != "") > 0
  message("[", tag, "] keep annotated TE rows: ", sum(keep), "/", nrow(map_df))
  if (sum(keep) == 0) stop("[", tag, "] 0 个 TE 行被注释到")
  out <- count_mat[map_df$feature_id[keep], , drop = FALSE]
  attr(out, "feature_map") <- map_df[keep, , drop = FALSE]
  out
}

aggregate_by_annotation <- function(count_mat, anno_df, key_col, preferred_match_cols = NULL) {
  if (!key_col %in% colnames(anno_df)) stop("annotation 中缺少列: ", key_col)

  map_df <- match_annotation_table(rownames(count_mat), anno_df, preferred_match_cols = preferred_match_cols)
  cnt_df <- data.frame(feature_id = rownames(count_mat), count_mat, check.names = FALSE, stringsAsFactors = FALSE)
  df <- cnt_df %>%
    dplyr::left_join(map_df[, c("feature_id", key_col), drop = FALSE], by = "feature_id") %>%
    dplyr::filter(!is.na(.data[[key_col]]), .data[[key_col]] != "")

  if (nrow(df) == 0) stop("无法按 ", key_col, " 聚合：匹配后 0 行")

  count_cols <- setdiff(colnames(df), c("feature_id", key_col))
  agg <- df %>%
    group_by(.data[[key_col]]) %>%
    summarise(across(all_of(count_cols), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  agg <- as.data.frame(agg, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(agg) <- agg[[key_col]]
  mat <- as.matrix(agg[, setdiff(colnames(agg), key_col), drop = FALSE])
  storage.mode(mat) <- "numeric"
  mat
}

read_te_annotation <- function(path) {
  dt <- fread(path, data.table = FALSE)
  if (!"locus_id" %in% colnames(dt)) stop("te_annotation.tsv 缺少 locus_id")
  for (nm in intersect(c("locus_id","transcript_id","repName","repFamily","repClass","gene_name"), colnames(dt))) {
    dt[[nm]] <- as.character(dt[[nm]])
  }
  if ("milliDiv" %in% colnames(dt)) dt$milliDiv <- suppressWarnings(as.numeric(dt$milliDiv))
  dt <- unique(dt)
  dt <- drop_short_repeat_annotation(dt, tag = "read_te_annotation")
  message("[read_te_annotation] rows=", nrow(dt))
  dt
}

read_salmonte_clades <- function(files) {
  if (length(files) == 0) return(NULL)

  dt <- fread(files[1], data.table = FALSE)
  if (nrow(dt) == 0 || ncol(dt) == 0) return(NULL)

  raw_cols <- colnames(dt)
  std_cols <- tolower(gsub("[^a-z0-9]+", "_", raw_cols))
  colnames(dt) <- std_cols

  pick_first <- function(cands) {
    hit <- intersect(cands, colnames(dt))
    if (length(hit) == 0) return(NA_character_)
    hit[1]
  }

  feature_col <- pick_first(c("te","name","subfamily","element","repname"))
  clade_col   <- pick_first(c("clade","group"))
  class_col   <- pick_first(c("class","superfamily","type"))
  if (is.na(feature_col)) {
    message("[read_salmonte_clades] 未识别 feature 列，跳过 clades.csv 注释")
    return(NULL)
  }

  out <- data.frame(
    feature_id = as.character(dt[[feature_col]]),
    repName = as.character(dt[[feature_col]]),
    repFamily = if (!is.na(clade_col)) as.character(dt[[clade_col]]) else if (!is.na(class_col)) as.character(dt[[class_col]]) else NA_character_,
    repClass = if (!is.na(class_col)) as.character(dt[[class_col]]) else if (!is.na(clade_col)) as.character(dt[[clade_col]]) else NA_character_,
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$feature_id) & out$feature_id != "", , drop = FALSE]
  out <- unique(out)
  message("[read_salmonte_clades] rows=", nrow(out))
  out
}


pick_col_or_na <- function(df, nm) {
  if (nm %in% colnames(df)) as.character(df[[nm]]) else rep(NA_character_, nrow(df))
}

parse_telocal_compound_id <- function(x) {
  x <- as.character(x)
  x <- gsub('^"|"$', "", x)
  spl <- strsplit(x, ":", fixed = TRUE)

  data.frame(
    feature_id = x,
    locus_id = vapply(spl, function(z) if (length(z) >= 1) z[1] else NA_character_, character(1)),
    repName = vapply(spl, function(z) if (length(z) >= 2) z[2] else NA_character_, character(1)),
    repFamily = vapply(spl, function(z) if (length(z) >= 3) z[3] else NA_character_, character(1)),
    repClass = vapply(spl, function(z) if (length(z) >= 4) paste(z[4:length(z)], collapse = ":") else NA_character_, character(1)),
    stringsAsFactors = FALSE
  )
}

annotate_te_res <- function(res_df, te_anno, preferred_match_cols = NULL, label_level = NULL, color_level = NULL) {
  map_df <- match_annotation_table(res_df$feature_id, te_anno, preferred_match_cols = preferred_match_cols)
  keep_cols <- intersect(c("feature_id","match_col","locus_id","repName","repFamily","repClass","gene_name","milliDiv"), colnames(map_df))
  out <- res_df %>% left_join(map_df[, keep_cols, drop = FALSE], by = "feature_id")
    # 对 TElocal locus 这种复合 ID 做兜底拆分
  parsed_telocal <- parse_telocal_compound_id(out$feature_id)

  for (nm in c("locus_id", "repName", "repFamily", "repClass")) {
    if (!nm %in% colnames(out)) out[[nm]] <- NA_character_
    miss <- is.na(out[[nm]]) | out[[nm]] == ""
    out[[nm]][miss] <- parsed_telocal[[nm]][miss]
  }

  default_label <- dplyr::coalesce(
    pick_col_or_na(out, "repName"),
    pick_col_or_na(out, "gene_name"),
    pick_col_or_na(out, "locus_id"),
    pick_col_or_na(out, "repFamily"),
    pick_col_or_na(out, "repClass"),
    gsub("\\|", ":", as.character(out$feature_id))
  )

  if (!is.null(label_level) && label_level %in% colnames(out)) {
    out$te_label_plot <- dplyr::coalesce(pick_col_or_na(out, label_level), default_label)
  } else {
    out$te_label_plot <- default_label
  }

  if (!is.null(color_level) && color_level %in% colnames(out)) {
    out$te_color_plot <- pick_col_or_na(out, color_level)
  } else {
    out$te_color_plot <- dplyr::coalesce(
      pick_col_or_na(out, "repFamily"),
      pick_col_or_na(out, "repClass"),
      pick_col_or_na(out, "repName")
    )
  }
  out$te_color_plot[is.na(out$te_color_plot) | out$te_color_plot == ""] <- "Unknown"
  out$te_heatmap_group <- out$te_color_plot
  out
}

extract_rds_matrix <- function(path) {
  obj <- readRDS(path)

  # 直接就是 matrix/data.frame
  if (is.matrix(obj)) return(obj)
  if (is.data.frame(obj)) return(as.matrix(obj))

  # SummarizedExperiment
  if (inherits(obj, "SummarizedExperiment")) {
    if (length(assays(obj)) == 0) stop("RDS 中 SummarizedExperiment 无 assay")
    return(as.matrix(assay(obj)))
  }

  # list：优先找 raw counts / counts
  if (is.list(obj)) {
    cand_names <- names(obj)

    hit <- c(
      grep("raw[_]?counts|counts", cand_names, ignore.case = TRUE, value = TRUE),
      grep("expr|abundance|matrix|mat", cand_names, ignore.case = TRUE, value = TRUE)
    )
    hit <- unique(hit)

    for (nm in hit) {
      x <- obj[[nm]]
      if (is.matrix(x) || is.data.frame(x)) return(as.matrix(x))
      if (inherits(x, "SummarizedExperiment")) return(as.matrix(assay(x)))
    }
  }

  stop("无法从 RDS 提取表达矩阵: ", basename(path))
}


run_and_plot_deseq <- function(count_mat, sample_table, contrast_list, out_prefix, do_go = FALSE, species = "hg38", label_col = "feature_id", gene_anno = NULL, plot_color_col = NULL, heatmap_row_anno_col = NULL, te_label_level = NULL, te_color_level = NULL) {
  tool_outdir <- file.path(outdir, out_prefix)
  dir.create(tool_outdir, recursive = TRUE, showWarnings = FALSE)
  figure_dir <- file.path(tool_outdir, "figures")
  de_dir <- file.path(tool_outdir, "de_matrices")
  pathway_dir <- file.path(tool_outdir, "pathway")
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pathway_dir, recursive = TRUE, showWarnings = FALSE)

  safe_plot <- function(plot_tag, expr_fun, outfile) {
    tryCatch(
      expr_fun(),
      error = function(e) {
        message("[", plot_tag, "] ERROR: ", conditionMessage(e))
        writeLines(paste0("[", plot_tag, "] ERROR: ", conditionMessage(e)), paste0(outfile, ".txt"))
        save_placeholder_plot(outfile, plot_tag, paste("plot failed:", conditionMessage(e)))
        NULL
      }
    )
  }

  # 每个工具只画一次全局 PCA / Pearson
  safe_run(paste0(out_prefix, "_global_qc"), {
    x_all <- ensure_sample_order(count_mat, sample_table)
    count_mat_all <- as.matrix(x_all$mat)
    sample_table_all <- x_all$sample_table

    if (ncol(count_mat_all) < 2) {
      message("[", out_prefix, "_global_qc] too few samples, skip global PCA/Pearson")
    } else {
      libsize <- colSums(count_mat_all, na.rm = TRUE)
      libsize[libsize == 0] <- 1
      norm_mat_all <- t(t(count_mat_all) / libsize * 1e6)
      qc_mat_all <- log2(norm_mat_all + 1)

      safe_plot(
        paste0(out_prefix, "_PCA_global"),
        function() plot_pca_simple(
          qc_mat_all,
          sample_table_all,
          file.path(figure_dir, paste0(out_prefix, "_PCA_global.pdf")),
          paste0(out_prefix, "_global")
        ),
        file.path(figure_dir, paste0(out_prefix, "_PCA_global.pdf"))
      )

      safe_plot(
        paste0(out_prefix, "_Pearson_global"),
        function() plot_corr_simple(
          qc_mat_all,
          file.path(figure_dir, paste0(out_prefix, "_Pearson_global.pdf")),
          paste0(out_prefix, "_global")
        ),
        file.path(figure_dir, paste0(out_prefix, "_Pearson_global.pdf"))
      )
    }
  })

  worker_fun <- function(ct) {
    tag <- sanitize_name(paste(out_prefix, ct[2], "vs", ct[3], sep = "_"))
    gsea_df_out <- NULL

    safe_run(tag, {
      fit <- run_deseq_simple(count_mat, sample_table, ct)
      res <- fit$res

      this_label_col <- label_col
      this_color_col <- plot_color_col
      this_heatmap_row_anno_col <- heatmap_row_anno_col

      if (!is.null(gene_anno)) {
        if (all(c("gene_id", "gene_name") %in% colnames(gene_anno))) {
          res <- annotate_gene_res(res, gene_anno)
          if ("gene_name_plot" %in% colnames(res)) this_label_col <- "gene_name_plot"
        } else if (any(c("repName","repFamily","repClass","locus_id") %in% colnames(gene_anno))) {
          pref_cols <- unique(na.omit(c(
            te_label_level,
            te_color_level,
            "repName","repFamily","repClass","locus_id","gene_name","feature_id"
          )))
          res <- annotate_te_res(
            res, gene_anno,
            preferred_match_cols = pref_cols,
            label_level = te_label_level,
            color_level = te_color_level
          )
          if ("te_label_plot" %in% colnames(res)) this_label_col <- "te_label_plot"
          if ("te_color_plot" %in% colnames(res)) this_color_col <- "te_color_plot"
          if ("te_heatmap_group" %in% colnames(res)) this_heatmap_row_anno_col <- "te_heatmap_group"
        }
      }

      res_export <- clean_res_for_export(res, tag = tag)
      write.csv(res_export, file.path(de_dir, paste0(tag, "_DE.csv")), row.names = FALSE)

      volcano_metric <- plot_metric_for_res(res_export, "padj")
      volcano_pvalue_metric <- plot_metric_for_res(res_export, "pvalue")
      plot_tasks <- list(
        list(
          tag = paste0(tag, "_volcano"),
          outfile = file.path(figure_dir, paste0(tag, "_volcano.pdf")),
          fun = function(outfile) plot_volcano_simple(
            res_export, outfile, tag,
            label_col = this_label_col,
            top_n = label_top_n,
            color_col = this_color_col,
            sig_metric = volcano_metric
          )
        ),
        list(
          tag = paste0(tag, "_volcano_pvalue"),
          outfile = file.path(figure_dir, paste0(tag, "_volcano_pvalue.pdf")),
          fun = function(outfile) plot_volcano_simple(
            res_export, outfile, tag,
            label_col = this_label_col,
            top_n = label_top_n,
            color_col = this_color_col,
            sig_metric = volcano_pvalue_metric
          )
        ),
        list(
          tag = paste0(tag, "_MA"),
          outfile = file.path(figure_dir, paste0(tag, "_MA.pdf")),
          fun = function(outfile) plot_ma_simple(
            res_export, outfile, tag,
            label_col = this_label_col,
            top_n = label_top_n,
            color_col = this_color_col,
            sig_metric = volcano_metric
          )
        ),
        list(
          tag = paste0(tag, "_MA_pvalue"),
          outfile = file.path(figure_dir, paste0(tag, "_MA_pvalue.pdf")),
          fun = function(outfile) plot_ma_simple(
            res_export, outfile, tag,
            label_col = this_label_col,
            top_n = label_top_n,
            color_col = this_color_col,
            sig_metric = volcano_pvalue_metric
          )
        ),
        list(
          tag = paste0(tag, "_heatmap"),
          outfile = file.path(figure_dir, paste0(tag, "_heatmap.pdf")),
          fun = function(outfile) plot_top_heatmap_simple(
            fit$vst_mat,
            res_export,
            fit$sample_table,
            outfile,
            tag,
            top_n = heatmap_top_n,
            label_col = this_label_col,
            row_anno_col = this_heatmap_row_anno_col
          )
        )
      )
      run_one_plot <- function(task) safe_plot(task$tag, function() task$fun(task$outfile), task$outfile)
      if (length(contrast_list) == 1 && .Platform$OS.type == "unix" && plot_n_cores > 1) {
        parallel::mclapply(plot_tasks, run_one_plot, mc.cores = min(plot_n_cores, length(plot_tasks)), mc.preschedule = FALSE)
      } else {
        lapply(plot_tasks, run_one_plot)
      }

      if (do_go) {
        run_go_simple(res_export, file.path(pathway_dir, tag), species = species)

        if (!is.null(gene_anno) && all(c("gene_id", "gene_name") %in% colnames(gene_anno))) {
          gsea_df_out <<- run_gsea_simple(
            res_df = res_export,
            prefix = file.path(pathway_dir, tag),
            species = species,
            fit = fit,
            contrast = ct,
            gene_anno = gene_anno
          )
        } else {
          message("[", tag, "] skip GSEA: current matrix is not gene-level or gene annotation unavailable")
        }
      }
    })

    list(tag = tag, gsea_df = gsea_df_out)
  }

  res_list <- if (.Platform$OS.type == "unix" && length(contrast_list) > 1 && plot_n_cores > 1) {
    parallel::mclapply(
      contrast_list,
      worker_fun,
      mc.cores = min(plot_n_cores, length(contrast_list)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(contrast_list, worker_fun)
  }

  gsea_collect <- lapply(res_list, function(x) {
    if (is.null(x) || is.null(x$gsea_df) || !is.data.frame(x$gsea_df) || nrow(x$gsea_df) == 0) return(NULL)
    x$gsea_df
  })
  gsea_collect <- gsea_collect[!vapply(gsea_collect, is.null, logical(1))]

  if (length(gsea_collect) > 0) {
    gsea_all <- bind_rows(gsea_collect)
    write.csv(gsea_all, file.path(pathway_dir, paste0(out_prefix, "_GSEA_all_contrasts_all_categories.csv")), row.names = FALSE)
  }
}

get_symbol_fallback <- function(gene_ids, species = "hg38") {
  gene_ids <- unique(sub("\\..*$", "", as.character(gene_ids)))
  gene_ids <- gene_ids[!is.na(gene_ids) & gene_ids != ""]
  if (length(gene_ids) == 0) return(setNames(character(0), character(0)))
  species <- as.character(species)[1]

  db_obj <- if (grepl("^hg|human", species, ignore.case = TRUE)) org.Hs.eg.db else org.Mm.eg.db
  sym_map <- AnnotationDbi::mapIds(
    db_obj,
    keys = gene_ids,
    keytype = "ENSEMBL",
    column = "SYMBOL",
    multiVals = "first"
  )
  as.character(sym_map)
}

load_tx2gene <- function(tx2gene_path) {
  if (all(is.na(tx2gene_path)) || !file.exists(tx2gene_path)) {
    stop("tx2gene_path 不存在：", tx2gene_path)
  }

  # 允许直接给 GTF/GFF，而不是必须 2 列 tsv
  if (grepl("\\.(gtf|gff)(\\.gz)?$", tx2gene_path, ignore.case = TRUE)) {
    gr <- rtracklayer::import(tx2gene_path)
    tx2gene <- data.frame(
      tx   = as.character(mcols(gr)$transcript_id),
      gene = as.character(mcols(gr)$gene_id),
      type = as.character(mcols(gr)$type),
      stringsAsFactors = FALSE
    )
    tx2gene <- tx2gene[
      tx2gene$type == "transcript" &
      !is.na(tx2gene$tx) &
      !is.na(tx2gene$gene),
      c("tx", "gene"),
      drop = FALSE
    ]
  } else {
    tx2gene <- fread(tx2gene_path, data.table = FALSE)
    if (ncol(tx2gene) < 2) stop("tx2gene 至少需要2列")
    colnames(tx2gene)[1:2] <- c("tx", "gene")
    tx2gene <- tx2gene[, c("tx", "gene"), drop = FALSE]
  }

  # 去掉版本号，避免 ENSTxxxx.1 / ENSTxxxx 不匹配
  tx2gene$tx   <- sub("\\..*$", "", tx2gene$tx)
  tx2gene$gene <- sub("\\..*$", "", tx2gene$gene)

  unique(tx2gene)
}


load_salmon_header_gene_anno <- function(files) {
  if (length(files) == 0) return(data.frame(gene_id = character(0), gene_name = character(0), stringsAsFactors = FALSE))
  out_list <- lapply(files, function(f) {
    dt <- fread(f, select = "Name", data.table = FALSE)
    raw_name <- as.character(dt$Name)
    parts <- strsplit(raw_name, "\\|", perl = TRUE)
    tx <- vapply(parts, function(x) if (length(x) >= 1) x[1] else NA_character_, character(1))
    gene_id <- vapply(parts, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    gene_name <- vapply(parts, function(x) if (length(x) >= 3) x[3] else NA_character_, character(1))
    data.frame(
      tx = sub("\\..*$", "", tx),
      gene_id = sub("\\..*$", "", gene_id),
      gene_name = gene_name,
      stringsAsFactors = FALSE
    )
  })
  df <- dplyr::bind_rows(out_list)
  df$gene_name[df$gene_name == ""] <- NA_character_
  df <- df %>%
    dplyr::filter(!is.na(gene_id), gene_id != "") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      gene_name = {
        z <- unique(gene_name[!is.na(gene_name) & gene_name != "" & !grepl("^ENS", gene_name)])
        if (length(z) > 0) z[1] else NA_character_
      },
      .groups = "drop"
    )
  as.data.frame(df, stringsAsFactors = FALSE)
}

merge_gene_anno <- function(primary_df, extra_df = NULL) {
  if (is.null(extra_df) || nrow(extra_df) == 0) return(primary_df)
  x <- primary_df %>%
    dplyr::full_join(extra_df, by = "gene_id", suffix = c("", "_extra"))
  x$gene_name <- dplyr::coalesce(
    ifelse(!is.na(x$gene_name) & x$gene_name != "" & !grepl("^ENS", x$gene_name), x$gene_name, NA_character_),
    ifelse(!is.na(x$gene_name_extra) & x$gene_name_extra != "" & !grepl("^ENS", x$gene_name_extra), x$gene_name_extra, NA_character_),
    x$gene_name,
    x$gene_name_extra,
    x$gene_id
  )
  x <- x[, c("gene_id", "gene_name"), drop = FALSE]
  x <- x[!duplicated(x$gene_id), , drop = FALSE]
  x
}

load_gene_anno_from_gtf <- function(gtf_path, salmon_files = character(0), species = "hg38") {
  species <- as.character(species)[1]
  gr <- rtracklayer::import(gtf_path)
  get_mcol_chr <- function(nm) {
    x <- mcols(gr)[[nm]]
    if (is.null(x)) return(rep(NA_character_, length(gr)))
    as.character(x)
  }
  df <- data.frame(
    gene_id = get_mcol_chr("gene_id"),
    gene_name = get_mcol_chr("gene_name"),
    transcript_id = get_mcol_chr("transcript_id"),
    gene_type = dplyr::coalesce(
      get_mcol_chr("gene_type"),
      get_mcol_chr("gene_biotype"),
      get_mcol_chr("transcript_type"),
      get_mcol_chr("transcript_biotype")
    ),
    stringsAsFactors = FALSE
  )

  df$gene_id <- sub("\\..*$", "", df$gene_id)
  df$transcript_id <- sub("\\..*$", "", df$transcript_id)
  df$gene_name[df$gene_name == ""] <- NA_character_

  df <- df %>%
    dplyr::filter(!is.na(gene_id), gene_id != "") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarise(
      gene_name = {
        z <- unique(gene_name[!is.na(gene_name) & gene_name != "" & !grepl("^ENS", gene_name)])
        if (length(z) == 0) z <- unique(gene_name[!is.na(gene_name) & gene_name != ""])
        if (length(z) > 0) z[1] else NA_character_
      },
      gene_type = {
        z <- unique(gene_type[!is.na(gene_type) & gene_type != ""])
        if (length(z) > 0) z[1] else NA_character_
      },
      .groups = "drop"
    )

  extra_df <- load_salmon_header_gene_anno(salmon_files)
  if (nrow(extra_df) > 0) {
    message("[load_gene_anno_from_gtf] merge salmon header gene_name rows: ", nrow(extra_df))
    df <- merge_gene_anno(df, extra_df)
  }

  fallback_map <- suppressMessages(
    AnnotationDbi::mapIds(
      if (grepl("^hg|human", species, ignore.case = TRUE)) org.Hs.eg.db else org.Mm.eg.db,
      keys = df$gene_id,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
  )

  idx <- match(df$gene_id, names(fallback_map))
  miss <- is.na(df$gene_name) | df$gene_name == "" | grepl("^ENS", df$gene_name)
  df$gene_name[miss] <- ifelse(!is.na(idx[miss]) & !is.na(fallback_map[idx[miss]]) & fallback_map[idx[miss]] != "", fallback_map[idx[miss]], df$gene_name[miss])
  df$gene_name[is.na(df$gene_name) | df$gene_name == ""] <- df$gene_id[is.na(df$gene_name) | df$gene_name == ""]

  unresolved_n <- sum(grepl("^ENS", df$gene_name))
  message("[load_gene_anno_from_gtf] annotation rows: ", nrow(df), "; unresolved_symbol=", unresolved_n)
  df
}

annotate_gene_res <- function(res_df, gene_anno) {
  if (!"feature_id" %in% colnames(res_df)) {
    id_col <- intersect(c("gene_id", "Geneid", "id", "ID", "X"), colnames(res_df))[1]
    if (is.na(id_col) || !nzchar(id_col)) id_col <- colnames(res_df)[1]
    res_df$feature_id <- as.character(res_df[[id_col]])
  }
  res_df$feature_id <- as.character(res_df$feature_id)
  res_df$gene_id_clean <- sub("\\..*$", "", res_df$feature_id)

  anno_cols <- intersect(c("gene_id", "gene_name", "gene_type"), colnames(gene_anno))
  anno <- unique(gene_anno[, anno_cols, drop = FALSE])
  if ("gene_name" %in% colnames(anno)) colnames(anno)[colnames(anno) == "gene_name"] <- "gene_name_anno"
  if ("gene_type" %in% colnames(anno)) colnames(anno)[colnames(anno) == "gene_type"] <- "gene_type_anno"
  out <- res_df %>%
    left_join(anno, by = c("gene_id_clean" = "gene_id"))

  existing_gene_name <- if ("gene_name" %in% colnames(out)) as.character(out$gene_name) else rep(NA_character_, nrow(out))
  anno_gene_name <- if ("gene_name_anno" %in% colnames(out)) as.character(out$gene_name_anno) else rep(NA_character_, nrow(out))
  out$gene_name <- dplyr::coalesce(
    ifelse(!is.na(existing_gene_name) & existing_gene_name != "", existing_gene_name, NA_character_),
    ifelse(!is.na(anno_gene_name) & anno_gene_name != "", anno_gene_name, NA_character_),
    out$feature_id
  )
  if ("gene_type_anno" %in% colnames(out)) {
    existing_gene_type <- if ("gene_type" %in% colnames(out)) as.character(out$gene_type) else rep(NA_character_, nrow(out))
    out$gene_type <- dplyr::coalesce(
      ifelse(!is.na(existing_gene_type) & existing_gene_type != "", existing_gene_type, NA_character_),
      ifelse(!is.na(out$gene_type_anno) & out$gene_type_anno != "", as.character(out$gene_type_anno), NA_character_)
    )
  }

  existing_plot <- if ("gene_name_plot" %in% colnames(out)) as.character(out$gene_name_plot) else rep(NA_character_, nrow(out))
  out$gene_name_plot <- dplyr::coalesce(
    ifelse(!is.na(existing_plot) & existing_plot != "", existing_plot, NA_character_),
    ifelse(!is.na(out$gene_name) & out$gene_name != "", out$gene_name, NA_character_),
    out$feature_id
  )
  bad <- grepl("^OTTHUM(G|T|P)", out$gene_name_plot)
  out$gene_name_plot[bad] <- out$feature_id[bad]
  unresolved <- grepl("^ENS", out$gene_name_plot)
  if (any(unresolved, na.rm = TRUE)) {
    message("[annotate_gene_res] unresolved labels kept as ENS ids: ", sum(unresolved, na.rm = TRUE))
  }
  out
}

read_salmon_quant <- function(files, tx2gene_path = NA_character_) {
  tx2gene <- load_tx2gene(tx2gene_path)

  mats <- lapply(files, function(f) {
    dt <- fread(f, data.table = FALSE)
    need <- c("Name", "NumReads")
    if (!all(need %in% colnames(dt))) stop("Salmon quant.sf 缺少列: ", basename(f))

    raw_name <- as.character(dt$Name)
    tx_from_name <- sub("\\|.*$", "", raw_name)
    gene_from_name <- vapply(strsplit(raw_name, "\\|", fixed = FALSE), function(x) {
      hit <- grep("^ENS(G|MUSG)[0-9]+(\\.[0-9]+)?$", x, value = TRUE)
      if (length(hit) > 0) hit[1] else NA_character_
    }, character(1))

    tx_from_name <- sub("\\..*$", "", tx_from_name)
    gene_from_name <- sub("\\..*$", "", gene_from_name)

    dt$tx <- tx_from_name
    dt$gene_inline <- gene_from_name
    dt$count <- suppressWarnings(as.numeric(dt$NumReads))

    dt2 <- dt %>%
      dplyr::select(tx, gene_inline, count) %>%
      dplyr::left_join(tx2gene, by = "tx")

    dt2$gene_final <- ifelse(is.na(dt2$gene) | dt2$gene == "", dt2$gene_inline, dt2$gene)
    dt2 <- dt2 %>%
      dplyr::filter(!is.na(gene_final), gene_final != "") %>%
      dplyr::group_by(gene_final) %>%
      dplyr::summarise(count = sum(count, na.rm = TRUE), .groups = "drop")

    sample_name <- sub("_salmon$", "", basename(dirname(f)))
    colnames(dt2) <- c("feature_id", sample_name)
    dt2
  })

  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  finalize_feature_matrix(merged, tag = "Salmon")
}

# ==== Telescope reader ====
read_telescope_report <- function(files) {
  mats <- lapply(files, function(f) {
    dt <- tryCatch(
      utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE),
      error = function(e) { message("[read_telescope] skip ", basename(f), ": ", conditionMessage(e)); return(NULL) }
    )
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    need_cols <- c("transcript", "final_count")
    miss <- setdiff(need_cols, colnames(dt))
    if (length(miss) > 0) {
      message("[read_telescope] skip ", basename(f), ": missing columns ", paste(miss, collapse = ", "))
      return(NULL)
    }
    sample_name <- sub("[-_]telescope_report.tsv$", "", basename(f))
    out <- dt[, c("transcript", "final_count"), drop = FALSE]
    colnames(out) <- c("feature_id", sample_name)
    out <- out[!is.na(out[[sample_name]]) & out[[sample_name]] > 0, , drop = FALSE]
    out
  })
  mats <- Filter(Negate(is.null), mats)
  if (length(mats) == 0) stop("[read_telescope] 无有效 Telescope 报告")
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  num_cols <- setdiff(colnames(merged), "feature_id"); merged[num_cols][is.na(merged[num_cols])] <- 0; merged[["feature_id"]][is.na(merged[["feature_id"]])] <- "unknown"
  rownames(merged) <- merged[["feature_id"]]
  merged["feature_id"] <- NULL
  message("[read_telescope] matrix dim = ", paste(dim(merged), collapse = " x "))
  merged
}

parse_telescope_feature_map <- function(feature_ids, te_anno) {
  map_df <- match_annotation_table(feature_ids, te_anno, preferred_match_cols = c("repName", "gene_name", "locus_id"))
  message("[parse_telescope_feature_map] matched ", sum(!is.na(map_df)), " / ", nrow(map_df), " features to repName")
  map_df
}

read_telescope_conf_counts <- function(files) {
  mats <- lapply(files, function(f) {
    dt <- tryCatch(
      utils::read.table(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "#", check.names = FALSE),
      error = function(e) { message("[read_telescope_conf] skip ", basename(f), ": ", conditionMessage(e)); return(NULL) }
    )
    if (is.null(dt) || nrow(dt) == 0 || !"final_conf" %in% colnames(dt)) return(NULL)
    sample_name <- sub("[-_]telescope_report.tsv$", "", basename(f))
    out <- dt[, c("transcript", "final_conf"), drop = FALSE]
    colnames(out) <- c("feature_id", sample_name)
    out
  })
  mats <- Filter(Negate(is.null), mats)
  if (length(mats) == 0) return(NULL)
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), mats)
  num_cols <- setdiff(colnames(merged), "feature_id"); merged[num_cols][is.na(merged[num_cols])] <- 0; merged[["feature_id"]][is.na(merged[["feature_id"]])] <- "unknown"
  rownames(merged) <- merged[["feature_id"]]
  merged["feature_id"] <- NULL
  merged
}

message("[init] 通用函数加载完成")
