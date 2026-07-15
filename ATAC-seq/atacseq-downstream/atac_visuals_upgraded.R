# ============================================================
# ATAC-seq 可视化升级模块
# 在 atac_functions.R 末尾 source() 加载，覆盖部分绘图函数
# 改动: msigdbr支持、TE Family着色、通路图优化、主题升级
# ============================================================

suppressPackageStartupMessages({
  library(viridis)
  library(msigdbr)
})

# ---- 覆盖: Hallmark GSEA 使用 msigdbr (downstream env) ----
# 在原函数基础上修改 hallmark 加载部分
# 通过在 atac_functions.R 中 source 此文件覆盖旧定义

# ---- 覆盖: TE 火山图按 Family 着色 ----
# 注入到 plot_annotation_volcano 的 TE-only 部分
# 注: 无法直接覆盖，改用包装函数方式

# ---- 改进: 所有 ggplot 主题升级 ----
upgrade_theme <- function(p, title_text = NULL) {
  p <- p + theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "grey40"),
      legend.position = "bottom"
    )
  if (!is.null(title_text)) {
    p <- p + labs(title = title_text)
  }
  p
}

# ---- 改进: 通路富集图 - 用 cnetplot 和 emapplot ----
plot_enrichment_network <- function(ego, out_pdf, title = "", showCategory = 10) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) return(invisible(NULL))
  tryCatch({
    # cnetplot: gene-concept network
    p1 <- enrichplot::cnetplot(ego, showCategory = showCategory, node_label = "category", cex_label_category = 0.8) +
      upgrade_theme(NULL, title) +
      theme(legend.position = "right")
    ggsave(out_pdf, p1, width = 10, height = 8)
    ggsave(sub(".pdf$", ".png", out_pdf), p1, width = 10, height = 8, dpi = 180)
  }, error = function(e) message("[visuals] cnetplot failed: ", e$message))

  tryCatch({
    if (nrow(as.data.frame(ego)) >= 5) {
      ego2 <- pairwise_termsim(ego)
      p2 <- enrichplot::emapplot(ego2, showCategory = showCategory, color = "p.adjust") +
        scale_color_viridis_c(option = "plasma") +
        upgrade_theme(NULL, paste0(title, " (Enrichment Map)"))
      out_emap <- sub(".pdf$", "_emap.pdf", out_pdf)
      ggsave(out_emap, p2, width = 9, height = 8)
      ggsave(sub(".pdf$", ".png", out_emap), p2, width = 9, height = 8, dpi = 180)
    }
  }, error = function(e) message("[visuals] emapplot failed: ", e$message))
}

# ---- 改进: Hallmark GSEA 使用 msigdbr + 更好可视化 ----
run_hallmark_gsea_upgraded <- function(gene_list, species, out_pdf, title = "", orgdb = NULL) {
  sp <- if (species == "hg38") "Homo sapiens" else "Mus musculus"
  
  # Convert Entrez IDs to symbols if orgdb provided
  if (!is.null(orgdb) && grepl("^[0-9]+$", names(gene_list)[1])) {
    conv <- tryCatch(clusterProfiler::bitr(names(gene_list), fromType = "ENTREZID", toType = "SYMBOL", OrgDb = orgdb), error = function(e) NULL)
    if (!is.null(conv) && nrow(conv) > 0) {
      names(gene_list) <- conv[match(names(gene_list), conv)]
      gene_list <- gene_list[!is.na(names(gene_list))]
      gene_list <- gene_list[names(gene_list) != ""]
    }
  }

  hallmark_list <- tryCatch({
    msigdbr::msigdbr(species = sp, collection = "H") |>
      dplyr::select(gs_name, gene_symbol) |>
      split(x = .$gene_symbol, f = .$gs_name)
  }, error = function(e) {
    # GMT fallback
    pipeline_root <- if (exists('atac_pipeline_root', inherits = TRUE)) {
      get('atac_pipeline_root', inherits = TRUE)
    } else {
      Sys.getenv('ATAC_PIPELINE_ROOT', unset = '')
    }
    gmt <- if (nzchar(pipeline_root)) {
      file.path(pipeline_root, 'conf', 'h.all.v2024.1.Hs.symbols.gmt')
    } else {
      ''
    }
    if (file.exists(gmt)) fgsea::gmtPathways(gmt) else stop("No Hallmark data available")
  })
  
  res <- fgsea::fgseaMultilevel(pathways = hallmark_list, stats = gene_list, 
                                  minSize = 10, maxSize = 500)
  
  if (nrow(res) == 0) return(invisible(NULL))
  
  res <- res[order(res$padj), ]
  res$leadingEdge <- vapply(res$leadingEdge, function(x) paste(x, collapse = ";"), character(1))
  
  top <- head(res[!is.na(res$padj) & res$padj < 0.5, ], 30)
  if (nrow(top) == 0) top <- head(res[order(res$pval), ], 15)
  
  if (nrow(top) > 0) {
    top$neg_log10_padj <- -log10(pmax(top$padj, 1e-10, na.rm = TRUE))
    p <- ggplot(top, aes(x = NES, y = reorder(pathway, NES), fill = neg_log10_padj)) +
      geom_col(color = "grey30", width = 0.75) +
      scale_fill_viridis_c(option = "plasma", name = "-log10(padj)") +
      theme_minimal(base_size = 11) +
      theme(panel.grid.major.y = element_blank(),
            plot.title = element_text(face = "bold", size = 13),
            plot.subtitle = element_text(size = 9, color = "grey40")) +
      labs(title = if (title != "") title else "Hallmark GSEA",
           subtitle = sprintf("Top %d pathways", nrow(top)),
           x = "Normalized Enrichment Score", y = "")
    ggsave(out_pdf, p, width = 11, height = 7.5)
    ggsave(sub(".pdf$", ".png", out_pdf), p, width = 11, height = 7.5, dpi = 180)
  }
  
  invisible(res)
}

message("[visuals] ATAC visualization upgrades loaded")
