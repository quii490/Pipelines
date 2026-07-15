# ==============================================================================
#                      模块3: 生成所有可视化图表
# ==============================================================================
message("\n--- [模块3] 正在执行：生成所有可视化图表 ---")
message("--> [运行] 无论是否存在缓存，此模块总是重新生成图表。")

# 加载所有库
#suppressPackageStartupMessages({
#  library(tidyverse); library(ggplot2); library(ggrepel)
#  library(viridis); library(ggpubr); library(circlize)
#})

suppressPackageStartupMessages({
  library(tidyverse); library(ggplot2); library(ggrepel)
  library(viridis); library(circlize); library(patchwork)
})

# --- 加载上一步的结果 ---
plot_data <- readRDS(cache_files$plot_data)

# Keep analytical tables complete, but cap only display-layer data.  TE count
# matrices commonly contain millions of copies; sending all rows through
# repeated ggplot/loess operations makes RUN_ANALYSIS unnecessarily slow.
MAX_TE_PLOT_ROWS <- if (exists("MAX_TE_PLOT_ROWS")) MAX_TE_PLOT_ROWS else 200000
MAX_TE_SCATTER_POINTS <- if (exists("MAX_TE_SCATTER_POINTS")) MAX_TE_SCATTER_POINTS else 100000
MAX_TE_FACET_LEVELS <- if (exists("MAX_TE_FACET_LEVELS")) MAX_TE_FACET_LEVELS else 100
sample_for_plot <- function(df, n, seed = 20240713L) {
  if (nrow(df) <= n) return(df)
  set.seed(seed)
  dplyr::slice_sample(df, n = n)
}

# --- 创建输出目录 ---
#dir.create(file.path(figures_dir, "Correlations"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_Family_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_repName_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Correlation_Gene_Level"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Boxplots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(figures_dir, "Heatmaps"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# ## 第 6 节：定义核心绘图函数 ##
# ==============================================================================
message("=== 第 6 节：定义核心绘图函数 `generate_correlation_plot` ===")

#' @title 生成富集相关性散点图
#' @description 比较两个样本在基因或TE上的富集程度，并标注高差异目标
#' @param ... (参数与前一版相同)
#' @return 返回一个包含图中被标注数据点信息的数据框
generate_correlation_plot <- function(data, x_name, y_name, feature_type, te_classes = NULL,
                                      grouping_vars = NULL, label_var = NULL,
                                      n_labels, output_dir_path) {
  
  message(paste("    -> 正在绘制:", x_name, "vs", y_name, "| 类型:", feature_type))
  
  x_col <- sym(x_name)
  y_col <- sym(y_name)
  
  if (!feature_type %in% unique(data$featureType)) {
    message("  [SKIP] No ", feature_type, " data available")
    return(NULL)
  }
  base_data <- data %>% filter(featureType == feature_type, sample %in% c(x_name, y_name))
  
  # --- TE 和 Gene 的数据处理逻辑分支 ---
  if (feature_type == "TE") {
    plot_df_base <- base_data %>%
      group_by(across(all_of(grouping_vars))) %>%
      summarise(
        !!x_col := mean(lfc_over_igg[sample == x_name], na.rm = TRUE),
        !!y_col := mean(lfc_over_igg[sample == y_name], na.rm = TRUE),
        n_elements = n() / 2, mean_div = mean(milliDiv, na.rm = TRUE), .groups = "drop"
      ) 
      #%>%
      #filter(Class %in% te_classes, n_elements > 10)
    if (!identical(te_classes, "all")) {
      plot_df_base <- plot_df_base %>% filter(Class %in% te_classes)
    }
    
    plot_df_base <- plot_df_base %>% filter(n_elements > 10)
    
    if(nrow(plot_df_base) < n_labels) {
      message(paste("       警告: TE数据点不足 (", nrow(plot_df_base), ")，跳过绘图。"))
      return(NULL)
    }
    
    plot_df <- plot_df_base %>%
      mutate(
        age_group = case_when(
          is.na(mean_div) ~ "unknown",
          mean_div < 50 ~ "young (<5%)",
          mean_div < 150 ~ "middle (5-15%)",
          TRUE ~ "old (>=15%)"
        ),
        age_group = factor(age_group, levels = c("young (<5%)", "middle (5-15%)", "old (>=15%)", "unknown")),
        point_size = log10(pmax(n_elements, 1))
      )
    point_aes <- aes(color = Family, size = point_size, shape = age_group)
    size_lab <- "log10(TE element count)"; shape_lab <- "TE age (milliDiv)"; color_lab <- "TE family"
    
  } else { # Gene
    plot_df <- base_data %>%
      select(sample, repName, Family, lfc_over_igg) %>%
      pivot_wider(names_from = sample, values_from = lfc_over_igg, values_fn = mean) %>%
      filter(!is.na(.data[[x_name]]) & !is.na(.data[[y_name]]))
    
    if(nrow(plot_df) < n_labels) {
      message(paste("       警告: Gene数据点不足 (", nrow(plot_df), ")，跳过绘图。"))
      return(NULL)
    }
    point_aes <- aes(color = Family)
    label_var <- "repName"
    size_lab <- NULL; shape_lab <- NULL; color_lab <- "gene biotype"
  }
  
  # --- 筛选需要标注的数据点 (离原点最远的点) ---
  label_data <- plot_df %>%
    mutate(dist_from_origin = sqrt((!!x_col)^2 + (!!y_col)^2)) %>%
    slice_max(order_by = dist_from_origin, n = n_labels)
  
  pearson_value <- cor(plot_df[[x_name]], plot_df[[y_name]], use = "complete.obs", method = "pearson")
  spearman_value <- cor(plot_df[[x_name]], plot_df[[y_name]], use = "complete.obs", method = "spearman")
  p_value <- tryCatch(
    cor.test(plot_df[[x_name]], plot_df[[y_name]], method = "spearman", exact = FALSE)$p.value,
    error = function(e) NA_real_
  )
  label_text <- paste0(
    "Pearson r = ", round(pearson_value, 3),
    "\nSpearman ρ = ", round(spearman_value, 3),
    "\nP = ", signif(p_value, 3)
  )

  # --- 核心绘图代码 ---
  p <- ggplot(data = plot_df, aes(x = !!x_col, y = !!y_col)) +
    geom_point(point_aes, alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, linetype = "longdash", color = "grey65", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    #stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 5) +
    annotate("text",x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1, label = label_text, size = 5) +
    geom_text_repel(
      data = label_data, aes(label = .data[[label_var]]),
      size = 3.5, box.padding = 0.5, min.segment.length = 0,
      segment.color = 'grey50', max.overlaps = Inf, force = 3
    ) +
    scale_color_viridis_d(option = "D", end = 0.9) +
    theme_bw(base_size = 14) +
    labs(
      title = paste(feature_type, "enrichment correlation:", x_name, "vs", y_name),
      subtitle = if(feature_type == "TE") paste("Points are grouped by", label_var) else "Each point represents one gene",
      x = paste(x_name, " enrichment (log2FC vs IgG)"),
      y = paste(y_name, " enrichment (log2FC vs IgG)"),
      size = size_lab, shape = shape_lab, color = color_lab
    )
  
  if (feature_type == "TE") { p <- p + scale_size_continuous(range = c(0.5, 4.5), name = "log10(TE element count)") }
  
  # --- 保存图像 ---
  file_base <- paste0("Correlation_", feature_type, "_", x_name, "_vs_", y_name)
  ggsave(file.path(output_dir_path, paste0(file_base, ".pdf")), p, width = 12, height = 10, bg = "white")
  ggsave(file.path(output_dir_path, paste0(file_base, ".png")), p, width = 12, height = 10, dpi = 300, bg = "white")
  
  return(label_data)
  #if (for_matrix) {
  #  p <- p +
  #    theme(
  #      legend.position = "none",
  #      plot.title = element_text(size = 10),
  #      axis.title = element_blank()
  #    )
  #}
  #return(list(data = label_data, plot = p))
  
}





generate_correlation_matrix_plot <- function(data, feature_type,
                                             te_classes=NULL, grouping_vars=NULL,
                                             label_var=NULL, n_labels=10,
                                             output_file, matrix_samples=NULL) {
  
  message("生成整体矩阵图...")
  
  if (is.null(matrix_samples)) matrix_samples <- unique(data$sample)
  sample_pairs <- combn(matrix_samples, 2, simplify = FALSE)
  
  all_plot_data <- list()
  all_label_data <- list()
  
  for (pair in sample_pairs) {
    x_name <- pair[1]; y_name <- pair[2]
    x_col <- sym(x_name); y_col <- sym(y_name)
    
    if (!feature_type %in% unique(data$featureType)) {
      message("  [SKIP] No ", feature_type, " data available")
      return(NULL)
    }
    base_data <- data %>% filter(featureType == feature_type, sample %in% c(x_name, y_name))
    
    if (feature_type=="TE") {
      plot_df <- base_data %>%
        group_by(across(all_of(grouping_vars))) %>%
        summarise(
          x_val = mean(lfc_over_igg[sample==x_name], na.rm=TRUE),
          y_val = mean(lfc_over_igg[sample==y_name], na.rm=TRUE),
          n_elements = n()/2,
          mean_div = mean(milliDiv, na.rm=TRUE),
          .groups="drop"
        )
      if (!identical(te_classes,"all")) plot_df <- plot_df %>% filter(Class %in% te_classes)
      plot_df <- plot_df %>% filter(n_elements>10)
      
      if (nrow(plot_df)<1) next
      
      plot_df <- plot_df %>%
        mutate(
          age_group = case_when(
            is.na(mean_div) ~ "unknown",
            mean_div < 50 ~ "young (<5%)",
            mean_div < 150 ~ "middle (5-15%)",
            TRUE ~ "old (>=15%)"
          ),
          age_group = factor(age_group, levels = c("young (<5%)", "middle (5-15%)", "old (>=15%)", "unknown")),
          point_size = log10(pmax(n_elements, 1))
        )
    } else { # Gene 
      plot_df <- base_data %>%
      select(sample, repName, Family, lfc_over_igg) %>%
      pivot_wider(names_from = sample, values_from = lfc_over_igg, values_fn = mean) %>%
      filter(!is.na(.data[[x_name]]) & !is.na(.data[[y_name]])) %>%
      rename_with(~c("x_val","y_val"), .cols = c(x_name, y_name))
      if (nrow(plot_df)<1) next
    }
    
    
    plot_df <- plot_df %>% mutate(pair=paste(x_name, y_name, sep="_"))
    
    label_data <- plot_df %>%
      mutate(dist_from_origin=sqrt(x_val^2+y_val^2)) %>%
      slice_max(order_by=dist_from_origin, n=n_labels) %>%
      mutate(pair=paste(x_name,y_name,sep="_"))
    
    all_plot_data[[length(all_plot_data)+1]] <- plot_df
    all_label_data[[length(all_label_data)+1]] <- label_data
  }
  
  combined_data <- bind_rows(all_plot_data)
  combined_labels <- bind_rows(all_label_data)
  
  # 如果没有有效数据，直接退出
  if (nrow(combined_data)==0) {
    message("无有效数据，跳过绘图: ", output_file)
    return(NULL)
  }
  
  # 为每行生成 row/col
  combined_data <- combined_data %>%
    rowwise() %>%
    mutate(row = factor(strsplit(pair,"_")[[1]][1], levels=matrix_samples),
           col = factor(strsplit(pair,"_")[[1]][2], levels=matrix_samples)) %>%
    ungroup()
  
  combined_labels <- combined_labels %>%
    rowwise() %>%
    mutate(row = factor(strsplit(pair,"_")[[1]][1], levels=matrix_samples),
           col = factor(strsplit(pair,"_")[[1]][2], levels=matrix_samples)) %>%
    ungroup()
  
  # --- 设置 aes ---
  if (feature_type=="TE") {
    geom_point_aes <- aes(color=Family, size=point_size, shape=age_group)
  } else {
    geom_point_aes <- aes(color=Family)
  }
  
  p <- ggplot(combined_data, aes(x=x_val, y=y_val)) +
    geom_point(geom_point_aes, alpha=0.7) +
    geom_vline(xintercept=0, linetype="dotted", color="grey50") +
    geom_hline(yintercept=0, linetype="dotted", color="grey50") +
    geom_text_repel(data=combined_labels,
                    aes(x=x_val, y=y_val, label=.data[[label_var]]), size=2.5, max.overlaps=Inf) +
    facet_grid(row~col, scales="fixed") +
    scale_color_viridis_d(option="D", end=0.9) +
    {if (feature_type == "TE") scale_size_continuous(range = c(0.5, 4), name = "log10(TE element count)") else NULL} +
    theme_bw(base_size=14) +
    theme(strip.text=element_text(size=8),
          axis.title=element_blank(),
          legend.position="bottom")
  
  output_base <- tools::file_path_sans_ext(output_file)
  ggsave(paste0(output_base, ".pdf"), p, width=20, height=20, bg="white")
  ggsave(paste0(output_base, ".png"), p, width=20, height=20, dpi=300, bg="white")
  message("矩阵图保存完成: ", output_base, ".pdf / .png")
}





message("--- 绘图函数定义完成 ---\n")


# ==============================================================================
# ## 第 7 节：执行批量相关性分析 ##
# ==============================================================================
message("=== 第 7 节：开始批量生成相关性分析图表 ===")

sample_list <- unique(plot_data$sample)
sample_pairs <- combn(sample_list, 2, simplify = FALSE)
all_significant_hits <- list()
all_plots_te_family <- list()

for (pair in sample_pairs) {
  x_name <- pair[1]
  y_name <- pair[2]
  
  # --- a. TE 家族 (Family) 水平 ---
  family_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "TE",
    te_classes = TE_CLASSES_OF_INTEREST, grouping_vars = c("Family", "Class"),
    label_var = "Family", n_labels = N_LABELS_TE_FAMILY,
    output_dir_path = file.path(figures_dir, "Correlation_Family_Level")
  )
  if (!is.null(family_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      family_hits %>% mutate(feature_type = "TE", analysis_level = "Family", comparison = paste(x_name, "vs", y_name))
    ))
  }
  
  # --- b. TE 元件 (repName) 水平 ---
  repname_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "TE",
    te_classes = TE_CLASSES_OF_INTEREST, grouping_vars = c("repName", "Family", "Class"),
    label_var = "repName", n_labels = N_LABELS_TE_REPNAME,
    output_dir_path = file.path(figures_dir, "Correlation_repName_Level")
  )
  if (!is.null(repname_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      repname_hits %>% mutate(feature_type = "TE", analysis_level = "repName", comparison = paste(x_name, "vs", y_name))
    ))
  }
  
  # --- c. 基因 (Gene) 水平 ---
  gene_hits <- generate_correlation_plot(
    data = plot_data, x_name = x_name, y_name = y_name, feature_type = "Gene",
    n_labels = N_LABELS_GENES,
    output_dir_path = file.path(figures_dir, "Correlation_Gene_Level")
  )
  if (!is.null(gene_hits)) {
    all_significant_hits <- append(all_significant_hits, list(
      gene_hits %>% mutate(feature_type = "Gene", analysis_level = "Gene", comparison = paste(x_name, "vs", y_name))
    ))
  }
}

# --- 汇总并保存所有图中被标注的"显著"目标 ---
if (length(all_significant_hits) > 0) {
  final_summary_df <- bind_rows(all_significant_hits)
  write_csv(final_summary_df, file.path(figures_dir, "significant_hits_summary.csv"))
  message("--- 所有相关性图表生成完毕，并已保存 `significant_hits_summary.csv` ---\n")
} else {
  message("--- 未能生成任何相关性图表或显著目标 ---\n")
}


# 调用函数
generate_correlation_matrix_plot(
  data = plot_data,              
  feature_type = "TE",          
  te_classes = TE_CLASSES_OF_INTEREST, 
  grouping_vars = c("Family", "Class"), 
  label_var = "Family",        
  n_labels = N_LABELS_TE_FAMILY,
  output_file =  file.path(figures_dir, "Correlation_Family_Level/TE_Family_Correlation_Matrix.pdf"),    
  matrix_samples = unique(plot_data$sample) 
)

generate_correlation_matrix_plot(
  data = plot_data,
  feature_type = "TE",
  te_classes = TE_CLASSES_OF_INTEREST,
  grouping_vars = c("repName", "Family", "Class"),
  label_var = "repName",
  n_labels = N_LABELS_TE_REPNAME,
  output_file = file.path(figures_dir, "Correlation_repName_Level/TE_repName_Correlation_Matrix.pdf"),
  matrix_samples = unique(plot_data$sample)
)

generate_correlation_matrix_plot(
  data = plot_data,
  feature_type = "Gene",
  grouping_vars = c("repName"),
  label_var = "repName",
  n_labels = N_LABELS_GENES,
  output_file = file.path(figures_dir, "Correlation_Gene_Level/Gene_Correlation_Matrix.pdf"),
  matrix_samples = unique(plot_data$sample)
)



# ==============================================================================
# ## 第 8 节：特定目标筛选与箱线图可视化 ##
# ==============================================================================
#message("=== 第 8 节：绘制箱线图 ===")
message("=== 第 8 节：TE & Gene Heatmap 和 Boxplot 可视化 ===")

if(!exists("TE_repname_OI")) TE_repname_OI <- NULL

# --- TE 数据 ---
te_data_full <- plot_data %>%
  filter(featureType=="TE", Class %in% TE_CLASSES_OF_INTEREST)
te_data <- sample_for_plot(te_data_full, MAX_TE_PLOT_ROWS)
if(!is.null(TE_repname_OI) && length(TE_repname_OI)>0 && any(TE_repname_OI != "")){
  te_data_rep <- te_data %>% filter(repName %in% TE_repname_OI)
} else {
  # matrixStats (attached by DESeq2) also exports count(); qualify dplyr here
  # so the list-column error cannot recur in a clean pipeline process.
  top_repnames <- te_data %>% dplyr::count(repName, sort = TRUE) %>% slice_head(n = MAX_TE_FACET_LEVELS) %>% pull(repName)
  te_data_rep <- te_data %>% filter(repName %in% top_repnames)
}
message(TE_repname_OI)

# --- Gene 数据 ---
gene_data <- plot_data %>%
  filter(featureType=="Gene")

# --- 函数：计算相关性矩阵 ---
compute_cor_matrix <- function(df, value_col="lfc_over_igg", feature_key=NULL){
  samples <- unique(df$sample)
  if (length(samples) == 0) return(tibble())
  if (is.null(feature_key)) {
    feature_key <- if ("geneID" %in% colnames(df)) "geneID" else if ("repName" %in% colnames(df)) "repName" else "Family"
  }
  cor_wide <- df %>%
    select(all_of(feature_key), sample, value = all_of(value_col)) %>%
    group_by(across(all_of(c(feature_key, "sample")))) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = value)
  cor_mat <- cor(
    as.matrix(cor_wide %>% select(any_of(samples))),
    use = "pairwise.complete.obs",
    method = "spearman"
  )
  as.data.frame(as.table(cor_mat)) %>%
    setNames(c("Sample1","Sample2","Correlation"))
}

# ================= TE 图 =================
if(nrow(te_data)>0){

  samples <- unique(te_data$sample)

  # --- 总体 Heatmap ---
  cor_long <- compute_cor_matrix(te_data)
  p_heatmap <- ggplot(cor_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=3) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    theme_minimal() +
    labs(title="TE enrichment correlation", subtitle="Spearman correlation across TE copies") + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
  ggsave(file.path(figures_dir,"Heatmaps","TE_All_Heatmap.pdf"), p_heatmap, width=8, height=6)

  # --- 总体 Boxplot ---
  p_box <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    theme_bw(base_size=14) +
    labs(title="TE log2FC distribution across samples", x="Sample", y="log2FC vs IgG") +
    theme(
      legend.position="none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
      strip.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9)
    )
  ggsave(file.path(figures_dir,"Boxplots","TE_All_Boxplot.pdf"), p_box, width=8, height=6)

  # --- 按 Family (Class facet) Boxplot ---
  p_box_family <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~Family, scales="free_y") +
    theme_bw(base_size=12) +
    labs(title="TE log2FC by Family", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none", 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_Family.pdf"), p_box_family, width=12, height=8)
  
  p_box_class <- ggplot(te_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~Class, scales="free_y") +
    theme_bw(base_size=12) +
    labs(title="TE log2FC by Class", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none", 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_Class.pdf"), p_box_class, width=12, height=8)

  # --- 按 repName Boxplot ---
  p_box_repname <- ggplot(te_data_rep, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    facet_wrap(~repName, scales="free_y") +
    theme_bw(base_size=10) +
    labs(title="TE log2FC by repName", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none")
  ggsave(file.path(figures_dir,"Boxplots","TE_Boxplot_by_repName.pdf"), p_box_repname, width=16, height=12)

  # Aggregate once to family/sample before class/family heatmaps.  The old
  # implementation repeatedly regrouped millions of TE-copy rows inside each
  # facet loop.
  te_family <- te_data %>%
    group_by(Class, Family, sample) %>%
    summarise(value = mean(lfc_over_igg, na.rm = TRUE), .groups = "drop")

  # --- 按 Class facet Heatmap ---
  cor_class_list <- list()
  for(cl in unique(te_family$Class)){
    df_cl <- te_family %>% filter(Class==cl)
    cor_cl <- compute_cor_matrix(df_cl, value_col = "value", feature_key = "Family")
    cor_cl$Class <- cl
    cor_class_list[[length(cor_class_list)+1]] <- cor_cl
  }
  cor_class_long <- bind_rows(cor_class_list)
  p_heatmap_class <- ggplot(cor_class_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=2.5) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    facet_wrap(~Class) +
    theme_minimal() +
    labs(title="TE Correlation Heatmap by Class") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Heatmaps","TE_Heatmap_by_Class.pdf"), p_heatmap_class, width=12, height=10)
  
  cor_family_list <- list()
  for(fa in unique(te_family$Family)){
    df_fa <- te_family %>% filter(Family == fa)
    cor_fa <- compute_cor_matrix(df_fa, value_col = "value", feature_key = "Family")
    cor_fa$Family <- fa
    cor_family_list[[length(cor_family_list)+1]] <- cor_fa
  }
  cor_family_long <- bind_rows(cor_family_list)
  p_heatmap_family <- ggplot(cor_family_long, aes(x = Sample1, y = Sample2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), size = 2.5) +
    scale_fill_viridis_c(option = "C", limits = c(-1, 1)) +
    facet_wrap(~Family) +
    theme_minimal() +
    labs(title = "TE Correlation Heatmap by Family") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir, "Heatmaps", "TE_Heatmap_by_Family.pdf"), p_heatmap_family, width = 12, height = 10)

  # --- 按 repName facet Heatmap ---
  te_rep <- te_data %>%
    group_by(Class, Family, repName, sample) %>%
    summarise(value = mean(lfc_over_igg, na.rm = TRUE), .groups = "drop")
  cor_rep_list <- list()
  for(rep in unique(te_data_rep$repName)){
    df_rep <- te_rep %>% filter(repName==rep)
    cor_rep <- compute_cor_matrix(df_rep, value_col = "value", feature_key = "repName")
    cor_rep$repName <- rep
    cor_rep_list[[length(cor_rep_list)+1]] <- cor_rep
  }
  cor_rep_long <- bind_rows(cor_rep_list)
  p_heatmap_rep <- ggplot(cor_rep_long, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=2.5) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    facet_wrap(~repName) +
    theme_minimal() +
    labs(title="TE Correlation Heatmap by repName") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Heatmaps","TE_Heatmap_by_repName.pdf"), p_heatmap_rep, width=16, height=12)

  message("--- TE Heatmap 和 Boxplot 已完成 ---")
}

# ================= Gene 图 =================
if(nrow(gene_data)>0){

  samples <- unique(gene_data$sample)

  # --- Gene 总体 Heatmap ---
  cor_long_gene <- compute_cor_matrix(gene_data)
  p_heatmap_gene <- ggplot(cor_long_gene, aes(x=Sample1, y=Sample2, fill=Correlation)) +
    geom_tile(color="white") +
    geom_text(aes(label=round(Correlation,2)), size=3) +
    scale_fill_viridis_c(option="C", limits=c(-1,1)) +
    theme_minimal() +
    labs(title="Gene enrichment correlation", subtitle="Spearman correlation across genes") +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Heatmaps","Gene_All_Heatmap.pdf"), p_heatmap_gene, width=8, height=6)

  # --- Gene 总体 Boxplot ---
  p_box_gene <- ggplot(gene_data, aes(x=sample, y=lfc_over_igg, fill=sample)) +
    geom_boxplot() +
    theme_bw(base_size=14) +
    labs(title="Gene log2FC distribution across samples", x="Sample", y="log2FC vs IgG") +
    theme(legend.position="none", 
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8)
    )
  ggsave(file.path(figures_dir,"Boxplots","Gene_All_Boxplot.pdf"), p_box_gene, width=8, height=6)

  message("--- Gene Heatmap 和 Boxplot 已完成 ---")
}


message("--- 第 8 节分析完成 ---\n")

# ==============================================================================
# ## 第 9 节：增强 TE 可视化
#
# These plots retain count_draw's copy/family/subfamily data model while adding
# publication-style overviews. They are exploratory summaries, not tests.
# ==============================================================================
message("=== 第 9 节：增强 TE 可视化 ===")
enhanced_te_dir <- file.path(figures_dir, "Enhanced_TE")
dir.create(enhanced_te_dir, recursive = TRUE, showWarnings = FALSE)

save_dual <- function(plot, stem, width, height) {
  ggsave(file.path(enhanced_te_dir, paste0(stem, ".pdf")), plot, width = width, height = height, bg = "white")
  ggsave(file.path(enhanced_te_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
}

te_enhanced <- plot_data %>%
  filter(featureType == "TE", !is.na(Family), !is.na(sample))

if (nrow(te_enhanced) > 0) {
  # 9a. Sample-level PCA from TE subfamily enrichment.
  pca_wide <- te_enhanced %>%
    group_by(repName, sample) %>%
    summarise(value = mean(lfc_over_igg, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = sample, values_from = value)
  if (ncol(pca_wide) >= 4 && nrow(pca_wide) >= 3) {
    pca_mat <- pca_wide %>% select(-repName) %>% as.data.frame()
    pca_mat <- pca_mat[, vapply(pca_mat, function(x) sum(is.finite(x)) >= 3, logical(1)), drop = FALSE]
    pca_mat <- pca_mat[, vapply(pca_mat, function(x) sd(x, na.rm = TRUE) > 0, logical(1)), drop = FALSE]
    if (ncol(pca_mat) >= 2) {
      pca_mat[!is.finite(as.matrix(pca_mat))] <- 0
      pca_fit <- prcomp(t(pca_mat), center = TRUE, scale. = TRUE)
      pca_df <- as.data.frame(pca_fit$x) %>% rownames_to_column("sample")
      pca_df$sample <- factor(pca_df$sample, levels = unique(te_enhanced$sample))
      variance <- (pca_fit$sdev^2) / sum(pca_fit$sdev^2)
      p_pca <- ggplot(pca_df, aes(PC1, PC2, color = sample, label = sample)) +
        geom_hline(yintercept = 0, color = "grey85") +
        geom_vline(xintercept = 0, color = "grey85") +
        geom_point(size = 3.5) +
        ggrepel::geom_text_repel(show.legend = FALSE, max.overlaps = Inf) +
        scale_color_viridis_d(option = "D", end = 0.9) +
        theme_bw(base_size = 13) +
        labs(title = "TE enrichment PCA", subtitle = "Exploratory PCA of subfamily-level log2 enrichment",
             x = paste0("PC1 (", round(100 * variance[1], 1), "%)"),
             y = paste0("PC2 (", round(100 * variance[2], 1), "%)"), color = "Sample")
      save_dual(p_pca, "TE_PCA", 8, 6)
      write_csv(pca_df, file.path(enhanced_te_dir, "TE_PCA_coordinates.csv"))
    }
  }

  # 9b. Family dot plot: effect size, copy support and class.
  family_stats <- te_enhanced %>%
    group_by(sample, Class, Family) %>%
    summarise(mean_lfc = mean(lfc_over_igg, na.rm = TRUE),
              median_lfc = median(lfc_over_igg, na.rm = TRUE),
              n_elements = n(), mean_milliDiv = mean(milliDiv, na.rm = TRUE),
              .groups = "drop") %>%
    mutate(abs_lfc = abs(mean_lfc))
  dot_data <- family_stats %>%
    group_by(Class) %>%
    slice_max(abs_lfc, n = 30, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Family = forcats::fct_reorder(Family, mean_lfc))
  if (nrow(dot_data) > 0) {
    p_dot <- ggplot(dot_data, aes(x = mean_lfc, y = Family, color = mean_lfc, size = log10(n_elements + 1))) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
      geom_point(alpha = 0.85) +
      facet_wrap(~ Class, scales = "free_y", ncol = 1) +
      scale_color_gradient2(low = "#3267a8", mid = "white", high = "#c43c39", midpoint = 0) +
      scale_size_continuous(range = c(1.5, 6), name = "log10(copy count + 1)") +
      theme_bw(base_size = 11) +
      theme(panel.grid.major.y = element_blank(), legend.position = "bottom") +
      labs(title = "TE family enrichment", subtitle = "Top families per TE class",
           x = "Mean log2 enrichment over control", y = NULL, color = "Mean log2 enrichment")
    save_dual(p_dot, "TE_Family_Dotplot", 10, max(7, 0.22 * nrow(dot_data)))
    write_csv(family_stats, file.path(enhanced_te_dir, "TE_family_summary.csv"))
  }

  # 9c. Age-stratified distributions and continuous milliDiv relationship.
  age_data <- te_enhanced %>%
    mutate(age_group = case_when(is.na(milliDiv) ~ "unknown",
      milliDiv < 50 ~ "young (<5%)", milliDiv < 150 ~ "middle (5–15%)",
      TRUE ~ "old (≥15%)"))
  age_plot <- sample_for_plot(age_data, MAX_TE_PLOT_ROWS, seed = 20240714L)
  p_age <- ggplot(age_plot, aes(x = sample, y = lfc_over_igg, fill = age_group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_boxplot(outlier.alpha = 0.12, position = position_dodge(width = 0.8)) +
    scale_fill_viridis_d(option = "C", end = 0.9) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom") +
    labs(title = "TE enrichment by repeat age", x = NULL,
         y = "log2 enrichment over control", fill = "milliDiv group")
  save_dual(p_age, "TE_Age_Distribution", 10, 6)

  age_scatter <- age_data %>%
    filter(is.finite(milliDiv), is.finite(lfc_over_igg)) %>%
    sample_for_plot(MAX_TE_SCATTER_POINTS, seed = 20240715L)
  if (nrow(age_scatter) >= 10) {
    p_age_scatter <- ggplot(age_scatter, aes(x = milliDiv, y = lfc_over_igg, color = Class)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
      geom_point(alpha = 0.18, size = 0.7) +
      # A capped display sample plus linear smooth avoids the multi-hour loess
      # fit observed on the 5.6M-row LZH TE table.
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.7) +
      scale_color_viridis_d(option = "D", end = 0.9) +
      facet_wrap(~ sample) + theme_bw(base_size = 11) +
      labs(title = "TE enrichment versus milliDiv", x = "RepeatMasker milliDiv",
           y = "log2 enrichment over control", color = "Class")
    save_dual(p_age_scatter, "TE_Enrichment_vs_milliDiv", 12, 8)
  }
}

message("--- 增强 TE 可视化完成 ---\n")

# ==============================================================================
# ## 第 10 节：分析结束 ##
# ==============================================================================
message("==========================================================")
message("          分析流程全部顺利完成！")
message(paste("          所有结果已保存至:", file.path(getwd(), output_dir)))
message("==========================================================")
