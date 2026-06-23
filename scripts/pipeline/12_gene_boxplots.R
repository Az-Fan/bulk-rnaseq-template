# ==============================================================================
# Step12 指定基因箱线图
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())
# !!! 在这里输入你想要绘制 Boxplot 的基因名 (SYMBOL) !!!
# 示例：
# genes_to_plot_list <- c("GAPDH", "ACTB", "TP53", "MYC")
genes_to_plot_raw <- Sys.getenv("SELECTED_GENES", unset = "")
genes_to_plot_list <- trimws(unlist(strsplit(genes_to_plot_raw, ",")))
genes_to_plot_list <- genes_to_plot_list[genes_to_plot_list != ""]


# 0. 设置输出文件夹结构 (根据项目根目录的 'results' 结构)
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
dir.create(output_dir, showWarnings = FALSE)

# 数据处理相关的输出
data_processed_dir <- file.path(output_dir, "data_processed")
dir.create(data_processed_dir, showWarnings = FALSE)

# 差异分析结果输出 (仅用于加载数据)
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
dir.create(diff_analysis_dir, showWarnings = FALSE) # 确保存在，即使只读

# 图表输出
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# 质量控制图子文件夹 (本脚本不生成，但确保存在)
qc_plots_dir <- file.path(plots_dir, "qc_plots")
dir.create(qc_plots_dir, showWarnings = FALSE)

# 单个基因 Boxplot 的子文件夹
individual_boxplot_dir <- file.path(plots_dir, "individual_gene_boxplots")
dir.create(individual_boxplot_dir, showWarnings = FALSE)


library(dplyr)
library(tibble)
library(ggplot2) # 用于绘制boxplot
library(tidyr)   # 用于数据整理

# 1. 读入原始数据 (确保这些文件存在于正确路径)
# 注意：第一个脚本在移除离群样本后会保存 exprSet_vst_filtered.Rdata
# 如果没有移除，则保存为 exprSet_vst_unfiltered.Rdata
# 这里应该优先加载过滤后的数据
if (file.exists(file.path(data_processed_dir, "exprSet_vst_filtered.Rdata"))) {
  load(file.path(data_processed_dir, "exprSet_vst_filtered.Rdata"))
  cat("已加载过滤后的VST表达矩阵: exprSet_vst_filtered.Rdata\n")
} else if (file.exists(file.path(data_processed_dir, "exprSet_vst.Rdata"))) {
  load(file.path(data_processed_dir, "exprSet_vst.Rdata"))
  cat("已加载未过滤的VST表达矩阵: exprSet_vst.Rdata\n")
} else {
  stop("Error: VST expression data (exprSet_vst_filtered.Rdata or exprSet_vst.Rdata) not found.")
}

load(file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata")) # 加载差异分析结果
load(file.path(data_processed_dir, "metadata.Rdata"))              # 加载metadata

# 将加载的res_annotated重新命名为res，以匹配原始脚本逻辑
res <- res_annotated

# 2. 提取 ID 与 Symbol 的对应关系
if(!"gene_id" %in% colnames(res) | !"gene" %in% colnames(res)){
  stop("Error: 'res' data frame must contain 'gene_id' and 'gene' columns for mapping.")
}

gene_map <- res %>%
  dplyr::select(gene_id, gene) %>%
  filter(gene != "" & !is.na(gene)) %>%
  distinct(gene_id, .keep_all = TRUE)

# 3. 处理表达矩阵
# 确保 exprSet_vst 的行名是 gene_id
expr_clean <- as.data.frame(exprSet_vst) %>%
  rownames_to_column("gene_id") %>%
  inner_join(gene_map, by = "gene_id") %>%
  dplyr::select(-gene_id)

# 4. 基因名去重 (取平均表达量最高的)
cat("去重前基因数量：", nrow(expr_clean), "\n")
expr_clean$avg_exp <- rowMeans(dplyr::select(expr_clean, -gene))

expr_clean <- expr_clean %>%
  arrange(desc(avg_exp)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(-avg_exp)
cat("去重后基因数量：", nrow(expr_clean), "\n")


# 5. 格式整理：转换为长格式并添加分组信息
rownames(expr_clean) <- expr_clean$gene
expr_clean <- dplyr::select(expr_clean, -gene)

# 匹配 metadata 的样本顺序
# 确保 expr_clean 的列名与 metadata 的样本名一致
if (!all(colnames(expr_clean) %in% rownames(metadata))) {
  stop("Error: Sample names in VST expression matrix do not match metadata sample names.")
}

metadata_ordered <- metadata[colnames(expr_clean), , drop = FALSE]

metadata_for_join <- metadata_ordered
if (!"sample" %in% colnames(metadata_for_join)) {
  # 如果 metadata_for_join 还没有 'sample' 列，那么将其行名作为 'sample' 列
  metadata_for_join <- metadata_for_join %>% rownames_to_column("sample")
}

# 转换为长格式数据，以便ggplot2绘图
plot_data_long <- expr_clean %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  left_join(metadata_for_join, by = "sample") # 使用转换后的 metadata_for_join

# 6. 保存用于后续绘图的整合数据 (如果需要的话)
save(plot_data_long, file = file.path(data_processed_dir, "plot_data_long_for_selected_boxplots.Rdata"))
cat("恭喜！用于 Boxplot 绘图的整合长格式数据文件已生成并保存至：", file.path(data_processed_dir, "plot_data_long_for_selected_boxplots.Rdata"), "\n")

# 7. 为特定基因绘制并保存 Boxplot
cat("\n开始为指定基因绘制 Boxplot...\n")



# 检查指定基因是否在数据中
available_genes <- unique(plot_data_long$gene)
genes_found <- genes_to_plot_list[genes_to_plot_list %in% available_genes]
genes_not_found <- genes_to_plot_list[!genes_to_plot_list %in% available_genes]

if (length(genes_not_found) > 0) {
  warning("以下基因未在数据中找到，将跳过绘制: ", paste(genes_not_found, collapse = ", "), "\n")
}

if (length(genes_found) == 0) {
  cat("警告：没有找到任何指定的基因，未生成任何 Boxplot。\n")
} else {
  count_plots <- 0
  for (gene_name in genes_found) {
    # 获取当前基因的数据
    current_gene_data <- plot_data_long %>% filter(gene == gene_name)
    
    # 绘制 Boxplot
    p <- ggplot(current_gene_data, aes(x = group, y = expression, fill = group)) +
      geom_boxplot(outlier.shape = NA) + # 不显示离群点
      geom_jitter(width = 0.2, alpha = 0.6, size = 2) + # 添加散点以显示每个样本
      labs(title = paste0("Expression of ", gene_name),
           x = "Group",
           y = "VST Expression") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.title = element_text(size = 12),
            axis.text = element_text(size = 10),
            legend.position = "none") +
      scale_fill_manual(values = setNames(c("#6495ED", "#FF6347"), levels(metadata_for_join$group)))
    
    # 保存图片，以基因名命名
    safe_gene_name <- gsub("[^A-Za-z0-9_.-]", "_", gene_name)
    file_path <- file.path(individual_boxplot_dir, paste0("Boxplot_", safe_gene_name, ".png"))
    ggsave(file_path, plot = p, width = 6, height = 5, dpi = 300)
    
    count_plots <- count_plots + 1
  }
  cat(paste0("已绘制 ", count_plots, " 个指定基因的 Boxplot，并保存到 '", individual_boxplot_dir, "' 文件夹中。\n"))
}
