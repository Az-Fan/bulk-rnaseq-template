# ==============================================================================
# Step06 TF数据库GSEA分析
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. 设置输出文件夹结构
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
dir.create(output_dir, showWarnings = FALSE)

# 功能富集分析根目录
functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, showWarnings = FALSE)

# TF 分析主目录
tf_analysis_dir <- file.path(functional_enrichment_dir, "tf_analysis")
dir.create(tf_analysis_dir, showWarnings = FALSE)

# ==============================================================================
# PART 1: 全局参数设置
# ==============================================================================

# --- 1. TF 基因集列表 ---
# 定义你要分析的两个基因集（文件位于 resources/tf_databases/）
# 格式：list(name = "自定义名字", file = "文件路径")
tf_datasets <- list(
  list(name = "ENCODE_2015", file = "resources/tf_databases/ENCODE_TF_ChIP-seq_2015.txt"),
  list(name = "ChEA_Consensus", file = "resources/tf_databases/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.txt")
)

# --- 2. 阈值设置 ---
calc_p_cutoff   <- 1      # 计算阈值 (保留所有)
plot_p_cutoff   <- 0.05   # 绘图显著性阈值
plot_top_n      <- 5      # 气泡图和单图绘制数量 (Top N)

# --- 3. 绘图参数 ---
single_plot_width <- 8
single_plot_height <- 6

# ==============================================================================
# PART 2: 数据加载与准备
# ==============================================================================
# 加载差异分析数据
diff_rdata_path <- file.path(output_dir, "diff_analysis/DEseq2_Diff_Annotated.Rdata")
if(file.exists(diff_rdata_path)){
  load(file = diff_rdata_path)
} else {
  stop(paste0("未找到差异分析结果文件: ", diff_rdata_path))
}

library(dplyr)
library(stringr)
library(clusterProfiler)
library(ggplot2)
library(enrichplot)
library(BiocParallel)

# 准备 GeneList (logFC 排序)
gene_df <- as.data.frame(res_annotated) %>% 
  dplyr::select(gene, logFC) %>% 
  filter(gene != "" & !is.na(gene)) %>% 
  distinct(gene, .keep_all = TRUE)

geneList <- gene_df$logFC
names(geneList) <- gene_df$gene
geneList <- sort(geneList, decreasing = TRUE)

message(paste0("GeneList 准备完成，共 ", length(geneList), " 个基因。"))


# ==============================================================================
# PART 3: 循环进行 TF GSEA 分析
# ==============================================================================

for (tf_set in tf_datasets) {
  
  dataset_name <- tf_set$name
  dataset_file <- tf_set$file
  
  message(paste0("\n=================================================="))
  message(paste0("正在分析数据集: ", dataset_name))
  message(paste0("文件路径: ", dataset_file))
  
  # 1. 创建该数据集的专属输出目录
  # 例如: results/functional_enrichment/tf_analysis/ENCODE_2015/
  current_out_dir <- file.path(tf_analysis_dir, dataset_name)
  dir.create(current_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 2. 读取 GMT 文件
  if (!file.exists(dataset_file)) {
    message(paste0("警告：文件不存在，跳过: ", dataset_file))
    next
  }
  
  # read.gmt 来自 clusterProfiler
  tf_term2gene <- read.gmt(dataset_file)
  
  # 3. 运行 GSEA
  message("正在运行 GSEA...")
  gsea_res <- GSEA(geneList, 
                   TERM2GENE = tf_term2gene,
                   pvalueCutoff = calc_p_cutoff, 
                   pAdjustMethod = "BH",
                   seed = 123,
                   BPPARAM = BiocParallel::SerialParam()) # 固定随机种子
  
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    message("  -> 未发现显著富集结果，跳过后续绘图。")
    next
  }
  
  # 4. 数据清洗 (提取 TF 名字)
  # 逻辑：去除后缀 (通常是 _hg19, _human 等) 并转为首字母大写
  # 假设格式如 "TP53_Chea_..." 或 "TP53"
  clean_desc <- gsea_res@result$Description
  # 尝试只保留第一个下划线前的部分作为 TF 名 (针对常见数据库格式)
  clean_desc <- gsub("_.*", "", clean_desc)
  clean_desc <- str_to_title(clean_desc) # 首字母大写 (Tp53) 或 str_to_upper 全大写 (TP53)
  
  gsea_res@result$Description <- clean_desc
  
  # 5. 保存结果表格
  write.csv(gsea_res@result, 
            file = file.path(current_out_dir, paste0("GSEA_Full_Table_", dataset_name, ".csv")), 
            row.names = FALSE)
  message(paste0("  -> 结果表格已保存。"))
  
  
  # ==========================================================================
  # PART 4: 绘制气泡图 (Dotplot)
  # ==========================================================================
  
  # 过滤显著结果用于气泡图
  gsea_plot_obj <- gsea_res
  gsea_plot_obj@result <- gsea_res@result %>% filter(p.adjust < plot_p_cutoff)
  
  if (nrow(gsea_plot_obj@result) > 0) {
    
    # 动态调整高度
    n_show <- min(plot_top_n, nrow(gsea_plot_obj@result))
    # 如果激活和抑制都有，split后行数翻倍，简单估算高度
    dynamic_height <- max(8, n_show * 0.4 + 2)
    
    p_dot <- dotplot(gsea_plot_obj,
                     showCategory = plot_top_n, 
                     split = ".sign",            
                     label_format = 50,
                     title = paste0("TF Enrichment: ", dataset_name)) +
      facet_grid(.~.sign) +                      
      scale_color_gradient(low = "blue", high = "red") + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    # 保存气泡图
    ggsave(filename = file.path(current_out_dir, paste0("Summary_Dotplot_", dataset_name, ".pdf")), 
           plot = p_dot, width = 10, height = dynamic_height)
    ggsave(filename = file.path(current_out_dir, paste0("Summary_Dotplot_", dataset_name, ".png")), 
           plot = p_dot, width = 10, height = dynamic_height, dpi = 300)
    
    message("  -> 汇总气泡图已保存。")
  } else {
    message("  -> 无显著通路满足绘图阈值，跳过气泡图。")
  }
  
  
  # ==========================================================================
  # PART 5: 批量绘制 GSEA 曲线图 (Single Plots)
  # ==========================================================================
  
  # 创建图片子文件夹
  plots_sub_dir <- file.path(current_out_dir, "GSEA_Curve_Plots")
  if(!dir.exists(plots_sub_dir)) dir.create(plots_sub_dir)
  old_png <- list.files(plots_sub_dir, pattern = "\\.png$", full.names = TRUE)
  if (length(old_png) > 0) file.remove(old_png)
  
  # 筛选 Top N 显著通路 (按 p.adjust < 0.05 且 NES 绝对值排序)
  top_tfs <- gsea_res@result %>% 
    filter(p.adjust < plot_p_cutoff) %>% 
    arrange(desc(abs(NES))) %>% 
    head(plot_top_n) 
  
  if (nrow(top_tfs) > 0) {
    message(paste0("  -> 正在绘制 Top ", nrow(top_tfs), " 条 GSEA 曲线图..."))
    
    for (i in 1:nrow(top_tfs)) {
      pid     <- top_tfs$ID[i]
      desc    <- top_tfs$Description[i]
      nes_val <- top_tfs$NES[i]
      
      line_col <- if(nes_val > 0) "#B31B21" else "#1465AC"
      direction <- if(nes_val > 0) "UP" else "DOWN"
      
      # 文件名
      safe_filename <- gsub("[[:punct:]]", "_", desc)
      file_path <- file.path(plots_sub_dir, paste0(i, "_", direction, "_", safe_filename, ".png"))
      
      p <- gseaplot2(gsea_res, 
                     geneSetID = pid, 
                     title = paste0(desc, " (", dataset_name, ")"), 
                     color = line_col,
                     pvalue_table = TRUE,
                     ES_geom = "line",
                     base_size = 14)
      
      ggsave(filename = file_path, plot = p, 
             width = single_plot_width, height = single_plot_height, 
             dpi = 300)
    }
    message("  -> 单图绘制完成。")
  }
}

message("\n所有 TF 分析流程结束！结果存放在 results/functional_enrichment/tf_analysis/")
