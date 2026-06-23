# ==============================================================================
# Step07 指定GSEA通路热图
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

# 输入数据路径 (加载之前处理好的数据)
data_processed_dir <- file.path(output_dir, "data_processed")
gsea_tables_dir <- file.path(output_dir, "functional_enrichment/gsea_analysis/tables")

# 输出路径 (新建一个文件夹专门存通路热图)
gsea_analysis_dir <- file.path(output_dir, "functional_enrichment/gsea_analysis")
heatmap_out_dir <- file.path(gsea_analysis_dir, "pathway_heatmaps")
dir.create(heatmap_out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. 参数设置 (Control Panel - 只需要改这里)
# ==============================================================================
# 输入该通路所在的 GSEA 结果表格 ---
# (请确保文件名正确，位于 results/functional_enrichment/gsea_analysis/tables/ 下)
# 常用文件示例: "GSEA_Full_Table_KEGG.csv" 或 "GSEA_Full_Table_Hallmark.csv"
target_database <- Sys.getenv("GSEA_HEATMAP_DATABASE", unset = "Hallmark")
target_csv_name <- paste0("GSEA_Full_Table_", target_database, ".csv")


# ==============================================================================
# 2. 数据加载与准备
# ==============================================================================
library(pheatmap)
library(dplyr)
library(stringr)
library(RColorBrewer)

# 1. 加载表达矩阵和分组信息
# 注意：这里加载的是之前生成的 exprSet_for_heatmap.Rdata
load(file = file.path(data_processed_dir, "exprSet_for_heatmap.Rdata"))
load(file = file.path(data_processed_dir, "metadata.Rdata"))

# 变量名兼容 (之前的脚本保存的变量名是 exprSet_heatmap)
exprSet <- exprSet_heatmap 

# 2. 加载 GSEA 结果表格
csv_file_path <- file.path(gsea_tables_dir, target_csv_name)
if(!file.exists(csv_file_path)){
  stop(paste0("错误：找不到 GSEA 结果文件: ", csv_file_path))
}
gsea_df <- read.csv(csv_file_path)


# ---  输入你想绘制的通路名称 ---
# (需与 GSEA 表格中 Description 列一致，支持部分匹配，不区分大小写)
# target_pathway <- "Mtorc1 signaling" 
target_pathway <- Sys.getenv("GSEA_HEATMAP_PATHWAY", unset = "")
if (target_pathway == "") {
  message("No GSEA heatmap pathway selected; Step07 skipped.")
  quit(save = "no", status = 0)
}

# 3. 准备分组注释
annotation_col <- data.frame(group = metadata$group)
rownames(annotation_col) <- metadata$sample

# ==============================================================================
# 3. 提取通路基因 (Core Enrichment)
# ==============================================================================

# 查找通路
match_index <- grep(pattern = target_pathway, x = gsea_df$Description, ignore.case = TRUE)

if(length(match_index) == 0){
  warning(paste0("未找到目标通路 '", target_pathway, "'，本脚本跳过绘图。可手动修改 target_pathway 后重跑。"))
  quit(save = "no", status = 0)
} else if (length(match_index) > 1){
  message("提示：匹配到多个通路，默认使用第一个匹配项：")
  print(gsea_df$Description[match_index])
  match_index <- match_index[1]
}

# 获取完整通路名称 (用于文件名)
full_pathway_name <- gsea_df$Description[match_index]

# 提取核心富集基因
# 注意：之前的 GSEA 分析使用的是 Symbol，所以这里不需要 ID 转换
core_gene_str <- gsea_df$core_enrichment[match_index]
core_genes <- unlist(strsplit(core_gene_str, "/"))

message(paste0("正在处理通路：", full_pathway_name))
message(paste0("提取到核心基因数量：", length(core_genes)))

# ==============================================================================
# 4. 匹配表达矩阵
# ==============================================================================

# 取交集：确保提取的基因在表达矩阵中存在
valid_genes <- intersect(core_genes, rownames(exprSet))

# 提取热图数据
heatdata <- exprSet[valid_genes, ]

# 检查是否提取到数据
if(nrow(heatdata) < 2){
  stop("错误：该通路匹配到的基因在表达矩阵中少于2个，无法绘制热图。")
}

# ==============================================================================
# 5. 绘制并保存热图
# ==============================================================================

# 文件名处理
safe_name <- gsub(" ", "_", full_pathway_name)
safe_name <- gsub("[^A-Za-z0-9_]", "", safe_name) # 去除特殊字符
save_path_pdf <- file.path(heatmap_out_dir, paste0("Heatmap_", safe_name, ".pdf"))
save_path_png <- file.path(heatmap_out_dir, paste0("Heatmap_", safe_name, ".png"))

# 动态调整参数
# 如果基因很少(比如<20)，格子可以高一点；如果很多(>50)，格子矮一点
dynamic_cellheight <- ifelse(nrow(heatdata) > 50, 8, 
                             ifelse(nrow(heatdata) < 20, 15, 10))

dynamic_fontsize_row <- ifelse(nrow(heatdata) > 50, 6, 
                               ifelse(nrow(heatdata) < 20, 10, 8))

# 动态计算 PDF 高度，防止基因名被切掉
# 基础高度 5 + 每行基因增加的高度
pdf_height_calc <- 5 + (nrow(heatdata) * dynamic_cellheight / 72) 

# 为了生成 PDF，我们使用 pdf() 函数包裹 pheatmap，这样可以精确控制画布大小
pdf(file = save_path_pdf, width = 8, height = pdf_height_calc)

pheatmap(heatdata, 
         cluster_rows = TRUE,    # 行聚类
         cluster_cols = TRUE,    # 列聚类
         annotation_col = annotation_col, 
         annotation_legend = TRUE, 
         
         # --- 显示设置 ---
         show_rownames = TRUE,   
         show_colnames = FALSE,  
         
         # --- 尺寸与字体 ---
         cellwidth = 25, 
         cellheight = dynamic_cellheight, 
         fontsize = 10,
         fontsize_row = dynamic_fontsize_row,
         
         # --- 颜色 ---
         scale = "row",          
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         border_color = NA,
         main = paste0("GSEA Core Enrichment: ", full_pathway_name)
)
dev.off() # 关闭 PDF 设备

# 再保存一份 PNG 用于快速查看
pheatmap(heatdata, 
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_col = annotation_col, annotation_legend = TRUE,
         show_rownames = TRUE, show_colnames = FALSE,
         cellwidth = 25, cellheight = dynamic_cellheight, 
         fontsize = 10, fontsize_row = dynamic_fontsize_row,
         scale = "row", color = colorRampPalette(c("navy", "white", "firebrick3"))(100), 
         border_color = NA,
         main = paste0("GSEA Core Enrichment: ", full_pathway_name),
         filename = save_path_png,
         width = 8, height = pdf_height_calc # 使用计算好的高度
)

message(paste0("热图绘制成功！\nPDF文件：", save_path_pdf, "\nPNG文件：", save_path_png))
