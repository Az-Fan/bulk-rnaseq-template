# ==============================================================================
# Step11 GSVA与TF-Pathway相关性
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

functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, showWarnings = FALSE)

# GSVA 分析主目录
gsva_analysis_dir <- file.path(functional_enrichment_dir, "gsva_analysis")
dir.create(gsva_analysis_dir, showWarnings = FALSE)

# TF 分析结果目录 (用于加载)
tf_results_dir <- file.path(functional_enrichment_dir, "tf_analysis/decoupler_inference")

# ==============================================================================
# 1. 参数设置 (Control Panel)
# ==============================================================================
# 目标设置
target_tf <- "FOSL1"       # 重点关注的 TF

# 数据库设置 (Hallmark)
gs_cat    <- "H"           
gs_subcat <- NULL          

# 绘图筛选设置
n_top_pathways_heatmap <- 30  # GSVA差异热图展示数量
n_top_tfs_cor          <- 30  # 相关性热图展示前 N 个 TF

# ==============================================================================
# 2. 环境准备与数据加载
# ==============================================================================
library(GSVA)
library(msigdbr)
library(pheatmap)
library(tidyverse)
library(limma)

message(">>> [Step 1] 正在加载原始数据...")

# 加载数据路径
diff_rdata_path <- file.path(output_dir, "diff_analysis/DEseq2_Diff_Annotated.Rdata")
# 智能寻找 VST 数据
if(file.exists(file.path(output_dir, "data_processed/exprSet_vst_filtered.Rdata"))){
  vst_rdata_path <- file.path(output_dir, "data_processed/exprSet_vst_filtered.Rdata")
} else {
  vst_rdata_path <- file.path(output_dir, "data_processed/exprSet_vst.Rdata")
}
metadata_rdata_path <- file.path(output_dir, "data_processed/metadata.Rdata")

load(diff_rdata_path) # res_annotated
load(vst_rdata_path)  # exprSet_vst (或 filtered/unfiltered)
# 统一 exprSet_vst 变量名
if(!exists("exprSet_vst")) {
  if(exists("exprSet_vst_filtered")) exprSet_vst <- exprSet_vst_filtered
  if(exists("exprSet_vst_unfiltered")) exprSet_vst <- exprSet_vst_unfiltered
}
load(metadata_rdata_path)

# 加载 TF 活性矩阵
tf_file <- file.path(tf_results_dir, "TF_Activity_Sample_Matrix.csv")
if(file.exists(tf_file)){
  # check.names = FALSE 防止样本名被 R 自动修改 (例如 - 变成 .)
  tf_mat <- read.csv(tf_file, row.names = 1, check.names = FALSE)
  tf_mat <- as.matrix(tf_mat)
  message(paste(">>> TF 活性矩阵加载成功，维度:", paste(dim(tf_mat), collapse = " x ")))
} else {
  msg <- paste0("未找到 TF 活性矩阵文件: ", tf_file, "。Step11 跳过（通常是 Step10 网络资源不可用导致）。")
  message(msg)
  writeLines(msg, con = file.path(gsva_analysis_dir, "SKIPPED_gsva_due_to_missing_tf_matrix.txt"))
  quit(save = "no", status = 0)
}

# ==============================================================================
# 3. 数据预处理 (ID 转换)
# ==============================================================================
message(">>> [Step 2] 准备 GSVA 输入矩阵 (转换为 Symbol)...")

id_map <- res_annotated %>% 
  dplyr::select(gene_id, gene) %>% 
  filter(gene != "" & !is.na(gene)) %>% 
  distinct(gene_id, .keep_all = TRUE)

expr_df <- as.data.frame(exprSet_vst)
expr_df$gene_id <- rownames(expr_df)

expr_mat_symbol <- expr_df %>%
  inner_join(id_map, by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  group_by(gene) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("gene") %>%
  as.matrix()

message(paste(">>> 表达矩阵准备就绪，维度:", paste(dim(expr_mat_symbol), collapse = " x ")))

# ==============================================================================
# 4. 运行 GSVA 分析
# ==============================================================================
message(paste0(">>> [Step 3] 获取基因集 (", gs_cat, ")..."))

if (is.null(gs_subcat)) {
  m_df <- msigdbr(species = "Homo sapiens", category = gs_cat)
} else {
  m_df <- msigdbr(species = "Homo sapiens", category = gs_cat, subcategory = gs_subcat)
}
gs_list <- split(m_df$gene_symbol, m_df$gs_name)

message(">>> [Step 4] 正在运行 GSVA 计算...")
gsva_par <- gsvaParam(exprData = expr_mat_symbol, geneSets = gs_list, kcdf = "Gaussian")
gsva_res <- gsva(gsva_par)

write.csv(gsva_res, file.path(gsva_analysis_dir, "GSVA_Score_Matrix.csv"))
message(">>> GSVA 计算完成并保存。")

# ==============================================================================
# 5. 绘制 GSVA 通路活性热图
# ==============================================================================
message(">>> [Step 5] 绘制 GSVA 通路差异热图...")

pathway_sd <- apply(gsva_res, 1, sd)
n_plot <- min(n_top_pathways_heatmap, nrow(gsva_res))
top_pathways <- names(sort(pathway_sd, decreasing = TRUE))[1:n_plot]
heatmap_data_gsva <- gsva_res[top_pathways, , drop=FALSE]

annotation_col <- data.frame(Group = metadata$group)
rownames(annotation_col) <- rownames(metadata)

pdf_gsva <- file.path(gsva_analysis_dir, "GSVA_Pathway_Heatmap.pdf")
pheatmap(heatmap_data_gsva,
         annotation_col = annotation_col,
         cluster_rows = TRUE, cluster_cols = TRUE,
         scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         border_color = NA, fontsize_row = 8,
         main = paste0("Top ", n_plot, " Variable Pathways (GSVA)"),
         filename = pdf_gsva, width = 10, height = 8)

message(paste(">>> GSVA 热图已保存:", pdf_gsva))

# ==============================================================================
# 6. TF - Pathway 相关性分析
# ==============================================================================
message(">>> [Step 6] 计算 TF 与 通路的相关性...")

# 检查样本名匹配
common_samples <- intersect(colnames(tf_mat), colnames(gsva_res))
if(length(common_samples) == 0) {
  stop("❌ 错误：TF 矩阵和 GSVA 矩阵没有共同的样本名！请检查样本命名是否一致。")
}
message(paste(">>> 成功匹配到", length(common_samples), "个共同样本。"))

sub_tf_mat <- tf_mat[, common_samples, drop=FALSE]
sub_gsva_res <- gsva_res[, common_samples, drop=FALSE]

# 筛选 Top TF (SD 排序)
tf_sd <- apply(sub_tf_mat, 1, sd)
tf_sd_sorted <- sort(tf_sd, decreasing = TRUE)
n_limit <- min(n_top_tfs_cor, length(tf_sd_sorted))
top_vars_tf <- names(tf_sd_sorted)[1:n_limit]

# 强制添加感兴趣的 TF
interest_tfs <- c(target_tf) # 可以添加更多
valid_interest_tfs <- intersect(interest_tfs, rownames(sub_tf_mat))
use_tfs <- unique(c(top_vars_tf, valid_interest_tfs))
use_tfs <- use_tfs[!is.na(use_tfs)]

sub_tf_mat <- sub_tf_mat[use_tfs, , drop=FALSE]

# 筛选 Pathway (如果太多只取 Top 50)
if (nrow(sub_gsva_res) > 60) {
  use_pathways <- names(sort(apply(sub_gsva_res, 1, sd), decreasing = TRUE))[1:50]
  sub_gsva_res <- sub_gsva_res[use_pathways, , drop=FALSE]
}

# 计算相关性
cor_mat <- cor(t(sub_tf_mat), t(sub_gsva_res), method = "pearson")
write.csv(cor_mat, file.path(gsva_analysis_dir, "TF_Pathway_Correlation_Matrix.csv"))

# ==============================================================================
# 7. 绘制相关性热图
# ==============================================================================
message(">>> [Step 7] 绘制相关性热图...")

pdf_cor <- file.path(gsva_analysis_dir, "TF_Pathway_Correlation_Heatmap.pdf")
pheatmap(cor_mat,
         main = paste("Correlation: TF Activity vs", gs_cat, "Pathways"),
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = seq(-1, 1, length.out = 100),
         border_color = "grey95",
         treeheight_row = 20, treeheight_col = 20,
         fontsize_row = 8, fontsize_col = 8,
         filename = pdf_cor, width = 12, height = 10)

message(paste(">>> 相关性热图已保存:", pdf_cor))

# ==============================================================================
# 8. 打印目标 TF 的关联通路
# ==============================================================================
if (target_tf %in% rownames(cor_mat)) {
  target_vals <- cor_mat[target_tf, ]
  message(paste0("\n=== [", target_tf, "] 调控预测分析 ==="))
  message(">>> 最强正相关通路 (可能促进):")
  print(round(sort(target_vals, decreasing = TRUE)[1:5], 3))
  message("\n>>> 最强负相关通路 (可能抑制):")
  print(round(sort(target_vals, decreasing = FALSE)[1:5], 3))
}

message(paste0("\n>>> 全部分析完成！"))
