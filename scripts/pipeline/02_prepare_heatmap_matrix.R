# ==============================================================================
# Step02 准备热图输入矩阵
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# 0. 设置输出文件夹结构 (根据项目根目录的 'results' 结构)
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results_01")
dir.create(output_dir, showWarnings = FALSE)

# 数据处理相关的输出
data_processed_dir <- file.path(output_dir, "data_processed")
dir.create(data_processed_dir, showWarnings = FALSE)

# 差异分析结果输出 (仅用于加载数据)
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
dir.create(diff_analysis_dir, showWarnings = FALSE) # 确保存在，即使只读

# 这个脚本不直接生成图表，但为了完整性，可以定义图表目录
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)



library(dplyr)
library(tibble)

# 1. 读入原始数据 (确保这些文件存在于正确路径)
# 注意：这里加载exprSet_vst时，需要考虑第一个脚本是否进行了离群样本过滤
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

########################################
### 获取清洁数据
### 1.获取基因 ID 与 Symbol 的映射关系
# 确保 'gene_id' 和 'gene' 列存在
if(!"gene_id" %in% colnames(res) | !"gene" %in% colnames(res)){
  stop("Error: 'res' data frame must contain 'gene_id' and 'gene' columns for mapping.")
}
gene_map <- res %>%
  dplyr::select(gene_id, gene) %>%
  filter(gene != "" & !is.na(gene)) %>% # 过滤掉空的或NA的基因Symbol
  distinct(gene_id, .keep_all = TRUE) # 如果有重复的gene_id，保留第一个

### 2.准备表达量数据：将 VST 表达矩阵的行名转换为 'gene_id' 列
exprSet_vst_df <- as.data.frame(exprSet_vst) %>%
  rownames_to_column("gene_id")

### 3.交叉合并：将基因 Symbol 添加到表达矩阵
# 使用 inner_join 确保只保留有 Symbol 的基因，且 gene_id 匹配
exprSet_merged <- inner_join(gene_map, exprSet_vst_df, by = "gene_id")

### 4.基因名称去重 (保留平均表达量最高的方法)
# 列转行名之前一定要去重，因为行名不支持重复
cat("去重前基因数量：", nrow(exprSet_merged), "\n")
exprSet_heatmap <- exprSet_merged %>%
  dplyr::select(-gene_id) %>% # 移除 gene_id 列，因为我们最终要用 gene 作为行名
  mutate(avg_expression = rowMeans(dplyr::select(., -gene))) %>% # 计算每行（基因）的平均表达量
  arrange(desc(avg_expression)) %>% # 按照平均表达量降序排列
  distinct(gene, .keep_all = TRUE) %>% # 根据 gene 列去重，保留排在前面的（即平均表达量最高的）
  dplyr::select(-avg_expression) # 移除 avg_expression 辅助列
cat("去重后基因数量：", nrow(exprSet_heatmap), "\n")


### 5.列变成行名
rownames(exprSet_heatmap) <- exprSet_heatmap$gene
exprSet_heatmap <- exprSet_heatmap %>% dplyr::select(-gene) # 移除 gene 列，因为它已经是行名了

### 保存数据，用于后续画热图
save(exprSet_heatmap, file = file.path(data_processed_dir, "exprSet_for_heatmap.Rdata"))
cat("\n恭喜！用于热图的表达数据已生成并保存至：", file.path(data_processed_dir, "exprSet_for_heatmap.Rdata"), "\n")

# 如果需要，这里可以打印前几行数据进行检查
# print(head(exprSet_heatmap))
