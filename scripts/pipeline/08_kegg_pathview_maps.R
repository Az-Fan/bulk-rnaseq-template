# ==============================================================================
# Step08 KEGG Pathview绘图
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch",
#    SHARED_CACHE_DIR = "shared_cache")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. 设置输出文件夹结构
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# 共享缓存目录（跨项目复用，避免重复下载）
shared_cache_dir <- Sys.getenv("SHARED_CACHE_DIR", unset = file.path(getwd(), "shared_cache"))
pathview_cache_dir <- file.path(shared_cache_dir, "pathview_kegg")
dir.create(pathview_cache_dir, recursive = TRUE, showWarnings = FALSE)

# 功能富集分析相关目录
functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, recursive = TRUE, showWarnings = FALSE)

kegg_analysis_dir <- file.path(functional_enrichment_dir, "kegg_analysis")
dir.create(kegg_analysis_dir, recursive = TRUE, showWarnings = FALSE)

# 新增：Pathview 地图专用文件夹
pathview_maps_dir <- file.path(kegg_analysis_dir, "pathview_maps")
dir.create(pathview_maps_dir, recursive = TRUE, showWarnings = FALSE)

# 差异分析结果路径
diff_rdata_path <- file.path(output_dir, "diff_analysis/DEseq2_Diff_Annotated.Rdata")


# ==============================================================================
# 1. 环境准备与数据加载
# ==============================================================================
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# 加载差异分析结果
if(file.exists(diff_rdata_path)){
  load(file = diff_rdata_path)
  # 兼容性处理：如果保存的对象名是 res_annotated，将其赋值给 res 以匹配你原代码逻辑
  if(exists("res_annotated") && !exists("res")) {
    res <- res_annotated
  }
} else {
  stop(paste0("请先运行之前的差异分析脚本，生成 ", diff_rdata_path))
}

# ==============================================================================
# 2. 制作 pathview 专属基因列表 (Entrez ID + logFC)
# ==============================================================================
print("正在构建基因列表...")

res_df <- as.data.frame(res)
if(!"gene" %in% colnames(res_df)) res_df$gene <- rownames(res_df)

# 安全 ID 映射：规避部分环境中 bitr() 的 sqlite 键校验异常
map_symbol_to_entrez <- function(symbol_vec) {
  symbol_vec <- unique(symbol_vec[!is.na(symbol_vec) & symbol_vec != ""])
  if (length(symbol_vec) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE))
  }
  ids <- AnnotationDbi::select(
    x = org.Hs.eg.db,
    keys = symbol_vec,
    columns = c("SYMBOL", "ENTREZID"),
    keytype = "SYMBOL",
    skipValidKeysTest = TRUE
  )
  ids <- ids[!is.na(ids$ENTREZID) & ids$ENTREZID != "", c("SYMBOL", "ENTREZID")]
  ids <- ids[!duplicated(ids$SYMBOL), , drop = FALSE]
  ids
}

# 关键步骤：KEGG 作图必须使用 Entrez ID
# 我们将 Gene Symbol 转换为 Entrez ID
if("ENTREZID" %in% colnames(res_df)) {
  # 如果之前已经有了 ENTREZID 列，直接使用，防止 bitr 重复查询
  merged_data <- res_df %>% filter(!is.na(ENTREZID))
} else {
  ids <- map_symbol_to_entrez(res_df$gene)
  # 合并 ID 并去重
  merged_data <- merge(res_df, ids, by.x = "gene", by.y = "SYMBOL")
}

# 保留 logFC 绝对值最大的那个，防止多对一问题
merged_data <- merged_data %>% 
  arrange(desc(abs(logFC))) %>% 
  distinct(ENTREZID, .keep_all = TRUE)

# 制作 named vector (值是 logFC，名字是 EntrezID)
geneList <- merged_data$logFC
names(geneList) <- merged_data$ENTREZID

print(paste0("基因列表构建完成，包含 ", length(geneList), " 个基因。"))

# ==============================================================================
# 3. 封装绘图函数
# ==============================================================================
# 参数说明：
# pathway_id: KEGG 通路 ID (例如 "hsa04151")
# pathway_name: (可选) 图片文件后缀名，如果不填默认用 ID

plot_kegg_pathway <- function(pathway_id, output_suffix = NULL) {
  
  # 获取当前工作目录，以便稍后恢复
  current_dir <- getwd()
  
  # 设置工作目录到我们的 pathview 专用文件夹
  # pathview 会将文件下载并生成到当前工作目录
  setwd(pathview_maps_dir)
  
  tryCatch({
    print(paste0("正在下载并绘制通路: ", pathway_id, "..."))
    
    # 核心绘图函数 (完全保持你的原样)
    pathview(gene.data  = geneList,
             pathway.id = pathway_id,
             species    = "hsa",       # 人类
             kegg.dir   = pathview_cache_dir, # 共享 KEGG 缓存目录
             limit      = list(gene=max(abs(geneList)), cpd=1), # 颜色范围自动适应最大值
             low = list(gene = "blue"),  # 下调颜色
             mid = list(gene = "white"), # 无差异颜色
             high = list(gene = "red"),  # 上调颜色
             kegg.native = TRUE,         # TRUE=生成png图(原图风格), FALSE=生成pdf
             out.suffix = output_suffix  # 文件名后缀
    )
    
    print(paste0("成功！图片已保存在: ", pathview_maps_dir, "/", pathway_id, "xxx.png"))
    
  }, error = function(e) {
    print(paste0("绘图失败: ", e$message))
    print("可能原因: 网络连接问题(无法连接KEGG官网) 或 ID 不存在。")
  }, finally = {
    # 恢复工作目录
    setwd(current_dir)
  })
}


# ==============================================================================
# 4. 在这里输入名称进行绘制！
# ==============================================================================

# 提示：请查看 results/functional_enrichment/kegg_analysis/Enrichment_Full_KEGG.csv 
# 获取你感兴趣的通路 ID

# 示例 1: PI3K-Akt signaling pathway
plot_kegg_pathway("hsa04151")

# 示例 2: Cell cycle
plot_kegg_pathway("hsa04110")

# 示例 3: TNF signaling pathway
# plot_kegg_pathway("hsa04668")

# 在下面输入你想画的 ID 即可:
# plot_kegg_pathway("hsaXXXXX")

print("Pathview 任务结束。")
