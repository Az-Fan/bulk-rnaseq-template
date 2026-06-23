# ==============================================================================
# Step13 PROGENy decouple分析（复用 Step01 结果）
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch",
#    SHARED_CACHE_DIR = "shared_cache")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

# ==============================================================================
# 步骤1：加载所有必需包（默认不在脚本内安装，避免联网失败）
# ==============================================================================
# 共享缓存目录（跨项目复用）
shared_cache_dir <- Sys.getenv("SHARED_CACHE_DIR", unset = file.path(getwd(), "shared_cache"))
dir.create(shared_cache_dir, recursive = TRUE, showWarnings = FALSE)

# 本地资源缓存目录（优先使用）
resource_cache_dir <- file.path(getwd(), "resources", "cache")
dir.create(resource_cache_dir, recursive = TRUE, showWarnings = FALSE)
progeny_cache_file <- file.path(resource_cache_dir, "progeny_human_top500.rds")

Sys.setenv(XDG_CONFIG_HOME = file.path(shared_cache_dir, "xdg_config"))
dir.create(Sys.getenv("XDG_CONFIG_HOME"), recursive = TRUE, showWarnings = FALSE)

required_pkgs <- c("decoupleR", "dplyr", "tibble", "tidyr", "ggplot2", "pheatmap", "ggrepel", "RColorBrewer", "openxlsx")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0("❌ 缺少必要包：", paste(missing_pkgs, collapse = ", "), "。请先安装后再运行 decouple.R"))
}

# 加载所有包
library(decoupleR)       # 通路活性计算
library(dplyr)           # 数据处理
library(tibble)          # 数据格式转换
library(tidyr)           # 数据重塑
library(ggplot2)         # 可视化
library(pheatmap)        # 热图
library(ggrepel)         # 基因标注
library(RColorBrewer)    # 配色
library(openxlsx)        # 结果保存为Excel

# ==============================================================================
# 步骤2：自定义参数（仅需修改这部分！）
# ==============================================================================
# -------------------------- 必改参数 --------------------------
# 1. 数据路径：读取 Step01 已产出的 VST、metadata 和差异分析结果
output_root <- Sys.getenv("OUTPUT_DIR", unset = "results/final/nobatch")
vst_path <- file.path(output_root, "data_processed", "exprSet_vst.Rdata")
metadata_path <- file.path(output_root, "data_processed", "metadata.csv")
deg_path <- file.path(output_root, "diff_analysis", "DEG_results_annotated.csv")

# 2. 物种：human/mouse（根据你的研究物种选择）
organism <- "human"

# 3. 差异分析分组：处理组vs对照组（需和 design 中的 condition 一致）
treat_group <- Sys.getenv("TREAT_GROUP", unset = "Proliferating EC")
control_group <- Sys.getenv("CONTROL_GROUP", unset = "Static EC")

# 4. 结果保存路径（自动创建文件夹，无需手动建）
output_dir <- file.path(output_root, "functional_enrichment", "tf_analysis", "decouple_progeny")

# -------------------------- 可选参数（默认即可） --------------------------
progeny_top <- 500      # PROGENy通路靶基因数量（500是教程推荐值）
minsize <- 5            # 过滤靶基因数<5的通路
plot_width <- 10        # 图片宽度（英寸）
plot_height <- 8        # 图片高度（英寸）

# ==============================================================================
# 步骤3：创建结果文件夹 + 数据加载与校验（防错）
# ==============================================================================
# 创建结果文件夹
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("✅ 结果文件夹创建成功：", output_dir, "\n")
}

# 加载 Step01 输出
if (!file.exists(vst_path)) stop("❌ 未找到 Step01 VST 矩阵：", vst_path)
if (!file.exists(metadata_path)) stop("❌ 未找到 Step01 metadata：", metadata_path)
if (!file.exists(deg_path)) stop("❌ 未找到 Step01 DEG 结果：", deg_path)

load(vst_path)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
deg_ttop <- read.csv(deg_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!exists("exprSet_vst")) {
  stop("❌ VST Rdata 中未找到 exprSet_vst 对象。")
}
if (!all(c("sample", "group") %in% colnames(metadata))) {
  stop("❌ metadata.csv 必须包含 sample 和 group 列。")
}
if (!all(c("gene_id", "gene", "stat") %in% colnames(deg_ttop))) {
  stop("❌ DEG_results_annotated.csv 必须包含 gene_id、gene 和 stat 列。")
}

metadata <- metadata[metadata$sample %in% colnames(exprSet_vst), , drop = FALSE]
metadata <- metadata[match(colnames(exprSet_vst), metadata$sample), , drop = FALSE]
if (!all(metadata$sample == colnames(exprSet_vst))) {
  stop("❌ VST 表达矩阵列名与 metadata$sample 不匹配。")
}

# 样本级 PROGENy 使用 Step01 的 VST 表达矩阵，避免直接使用 raw counts。
id_map <- deg_ttop %>%
  dplyr::select(gene_id, gene) %>%
  filter(!is.na(gene), gene != "") %>%
  distinct(gene_id, .keep_all = TRUE)

expr_df <- as.data.frame(exprSet_vst, check.names = FALSE)
expr_df$gene_id <- rownames(expr_df)
expr_mat_symbol <- expr_df %>%
  inner_join(id_map, by = "gene_id") %>%
  dplyr::select(-gene_id) %>%
  group_by(gene) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  column_to_rownames("gene") %>%
  as.matrix()
expr_mat_symbol[is.na(expr_mat_symbol)] <- 0
cat("✅ VST 表达矩阵加载成功，维度：", nrow(expr_mat_symbol), "个基因 ×", ncol(expr_mat_symbol), "个样本\n")

# contrast-level PROGENy 使用主 DESeq2 模型输出的 Wald statistic。
deg_ttop <- deg_ttop %>%
  filter(!is.na(gene), gene != "", !is.na(stat), is.finite(stat)) %>%
  arrange(desc(abs(stat))) %>%
  distinct(gene, .keep_all = TRUE)

deg_matrix <- deg_ttop %>%
  transmute(ID = gene, stat = stat) %>%
  column_to_rownames(var = "ID") %>%
  as.matrix()
cat("✅ 已加载 Step01 DESeq2 Wald statistic：", nrow(deg_matrix), "个有效基因\n")

# ==============================================================================
# 步骤6：加载PROGENy通路-靶基因网络
# ==============================================================================
net <- NULL
if (file.exists(progeny_cache_file)) {
  net <- tryCatch(readRDS(progeny_cache_file), error = function(e) NULL)
  if (!is.null(net)) cat("✅ 已加载本地 PROGENy 缓存：", progeny_cache_file, "\n")
}

if (is.null(net)) {
  net <- tryCatch(
    get_progeny(
      organism = organism,
      top = progeny_top
    ),
    error = function(e) NULL
  )
  if (!is.null(net) && nrow(net) > 0) {
    saveRDS(net, progeny_cache_file)
    cat("✅ 已写入 PROGENy 缓存：", progeny_cache_file, "\n")
  }
}
if (is.null(net) || nrow(net) == 0) {
  msg <- "未能获取 PROGENy 网络（可能网络不可用），Step13 已跳过。"
  cat("⚠️ ", msg, "\n", sep = "")
  writeLines(msg, con = file.path(output_dir, "SKIPPED_decouple_progeny_due_to_network.txt"))
  quit(save = "no", status = 0)
}
cat("✅ 加载PROGENy通路网络完成：包含", length(unique(net$source)), "个通路\n")

# ==============================================================================
# 步骤7：计算样本水平通路活性 + 可视化（热图）
# ==============================================================================
# 计算样本通路活性
sample_acts <- run_mlm(
  mat = expr_mat_symbol,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'weight',
  minsize = minsize
)

# 转换为热图矩阵（标准化z-score）
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'condition', names_from = 'source', values_from = 'score') %>%
  column_to_rownames('condition') %>%
  as.matrix() %>%
  scale() # z-score标准化

# 绘制并保存热图
png(file.path(output_dir, "pathway_activity_heatmap.png"), 
    width = plot_width, height = plot_height, units = "in", res = 300)
pheatmap(
  mat = sample_acts_mat,
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  border_color = "white",
  breaks = c(seq(-2, 0, length.out = 51), seq(0.05, 2, length.out = 49)),
  cellwidth = 20,
  cellheight = 20,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "Pathway Activity Heatmap (Z-score)"
)
dev.off()
cat("✅ 通路活性热图已保存至：", file.path(output_dir, "pathway_activity_heatmap.png"), "\n")

# ==============================================================================
# 步骤8：计算差异驱动的通路活性 + 可视化（条形图）
# ==============================================================================
# 基于stat值计算通路活性
contrast_acts <- run_mlm(
  mat = deg_matrix,
  net = net,
  .source = 'source',
  .target = 'target',
  .mor = 'weight',
  minsize = minsize
)

write.csv(
  sample_acts,
  file.path(output_dir, "PROGENy_Sample_Activity.csv"),
  row.names = FALSE
)
write.csv(
  contrast_acts,
  file.path(output_dir, "PROGENy_Contrast_Activity.csv"),
  row.names = FALSE
)

# 绘制并保存条形图
png(file.path(output_dir, "pathway_activity_barplot.png"), 
    width = plot_width, height = plot_height, units = "in", res = 300)
ggplot(contrast_acts, aes(x = reorder(source, score), y = score)) +
  geom_bar(aes(fill = score), color = "black", stat = "identity") +
  scale_fill_gradient2(low = rev(brewer.pal(11, "RdBu"))[2], 
                       mid = "whitesmoke", 
                       high = rev(brewer.pal(11, "RdBu"))[10], 
                       midpoint = 0) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid = element_blank()
  ) +
  labs(x = "Pathways", y = "Activity Score", title = "Pathway Activity (Treat vs Control)")
dev.off()
cat("✅ 通路活性条形图已保存至：", file.path(output_dir, "pathway_activity_barplot.png"), "\n")

# ==============================================================================
# 步骤9：关键通路靶基因可视化（散点图，示例：MAPK通路）
# ==============================================================================
# 选择关注的通路（可替换为你的显著通路，如"JAK-STAT"、"Hypoxia"）
focus_pathway <- "MAPK"
if (!focus_pathway %in% unique(net$source)) {
  focus_pathway <- unique(net$source)[1]
  cat("⚠️ 指定通路 MAPK 不在 PROGENy 网络中，自动改为：", focus_pathway, "\n")
}

# 提取该通路靶基因
pathway_genes <- net %>%
  filter(source == focus_pathway) %>%
  arrange(target) %>%
  mutate(ID = target) %>%
  column_to_rownames('target')

# 取交集
inter_genes <- intersect(rownames(deg_matrix), rownames(pathway_genes))
pathway_genes <- pathway_genes[inter_genes, ]
pathway_genes$t_value <- deg_matrix[inter_genes, "stat"]

# 分组标记（一致/不一致）
pathway_genes <- pathway_genes %>%
  mutate(
    color = case_when(
      weight > 0 & t_value > 0 ~ "1",
      weight > 0 & t_value < 0 ~ "2",
      weight < 0 & t_value > 0 ~ "2",
      weight < 0 & t_value < 0 ~ "1",
      TRUE ~ "3"
    )
  )

# 绘制并保存散点图
png(file.path(output_dir, paste0(focus_pathway, "_target_genes.png")), 
    width = plot_width, height = plot_height, units = "in", res = 300)
ggplot(pathway_genes, aes(x = weight, y = t_value)) +
  geom_point(size = 2.5, color = "black") +
  geom_point(aes(color = color), size = 1.5) +
  scale_color_manual(values = c(rev(brewer.pal(11, "RdBu"))[10], 
                                rev(brewer.pal(11, "RdBu"))[2], 
                                "grey")) +
  geom_label_repel(aes(label = ID), max.overlaps = 20) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_hline(yintercept = 0, linetype = 'dotted') +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(x = "Pathway Weight", y = "DESeq2 stat value", title = paste0(focus_pathway, " Target Genes"))
dev.off()
cat("✅", focus_pathway, "通路靶基因散点图已保存至：", 
    file.path(output_dir, paste0(focus_pathway, "_target_genes.png")), "\n")

# ==============================================================================
# 步骤10：保存所有数值结果（Excel）
# ==============================================================================
# 创建结果列表
results_list <- list(
  "Step01_DESeq2_Wald_statistic" = deg_ttop,
  "样本水平通路活性" = sample_acts,
  "差异驱动通路活性" = contrast_acts
)
results_list[[paste0(focus_pathway, "通路靶基因")]] <- pathway_genes

# 保存为Excel
write.xlsx(results_list, file.path(output_dir, "all_analysis_results.xlsx"), rowNames = TRUE)
cat("✅ 所有数值结果已保存至：", file.path(output_dir, "all_analysis_results.xlsx"), "\n")

# 运行完成提示
cat("\n🎉 全部分析流程运行完成！所有结果已保存至：", output_dir, "\n")
