# ==============================================================================
# Step10 TF活性推断（decoupleR）
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch",
#    SHARED_CACHE_DIR = "shared_cache")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# 共享缓存目录（跨项目复用，避免重复下载 OmniPath 资源）
shared_cache_dir <- Sys.getenv("SHARED_CACHE_DIR", unset = file.path(getwd(), "shared_cache"))
dir.create(shared_cache_dir, recursive = TRUE, showWarnings = FALSE)

xdg_config_dir <- file.path(shared_cache_dir, "xdg_config")
Sys.setenv(XDG_CONFIG_HOME = xdg_config_dir)
dir.create(Sys.getenv("XDG_CONFIG_HOME"), recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0. 输出文件夹设置
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
dir.create(output_dir, showWarnings = FALSE)

functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, showWarnings = FALSE)

tf_analysis_dir <- file.path(functional_enrichment_dir, "tf_analysis")
dir.create(tf_analysis_dir, showWarnings = FALSE)

tf_boxplot_dir <- file.path(tf_analysis_dir, "boxplots")
dir.create(tf_boxplot_dir, recursive = TRUE, showWarnings = FALSE)

# decoupleR 专用文件夹
decoupler_out_dir <- file.path(tf_analysis_dir, "decoupler_inference")
dir.create(decoupler_out_dir, showWarnings = FALSE)

# 本地资源缓存目录（优先使用）
resource_cache_dir <- file.path(getwd(), "resources", "cache")
dir.create(resource_cache_dir, recursive = TRUE, showWarnings = FALSE)
collectri_cache_file <- file.path(resource_cache_dir, "collectri_human.rds")

# 特定 TF 分析结果文件夹
single_tf_out_dir <- file.path(decoupler_out_dir, "single_tf_plots")
dir.create(single_tf_out_dir, showWarnings = FALSE)

# ==============================================================================
# 1. 参数设置 (Control Panel)
# ==============================================================================
# 全局分析参数
n_top_bar       <- 25         # 差异活性条形图展示数量 (激活/抑制各Top N)
n_top_heatmap   <- 30         # 样本活性热图展示数量
n_top_jaccard   <- 30         # Jaccard 分析选取的 Top TF 数量

# --- 【在此处输入你想关注的 TF】 (PART 7 会用到) ---
target_tf_list <- c("FOSL1")  

# ==============================================================================
# 2. 环境准备与数据加载
# ==============================================================================
packages <- c("decoupleR", "tidyverse", "pheatmap", "ggrepel", "clusterProfiler", "enrichplot", "ggpubr", "grid", "BiocParallel")
missing_pkgs <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0("缺少必要包: ", paste(missing_pkgs, collapse = ", "), "。请先安装后再运行。"))
}

library(decoupleR)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(ggpubr)
library(grid)
library(BiocParallel)

# 加载数据路径
diff_rdata_path <- file.path(output_dir, "diff_analysis/DEseq2_Diff_Annotated.Rdata")
# 智能寻找 VST 数据
if(file.exists(file.path(output_dir, "data_processed/exprSet_vst_filtered.Rdata"))){
  vst_rdata_path <- file.path(output_dir, "data_processed/exprSet_vst_filtered.Rdata")
} else {
  vst_rdata_path <- file.path(output_dir, "data_processed/exprSet_vst.Rdata")
}
metadata_rdata_path <- file.path(output_dir, "data_processed/metadata.Rdata")

message(">>> 正在加载数据...")
load(diff_rdata_path) # 加载 res_annotated
load(vst_rdata_path)  # 加载 exprSet_vst
# 统一 exprSet_vst 变量名
if(!exists("exprSet_vst")) {
  if(exists("exprSet_vst_filtered")) exprSet_vst <- exprSet_vst_filtered
  if(exists("exprSet_vst_unfiltered")) exprSet_vst <- exprSet_vst_unfiltered
}
load(metadata_rdata_path)

# ==============================================================================
# 3. 数据清洗 (构建输入矩阵)
# ==============================================================================
message(">>> [Step 1] 构建差异矩阵与表达矩阵...")

# --- A. 构建差异统计量矩阵（DESeq2 Wald statistic） ---
df_base <- as.data.frame(res_annotated) %>%
  filter(!is.na(gene) & gene != "" & !is.na(stat) & is.finite(stat)) %>%
  arrange(desc(abs(stat))) %>%
  distinct(gene, .keep_all = TRUE) # 去重

deg_mat <- as.matrix(df_base[, "stat", drop = FALSE])
rownames(deg_mat) <- df_base$gene

# --- B. 构建表达量矩阵 (Gene Symbol) ---
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

# ==============================================================================
# 4. 获取网络与计算活性 (decoupleR) - [修复版]
# ==============================================================================
message(">>> [Step 2] 获取网络并计算 ulm 活性...")

# --- 修复 1: 解决日志权限报错 (Permission denied) - [兼容版] ---
options(omnipath.logfile = NULL)
options(omnipath.console_loglevel = 'warn')

# 设置缓存目录
if (requireNamespace("OmnipathR", quietly = TRUE) &&
    exists("omnipath_set_cachedir", where = asNamespace("OmnipathR"), mode = "function")) {
  cache_dir <- file.path(shared_cache_dir, "omnipath")
  if(!dir.exists(cache_dir)) dir.create(cache_dir, recursive = TRUE)
  OmnipathR::omnipath_set_cachedir(cache_dir)
}

# --- 修复 2: 健壮的网络获取 ---
net <- NULL
if (file.exists(collectri_cache_file)) {
  message("  -> 读取本地 CollecTRI 缓存: ", collectri_cache_file)
  net <- tryCatch(readRDS(collectri_cache_file), error = function(e) NULL)
}

if (is.null(net)) {
  net <- tryCatch({
    message("  -> 尝试从 OmniPath 服务器下载 CollecTRI 网络...")
    decoupleR::get_collectri(organism = 'human', split_complexes = FALSE)
  }, error = function(e) {
    message("  -> 标准下载失败。")
    if (requireNamespace("OmnipathR", quietly = TRUE)) {
      message("  -> 尝试使用 OmnipathR 备用数据源...")
      return(tryCatch({
        net_raw <- OmnipathR::import_transcriptional_interactions(resources = "CollecTRI", organism = 9606)
        net_raw %>%
          dplyr::filter(!is.na(is_stimulation)) %>% 
          dplyr::mutate(mor = ifelse(is_stimulation == 1, 1, -1)) %>%
          dplyr::select(source = source_genesymbol, target = target_genesymbol, mor) %>%
          dplyr::distinct()
      }, error = function(e2) NULL))
    }
    NULL
  })
  if (!is.null(net) && nrow(net) > 0) {
    saveRDS(net, collectri_cache_file)
    message("  -> 已缓存 CollecTRI 到: ", collectri_cache_file)
  }
}

if (is.null(net) || nrow(net) == 0) {
  msg <- "未能获取 CollecTRI 网络（可能是无外网或服务不可用），Step10 已跳过，不影响其余流程。"
  message(msg)
  writeLines(msg, con = file.path(decoupler_out_dir, "SKIPPED_tf_activity_due_to_network.txt"))
  quit(save = "no", status = 0)
}

if (!all(c("source", "target", "mor") %in% colnames(net))) {
  if ("weight" %in% colnames(net)) {
    net$mor <- net$weight
  } else {
    stop("CollecTRI 缓存缺少 source/target/mor 列。")
  }
}
net <- net %>% dplyr::select(source, target, mor) %>% dplyr::distinct()

# 4.1 差异活性 (ulm)
contrast_acts <- run_ulm(mat = deg_mat, net = net, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
contrast_df <- as.data.frame(contrast_acts) %>% arrange(desc(score))
write.csv(contrast_df, file = file.path(decoupler_out_dir, "TF_Activity_Contrast_Full.csv"), row.names = FALSE)

# 4.2 样本活性 (ulm)
sample_acts <- run_ulm(mat = expr_mat_symbol, net = net, .source = 'source', .target = 'target', .mor = 'mor', minsize = 5)
sample_acts_mat <- sample_acts %>%
  pivot_wider(id_cols = 'source', names_from = 'condition', values_from = 'score') %>%
  column_to_rownames('source') %>%
  as.matrix()
write.csv(sample_acts_mat, file = file.path(decoupler_out_dir, "TF_Activity_Sample_Matrix.csv"))

# 4.3 GSEA 分析 (用于后续单图验证；ranking metric = DESeq2 Wald statistic)
message(">>> [Step 3] 运行 GSEA (用于一致性验证)...")
gene_list <- df_base$stat
names(gene_list) <- df_base$gene
gene_list <- sort(gene_list, decreasing = TRUE)# 2. 【关键修复】GSEA 要求 geneList 必须严格降序排列 (按数值，而非绝对值)
term2gene <- net %>% select(source, target) %>% distinct()# 3. 准备基因集
set.seed(123)# 4. 运行 GSEA
gsea_res <- GSEA(geneList = gene_list, 
                 TERM2GENE = term2gene, 
                 pvalueCutoff = 1, 
                 minGSSize = 5, 
                 maxGSSize = 500, 
                 eps = 0, 
                 nPermSimple = 10000,
                 BPPARAM = BiocParallel::SerialParam(),
                 verbose = FALSE)


# ==============================================================================
# 5. 生成全局概览 PDF 报告
# ==============================================================================
message(">>> [Step 4] 生成全局概览 PDF 报告...")
pdf_file <- file.path(decoupler_out_dir, "Global_TF_Analysis_Report.pdf")
pdf(pdf_file, width = 10, height = 12, onefile = TRUE)

# 5.1 条形图
plot_df <- contrast_df
top_rows <- c(1:min(n_top_bar, nrow(plot_df)), max(1, nrow(plot_df)-n_top_bar+1):nrow(plot_df))
tf_plot_data <- plot_df[unique(top_rows), ]
tf_plot_data$type <- ifelse(tf_plot_data$score > 0, "Activated", "Inhibited")
tf_plot_data$source <- factor(tf_plot_data$source, levels = tf_plot_data$source[order(tf_plot_data$score)])

p_bar <- ggplot(tf_plot_data, aes(x = source, y = score, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c("Activated" = "firebrick3", "Inhibited" = "navy")) +
  coord_flip() + theme_bw() +
  labs(title = "Top Predicted TFs (Activity)", y = "ulm Score (t-value)", x = "TF")
print(p_bar)

# 5.2 样本热图
top_var_tfs <- names(sort(apply(sample_acts_mat, 1, sd), decreasing = TRUE))[1:min(n_top_heatmap, nrow(sample_acts_mat))]
pheatmap(sample_acts_mat[top_var_tfs, ], scale = "row", 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         main = "Top Variable TF Activities Across Samples", fontsize_row = 8)

# 5.3 一致性散点图
ulm_res <- contrast_df %>% select(ID = source, ulm_score = score, ulm_pval = p_value)
gsea_res_df <- gsea_res@result %>% select(ID, gsea_nes = NES, gsea_pval = p.adjust)
merged_res <- inner_join(ulm_res, gsea_res_df, by = "ID")
merged_res$Type <- "Not Sig"
merged_res$Type[merged_res$ulm_pval < 0.05 & merged_res$gsea_pval < 0.05 & sign(merged_res$ulm_score) == sign(merged_res$gsea_nes)] <- "Consistent"

# 5.3a Top TF expression overview
# Purpose: visualize the expression of the top TFs supported by TF-target GSEA
# and decoupleR TF activity inference. X-axis label colors summarize whether
# the TF's own expression direction agrees with predicted activity:
# green = consistent, yellow = weak expression change, red = opposite.
make_top_tf_expression_plot <- function(
    contrast_df, merged_res, expr_mat_symbol, metadata,
    tf_boxplot_dir, n_top = 20, weak_delta = 0.10) {
  if (!"group" %in% colnames(metadata)) {
    message("  -> metadata has no 'group' column; skip Top TF expression boxplot.")
    return(invisible(NULL))
  }

  control_group <- Sys.getenv("CONTROL_GROUP", unset = "")
  treat_group <- Sys.getenv("TREAT_GROUP", unset = "")
  group_levels <- unique(as.character(metadata$group))
  if (control_group %in% group_levels && treat_group %in% group_levels) {
    group_levels <- c(control_group, treat_group, setdiff(group_levels, c(control_group, treat_group)))
  }

  gsea_ranked <- merged_res %>%
    filter(!is.na(gsea_pval), gsea_pval < 0.05) %>%
    arrange(desc(abs(ulm_score))) %>%
    pull(ID)
  activity_ranked <- contrast_df %>%
    arrange(desc(abs(score))) %>%
    pull(source)
  selected_tfs <- unique(c(gsea_ranked, activity_ranked))
  selected_tfs <- selected_tfs[selected_tfs %in% rownames(expr_mat_symbol)]
  selected_tfs <- head(selected_tfs, n_top)

  if (length(selected_tfs) == 0) {
    message("  -> no top TFs are present in the expression matrix; skip Top TF expression boxplot.")
    return(invisible(NULL))
  }

  metadata_ordered <- metadata[colnames(expr_mat_symbol), , drop = FALSE]
  metadata_ordered$sample <- rownames(metadata_ordered)
  metadata_ordered$group <- factor(as.character(metadata_ordered$group), levels = group_levels)

  plot_data <- as.data.frame(t(expr_mat_symbol[selected_tfs, , drop = FALSE])) %>%
    rownames_to_column("sample") %>%
    pivot_longer(cols = -sample, names_to = "TF", values_to = "expression") %>%
    left_join(metadata_ordered[, c("sample", "group")], by = "sample")

  activity_lookup <- contrast_df %>%
    select(TF = source, activity_score = score, activity_p = p_value)
  gsea_lookup <- merged_res %>%
    select(TF = ID, gsea_nes, gsea_pval)

  expression_summary <- plot_data %>%
    group_by(TF, group) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = mean_expression)

  if (control_group %in% colnames(expression_summary) && treat_group %in% colnames(expression_summary)) {
    expression_summary$expression_delta <- expression_summary[[treat_group]] - expression_summary[[control_group]]
  } else {
    expression_summary$expression_delta <- NA_real_
  }
  expression_delta <- expression_summary %>% select(TF, expression_delta)

  annotation <- tibble(TF = selected_tfs) %>%
    left_join(activity_lookup, by = "TF") %>%
    left_join(gsea_lookup, by = "TF") %>%
    left_join(expression_delta, by = "TF") %>%
    mutate(
      expression_activity_relation = case_when(
        is.na(expression_delta) | is.na(activity_score) ~ "weak",
        abs(expression_delta) < weak_delta ~ "weak",
        sign(expression_delta) == sign(activity_score) ~ "consistent",
        TRUE ~ "opposite"
      ),
      label_color = case_when(
        expression_activity_relation == "consistent" ~ "#2E8B57",
        expression_activity_relation == "opposite" ~ "#C43C35",
        TRUE ~ "#C17C00"
      ),
      TF = factor(TF, levels = selected_tfs)
    )

  plot_data <- plot_data %>%
    mutate(TF = factor(TF, levels = selected_tfs))

  label_cols <- annotation$label_color
  names(label_cols) <- as.character(annotation$TF)

  base_group_palette <- c("#2C7FB8", "#F06B2B", "#7A7A7A", "#8E6BBE", "#2CA25F", "#D95F0E")
  group_palette <- setNames(rep(base_group_palette, length.out = length(group_levels)), group_levels)

  p_tf_expression <- ggplot(plot_data, aes(x = TF, y = expression, color = group)) +
    geom_boxplot(aes(fill = group), width = 0.38, alpha = 0.20, outlier.shape = NA,
                 position = position_dodge(width = 0.62), linewidth = 0.35) +
    geom_point(position = position_jitterdodge(jitter.width = 0.06, dodge.width = 0.62),
               size = 1.8, alpha = 0.95) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    scale_fill_manual(values = group_palette, drop = FALSE) +
    labs(
      title = paste0("Top ", length(selected_tfs), " TF Expression: ", treat_group, " vs ", control_group),
      subtitle = "Gene names: Green=consistent, Yellow=weak, Red=opposite (expression vs predicted TF activity)",
      x = NULL,
      y = "Expression (VST normalized)",
      color = "Group",
      fill = "Group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 13),
      plot.subtitle = element_text(size = 8.5, color = "grey25"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold",
                                 color = label_cols[levels(plot_data$TF)]),
      legend.position = "bottom"
    )

  write.csv(annotation,
            file = file.path(tf_boxplot_dir, "Top20_TF_Expression_Boxplot_Annotation.csv"),
            row.names = FALSE)
  ggsave(file.path(tf_boxplot_dir, "TF_Expression_Vertical_Colored_Boxplot.png"),
         p_tf_expression, width = 12, height = 6.2, dpi = 300)
  ggsave(file.path(tf_boxplot_dir, "TF_Expression_Vertical_Colored_Boxplot.pdf"),
         p_tf_expression, width = 12, height = 6.2)
  message("  -> Top TF expression boxplot saved to: ", tf_boxplot_dir)
  invisible(p_tf_expression)
}

make_top_tf_expression_plot(
  contrast_df = contrast_df,
  merged_res = merged_res,
  expr_mat_symbol = expr_mat_symbol,
  metadata = metadata,
  tf_boxplot_dir = tf_boxplot_dir,
  n_top = 20
)

p_concord <- ggplot(merged_res, aes(x = ulm_score, y = gsea_nes, color = Type)) +
  geom_point(alpha = 0.6) + geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE) +
  scale_color_manual(values = c("Consistent" = "firebrick3", "Not Sig" = "grey")) +
  labs(title = "Method Concordance (ulm vs GSEA)", x = "ulm Score", y = "GSEA NES") + theme_bw()
print(p_concord)

# 5.4 Jaccard 热图函数
plot_jaccard <- function(tf_list, title, color_pal) {
  tf_targets <- split(net$target[net$source %in% tf_list], net$source[net$source %in% tf_list])
  if(length(tf_targets) < 2) return()
  n <- length(tf_targets)
  mat <- matrix(0, n, n, dimnames = list(names(tf_targets), names(tf_targets)))
  for (i in 1:n) for (j in i:n) {
    inter <- length(intersect(tf_targets[[i]], tf_targets[[j]]))
    uni <- length(union(tf_targets[[i]], tf_targets[[j]]))
    if (uni > 0) mat[i, j] <- mat[j, i] <- inter / uni
  }
  diag(mat) <- NA; max_val <- max(mat, na.rm = TRUE); diag(mat) <- 1
  if(max_val == 0) max_val <- 0.1
  pheatmap(mat, main = title, color = color_pal, breaks = c(seq(0, max_val, length.out=100), 1), border_color = NA)
}

# 绘制抑制/激活组 Jaccard
inh_tfs <- contrast_df %>% filter(score < 0) %>% arrange(score) %>% head(n_top_jaccard) %>% pull(source)
act_tfs <- contrast_df %>% filter(score > 0) %>% arrange(desc(score)) %>% head(n_top_jaccard) %>% pull(source)
plot_jaccard(inh_tfs, "Jaccard Similarity: Top Inhibited TFs", colorRampPalette(c("white", "aliceblue", "navy"))(101))
plot_jaccard(act_tfs, "Jaccard Similarity: Top Activated TFs", colorRampPalette(c("white", "mistyrose", "firebrick3"))(101))

dev.off()
message("全局报告生成完毕。")


# ==============================================================================
# PART 7: 特定 TF 快速出图模块 (按需运行此部分)
# ==============================================================================
message("\n>>> [Step 7] 开始特定 TF 快速出图...")

for (tf_name in target_tf_list) {
  message(paste0("正在处理 TF: ", tf_name))
  
  # 检查 TF 是否在结果中
  if (!tf_name %in% contrast_df$source) {
    message(paste0("  -> 警告: ", tf_name, " 未在推断结果中找到，可能该 TF 不在 CollecTRI 网络中，跳过。"))
    next
  }
  
  # 准备数据
  safe_tf_name <- gsub("[^A-Za-z0-9_]", "", tf_name)
  tf_out_file <- file.path(single_tf_out_dir, paste0(safe_tf_name, "_Analysis_Report.pdf"))
  
  pdf(tf_out_file, width = 8, height = 6)
  
  # 1. 差异活性条形图 (高亮目标 TF)
  # 取 Top 20 + 目标 TF
  subset_df <- head(contrast_df, 20)
  if (!tf_name %in% subset_df$source) subset_df <- rbind(subset_df, contrast_df[contrast_df$source == tf_name, ])
  subset_df$Highlight <- ifelse(subset_df$source == tf_name, "Target", "Other")
  subset_df$source <- factor(subset_df$source, levels = subset_df$source[order(subset_df$score)])
  
  p1 <- ggplot(subset_df, aes(x = source, y = score, fill = Highlight)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("Target" = "gold", "Other" = "grey70")) +
    coord_flip() + theme_bw() +
    labs(title = paste0("Activity Rank: ", tf_name), x = "TF", y = "ulm Score")
  print(p1)
  
  # 2. GSEA Plot
  if (tf_name %in% gsea_res@result$ID) {
    p2 <- gseaplot2(gsea_res, geneSetID = tf_name, title = paste0("GSEA Enrichment: ", tf_name), pvalue_table = TRUE)
    print(p2)
  } else {
    plot.new()
    text(0.5, 0.5, "GSEA result not found for this TF", cex = 1.2)
  }
  
  # 3. 靶基因表达热图
  # 提取该 TF 的靶基因 (Top 权重)
  targets <- net %>% filter(source == tf_name) %>% arrange(desc(abs(mor))) %>% head(30) %>% pull(target)
  # 修正：使用 expr_mat_symbol (表达矩阵) 而不是 sample_acts_mat (TF活性矩阵)
  valid_targets_expr <- intersect(targets, rownames(expr_mat_symbol))
  
  if (length(valid_targets_expr) > 2) {
    heat_data <- expr_mat_symbol[valid_targets_expr, ]
    anno <- data.frame(Group = metadata$group)
    rownames(anno) <- rownames(metadata)
    pheatmap(heat_data, scale = "row", annotation_col = anno, 
             main = paste0("Expression of Top Targets: ", tf_name), 
             color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
             fontsize_row = 8)
  }
  
  dev.off()
  message(paste0("  -> 报告已保存: ", tf_out_file))
}

message("特定 TF 分析全部完成。")
