# ==============================================================================
# Step03 绘制火山图/MA图/DEG热图
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())
library(dplyr)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(RColorBrewer) # 用于生成更美观的颜色板

# 0. 设置输出文件夹结构
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results_01")
dir.create(output_dir, showWarnings = FALSE)

# 数据处理相关的输出 (用于加载数据)
data_processed_dir <- file.path(output_dir, "data_processed")
dir.create(data_processed_dir, showWarnings = FALSE)

# 差异分析结果输出 (用于加载数据)
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
dir.create(diff_analysis_dir, showWarnings = FALSE)

# 图表输出的父目录
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# 质量控制图子文件夹 (MA图和PCA图会放在这里)
qc_plots_dir <- file.path(plots_dir, "qc_plots")
dir.create(qc_plots_dir, showWarnings = FALSE)

# 火山图子文件夹
volcano_plots_dir <- file.path(plots_dir, "volcano_plots")
dir.create(volcano_plots_dir, showWarnings = FALSE)

# 热图子文件夹
heatmaps_dir <- file.path(plots_dir, "heatmaps")
dir.create(heatmaps_dir, showWarnings = FALSE)


# ==============================================================================
# PART 1: 全局参数设置
# ==============================================================================
# --- 1. 阈值设置 ---
param_fc_cutoff   <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
param_p_cutoff    <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
param_label_fc    <- 2      # 标签阈值 (log2FoldChange)
param_top_n       <- as.integer(Sys.getenv("VOLCANO_LABEL_COUNT", unset = "20"))
param_heatmap_top <- as.integer(Sys.getenv("HEATMAP_TOP_PER_DIRECTION", unset = "20"))

# --- 2. 颜色设置 ---
col_volcano_up    <- "firebrick3" # 与热图红端一致
col_volcano_down  <- "navy"       # 与热图蓝端一致
col_volcano_ns    <- "grey75"     # 中性灰
hm_color_grad     <- rev(brewer.pal(n = 7, name = "RdBu")) # 更专业的红蓝渐变色板

# --- 3. 热图尺寸设置 ---
hm_cell_w         <- 25
hm_cell_h_text    <- 11

# ==============================================================================
# PART 2: 数据处理
# ==============================================================================
# 1. 加载数据
# 注意：这里加载exprSet_heatmap时，使用上一个脚本保存的文件
load(file = file.path(data_processed_dir, "exprSet_for_heatmap.Rdata"))
load(file = file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata")) # 加载重命名后的res_annotated
load(file = file.path(data_processed_dir, "metadata.Rdata"))

# 将加载的res_annotated重新命名为res，以匹配原始脚本逻辑
res <- res_annotated

# 2. 准备热图数据
# 筛选差异基因，用于热图
diffgene <- res %>%
  filter(gene != "" & !is.na(gene)) %>% # 确保基因Symbol有效
  filter(adj.P.Val < param_p_cutoff) %>%
  filter(abs(logFC) > param_fc_cutoff)

if(nrow(diffgene) < 2) {
  warning("差异基因太少 (<2)，无法绘制热图。请检查阈值或数据。\n")
  # 可以选择在这里退出脚本，或者继续，但不生成热图
} else {
  cat(paste0("检测到 ", nrow(diffgene), " 个差异基因用于热图绘制。\n"))
}


# 确保 exprSet_heatmap 中的基因与 diffgene 匹配，并按 diffgene 的顺序
# 筛选出差异基因的表达数据
available_diffgene <- diffgene %>%
  filter(gene %in% rownames(exprSet_heatmap)) %>%
  distinct(gene, .keep_all = TRUE)

top_up <- available_diffgene %>%
  filter(logFC > 0) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice_head(n = param_heatmap_top) %>%
  mutate(heatmap_direction = "Up")

top_down <- available_diffgene %>%
  filter(logFC < 0) %>%
  arrange(adj.P.Val, desc(abs(logFC))) %>%
  slice_head(n = param_heatmap_top) %>%
  mutate(heatmap_direction = "Down")

heatmap_selection <- bind_rows(top_up, top_down)
heat_genes <- heatmap_selection$gene
heatdata <- exprSet_heatmap[heat_genes, , drop = FALSE]
write.csv(
  heatmap_selection,
  file.path(diff_analysis_dir, "Heatmap_Top_Genes.csv"),
  row.names = FALSE
)

# 确保样本顺序与 metadata 一致
# metadata 已经经过排序，所以这里只需确保 annotation_col 的行名与 heatdata 的列名一致
annotation_col <- data.frame(group = metadata$group)
rownames(annotation_col) <- metadata$sample
annotation_col <- annotation_col[colnames(heatdata), , drop = FALSE] # 确保样本顺序一致

# 3. 准备火山图/MA图数据
# 使用差异分析结果 'res' 作为基础数据
data_plot <- res
data_plot$group <- "NS"
data_plot$group[data_plot$adj.P.Val < param_p_cutoff & data_plot$logFC > param_fc_cutoff] <- "Up"
data_plot$group[data_plot$adj.P.Val < param_p_cutoff & data_plot$logFC < -param_fc_cutoff] <- "Down"
data_plot$group <- factor(data_plot$group, levels = c("Up", "NS", "Down"))

# 避免 adj.P.Val 为 0 导致 -log10(0)=Inf，造成点“顶在边界”
padj_for_plot <- data_plot$adj.P.Val
min_nonzero_p <- suppressWarnings(min(padj_for_plot[is.finite(padj_for_plot) & padj_for_plot > 0], na.rm = TRUE))
if (!is.finite(min_nonzero_p)) min_nonzero_p <- 1e-300
padj_floor <- max(min_nonzero_p / 10, 1e-300)
padj_for_plot[!is.finite(padj_for_plot) | is.na(padj_for_plot) | padj_for_plot <= 0] <- padj_floor

data_plot$neglog10_p <- -log10(padj_for_plot)
y_cap <- suppressWarnings(as.numeric(quantile(data_plot$neglog10_p, probs = 0.999, na.rm = TRUE)))
if (!is.finite(y_cap)) y_cap <- max(data_plot$neglog10_p, na.rm = TRUE)
if (!is.finite(y_cap)) y_cap <- 50
data_plot$is_y_capped <- data_plot$neglog10_p > y_cap
data_plot$neglog10_p_plot <- pmin(data_plot$neglog10_p, y_cap)
volcano_subtitle <- NULL
if (any(data_plot$is_y_capped, na.rm = TRUE)) {
  volcano_subtitle <- paste0("Top 0.1% y-values capped at ", round(y_cap, 2))
}

# 筛选用于标签的基因
top_genes <- data_plot %>%
  filter(abs(logFC) >= param_label_fc & adj.P.Val < param_p_cutoff) %>%
  arrange(adj.P.Val) %>%
  head(param_top_n)

# ==============================================================================
# PART 3: 绘图对象构建 (火山图和MA图)
# ==============================================================================
# 火山图对象
p_volcano <- ggplot(data_plot, aes(x = logFC, y = neglog10_p_plot, color = group)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = col_volcano_up, "NS" = col_volcano_ns, "Down" = col_volcano_down)) +
  geom_vline(xintercept = c(-param_fc_cutoff, param_fc_cutoff), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = -log10(param_p_cutoff), lty = 4, col = "black", lwd = 0.5) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15))) +
  labs(
    x = "log2 (Fold Change)",
    y = "-log10 (adj.P.Val)",
    title = "Volcano Plot",
    subtitle = volcano_subtitle
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        panel.border = element_blank(), axis.line = element_line(color = "black", linewidth = 0.6),
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
  coord_cartesian(clip = "off") +
  # 对被截顶点用三角形提示
  geom_point(
    data = subset(data_plot, is_y_capped),
    aes(x = logFC, y = neglog10_p_plot, fill = group),
    shape = 24, size = 2.2, stroke = 0.2, color = "black", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c("Up" = col_volcano_up, "NS" = col_volcano_ns, "Down" = col_volcano_down), guide = "none") +
  # 标签
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, color = "black", show.legend = FALSE,
                  min.segment.length = 0.2, box.padding = 0.6, max.overlaps = Inf,
                  segment.color = "grey50")

# MA图对象 (放在qc_plots_dir)
p_ma <- ggplot(data = data_plot, aes(x = log2(AveExpr), y = logFC, color = group)) +
  geom_point(alpha = 0.6, size = 1.2) +
  scale_color_manual(values = c("Up" = col_volcano_up, "NS" = col_volcano_ns, "Down" = col_volcano_down)) +
  labs(y = "log2 (Fold Change)", x = "log2 (Base Mean)", title = "MA Plot") +
  geom_hline(yintercept = c(param_fc_cutoff, -param_fc_cutoff), lty = 4, col = "black", lwd = 0.5) +
  geom_hline(yintercept = 0, col = "black", lwd = 0.5) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "top", plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14), axis.text = element_text(size = 12),
        panel.border = element_blank(), axis.line = element_line(color = "black", linewidth = 0.6),
        plot.margin = margin(t = 20, r = 20, b = 10, l = 10, unit = "pt")) +
  coord_cartesian(clip = "off") +
  # 标签
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, color = "black", show.legend = FALSE,
                  min.segment.length = 0.2, box.padding = 0.6, max.overlaps = Inf,
                  segment.color = "grey50")

# ==============================================================================
# PART 4: 输出图表
# ==============================================================================

# --- 文件1：火山图 (单独保存) ---
file_volcano_plot_png <- file.path(volcano_plots_dir, "Volcano_Plot_DEGs.png")
ggsave(filename = file_volcano_plot_png, plot = p_volcano, width = 8, height = 7, dpi = 300)
cat(paste0(">>> 火山图已生成: ", file_volcano_plot_png, "\n"))
file_volcano_plot_pdf <- file.path(volcano_plots_dir, "Volcano_Plot_DEGs.pdf")
ggsave(filename = file_volcano_plot_pdf, plot = p_volcano, width = 8, height = 7)
cat(paste0(">>> 火山图 (PDF) 已生成: ", file_volcano_plot_pdf, "\n"))

# --- 文件2：MA图 (单独保存) ---
file_ma_plot_png <- file.path(qc_plots_dir, "MA_Plot_DEGs.png")
ggsave(filename = file_ma_plot_png, plot = p_ma, width = 8, height = 7, dpi = 300)
cat(paste0(">>> MA图已生成: ", file_ma_plot_png, "\n"))
file_ma_plot_pdf <- file.path(qc_plots_dir, "MA_Plot_DEGs.pdf")
ggsave(filename = file_ma_plot_pdf, plot = p_ma, width = 8, height = 7)
cat(paste0(">>> MA图 (PDF) 已生成: ", file_ma_plot_pdf, "\n"))


# --- 文件3：差异基因热图（上调 Top N + 下调 Top N） ---
old_heatmaps <- c(
  "Heatmap_DEGs_Overview.png", "Heatmap_DEGs_Overview.md",
  "Heatmap_DEGs_Detail.pdf", "Heatmap_DEGs_Detail.md"
)
unlink(file.path(heatmaps_dir, old_heatmaps), force = TRUE)

if (nrow(heatdata) >= 2) {
  file_heatmap_png <- file.path(
    heatmaps_dir,
    paste0("Heatmap_DEGs_Top", param_heatmap_top, "_Up_Down.png")
  )
  dynamic_height <- max(7.5, min(13, 3.8 + nrow(heatdata) * 0.19))
  pheatmap(
    heatdata,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_col,
    annotation_legend = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    scale = "row",
    color = colorRampPalette(hm_color_grad)(100),
    border_color = NA,
    cellwidth = hm_cell_w,
    cellheight = hm_cell_h_text,
    fontsize = 10,
    fontsize_row = 8,
    fontsize_col = 9,
    main = paste0(
      "Differential-expression heatmap (Top ",
      nrow(top_up), " Up + Top ", nrow(top_down), " Down)"
    ),
    filename = file_heatmap_png,
    width = 9,
    height = dynamic_height
  )
  cat(paste0(">>> Top Up/Down DEG heatmap generated: ", file_heatmap_png, "\n"))
} else {
  cat("因差异基因数量不足，未生成热图。\n")
}

cat("\n所有绘图任务完成！\n")
