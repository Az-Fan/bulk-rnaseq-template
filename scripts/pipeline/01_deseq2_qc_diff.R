# ==============================================================================
# Step01 差异分析与样本 QC（DESeq2）
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch",
#    COUNT_MATRIX_PATH = "data/matrix_gene.count.annot.xls", SHARED_CACHE_DIR = "shared_cache")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================


rm(list = ls())

# 0. 设置输出文件夹结构
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "result_01")
dir.create(output_dir, showWarnings = FALSE)

analysis_mode_tag <- Sys.getenv("ANALYSIS_MODE_TAG", unset = "nobatch")
batch_column <- Sys.getenv("BATCH_COLUMN", unset = "")
deg_padj_threshold <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
deg_lfc_threshold <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
min_count_threshold <- as.integer(Sys.getenv("MIN_COUNT", unset = "10"))
min_samples_threshold <- as.integer(Sys.getenv("MIN_SAMPLES", unset = "3"))

# 可选：自定义输入 count 矩阵路径（用于多项目复用）
count_matrix_path <- Sys.getenv("COUNT_MATRIX_PATH", unset = "data/matrix_gene.count.annot.xls")
if (!file.exists(count_matrix_path)) {
  stop(paste0("Error: 未找到输入矩阵文件: ", count_matrix_path))
}

# 可选：样本分组文件（推荐，适配任意样本命名）
sample_metadata_path <- Sys.getenv("SAMPLE_METADATA_PATH", unset = file.path(dirname(count_matrix_path), "sample_metadata.csv"))

# 可选：指定比较方向（默认兼容旧项目）
control_group <- Sys.getenv("CONTROL_GROUP", unset = "Static EC")
treat_group <- Sys.getenv("TREAT_GROUP", unset = "Proliferating EC")

# 可选：PCA/QC 后人工决定要排除的样本，多个样本用英文逗号分隔
exclude_samples_raw <- Sys.getenv("EXCLUDE_SAMPLES", unset = "")
exclude_samples <- trimws(unlist(strsplit(exclude_samples_raw, ",")))
exclude_samples <- exclude_samples[!is.na(exclude_samples) & exclude_samples != ""]

# 数据处理相关的输出
data_processed_dir <- file.path(output_dir, "data_processed")
dir.create(data_processed_dir, showWarnings = FALSE)

# 差异分析结果输出
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
dir.create(diff_analysis_dir, showWarnings = FALSE)

# 图表输出
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

# 质量控制图子文件夹
qc_plots_dir <- file.path(plots_dir, "qc_plots")
dir.create(qc_plots_dir, showWarnings = FALSE)


library(readxl)
library(data.table)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tibble)
library(limma)
library(pheatmap)
library(AnnotationDbi)
has_orgdb <- requireNamespace("org.Hs.eg.db", quietly = TRUE)
if (has_orgdb) {
  library(org.Hs.eg.db)
}

### 1. 读入表达矩阵
# 使用原始 count 矩阵，兼容 gene_id/id 等列名
cat("读取 count 矩阵：", count_matrix_path, "\n")
expr_raw <- data.table::fread(count_matrix_path, data.table = FALSE)

gene_id_col <- if ("gene_id" %in% colnames(expr_raw)) {
  "gene_id"
} else if ("id" %in% colnames(expr_raw)) {
  "id"
} else {
  colnames(expr_raw)[1]
}
if (gene_id_col != "gene_id") {
  cat("提示：输入矩阵未检测到 gene_id 列，使用列 ", gene_id_col, " 作为基因 ID。\n", sep = "")
}

gene_name_col <- if ("gene_name" %in% colnames(expr_raw)) {
  "gene_name"
} else if ("gene" %in% colnames(expr_raw)) {
  "gene"
} else if ("Symbol" %in% colnames(expr_raw)) {
  "Symbol"
} else if ("symbol" %in% colnames(expr_raw)) {
  "symbol"
} else {
  gene_id_col
}

gene_symbol_map <- expr_raw[, c(gene_id_col, gene_name_col)]
colnames(gene_symbol_map) <- c("gene_id", "gene_symbol")
gene_symbol_map$gene_symbol[is.na(gene_symbol_map$gene_symbol) | gene_symbol_map$gene_symbol == ""] <- gene_symbol_map$gene_id[is.na(gene_symbol_map$gene_symbol) | gene_symbol_map$gene_symbol == ""]
gene_symbol_map <- gene_symbol_map[!duplicated(gene_symbol_map$gene_id), ]
gene_symbol_lookup <- setNames(gene_symbol_map$gene_symbol, gene_symbol_map$gene_id)

metadata_input <- NULL
if (file.exists(sample_metadata_path)) {
  cat("读取样本分组文件：", sample_metadata_path, "\n")
  metadata_input <- read.csv(sample_metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("sample", "group") %in% colnames(metadata_input))) {
    stop("Error: sample_metadata.csv 必须至少包含 sample 和 group 两列。")
  }
  if (length(exclude_samples) > 0) {
    matched_exclude <- intersect(exclude_samples, metadata_input$sample)
    missing_exclude <- setdiff(exclude_samples, metadata_input$sample)
    if (length(matched_exclude) > 0) {
      cat("根据 EXCLUDE_SAMPLES 排除样本：", paste(matched_exclude, collapse = ", "), "\n")
      metadata_input <- metadata_input[!metadata_input$sample %in% matched_exclude, , drop = FALSE]
    }
    if (length(missing_exclude) > 0) {
      cat("警告：EXCLUDE_SAMPLES 中以下样本未在 metadata 中找到：", paste(missing_exclude, collapse = ", "), "\n")
    }
    if (nrow(metadata_input) < 2) {
      stop("Error: 排除样本后 metadata 少于 2 个样本，无法继续分析。")
    }
  }
  count_cols <- unique(metadata_input$sample)
  missing_samples <- setdiff(count_cols, colnames(expr_raw))
  if (length(missing_samples) > 0) {
    stop(paste0("Error: sample_metadata.csv 中以下样本不在 count 矩阵中: ", paste(missing_samples, collapse = ", ")))
  }
} else {
  count_cols <- grep("-(sub|confwe|confwent)$", colnames(expr_raw), value = TRUE, ignore.case = TRUE)
  if (length(count_cols) == 0) {
    count_cols <- grep("_count$", colnames(expr_raw), value = TRUE, ignore.case = TRUE)
  }
  if (length(count_cols) == 0) {
    stop("Error: 未检测到 sample_metadata.csv，且输入矩阵中也未找到可识别的样本列（-sub/-confwe/-confwent 或 *_count）。")
  }
  cat("未检测到 sample_metadata.csv，使用样本列名后缀自动分组。\n")
  if (length(exclude_samples) > 0) {
    matched_exclude <- intersect(exclude_samples, count_cols)
    if (length(matched_exclude) > 0) {
      cat("根据 EXCLUDE_SAMPLES 排除样本：", paste(matched_exclude, collapse = ", "), "\n")
      count_cols <- setdiff(count_cols, matched_exclude)
    }
  }
  if (length(count_cols) < 2) {
    stop("Error: 排除样本后可用样本少于 2 个，无法继续分析。")
  }
}

exprSet <- expr_raw[, c(gene_id_col, count_cols)]
colnames(exprSet) <- sub("confwent$", "confwe", colnames(exprSet), ignore.case = TRUE)
write.csv(exprSet, file.path(data_processed_dir, "counts_sub_confwe.csv"), row.names = FALSE)

### 2. 准备metadata文件
sample_ids <- colnames(exprSet)[-1]

if (!is.null(metadata_input)) {
  metadata <- metadata_input[match(sample_ids, metadata_input$sample), , drop = FALSE]
} else {
  # 自动分组兼容旧模板(-sub/-confwe)和常见 *_count 命名
  if (all(grepl("-(sub|confwe|confwent)$", sample_ids, ignore.case = TRUE))) {
    group <- ifelse(grepl("-sub$", sample_ids, ignore.case = TRUE), "Proliferating EC", "Static EC")
  } else {
    group <- sub("[-_][0-9]+.*$", "", sample_ids)
    group <- sub("_count$", "", group, ignore.case = TRUE)
  }
  metadata <- data.frame(
    sample = sample_ids,
    group = group,
    stringsAsFactors = FALSE
  )
}

present_groups <- unique(as.character(metadata$group))
if (!all(c(control_group, treat_group) %in% present_groups)) {
  if (length(present_groups) == 2) {
    control_group <- present_groups[1]
    treat_group <- present_groups[2]
    cat("警告：未匹配到 CONTROL_GROUP/TREAT_GROUP，自动使用：\n")
    cat("  control_group =", control_group, "\n")
    cat("  treat_group   =", treat_group, "\n")
  } else {
    stop("Error: 当前 group 不是两组，且 CONTROL_GROUP/TREAT_GROUP 未正确指定。")
  }
}

group_levels <- unique(c(control_group, treat_group, present_groups))
metadata$group <- factor(metadata$group, levels = group_levels)
rownames(metadata) <- sample_ids

group_counts <- table(metadata$group)
if (any(group_counts[c(control_group, treat_group)] < 2)) {
  stop("Each comparison group must contain at least two biological replicates.")
}
if (any(group_counts[c(control_group, treat_group)] < 3)) {
  warning("One or more groups contain fewer than three replicates; inference will be fragile.")
}

print("Metadata:")
print(metadata)
save(metadata, file = file.path(data_processed_dir, "metadata.Rdata"))
write.csv(metadata, file = file.path(data_processed_dir, "metadata.csv"), row.names = FALSE)
if (length(exclude_samples) > 0) {
  writeLines(exclude_samples, con = file.path(data_processed_dir, "Excluded_Samples.txt"))
}

# 当前分析模式，供后续批次矫正 PCA 和差异模型共用
batch_mode <- analysis_mode_tag

### 3. 构建dds对象
cat("exprSet样本列：", colnames(exprSet)[-1], "\n")
cat("metadata行名：", rownames(metadata), "\n")

# 确保exprSet的行名为基因ID，并与metadata的样本名匹配
# 假设exprSet的第一列是基因ID
rownames(exprSet) <- exprSet[, 1]
exprSet <- exprSet[, -1]

# 检查样本名一致性
if(!all(colnames(exprSet) == rownames(metadata))) {
  stop("Error: Sample names in exprSet columns and metadata rows do not match or are not in the same order!")
}

if (batch_column == "") {
  design_formula <- ~ group
  model_tag <- analysis_mode_tag
  model_label <- "No Batch Covariate"
} else {
  if (!batch_column %in% colnames(metadata)) {
    stop(paste0("Error: BATCH_COLUMN not found in metadata: ", batch_column))
  }

  batch_values <- metadata[[batch_column]]
  batch_values <- as.character(batch_values)

  if (all(is.na(batch_values) | batch_values == "" | batch_values == "None")) {
    stop(paste0("Error: batch column has no valid values: ", batch_column))
  }

  metadata[[batch_column]] <- factor(batch_values)

  if (nlevels(metadata[[batch_column]]) < 2) {
    stop(paste0("Error: batch column has fewer than 2 levels: ", batch_column))
  }

  design_formula <- reformulate(c(batch_column, "group"))
  model_tag <- analysis_mode_tag
  model_label <- paste0("Batch Covariate: ", batch_column)

  mm <- model.matrix(design_formula, data = metadata)
  if (qr(mm)$rank < ncol(mm)) {
    stop("Error: design matrix is not full rank. The selected batch variable is confounded with group.")
  }
}

save(metadata, file = file.path(data_processed_dir, "metadata.Rdata"))
write.csv(metadata, file = file.path(data_processed_dir, "metadata.csv"), row.names = FALSE)

dds_template <- DESeqDataSetFromMatrix(
  countData = round(exprSet), # DESeq2需要整数计数
  colData = metadata,
  design = design_formula
)

cat("原始dds对象基因数量：", nrow(dds_template), "\n")

# 低表达过滤
filter_before_n <- nrow(dds_template)
filter_rule <- paste0(
  "rowSums(counts(dds_template) >= ", min_count_threshold,
  ") >= ", min_samples_threshold
)
keep <- rowSums(counts(dds_template) >= min_count_threshold) >= min_samples_threshold
dds_template <- dds_template[keep, ]
filter_after_n <- nrow(dds_template)
cat("过滤后dds对象基因数量：", filter_after_n, "\n")
write.csv(
  data.frame(
    genes_before_filtering = filter_before_n,
    genes_after_filtering = filter_after_n,
    filtering_rule = filter_rule,
    stringsAsFactors = FALSE
  ),
  file.path(data_processed_dir, "Gene_Filtering_Summary.csv"),
  row.names = FALSE
)

writeLines(
  c(
    "This workflow starts from a gene-level count matrix.",
    "It does not replace FASTQ-level QC such as per-base quality, adapter content, duplication, alignment/quantification rate, strandedness, rRNA contamination, or sample identity checks.",
    "TF activity, PROGENy, GSVA correlations and PPI analyses are exploratory mechanism-generating analyses and require independent validation.",
    "With small sample sizes, sample-level correlations and network-derived rankings should not be interpreted as confirmatory evidence.",
    paste0("Primary DEG rule: adjusted P < ", deg_padj_threshold,
           " and |ashr-shrunken log2FC| > ", deg_lfc_threshold, ".")
  ),
  file.path(data_processed_dir, "Analysis_Limitations.txt")
)

### 4. 数据质量判断 (PCA)
vsd <- vst(dds_template, blind = FALSE)
exprSet_vst <- as.data.frame(assay(vsd))
save(exprSet_vst, file = file.path(data_processed_dir, "exprSet_vst.Rdata"))

if (batch_column != "") {
  if (!batch_column %in% colnames(metadata)) {
    stop(paste0("Error: BATCH_COLUMN not found in metadata: ", batch_column))
  }
  pca_intgroup <- c("group", batch_column)
} else {
  pca_intgroup <- "group"
}

pcaData <- plotPCA(vsd, intgroup = pca_intgroup, returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData$shape_group <- if (batch_column != "") pcaData[[batch_column]] else pcaData$group

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color = group, shape = shape_group, label = name)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, check_overlap = TRUE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot with Sample Names")

print(pca_plot)
ggsave(file.path(qc_plots_dir, "PCA_plot_initial.png"), pca_plot, width = 8, height = 6)
ggsave(file.path(qc_plots_dir, "PCA_plot_initial.pdf"), pca_plot, width = 8, height = 6)

# 4.1 批次矫正后 PCA（保留 group 生物学差异）
batch_var <- batch_column
batch_label <- if (batch_column != "") batch_column else "No Batch"
shape_col <- if (batch_column != "") metadata[[batch_column]] else metadata$group

if (batch_column != "") {
  design_keep_group <- model.matrix(~ group, data = metadata)
  exprSet_vst_batch_corrected <- limma::removeBatchEffect(
    as.matrix(exprSet_vst),
    batch = metadata[[batch_var]],
    design = design_keep_group
  )
  corrected_title <- paste0("PCA After Batch Correction (", batch_label, ")")
} else {
  cat("提示：当前为 no-batch 模式，跳过批次矫正并直接输出原始 VST PCA。\n")
  exprSet_vst_batch_corrected <- as.matrix(exprSet_vst)
  corrected_title <- "PCA (No Batch Correction Applied)"
}
save(
  exprSet_vst_batch_corrected,
  file = file.path(data_processed_dir, paste0("exprSet_vst_batch_corrected_", batch_mode, ".Rdata"))
)

pca_corrected <- prcomp(t(exprSet_vst_batch_corrected), center = TRUE, scale. = FALSE)
pc_var <- (pca_corrected$sdev^2) / sum(pca_corrected$sdev^2)
pca_corrected_df <- data.frame(
  sample = rownames(metadata),
  group = metadata$group,
  shape_group = shape_col,
  PC1 = pca_corrected$x[, 1],
  PC2 = pca_corrected$x[, 2],
  stringsAsFactors = FALSE
)

pca_corrected_plot <- ggplot(
  pca_corrected_df,
  aes(x = PC1, y = PC2, color = group, shape = shape_group, label = sample)
) +
  geom_point(size = 3) +
  geom_text(vjust = -1, check_overlap = TRUE) +
  xlab(paste0("PC1: ", round(pc_var[1] * 100), "% variance")) +
  ylab(paste0("PC2: ", round(pc_var[2] * 100), "% variance")) +
  theme_bw() +
  ggtitle(corrected_title)

print(pca_corrected_plot)
ggsave(
  file.path(qc_plots_dir, paste0("PCA_plot_batch_corrected_", batch_mode, ".png")),
  pca_corrected_plot,
  width = 8,
  height = 6
)
ggsave(
  file.path(qc_plots_dir, paste0("PCA_plot_batch_corrected_", batch_mode, ".pdf")),
  pca_corrected_plot,
  width = 8,
  height = 6
)

### 4.2 样本层面 QC（聚类、相关性、测序深度）
sample_qc_df <- data.frame(
  sample = colnames(exprSet),
  group = as.character(metadata$group),
  batch_column = if (batch_column != "") batch_column else "",
  batch_value = if (batch_column != "") as.character(metadata[[batch_column]]) else "",
  library_size = colSums(exprSet),
  detected_genes = colSums(exprSet > 0),
  zero_fraction = round(colMeans(exprSet == 0), 4),
  stringsAsFactors = FALSE
)

sample_dist_mat <- as.matrix(dist(t(exprSet_vst_batch_corrected), method = "euclidean"))
corr_mat <- cor(exprSet_vst_batch_corrected, method = "pearson")

corr_no_diag <- corr_mat
diag(corr_no_diag) <- NA_real_
sample_qc_df$mean_pearson_corr <- round(colMeans(corr_no_diag, na.rm = TRUE), 4)

corr_mean <- mean(sample_qc_df$mean_pearson_corr, na.rm = TRUE)
corr_sd <- sd(sample_qc_df$mean_pearson_corr, na.rm = TRUE)
corr_threshold <- corr_mean - 2 * corr_sd
if (!is.finite(corr_threshold)) {
  corr_threshold <- NA_real_
}
sample_qc_df$low_corr_flag <- if (is.finite(corr_threshold)) sample_qc_df$mean_pearson_corr < corr_threshold else FALSE
sample_qc_df$qc_flag <- ifelse(sample_qc_df$low_corr_flag, "CHECK", "PASS")

write.csv(sample_qc_df, file.path(data_processed_dir, "Sample_QC_Metrics.csv"), row.names = FALSE)
write.csv(sample_dist_mat, file.path(data_processed_dir, "Sample_Distance_Euclidean.csv"))
write.csv(corr_mat, file.path(data_processed_dir, "Sample_Correlation_Pearson.csv"))

if (any(sample_qc_df$low_corr_flag)) {
  flagged_samples <- sample_qc_df$sample[sample_qc_df$low_corr_flag]
  cat("警告：以下样本与其他样本平均相关性偏低，请重点检查：", paste(flagged_samples, collapse = ", "), "\n")
}

p_libsize <- ggplot(sample_qc_df, aes(x = reorder(sample, library_size), y = library_size, fill = group)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(title = "Library Size by Sample", x = NULL, y = "Raw Count Sum")

ggsave(file.path(qc_plots_dir, "Library_Size_Barplot.png"), p_libsize, width = 8, height = 5)
ggsave(file.path(qc_plots_dir, "Library_Size_Barplot.pdf"), p_libsize, width = 8, height = 5)

p_corr <- ggplot(sample_qc_df, aes(x = reorder(sample, mean_pearson_corr), y = mean_pearson_corr, fill = qc_flag)) +
  geom_col(width = 0.7) +
  coord_flip() +
  theme_bw() +
  labs(title = "Mean Pearson Correlation by Sample", x = NULL, y = "Mean Pearson Correlation") +
  scale_fill_manual(values = c("PASS" = "#4DAF4A", "CHECK" = "#E41A1C"))

if (is.finite(corr_threshold)) {
  p_corr <- p_corr + geom_hline(yintercept = corr_threshold, linetype = 2, color = "grey25")
}

ggsave(file.path(qc_plots_dir, "Sample_Correlation_QC_Barplot.png"), p_corr, width = 8, height = 5)
ggsave(file.path(qc_plots_dir, "Sample_Correlation_QC_Barplot.pdf"), p_corr, width = 8, height = 5)

if (batch_column != "") {
  sample_annotation <- data.frame(
    group = metadata$group,
    batch = metadata[[batch_var]],
    row.names = rownames(metadata),
    stringsAsFactors = FALSE
  )
  colnames(sample_annotation)[colnames(sample_annotation) == "batch"] <- batch_column
} else {
  sample_annotation <- data.frame(
    group = metadata$group,
    row.names = rownames(metadata),
    stringsAsFactors = FALSE
  )
}

dist_palette <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(100)
corr_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101)

pheatmap::pheatmap(
  sample_dist_mat,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  color = dist_palette,
  main = "Sample Distance Heatmap (Euclidean, VST)",
  filename = file.path(qc_plots_dir, "Sample_Distance_Heatmap.png"),
  width = 8,
  height = 7
)
pheatmap::pheatmap(
  sample_dist_mat,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  color = dist_palette,
  main = "Sample Distance Heatmap (Euclidean, VST)",
  filename = file.path(qc_plots_dir, "Sample_Distance_Heatmap.pdf"),
  width = 8,
  height = 7
)

corr_dist <- as.dist(1 - corr_mat)
pheatmap::pheatmap(
  corr_mat,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  color = corr_palette,
  breaks = seq(-1, 1, length.out = 102),
  clustering_distance_rows = corr_dist,
  clustering_distance_cols = corr_dist,
  main = "Sample Correlation Heatmap (Pearson, VST)",
  filename = file.path(qc_plots_dir, "Sample_Correlation_Heatmap.png"),
  width = 8,
  height = 7
)
pheatmap::pheatmap(
  corr_mat,
  annotation_col = sample_annotation,
  annotation_row = sample_annotation,
  color = corr_palette,
  breaks = seq(-1, 1, length.out = 102),
  clustering_distance_rows = corr_dist,
  clustering_distance_cols = corr_dist,
  main = "Sample Correlation Heatmap (Pearson, VST)",
  filename = file.path(qc_plots_dir, "Sample_Correlation_Heatmap.pdf"),
  width = 8,
  height = 7
)

sample_hclust <- hclust(as.dist(sample_dist_mat), method = "average")
png(file.path(qc_plots_dir, "Sample_Hierarchical_Clustering.png"), width = 1100, height = 700, res = 130)
plot(sample_hclust, main = "Sample Hierarchical Clustering (Euclidean, VST)", xlab = "", sub = "")
dev.off()
pdf(file.path(qc_plots_dir, "Sample_Hierarchical_Clustering.pdf"), width = 10, height = 6)
plot(sample_hclust, main = "Sample Hierarchical Clustering (Euclidean, VST)", xlab = "", sub = "")
dev.off()

sample_cluster_order <- data.frame(
  order = seq_along(sample_hclust$order),
  sample = sample_hclust$labels[sample_hclust$order],
  stringsAsFactors = FALSE
)
write.csv(sample_cluster_order, file.path(data_processed_dir, "Sample_Clustering_Order.csv"), row.names = FALSE)

pipeline_phase <- Sys.getenv("PIPELINE_PHASE", unset = "differential")
if (identical(pipeline_phase, "qc")) {
  review_lines <- c(
    "STAGE 01 QC REVIEW REQUIRED",
    "",
    "Review PCA, sample correlation, sample distance, clustering and library size.",
    "Decide whether a batch column should be included in the final design.",
    "Decide whether any sample should be excluded as an outlier.",
    "",
    paste0(
      "QC-flagged samples: ",
      if (any(sample_qc_df$qc_flag == "CHECK")) {
        paste(sample_qc_df$sample[sample_qc_df$qc_flag == "CHECK"], collapse = ", ")
      } else {
        "none"
      }
    ),
    "",
    "Before Stage 02:",
    "1. Update BatchColumn and ExcludeSamples in the selected project config.",
    "2. Confirm the comparison direction.",
    "3. Create config/stage_checkpoint.json with stage=02_differential and qc_reviewed=true."
  )
  writeLines(review_lines, file.path(data_processed_dir, "QC_REVIEW_REQUIRED.txt"))
  cat("\nStage 01 QC completed. Differential-expression modelling was intentionally not run.\n")
  quit(save = "no", status = 0)
}

annotate_res <- function(res_df) {
  if (has_orgdb) {
    gene_id_sample <- head(res_df$gene_id, 10)
    cat("\n前10个基因ID：\n")
    print(gene_id_sample)

    keytype_used <- "UNKNOWN"
    if(all(grepl("^ENSG", gene_id_sample))) {
      keytype_used <- "ENSEMBL"
    } else if(all(grepl("^ENST", gene_id_sample))) {
      keytype_used <- "ENSEMBLTRANS"
    } else if(all(grepl("^[0-9]+$", gene_id_sample))) {
      keytype_used <- "ENTREZID"
    } else if (all(grepl("^[A-Z0-9-]+$", gene_id_sample))) {
      keytype_used <- "SYMBOL"
    } else {
      cat("警告：未能自动识别keytype，默认为'SYMBOL'。\n")
      keytype_used <- "SYMBOL"
    }

    cat("检测到的keytype：", keytype_used, "\n")
    res_df$symbol <- mapIds(
      org.Hs.eg.db,
      keys = res_df$gene_id,
      column = "SYMBOL",
      keytype = keytype_used,
      multiVals = "first"
    )
    res_df$entrez <- mapIds(
      org.Hs.eg.db,
      keys = res_df$gene_id,
      column = "ENTREZID",
      keytype = keytype_used,
      multiVals = "first"
    )
    cat("注释完成，主结果保留所有 tested genes：", nrow(res_df), "\n")
  } else {
    cat("警告：未安装org.Hs.eg.db，跳过注释。\n")
    res_df$symbol <- NA_character_
    res_df$entrez <- NA_character_
  }

  mapped_symbol <- gene_symbol_lookup[res_df$gene_id]
  fallback_symbol <- unname(mapped_symbol)
  fallback_symbol[is.na(fallback_symbol) | fallback_symbol == ""] <- res_df$gene_id[is.na(fallback_symbol) | fallback_symbol == ""]
  res_df$symbol[is.na(res_df$symbol) | res_df$symbol == ""] <- fallback_symbol[is.na(res_df$symbol) | res_df$symbol == ""]
  res_df$ENTREZID <- res_df$entrez

  if ("baseMean" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "baseMean"] <- "AveExpr"
  }
  colnames(res_df)[colnames(res_df) == "pvalue"] <- "P.Value"
  colnames(res_df)[colnames(res_df) == "padj"] <- "adj.P.Val"
  colnames(res_df)[colnames(res_df) == "symbol"] <- "gene"
  res_df
}

run_deseq_model <- function(dds_in, design_formula, model_tag, model_label) {
  dds_model <- dds_in
  design(dds_model) <- design_formula
  dds_model <- DESeq(dds_model)

  contrast <- c("group", treat_group, control_group)
  dd1 <- results(dds_model, contrast = contrast, alpha = deg_padj_threshold)

  png(file.path(qc_plots_dir, paste0("MA_plot_raw_", model_tag, ".png")), width = 800, height = 600, res = 100)
  plotMA(dd1, ylim = c(-5, 5), main = paste0("MA Plot (Raw) - ", model_label))
  dev.off()

  dd2 <- lfcShrink(dds_model, contrast = contrast, res = dd1, type = "ashr")

  png(file.path(qc_plots_dir, paste0("MA_plot_shrunk_", model_tag, ".png")), width = 800, height = 600, res = 100)
  plotMA(dd2, ylim = c(-5, 5), main = paste0("MA Plot (LFC Shrunk) - ", model_label))
  dev.off()

  raw_res <- as.data.frame(dd1) %>%
    rownames_to_column("gene_id") %>%
    dplyr::transmute(
      gene_id,
      baseMean,
      logFC_raw = log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj
    )

  shrunk_res <- as.data.frame(dd2) %>%
    rownames_to_column("gene_id") %>%
    dplyr::transmute(
      gene_id,
      logFC = log2FoldChange
    )

  res <- dplyr::left_join(raw_res, shrunk_res, by = "gene_id")
  res_annotated <- annotate_res(res)
  save(res_annotated, file = file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata"))
  write.csv(res_annotated, file = file.path(diff_analysis_dir, "DEG_results_annotated.csv"), row.names = FALSE)

  significant_deg <- res_annotated[
    !is.na(res_annotated$adj.P.Val) &
      res_annotated$adj.P.Val < deg_padj_threshold &
      abs(res_annotated$logFC) > deg_lfc_threshold,
    ,
    drop = FALSE
  ]
  write.csv(
    significant_deg,
    file = file.path(diff_analysis_dir, "DEG_results_significant_ashr.csv"),
    row.names = FALSE
  )
  sig_n <- nrow(significant_deg)
  cat("\n", model_label, "分析完成！\n", sep = "")
  cat("差异基因数量 (adj.P.Val < 0.05 & |logFC| > 1)：", sig_n, "\n")

  data.frame(
    analysis_mode_tag = analysis_mode_tag,
    batch_column = batch_column,
    model = model_label,
    design = paste(deparse(design_formula), collapse = ""),
    design_formula = paste(deparse(design_formula), collapse = ""),
    lfc_shrinkage = "ashr",
    deg_rule = paste0(
      "adj.P.Val < ", deg_padj_threshold,
      " & abs(ashr_shrunken_logFC) > ", deg_lfc_threshold
    ),
    n_samples = nrow(metadata),
    n_control = sum(metadata$group == control_group),
    n_treat = sum(metadata$group == treat_group),
    genes_tested = nrow(res_annotated),
    n_sig_adjP05_abslogFC1 = sig_n,
    stringsAsFactors = FALSE
  )
}

### 5. 主模型运行
summary_df <- run_deseq_model(dds_template, design_formula, model_tag, model_label)

write.csv(summary_df, file.path(diff_analysis_dir, "DEG_model_summary.csv"), row.names = FALSE)

analysis_config <- data.frame(
  count_matrix_path = count_matrix_path,
  sample_metadata_path = sample_metadata_path,
  control_group = control_group,
  treat_group = treat_group,
  analysis_mode_tag = analysis_mode_tag,
  batch_column = batch_column,
  design_formula = paste(deparse(design_formula), collapse = ""),
  exclude_samples = paste(exclude_samples, collapse = ","),
  samples_used = paste(rownames(metadata), collapse = ","),
  n_samples = nrow(metadata),
  n_control = sum(metadata$group == control_group),
  n_treat = sum(metadata$group == treat_group),
  gene_filter_rule = filter_rule,
  genes_before_filtering = filter_before_n,
  genes_after_filtering = filter_after_n,
  stringsAsFactors = FALSE
)

write.csv(
  analysis_config,
  file.path(data_processed_dir, "Analysis_Config_Used.csv"),
  row.names = FALSE
)

cat("\n当前模型分析完成，汇总文件：", file.path(diff_analysis_dir, "DEG_model_summary.csv"), "\n")
