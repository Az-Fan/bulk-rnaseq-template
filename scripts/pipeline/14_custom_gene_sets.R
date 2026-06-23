# ==============================================================================
# Step14 自定义基因集（resources/gene_sets/custom_gene_sets.xlsx）热图 + GSEA
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. 输出目录
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
dir.create(output_dir, showWarnings = FALSE)

custom_set_dir <- file.path(output_dir, "functional_enrichment", "custom_gene_sets")
heatmap_dir <- file.path(custom_set_dir, "heatmaps")
gsea_dir <- file.path(custom_set_dir, "gsea")
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 1. 依赖与输入
# ==============================================================================
required_pkgs <- c("readxl", "dplyr", "stringr", "pheatmap", "clusterProfiler", "enrichplot", "ggplot2", "BiocParallel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0("缺少必要包: ", paste(missing_pkgs, collapse = ", "), "。请先安装后再运行。"))
}

library(readxl)
library(dplyr)
library(stringr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(BiocParallel)

gene_set_file <- file.path("resources", "gene_sets", "custom_gene_sets.xlsx")
if (!file.exists(gene_set_file)) {
  stop(paste0("未找到自定义基因集文件: ", gene_set_file))
}

load(file.path(output_dir, "data_processed", "exprSet_for_heatmap.Rdata"))  # exprSet_heatmap
load(file.path(output_dir, "data_processed", "metadata.Rdata"))              # metadata
load(file.path(output_dir, "diff_analysis", "DEseq2_Diff_Annotated.Rdata")) # res_annotated

if (!exists("exprSet_heatmap")) {
  stop("未找到 exprSet_heatmap，请先运行 Step02。")
}
if (!exists("metadata")) {
  stop("未找到 metadata，请先运行 Step01。")
}
if (!exists("res_annotated")) {
  stop("未找到 res_annotated，请先运行 Step01。")
}

# ==============================================================================
# 2. 读取并清洗 resources/gene_sets/custom_gene_sets.xlsx
# ==============================================================================
raw_tbl <- readxl::read_excel(gene_set_file, sheet = 1, col_names = FALSE)
colnames(raw_tbl) <- c("set_name", "genes")
raw_tbl <- raw_tbl %>%
  filter(!is.na(set_name), !is.na(genes), set_name != "", genes != "")

if (nrow(raw_tbl) == 0) {
  stop("gene_set.xlsx 第1个 sheet 为空或格式不正确（应为2列：基因集名 + 逗号分隔基因）。")
}

parse_genes <- function(x) {
  y <- unlist(strsplit(x, ","))
  y <- stringr::str_trim(y)
  y <- y[y != ""]
  unique(y)
}

custom_sets <- setNames(lapply(raw_tbl$genes, parse_genes), raw_tbl$set_name)
custom_sets <- custom_sets[vapply(custom_sets, length, integer(1)) >= 2]

if (length(custom_sets) == 0) {
  stop("gene_set.xlsx 中未解析出有效基因集（每个基因集至少需要2个基因）。")
}

# TERM2GENE 格式（供 GSEA 使用）
term2gene <- bind_rows(lapply(names(custom_sets), function(nm) {
  data.frame(term = nm, gene = custom_sets[[nm]], stringsAsFactors = FALSE)
}))

write.csv(term2gene, file.path(custom_set_dir, "custom_gene_sets_term2gene.csv"), row.names = FALSE)

# ==============================================================================
# 3. 基因集热图
# ==============================================================================
annotation_col <- data.frame(group = metadata$group)
rownames(annotation_col) <- metadata$sample

expr_gene_upper <- toupper(rownames(exprSet_heatmap))
row_lookup <- setNames(rownames(exprSet_heatmap), expr_gene_upper)

for (set_name in names(custom_sets)) {
  gs <- custom_sets[[set_name]]
  matched <- unname(row_lookup[toupper(gs)])
  matched <- unique(matched[!is.na(matched)])

  if (length(matched) < 2) {
    message("热图跳过: ", set_name, "（匹配基因 < 2）")
    next
  }

  heatdata <- exprSet_heatmap[matched, , drop = FALSE]
  annotation_use <- annotation_col[colnames(heatdata), , drop = FALSE]
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", set_name)
  png_path <- file.path(heatmap_dir, paste0("Heatmap_", safe_name, ".png"))
  pdf_path <- file.path(heatmap_dir, paste0("Heatmap_", safe_name, ".pdf"))
  pdf_h <- max(6, min(16, 4 + nrow(heatdata) * 0.18))

  pheatmap(
    heatdata,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_use,
    show_rownames = TRUE,
    show_colnames = FALSE,
    scale = "row",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    border_color = NA,
    main = paste0("Custom Set: ", set_name),
    filename = png_path,
    width = 9,
    height = pdf_h
  )

  pdf(pdf_path, width = 9, height = pdf_h)
  pheatmap(
    heatdata,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = annotation_use,
    show_rownames = TRUE,
    show_colnames = FALSE,
    scale = "row",
    color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    border_color = NA,
    main = paste0("Custom Set: ", set_name)
  )
  dev.off()
}

# ==============================================================================
# 4. 基因集 GSEA（针对所有 custom sets）
# ==============================================================================
rank_df <- as.data.frame(res_annotated) %>%
  dplyr::select(gene, logFC) %>%
  filter(!is.na(gene), gene != "", !is.na(logFC)) %>%
  distinct(gene, .keep_all = TRUE)

geneList <- rank_df$logFC
names(geneList) <- rank_df$gene
geneList <- sort(geneList, decreasing = TRUE)

gsea_res <- tryCatch(
  GSEA(
    geneList = geneList,
    TERM2GENE = term2gene,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    minGSSize = 3,
    maxGSSize = 500,
    BPPARAM = BiocParallel::SerialParam(),
    seed = 123
  ),
  error = function(e) {
    message("自定义基因集 GSEA 失败：", e$message)
    NULL
  }
)

if (!is.null(gsea_res) && nrow(gsea_res@result) > 0) {
  gsea_tbl <- gsea_res@result %>%
    arrange(p.adjust, desc(abs(NES)))
  write.csv(gsea_tbl, file.path(gsea_dir, "GSEA_CustomGeneSets_FullTable.csv"), row.names = FALSE)

  # 汇总气泡图
  dot_obj <- gsea_res
  dot_obj@result <- dot_obj@result %>% filter(p.adjust < 0.25)
  if (nrow(dot_obj@result) > 0) {
    p_dot <- dotplot(dot_obj, showCategory = min(20, nrow(dot_obj@result)), split = ".sign", label_format = 60) +
      facet_grid(. ~ .sign) +
      ggtitle("Custom Gene Sets GSEA") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(gsea_dir, "GSEA_CustomGeneSets_Dotplot.png"), p_dot, width = 10, height = 7, dpi = 300)
    ggsave(file.path(gsea_dir, "GSEA_CustomGeneSets_Dotplot.pdf"), p_dot, width = 10, height = 7)
  }

  # 每个基因集输出单条 GSEA 曲线图（若存在结果）
  for (set_name in names(custom_sets)) {
    hit <- gsea_tbl[gsea_tbl$ID == set_name, , drop = FALSE]
    if (nrow(hit) == 0) next
    safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", set_name)
    p <- gseaplot2(
      gsea_res,
      geneSetID = set_name,
      title = paste0(set_name, " (NES=", round(hit$NES[1], 3), ", padj=", signif(hit$p.adjust[1], 3), ")"),
      pvalue_table = TRUE,
      ES_geom = "line"
    )
    ggsave(file.path(gsea_dir, paste0("GSEA_", safe_name, ".png")), p, width = 10, height = 6, dpi = 300)
  }
}

# ==============================================================================
# 5. 个性化基因集分析（放在最后，按需改下面两个参数）
# ==============================================================================
custom_gene_set_name <- "MyPersonalSet"
custom_gene_set_genes <- c("GENE1", "GENE2", "GENE3")

is_placeholder <- identical(custom_gene_set_genes, c("GENE1", "GENE2", "GENE3"))
if (!is_placeholder) {
  personal_dir <- file.path(custom_set_dir, "personalized_set")
  dir.create(personal_dir, recursive = TRUE, showWarnings = FALSE)

  # 热图
  personal_match <- unname(row_lookup[toupper(custom_gene_set_genes)])
  personal_match <- unique(personal_match[!is.na(personal_match)])
  if (length(personal_match) >= 2) {
    heat_personal <- exprSet_heatmap[personal_match, , drop = FALSE]
    ann_personal <- annotation_col[colnames(heat_personal), , drop = FALSE]
    pheatmap(
      heat_personal,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      annotation_col = ann_personal,
      show_rownames = TRUE,
      show_colnames = FALSE,
      scale = "row",
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      border_color = NA,
      main = paste0("Personal Set: ", custom_gene_set_name),
      filename = file.path(personal_dir, paste0("Heatmap_", custom_gene_set_name, ".png")),
      width = 9,
      height = max(6, 4 + nrow(heat_personal) * 0.2)
    )
  }

  # GSEA（单一基因集）
  personal_term2gene <- data.frame(
    term = custom_gene_set_name,
    gene = unique(custom_gene_set_genes),
    stringsAsFactors = FALSE
  )
  personal_gsea <- tryCatch(
    GSEA(
      geneList = geneList,
      TERM2GENE = personal_term2gene,
      pvalueCutoff = 1,
      minGSSize = 3,
      maxGSSize = 500,
      BPPARAM = BiocParallel::SerialParam(),
      seed = 123
    ),
    error = function(e) NULL
  )
  if (!is.null(personal_gsea) && nrow(personal_gsea@result) > 0) {
    write.csv(personal_gsea@result, file.path(personal_dir, "Personal_Set_GSEA_Table.csv"), row.names = FALSE)
    p_personal <- gseaplot2(personal_gsea, geneSetID = custom_gene_set_name, pvalue_table = TRUE)
    ggsave(file.path(personal_dir, paste0("GSEA_", custom_gene_set_name, ".png")), p_personal, width = 10, height = 6, dpi = 300)
  }
}

cat("Step14 完成：已输出 custom_gene_sets.xlsx 的热图与 GSEA 结果。\n")
