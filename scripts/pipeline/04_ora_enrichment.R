# ==============================================================================
# Step04 ORA富集分析（GO/KEGG/Hallmark）
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. 设置输出文件夹结构
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results_01")
dir.create(output_dir, showWarnings = FALSE)

# 路径设置
data_processed_dir <- file.path(output_dir, "data_processed")
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, showWarnings = FALSE)

# 子文件夹
go_analysis_dir <- file.path(functional_enrichment_dir, "go_analysis")
kegg_analysis_dir <- file.path(functional_enrichment_dir, "kegg_analysis")
hallmark_analysis_dir <- file.path(functional_enrichment_dir, "hallmark_analysis")
adv_plots_dir <- file.path(functional_enrichment_dir, "advanced_plots") 

dir.create(go_analysis_dir, showWarnings = FALSE)
dir.create(kegg_analysis_dir, showWarnings = FALSE)
dir.create(hallmark_analysis_dir, showWarnings = FALSE)
dir.create(adv_plots_dir, showWarnings = FALSE)

# ==============================================================================
# PART 1: 参数设置
# ==============================================================================
param_fc_cutoff   <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
param_p_cutoff    <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
plot_p_cutoff     <- 0.05   
plot_q_cutoff     <- 0.05    

# 绘图尺寸
plot_show_n_max   <- as.integer(Sys.getenv("ORA_PLOT_TOP_N", unset = "10"))
plot_base_width   <- 14     # 气泡图宽度
plot_base_height  <- 8      # 基础高度

# ==============================================================================
# PART 2: 数据加载
# ==============================================================================
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(stringr)
library(ggplot2)
library(msigdbr)
library(enrichplot)
has_apear <- requireNamespace("aPEAR", quietly = TRUE)
has_plotly <- requireNamespace("plotly", quietly = TRUE)
has_htmlwidgets <- requireNamespace("htmlwidgets", quietly = TRUE)

# 加载差异结果
load(file = file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata"))
res_df <- as.data.frame(res_annotated)

# 安全 ID 映射：规避部分环境下 bitr() 内部 keys() 校验触发 sqlite 打开失败
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

# setReadable 的安全包装：若数据库映射异常则保留原对象继续后续流程
safe_set_readable <- function(res_obj, key_type = "ENTREZID") {
  if (is.null(res_obj)) return(NULL)
  tryCatch(
    setReadable(res_obj, OrgDb = org.Hs.eg.db, keyType = key_type),
    error = function(e) {
      p_hier <- tryCatch(make_apear("hier"), error = function(e2) NULL)
      if (!is.null(p_hier)) return(p_hier)
      message("setReadable 跳过：", e$message)
      res_obj
    }
  )
}

# ID 转换
if(!"ENTREZID" %in% colnames(res_df)){
  cat("正在进行 ID 转换...\n")
  ids <- map_symbol_to_entrez(res_df$gene)
  res_df <- left_join(res_df, ids, by = c("gene" = "SYMBOL"))
}
if (!"ENTREZID" %in% colnames(res_df) && "entrez" %in% colnames(res_df)) {
  res_df$ENTREZID <- res_df$entrez
}

universe_entrez <- unique(as.character(res_df$ENTREZID[!is.na(res_df$ENTREZID) & res_df$ENTREZID != ""]))
write.csv(
  data.frame(ENTREZID = universe_entrez, stringsAsFactors = FALSE),
  file.path(functional_enrichment_dir, "ORA_universe_entrez.csv"),
  row.names = FALSE
)
cat("ORA universe 基因数：", length(universe_entrez), "\n")

# 准备基因列表 (用于 cnetplot 染色)
genelist <- res_df$logFC
names(genelist) <- res_df$ENTREZID
genelist <- sort(genelist, decreasing = TRUE)

# 筛选差异基因
diff_df <- res_df %>% filter(!is.na(ENTREZID) & adj.P.Val < param_p_cutoff & abs(logFC) > param_fc_cutoff)
gene_cluster <- list(Up = diff_df$ENTREZID[diff_df$logFC > 0],
                     Down = diff_df$ENTREZID[diff_df$logFC < 0],
                     All = diff_df$ENTREZID)

# ==============================================================================
# PART 3: 富集分析 (分面气泡图专用 - Up/Down/All)
# ==============================================================================
cat("正在进行多组比较富集分析 (Up/Down/All)...\n")
m_go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(entrez_gene = as.character(entrez_gene))
m_go_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(entrez_gene = as.character(entrez_gene))
m_go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(entrez_gene = as.character(entrez_gene))
m_kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(entrez_gene = as.character(entrez_gene))
m_hall <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene) %>% mutate(entrez_gene = as.character(entrez_gene))

write.csv(
  data.frame(
    analysis = "ORA",
    source = "MSigDB via msigdbr",
    msigdbr_package_version = as.character(packageVersion("msigdbr")),
    collections = c("C5 GO:BP", "C5 GO:CC", "C5 GO:MF", "C2 CP:KEGG", "H Hallmark"),
    species = "Homo sapiens",
    retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  ),
  file.path(functional_enrichment_dir, "ORA_Gene_Set_Provenance.csv"),
  row.names = FALSE
)

# 全量计算，pvalueCutoff=1
res_go_bp <- compareCluster(gene_cluster, fun = "enricher", TERM2GENE = m_go_bp, universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
res_go_cc <- compareCluster(gene_cluster, fun = "enricher", TERM2GENE = m_go_cc, universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
res_go_mf <- compareCluster(gene_cluster, fun = "enricher", TERM2GENE = m_go_mf, universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
res_kegg <- compareCluster(gene_cluster, fun = "enricher", TERM2GENE = m_kegg, universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
res_hall <- compareCluster(gene_cluster, fun = "enricher", TERM2GENE = m_hall, universe = universe_entrez, pvalueCutoff = 1, qvalueCutoff = 1)

# ID 转回 Symbol
res_kegg <- safe_set_readable(res_kegg, "ENTREZID")
res_hall <- safe_set_readable(res_hall, "ENTREZID")

# ==============================================================================
# PART 3.5: 清洗通路名称
# ==============================================================================
clean_description <- function(res_obj, prefix_remove = "") {
  if (is.null(res_obj)) return(NULL)
  
  if (inherits(res_obj, "compareClusterResult")) {
    df <- res_obj@compareClusterResult
  } else if (inherits(res_obj, "enrichResult")) {
    df <- res_obj@result
  } else {
    warning("Unknown object class, skipping clean.")
    return(res_obj)
  }
  
  if (nrow(df) == 0) return(res_obj)
  
  if(prefix_remove != "") {
    df$Description <- gsub(paste0("^", prefix_remove, "_"), "", df$Description)
  }
  df$Description <- gsub("_", " ", df$Description)
  df$Description <- stringr::str_to_sentence(df$Description)
  
  if (inherits(res_obj, "compareClusterResult")) {
    res_obj@compareClusterResult <- df
  } else {
    res_obj@result <- df
  }
  return(res_obj)
}

cat("正在清洗 KEGG 和 Hallmark 通路名称...\n")
res_kegg <- clean_description(res_kegg, "KEGG")
res_hall <- clean_description(res_hall, "HALLMARK")
res_go_bp <- clean_description(res_go_bp, "")
res_go_cc <- clean_description(res_go_cc, "")
res_go_mf <- clean_description(res_go_mf, "")

# ==============================================================================
# PART 4: 保存全量表格
# ==============================================================================
save_table <- function(res, name, path) {
  if (is.null(res)) return()
  df <- as.data.frame(res)
  if (nrow(df) == 0) return()
  write.csv(df, file = file.path(path, paste0("Enrichment_Full_", name, ".csv")), row.names = FALSE)
}

save_table(res_go_bp, "GO_BP", go_analysis_dir)
save_table(res_go_cc, "GO_CC", go_analysis_dir)
save_table(res_go_mf, "GO_MF", go_analysis_dir)
save_table(res_kegg, "KEGG", kegg_analysis_dir)
save_table(res_hall, "Hallmark", hallmark_analysis_dir)

# ==============================================================================
# PART 5: 绘制分面气泡图
# ==============================================================================
plot_dot_faceted <- function(res_obj, title, path, filename) {
  if (is.null(res_obj)) return()
  
  res_filtered <- res_obj
  res_filtered@compareClusterResult <- res_filtered@compareClusterResult %>% 
    filter(p.adjust < plot_p_cutoff & qvalue < plot_q_cutoff)
  
  if (nrow(res_filtered@compareClusterResult) == 0) {
    cat(paste0("提示：", title, " 过滤后无显著通路，跳过绘图。\n"))
    return()
  }
  
  counts <- res_filtered@compareClusterResult %>% count(Cluster) %>% pull(n)
  actual_rows <- sum(pmin(counts, plot_show_n_max))
  h <- max(plot_base_height, actual_rows * 0.35 + 2)
  
  p <- dotplot(res_filtered, showCategory = plot_show_n_max) + 
    ggtitle(title) +
    scale_color_gradient(low = "red", high = "blue") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存 PNG 和 PDF
  ggsave(file.path(path, paste0(filename, ".png")), p, width = plot_base_width, height = h, dpi = 300, limitsize = FALSE)
  ggsave(file.path(path, paste0(filename, ".pdf")), p, width = plot_base_width, height = h, limitsize = FALSE)
  
  cat(paste0(">>> 图表已保存: ", filename, "\n"))
}

cat("正在绘制分面气泡图...\n")
plot_dot_by_cluster <- function(res_obj, title, path, filename) {
  if (is.null(res_obj)) return()
  filtered_df <- res_obj@compareClusterResult %>%
    filter(p.adjust < plot_p_cutoff & qvalue < plot_q_cutoff)
  if (nrow(filtered_df) == 0) return()

  for (cluster_name in intersect(c("Up", "Down", "All"), unique(filtered_df$Cluster))) {
    cluster_df <- filtered_df %>%
      filter(Cluster == cluster_name) %>%
      arrange(p.adjust) %>%
      slice_head(n = plot_show_n_max) %>%
      mutate(
        Description = sub("^Gobp\\s+", "", Description, ignore.case = TRUE),
        Description = factor(Description, levels = rev(Description)),
        GeneRatioNumeric = vapply(
          strsplit(as.character(GeneRatio), "/", fixed = TRUE),
          function(x) as.numeric(x[[1]]) / as.numeric(x[[2]]),
          numeric(1)
        )
      )
    if (nrow(cluster_df) == 0) next

    p <- ggplot(cluster_df, aes(x = GeneRatioNumeric, y = Description)) +
      geom_point(aes(size = Count, color = p.adjust), alpha = 0.9) +
      scale_color_gradient(low = "#D73027", high = "#4575B4", trans = "reverse") +
      labs(
        title = paste(title, "-", cluster_name),
        x = "Gene ratio", y = NULL, color = "Adjusted P", size = "Gene count"
      ) +
      scale_x_continuous(labels = scales::label_number(accuracy = 0.01)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold"),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 9)
      )
    h <- max(4.5, nrow(cluster_df) * 0.38 + 1.6)
    ggsave(
      file.path(path, paste0(filename, "_", cluster_name, ".png")),
      p, width = 8, height = h, dpi = 300
    )
  }
}

plot_dot_by_cluster(res_go_bp, "GO BP Enrichment", go_analysis_dir, "GO_BP_Dotplot")
plot_dot_by_cluster(res_go_cc, "GO CC Enrichment", go_analysis_dir, "GO_CC_Dotplot")
plot_dot_by_cluster(res_go_mf, "GO MF Enrichment", go_analysis_dir, "GO_MF_Dotplot")
plot_dot_by_cluster(res_kegg, "KEGG Enrichment", kegg_analysis_dir, "KEGG_Dotplot")
plot_dot_by_cluster(res_hall, "Hallmark Enrichment", hallmark_analysis_dir, "Hallmark_Dotplot")


# ==============================================================================
# PART 6: 绘制进阶图 (新增 PDF 保存)
# ==============================================================================
make_single_enrich_obj <- function(type) {
  if (type == "GO_BP") {
    res <- enricher(gene = gene_cluster$All, TERM2GENE = m_go_bp, universe = universe_entrez)
    return(clean_description(res, ""))
  } else if (type == "KEGG") {
    res <- enricher(gene = gene_cluster$All, TERM2GENE = m_kegg, universe = universe_entrez)
    if (is.null(res) || nrow(res@result) == 0) return(NULL)
    res <- safe_set_readable(res, "ENTREZID")
    return(clean_description(res, "KEGG"))
  }
}

plot_adv <- function(obj, title, filename) {
  if (is.null(obj) || nrow(obj) == 0) return()
  
  # 1. 网络图
  p1 <- cnetplot(obj, foldChange = genelist, showCategory = 5, circular = FALSE, colorEdge = TRUE) + 
    ggtitle(paste0(title, " - Gene-Concept Network"))
  ggsave(file.path(adv_plots_dir, paste0(filename, "_Cnet.png")), p1, width = 12, height = 10, dpi = 300)
  ggsave(file.path(adv_plots_dir, paste0(filename, "_Cnet.pdf")), p1, width = 12, height = 10) # 新增 PDF
  
  # 2. 基因热图
  p2 <- heatplot(obj, foldChange = genelist, showCategory = 10) + 
    ggtitle(paste0(title, " - Gene-Pathway Heatmap")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
  ggsave(file.path(adv_plots_dir, paste0(filename, "_Heat.png")), p2, width = 16, height = 8, dpi = 300)
  ggsave(file.path(adv_plots_dir, paste0(filename, "_Heat.pdf")), p2, width = 16, height = 8) # 新增 PDF
  
  # 3. 富集地图
  if (nrow(obj) > 1) {
    obj_sim <- pairwise_termsim(obj)
    p3 <- emapplot(obj_sim, showCategory = 20, layout = "kk") + 
      ggtitle(paste0(title, " - Enrichment Map"))
    ggsave(file.path(adv_plots_dir, paste0(filename, "_Emap.png")), p3, width = 10, height = 10, dpi = 300)
    ggsave(file.path(adv_plots_dir, paste0(filename, "_Emap.pdf")), p3, width = 10, height = 10) # 新增 PDF
  }
}

cat("正在绘制进阶图...\n")
go_bp_single <- make_single_enrich_obj("GO_BP")
plot_adv(go_bp_single, "GO BP", "GO_BP")

kegg_single <- make_single_enrich_obj("KEGG")
plot_adv(kegg_single, "KEGG", "KEGG")

# ==============================================================================
# PART 7: aPEAR 通路富集网络图（GO / KEGG）
# ==============================================================================
apear_dir <- file.path(functional_enrichment_dir, "apear_network")
dir.create(apear_dir, showWarnings = FALSE)

plot_apear_network <- function(
  enrich_obj,
  tag,
  p_cutoff = 0.05,
  min_terms = 8,
  inner_cutoff = 0.1,
  outer_cutoff = 0.5
) {
  if (!has_apear) {
    message("aPEAR 未安装，跳过 ", tag, " 网络图。")
    return(NULL)
  }
  if (is.null(enrich_obj) || nrow(enrich_obj) == 0) return(NULL)

  net_df <- enrich_obj@result %>% dplyr::filter(!is.na(p.adjust), p.adjust < p_cutoff)
  # 若条目太少则自动放宽阈值，避免网络图点位过少
  if (nrow(net_df) < min_terms) {
    fallback_cutoffs <- c(0.5, 1.0)
    for (pc in fallback_cutoffs) {
      tmp_df <- enrich_obj@result %>% dplyr::filter(!is.na(p.adjust), p.adjust < pc)
      if (nrow(tmp_df) >= min_terms || pc == 1.0) {
        net_df <- tmp_df
        message(
          "aPEAR ", tag, "：显著条目较少，p.adjust 阈值从 ", p_cutoff,
          " 放宽到 ", pc, "（当前条目数=", nrow(net_df), "）"
        )
        break
      }
    }
  }
  if (nrow(net_df) < 2) {
    message("aPEAR 跳过 ", tag, "：显著通路少于2个。")
    return(NULL)
  }

  # 节点越多时自动减小字体，避免标签重叠
  if (nrow(net_df) > 100) {
    net_df <- net_df %>%
      dplyr::arrange(.data$p.adjust, .data$pvalue) %>%
      dplyr::slice_head(n = 100)
    message("aPEAR ", tag, ": using the top 100 terms ranked by adjusted P value.")
  }
  label_font_size <- if (nrow(net_df) >= 80) 1.8 else if (nrow(net_df) >= 40) 2.2 else 2.6
  plot_w <- if (nrow(net_df) >= 60) 13 else 12
  plot_h <- if (nrow(net_df) >= 60) 10 else 9

  make_apear <- function(clust_method) {
    aPEAR::enrichmentNetwork(
      net_df,
      clustMethod = clust_method,
      drawEllipses = TRUE,
      fontSize = label_font_size,
      repelLabels = TRUE,
      innerCutoff = inner_cutoff,
      outerCutoff = outer_cutoff
    )
  }
  p_net <- tryCatch(
    make_apear("markov"),
    error = function(e) {
      message("aPEAR 跳过 ", tag, "：", e$message)
      NULL
    }
  )
  if (is.null(p_net)) return(NULL)
  png_path <- file.path(apear_dir, paste0(tag, "_aPEAR_network.png"))
  pdf_path <- file.path(apear_dir, paste0(tag, "_aPEAR_network.pdf"))
  ggsave(png_path, p_net, width = plot_w, height = plot_h, dpi = 300)
  ggsave(pdf_path, p_net, width = plot_w, height = plot_h)

  # 同款风格补充：viridis 配色版本
  if (requireNamespace("viridis", quietly = TRUE)) {
    p_viridis <- p_net + ggplot2::scale_color_viridis_c(option = "D")
    ggsave(
      file.path(apear_dir, paste0(tag, "_aPEAR_network_viridis.png")),
      p_viridis,
      width = plot_w,
      height = plot_h,
      dpi = 300
    )
  }

  # 可交互版本
  if (has_plotly && has_htmlwidgets) {
    p_int <- plotly::ggplotly(p_net, tooltip = c("ID", "Cluster", "Cluster size"))
    htmlwidgets::saveWidget(
      p_int,
      file = file.path(apear_dir, paste0(tag, "_aPEAR_network_interactive.html")),
      selfcontained = FALSE
    )
  }
  invisible(TRUE)
}

plot_apear_network(
  go_bp_single,
  "GO_BP",
  p_cutoff = plot_p_cutoff,
  min_terms = 12,
  inner_cutoff = 0.1,
  outer_cutoff = 0.5
)
# KEGG 通路通常更少，适当放宽阈值和聚类边界，尽量避免只剩少量节点
plot_apear_network(
  kegg_single,
  "KEGG",
  p_cutoff = max(plot_p_cutoff, 0.25),
  min_terms = 10,
  inner_cutoff = 0.05,
  outer_cutoff = 0.3
)

cat("\n所有分析与图表生成完毕！(包含 PNG 和 PDF)\n")
