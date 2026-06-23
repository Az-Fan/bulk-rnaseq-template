# ==============================================================================
# Step05 主GSEA分析
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

# 功能富集分析根目录
functional_enrichment_dir <- file.path(output_dir, "functional_enrichment")
dir.create(functional_enrichment_dir, showWarnings = FALSE)

# GSEA 分析主目录
gsea_analysis_dir <- file.path(functional_enrichment_dir, "gsea_analysis")
dir.create(gsea_analysis_dir, showWarnings = FALSE)

# 表格存放目录
gsea_tables_dir <- file.path(gsea_analysis_dir, "tables")
dir.create(gsea_tables_dir, showWarnings = FALSE)

# ==============================================================================
# PART 1: 全局参数设置
# ==============================================================================
# --- 1. GSEA 分析阈值 ---
calc_pvalue_cutoff <- 1      # 保留所有结果以供后续筛选

# --- 2. 绘图过滤阈值 ---
# 注意：这里设置了宽松的阈值，以便尽可能多地输出显著通路
# 你可以根据实际结果将其改为 0.05
plot_p_cutoff      <- 0.25   

# --- 3. 气泡图参数 (汇总图) ---
plot_show_n        <- 15     # 汇总气泡图中显示的通路数量

# --- 4. 单通路图参数 (gseaplot2) ---
single_plot_width  <- 10      # PNG 宽度
single_plot_height <- 6      # PNG 高度

# ==============================================================================
# PART 2: 数据处理
# ==============================================================================
# 加载差异分析结果
diff_rdata_path <- file.path(output_dir, "diff_analysis/DEseq2_Diff_Annotated.Rdata")
if(file.exists(diff_rdata_path)){
  load(file = diff_rdata_path)
} else {
  stop(paste0("未找到差异分析结果文件: ", diff_rdata_path))
}

library(dplyr)
library(stringr)
library(clusterProfiler)
library(msigdbr)
library(ggplot2)
library(enrichplot)
library(BiocParallel)

# 准备 GeneList
# 论文方法学约定：GSEA ranking metric = DESeq2 Wald statistic。
# 若旧结果中缺少 stat，则回退为 sign(logFC_raw) * -log10(P.Value)。
gene_df <- as.data.frame(res_annotated) %>%
  mutate(
    ranking_metric = dplyr::case_when(
      !is.na(stat) ~ stat,
      is.na(stat) & !is.na(logFC_raw) & !is.na(P.Value) ~ sign(logFC_raw) * -log10(pmax(P.Value, .Machine$double.xmin)),
      TRUE ~ NA_real_
    ),
    ranking_metric_source = dplyr::case_when(
      !is.na(stat) ~ "DESeq2 Wald statistic",
      is.na(stat) & !is.na(logFC_raw) & !is.na(P.Value) ~ "fallback sign(logFC_raw) * -log10(P.Value)",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::select(gene, stat, logFC_raw, P.Value, ranking_metric, ranking_metric_source) %>%
  filter(gene != "" & !is.na(gene) & !is.na(ranking_metric) & is.finite(ranking_metric)) %>%
  distinct(gene, .keep_all = TRUE)

write.csv(gene_df, file.path(gsea_tables_dir, "GSEA_ranked_gene_list.csv"), row.names = FALSE)

geneList <- gene_df$ranking_metric
names(geneList) <- gene_df$gene
geneList <- sort(geneList, decreasing = TRUE)

head(geneList)

# ==============================================================================
# PART 3: GSEA 分析核心流程
# ==============================================================================

# 定义需要分析的基因集列表
target_dbs <- list(
  list(name="Reactome", cat="C2", subcat="CP:REACTOME", prefix="REACTOME_"),
  list(name="Hallmark", cat="H",  subcat=NULL,       prefix="HALLMARK_"),
  list(name="KEGG",     cat="C2", subcat="CP:KEGG",  prefix="KEGG_"),
  list(name="GO_BP",    cat="C5", subcat="GO:BP",    prefix="GOBP_"),
  list(name="GTRD_TF",  cat="C3", subcat="TFT:GTRD", prefix="")
)

# 初始化列表用于存储结果对象
all_gsea_objects <- list() 
# 初始化列表用于存储汇总气泡图
summary_plots <- list()
gsea_provenance <- list()

# 循环进行分析
for (db in target_dbs) {
  
  message(paste0("正在进行 GSEA 分析: ", db$name, "..."))
  
  # 1. 获取基因集
  m_df <- msigdbr(species = "Homo sapiens", category = db$cat, subcategory = db$subcat)
  gsea_provenance[[db$name]] <- data.frame(
    analysis = "GSEA",
    database = db$name,
    source = "MSigDB via msigdbr",
    category = db$cat,
    subcategory = ifelse(is.null(db$subcat), "", db$subcat),
    msigdbr_package_version = as.character(packageVersion("msigdbr")),
    gene_sets = length(unique(m_df$gs_name)),
    species = "Homo sapiens",
    retrieved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stringsAsFactors = FALSE
  )
  # 使用 gene_symbol 进行富集 (与 geneList 对应)
  t2g <- m_df %>% dplyr::select(gs_name, gene_symbol)
  
  # 2. 运行 GSEA
  gsea_res <- GSEA(geneList, 
                   TERM2GENE = t2g,
                   pvalueCutoff = calc_pvalue_cutoff,
                   pAdjustMethod = "BH",
                   seed = 123,
                   BPPARAM = BiocParallel::SerialParam())
  
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    message(paste0("  -> ", db$name, " 未发现任何富集通路，跳过。"))
    next
  }
  
  # 3. 数据清洗 (美化 Description 列)
  res_df <- gsea_res@result
  if(db$prefix != "") {
    res_df$Description <- gsub(paste0("^", db$prefix), "", res_df$Description)
  }
  res_df$Description <- gsub("_", " ", res_df$Description)
  res_df$Description <- str_to_sentence(res_df$Description)
  
  # 将清洗后的描述写回 GSEA 对象 (绘图时会用到)
  gsea_res@result$Description <- res_df$Description
  
  # 保存对象到列表
  all_gsea_objects[[db$name]] <- gsea_res
  
  # 4. 保存完整结果表格
  # 路径: results/functional_enrichment/gsea_analysis/tables/
  write.csv(res_df, 
            file = file.path(gsea_tables_dir, paste0("GSEA_Full_Table_", db$name, ".csv")), 
            row.names = FALSE)
  
  # 5. 生成汇总气泡图对象 (仅作为记录，稍后可选择性保存)
  gsea_res_plot <- gsea_res
  # 过滤用于绘图
  gsea_res_plot@result <- gsea_res@result %>% filter(p.adjust < plot_p_cutoff)
  
  if (nrow(gsea_res_plot@result) > 0) {
    p <- dotplot(gsea_res_plot, 
                 showCategory = plot_show_n, 
                 split = ".sign", 
                 label_format = 60, 
                 title = paste0("GSEA Enrichment: ", db$name)) + 
      facet_grid(.~.sign) + 
      scale_color_gradient(low = "blue", high = "red") + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    summary_plots[[db$name]] <- p
  }
}


# ==============================================================================
# PART 4: 保存汇总气泡图 PDF (可选)
# ==============================================================================
if (length(gsea_provenance)) {
  write.csv(
    dplyr::bind_rows(gsea_provenance),
    file.path(gsea_tables_dir, "GSEA_Gene_Set_Provenance.csv"),
    row.names = FALSE
  )
}

if (length(summary_plots) > 0) {
  pdf_file <- file.path(gsea_analysis_dir, "GSEA_Summary_Dotplots.pdf")
  message(paste0("正在生成汇总气泡图 PDF: ", pdf_file))
  
  pdf(file = pdf_file, width = 12, height = 10)
  for (name in names(summary_plots)) {
    print(summary_plots[[name]])
  }
  dev.off()
}


# ==============================================================================
# PART 5: 批量绘制单通路 GSEA 图 (PNG) - 分文件夹存储 (Top 5)
# ==============================================================================

# 定义绘图函数 (新增 max_n 参数)
batch_save_gsea_plots <- function(gsea_list, root_dir, p_cut, max_n = 3) {
  
  message("\n--------------------------------------------------")
  message("开始批量绘制单通路 GSEA 曲线图...")
  message(paste0("根目录: ", root_dir))
  message(paste0("P.adjust 阈值: < ", p_cut))
  message(paste0("最大绘图数量: 每个数据库 Top ", max_n))
  
  # 遍历每个数据库 (Hallmark, KEGG...)
  for (db_name in names(gsea_list)) {
    
    obj <- gsea_list[[db_name]]
    
    # 1. 为每个数据库创建独立的子文件夹
    db_plot_dir <- file.path(root_dir, paste0(db_name, "_plots"))
    if(!dir.exists(db_plot_dir)) dir.create(db_plot_dir, recursive = TRUE)
    old_png <- list.files(db_plot_dir, pattern = "\\.png$", full.names = TRUE)
    if (length(old_png) > 0) file.remove(old_png)
    
    # 2. 筛选显著通路，并限制数量 (修改点)
    significant_pathways <- obj@result %>% 
      filter(p.adjust < p_cut) %>% 
      arrange(p.adjust) %>%  # 按显著性排序 (P值从小到大)
      head(max_n)            # 只取前 N 个
    
    if (nrow(significant_pathways) == 0) {
      message(paste0("  -> ", db_name, ": 无满足 P < ", p_cut, " 的通路，跳过。"))
      next
    }
    
    message(paste0("  -> ", db_name, ": 发现 ", nrow(significant_pathways), " 条通路 (已限制 Top ", max_n, ")，正在绘制..."))
    
    # 3. 遍历每一条显著通路进行绘图
    for (i in 1:nrow(significant_pathways)) {
      
      pid     <- significant_pathways$ID[i]          # 原始 ID
      desc    <- significant_pathways$Description[i] # 清洗后的描述
      nes_val <- significant_pathways$NES[i]         # NES 值
      
      # 智能配色
      line_color <- if(nes_val > 0) "#B31B21" else "#1465AC"
      
      # 文件名处理
      safe_filename <- gsub("[[:punct:]]", "_", desc) 
      safe_filename <- gsub("\\s+", "_", safe_filename)
      # 文件名包含数据库名和NES正负，方便排序查看
      direction <- if(nes_val > 0) "UP" else "DOWN"
      
      file_path <- file.path(db_plot_dir, paste0(i, "_", direction, "_", safe_filename, ".png"))
      
      # 绘制 GSEA图
      p <- gseaplot2(obj, 
                     geneSetID = pid, 
                     title = desc, 
                     color = line_color,
                     pvalue_table = TRUE,
                     ES_geom = "line",
                     base_size = 14)
      
      # 保存 PNG (假设 single_plot_width/height 在之前已定义)
      ggsave(filename = file_path, plot = p, 
             width = single_plot_width, height = single_plot_height, 
             dpi = 300)
    }
  }
  message("--------------------------------------------------")
  message("所有单图绘制完成！")
}


# 执行批量绘图
if (exists("all_gsea_objects") && length(all_gsea_objects) > 0) {
  # 这里设置 max_n = 5
  batch_save_gsea_plots(gsea_list = all_gsea_objects, 
                        root_dir    = gsea_analysis_dir, 
                        p_cut       = plot_p_cutoff,
                        max_n       = 5) 
} else {
  message("错误：未找到 GSEA 结果对象，请检查 Part 3。")
}

# ==============================================================================
# PART 6: GO-GSEA + aPEAR 网络可视化（按需求新增）
# ==============================================================================
apear_out_dir <- file.path(gsea_analysis_dir, "apear_network")
dir.create(apear_out_dir, recursive = TRUE, showWarnings = FALSE)

run_apear <- function() {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    message("跳过 aPEAR 部分：缺少 org.Hs.eg.db。")
    return(invisible(NULL))
  }

  res_df <- as.data.frame(res_annotated)
  if (!all(c("gene", "stat") %in% colnames(res_df))) {
    message("跳过 aPEAR 部分：res_annotated 缺少 gene/stat 列。")
    return(invisible(NULL))
  }

  # 优先使用 Step01 产出的 entrez 列；若缺失则现算
  if (!"entrez" %in% colnames(res_df)) {
    ids <- AnnotationDbi::select(
      x = org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(res_df$gene[!is.na(res_df$gene) & res_df$gene != ""]),
      columns = c("SYMBOL", "ENTREZID"),
      keytype = "SYMBOL",
      skipValidKeysTest = TRUE
    )
    ids <- ids[!is.na(ids$ENTREZID) & ids$ENTREZID != "", c("SYMBOL", "ENTREZID")]
    ids <- ids[!duplicated(ids$SYMBOL), , drop = FALSE]
    res_df <- dplyr::left_join(res_df, ids, by = c("gene" = "SYMBOL"))
    res_df$entrez <- res_df$ENTREZID
  }

  df2 <- res_df %>%
    dplyr::mutate(
      ranking_metric = dplyr::case_when(
        !is.na(stat) ~ stat,
        is.na(stat) & !is.na(logFC_raw) & !is.na(P.Value) ~ sign(logFC_raw) * -log10(pmax(P.Value, .Machine$double.xmin)),
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::filter(!is.na(gene), gene != "", !is.na(entrez), entrez != "", !is.na(ranking_metric), is.finite(ranking_metric)) %>%
    dplyr::select(gene, ranking_metric, entrez) %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  if (nrow(df2) < 20) {
    message("跳过 aPEAR 部分：有效 ENTREZ 基因数量过少。")
    return(invisible(NULL))
  }

  geneList_entrez <- df2$ranking_metric
  names(geneList_entrez) <- as.character(df2$entrez)
  geneList_entrez <- geneList_entrez[names(geneList_entrez) != ""]
  geneList_entrez <- geneList_entrez[!duplicated(names(geneList_entrez))]
  geneList_entrez <- sort(geneList_entrez, decreasing = TRUE)

  enrich_go <- gseGO(
    geneList = geneList_entrez,
    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
    ont = "ALL",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE
  )

  if (is.null(enrich_go) || nrow(enrich_go@result) == 0) {
    message("aPEAR 部分：gseGO 未得到显著结果。")
    return(invisible(NULL))
  }

  enrich_go <- clusterProfiler::simplify(enrich_go)
  write.csv(enrich_go@result, file.path(apear_out_dir, "GO_GSEA_Full_Table.csv"), row.names = FALSE)

  if (!requireNamespace("aPEAR", quietly = TRUE)) {
    message("aPEAR 未安装，已输出 GO_GSEA_Full_Table.csv；安装后可自动生成网络图。")
    return(invisible(NULL))
  }

  p_net <- aPEAR::enrichmentNetwork(enrich_go@result, drawEllipses = TRUE, fontSize = 2.5)
  ggsave(file.path(apear_out_dir, "GO_GSEA_EnrichmentNetwork.png"), p_net, width = 12, height = 8, dpi = 300)
  ggsave(file.path(apear_out_dir, "GO_GSEA_EnrichmentNetwork.pdf"), p_net, width = 12, height = 8)

  if (requireNamespace("plotly", quietly = TRUE) && requireNamespace("htmlwidgets", quietly = TRUE)) {
    p_interactive <- plotly::ggplotly(p_net, tooltip = c("ID", "Cluster", "Cluster size"))
    htmlwidgets::saveWidget(
      p_interactive,
      file = file.path(apear_out_dir, "GO_GSEA_EnrichmentNetwork_interactive.html"),
      selfcontained = FALSE
    )
  }

  message("aPEAR 网络分析已完成。")
  invisible(NULL)
}

tryCatch(
  run_apear(),
  error = function(e) message("aPEAR 部分执行失败（已跳过，不影响主流程）：", e$message)
)

cat("\n所有 GSEA 分析流程结束。\n")
