# ==============================================================================
# Step09 PPI网络分析（STRING）
# 运行建议（RStudio）：
# 1. 在 Console 先设置：Sys.setenv(OUTPUT_DIR = "results/final/nobatch",
#    SHARED_CACHE_DIR = "shared_cache")
# 2. 然后 Source 当前脚本（Ctrl+Shift+S）
# 说明：本文件为当前项目的正式顺序脚本，可直接修改参数后运行。
# ==============================================================================

rm(list = ls())

# ==============================================================================
# 0. 输出目录与参数
# ==============================================================================
output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
diff_analysis_dir <- file.path(output_dir, "diff_analysis")
ppi_dir <- file.path(output_dir, "functional_enrichment", "ppi_analysis")
cytoscape_dir <- file.path(ppi_dir, "cytoscape_input")
resource_cache_dir <- file.path(getwd(), "resources", "cache")
local_string_dir <- file.path(resource_cache_dir, "stringdb")
local_edges_file <- file.path(local_string_dir, "string_edges_project.csv")

# 共享缓存目录（跨项目复用，避免重复下载 STRING 数据）
shared_cache_dir <- Sys.getenv("SHARED_CACHE_DIR", unset = file.path(getwd(), "shared_cache"))
stringdb_cache_dir <- file.path(shared_cache_dir, "stringdb")

dir.create(ppi_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cytoscape_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(stringdb_cache_dir, recursive = TRUE, showWarnings = FALSE)

ppi_top_n <- as.integer(Sys.getenv("PPI_TOP_N", unset = "200"))
deg_padj_threshold <- as.numeric(Sys.getenv("DEG_PADJ_THRESHOLD", unset = "0.05"))
deg_lfc_threshold <- as.numeric(Sys.getenv("DEG_LFC_THRESHOLD", unset = "1"))
hub_show_n <- 30
conf_score <- 400
species_id <- 9606

skip_with_note <- function(msg) {
  message(msg)
  writeLines(msg, con = file.path(ppi_dir, "SKIPPED_ppi_due_to_network_or_data.txt"))
  quit(save = "no", status = 0)
}

# ==============================================================================
# 1. 依赖检查
# ==============================================================================
required_pkgs <- c("STRINGdb", "igraph", "ggraph", "dplyr", "ggplot2", "ggrepel")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop(paste0("缺少必要包: ", paste(missing_pkgs, collapse = ", "), "。请先安装后再运行 PPI 脚本。"))
}

library(STRINGdb)
library(igraph)
library(ggraph)
library(dplyr)
library(ggplot2)
library(ggrepel)

# STRING 初次下载别名文件体积较大，适当放宽超时避免中途断开
old_timeout <- getOption("timeout")
options(timeout = max(600, old_timeout))
on.exit(options(timeout = old_timeout), add = TRUE)

# ==============================================================================
# 2. 加载差异结果
# ==============================================================================
diff_rdata <- file.path(diff_analysis_dir, "DEseq2_Diff_Annotated.Rdata")
if (!file.exists(diff_rdata)) {
  stop(paste0("未找到差异结果文件: ", diff_rdata))
}
load(diff_rdata) # res_annotated

diff_df <- as.data.frame(res_annotated)
if (!"gene" %in% colnames(diff_df)) {
  diff_df$gene <- diff_df$gene_id
}

sig_genes <- diff_df %>%
  filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
  filter(adj.P.Val < deg_padj_threshold & abs(logFC) > deg_lfc_threshold) %>%
  arrange(desc(abs(logFC)))

if (nrow(sig_genes) == 0) {
  skip_with_note("没有满足阈值的差异基因，Step09 已跳过。")
}

n_limit <- min(nrow(sig_genes), ppi_top_n)
target_genes <- sig_genes[seq_len(n_limit), ]
cat("正在使用前", n_limit, "个差异基因构建 PPI 网络...\n")

# ==============================================================================
# 3. STRING 映射与网络构建
# ==============================================================================
g <- NULL

# 3.1 优先使用本地缓存边（离线模式）
if (file.exists(local_edges_file)) {
  message("检测到本地 STRING 缓存：", local_edges_file)
  local_edges <- tryCatch(read.csv(local_edges_file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (!is.null(local_edges) && all(c("from_symbol", "to_symbol", "combined_score") %in% colnames(local_edges))) {
    sig_symbols <- unique(target_genes$gene)
    use_edges <- local_edges %>%
      dplyr::filter(from_symbol %in% sig_symbols & to_symbol %in% sig_symbols)
    if (nrow(use_edges) > 0) {
      g <- igraph::graph_from_data_frame(
        d = use_edges[, c("from_symbol", "to_symbol", "combined_score")],
        directed = FALSE
      )
      g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(combined_score = "max"))
      V(g)$symbol <- V(g)$name
      V(g)$logFC <- target_genes$logFC[match(V(g)$symbol, target_genes$gene)]
      E(g)$combined_score[is.na(E(g)$combined_score)] <- conf_score
    }
  }
}

# 3.2 若本地缓存不可用，回退在线 STRINGdb 并顺带写入本地缓存
if (is.null(g) || length(V(g)) < 2) {
  # 3.2.1 优先尝试 STRING API（避免下载全量别名库）
  sig_symbols <- unique(target_genes$gene)
  sig_symbols <- sig_symbols[!is.na(sig_symbols) & sig_symbols != ""]
  api_edges <- tryCatch({
    ids_enc <- paste(utils::URLencode(sig_symbols, reserved = TRUE), collapse = "%0d")
    api_url <- paste0(
      "https://string-db.org/api/tsv/network?identifiers=",
      ids_enc,
      "&species=", species_id,
      "&required_score=", conf_score
    )
    x <- read.delim(api_url, stringsAsFactors = FALSE, check.names = FALSE)
    if (nrow(x) > 0 && all(c("preferredName_A", "preferredName_B", "score") %in% colnames(x))) {
      data.frame(
        from_symbol = x$preferredName_A,
        to_symbol = x$preferredName_B,
        combined_score = as.numeric(x$score) * 1000,
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  }, error = function(e) NULL)

  if (!is.null(api_edges) && nrow(api_edges) > 0) {
    g <- igraph::graph_from_data_frame(
      d = api_edges[, c("from_symbol", "to_symbol", "combined_score")],
      directed = FALSE
    )
    g <- igraph::simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(combined_score = "max"))
    V(g)$symbol <- V(g)$name
    V(g)$logFC <- target_genes$logFC[match(V(g)$symbol, target_genes$gene)]
    E(g)$combined_score[is.na(E(g)$combined_score)] <- conf_score

    dir.create(local_string_dir, recursive = TRUE, showWarnings = FALSE)
    write.csv(api_edges, local_edges_file, row.names = FALSE)
  }
}

# 3.2.2 API 失败再回退到 STRINGdb 全流程
if (is.null(g) || length(V(g)) < 2) {
  string_db <- tryCatch(
    STRINGdb$new(
      version = "11.5",
      species = species_id,
      score_threshold = conf_score,
      input_directory = stringdb_cache_dir
    ),
    error = function(e) NULL
  )
  if (is.null(string_db)) {
    skip_with_note("STRINGdb 初始化失败（通常是无外网或 DNS 不可用），Step09 已跳过。")
  }

  mapped_genes <- tryCatch(
    string_db$map(target_genes, "gene", removeUnmappedRows = TRUE),
    error = function(e) NULL
  )
  if (is.null(mapped_genes) || nrow(mapped_genes) < 2) {
    skip_with_note("STRING 映射结果不足（<2），Step09 已跳过。")
  }

  g <- tryCatch(
    string_db$get_subnetwork(mapped_genes$STRING_id),
    error = function(e) NULL
  )
  if (is.null(g) || length(V(g)) < 2) {
    skip_with_note("STRING 子网络节点不足（<2），Step09 已跳过。")
  }

  V(g)$symbol <- mapped_genes$gene[match(V(g)$name, mapped_genes$STRING_id)]
  V(g)$logFC <- mapped_genes$logFC[match(V(g)$name, mapped_genes$STRING_id)]

  # 保存为本地缓存，供后续离线直接调用
  dir.create(local_string_dir, recursive = TRUE, showWarnings = FALSE)
  all_inter <- tryCatch(string_db$get_interactions(unique(mapped_genes$STRING_id)), error = function(e) NULL)
  if (!is.null(all_inter) && nrow(all_inter) > 0) {
    all_inter$from_symbol <- mapped_genes$gene[match(all_inter$from, mapped_genes$STRING_id)]
    all_inter$to_symbol <- mapped_genes$gene[match(all_inter$to, mapped_genes$STRING_id)]
    all_inter <- all_inter[!is.na(all_inter$from_symbol) & !is.na(all_inter$to_symbol), c("from", "to", "from_symbol", "to_symbol", "combined_score")]
    write.csv(all_inter, local_edges_file, row.names = FALSE)
  }
}

if (length(V(g)) < 2) {
  skip_with_note("PPI 网络有效节点不足（<2），Step09 已跳过。")
}

V(g)$degree <- degree(g)

# ==============================================================================
# 4. Hub 网络图
# ==============================================================================
top_nodes <- names(sort(degree(g), decreasing = TRUE))[1:min(length(V(g)), hub_show_n)]
sub_g <- induced_subgraph(g, top_nodes)
V(sub_g)$color_type <- ifelse(V(sub_g)$logFC > 0, "Up", "Down")

pdf(file.path(ppi_dir, "PPI_Hub_Network.pdf"), width = 8, height = 7)
p <- ggraph(sub_g, layout = "fr") +
  geom_edge_link(alpha = 0.4, color = "grey80") +
  geom_node_point(aes(size = degree, color = color_type), alpha = 0.9) +
  scale_color_manual(values = c("Up" = "#B31B21", "Down" = "#1465AC")) +
  scale_size(range = c(3, 8)) +
  geom_node_text(aes(label = symbol), repel = TRUE, size = 3, fontface = "bold") +
  theme_void() +
  labs(
    title = paste0("PPI Network - Top ", hub_show_n, " Hub Genes"),
    subtitle = "Node size = Degree; Color = LogFC"
  ) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = "bold"))
print(p)
dev.off()

# ==============================================================================
# 5. 导出 Cytoscape 输入
# ==============================================================================
node_attr <- data.frame(
  STRING_id = V(g)$name,
  Symbol = V(g)$symbol,
  logFC = V(g)$logFC,
  Degree = V(g)$degree,
  stringsAsFactors = FALSE
)
write.csv(node_attr, file.path(cytoscape_dir, "Node_Attributes.csv"), row.names = FALSE)

edge_data <- igraph::as_data_frame(g, what = "edges")
edge_data$from_symbol <- node_attr$Symbol[match(edge_data$from, node_attr$STRING_id)]
edge_data$to_symbol <- node_attr$Symbol[match(edge_data$to, node_attr$STRING_id)]

final_edges <- edge_data %>%
  dplyr::select(Source = from_symbol, Target = to_symbol, Score = combined_score)
write.csv(final_edges, file.path(cytoscape_dir, "Network_Edges.csv"), row.names = FALSE)

cat("PPI 分析完成。输出目录：", ppi_dir, "\n")
