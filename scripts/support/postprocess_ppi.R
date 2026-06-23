options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ggrepel)
})
source("scripts/support/project_plot.R")

output_dir <- Sys.getenv("OUTPUT_DIR", unset = "results")
ppi_dir <- file.path(output_dir, "functional_enrichment", "ppi_analysis")
node_path <- file.path(ppi_dir, "cytoscape_input", "Node_Attributes.csv")
edge_path <- file.path(ppi_dir, "cytoscape_input", "Network_Edges.csv")
if (!file.exists(node_path) || !file.exists(edge_path)) {
  message("PPI postprocessing skipped: Cytoscape input files are absent.")
  quit(save = "no", status = 0)
}

nodes <- read.csv(node_path, check.names = FALSE)
edges <- read.csv(edge_path, check.names = FALSE)
if (!nrow(nodes) || !nrow(edges)) quit(save = "no", status = 0)

g <- graph_from_data_frame(
  edges[, c("Source", "Target", "Score")],
  directed = FALSE,
  vertices = data.frame(name = nodes$Symbol, stringsAsFactors = FALSE)
)
g <- simplify(g, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(Score = "max"))
weights <- E(g)$Score
weights[!is.finite(weights) | weights <= 0] <- 1
V(g)$logFC <- nodes$logFC[match(V(g)$name, nodes$Symbol)]

V(g)$Degree <- degree(g)
V(g)$Betweenness <- betweenness(g, directed = FALSE, weights = 1 / weights, normalized = TRUE)
V(g)$Closeness <- closeness(g, weights = 1 / weights, normalized = TRUE)
V(g)$PageRank <- page_rank(g, directed = FALSE, weights = weights)$vector
V(g)$Eigenvector <- eigen_centrality(g, directed = FALSE, weights = weights)$vector
V(g)$Module <- if (vcount(g) >= 3 && ecount(g) >= 2) {
  cluster_louvain(g, weights = weights)$membership
} else rep(1L, vcount(g))

zscore <- function(x) {
  if (length(unique(x[is.finite(x)])) <= 1) return(rep(0, length(x)))
  as.numeric(scale(x))
}
V(g)$HubScore <- rowMeans(
  cbind(
    zscore(log1p(V(g)$Degree)),
    zscore(log1p(V(g)$Betweenness)),
    zscore(V(g)$PageRank),
    zscore(V(g)$Eigenvector),
    zscore(abs(V(g)$logFC))
  ),
  na.rm = TRUE
)

metrics <- data.frame(
  Symbol = V(g)$name,
  logFC = V(g)$logFC,
  Degree = V(g)$Degree,
  Betweenness = V(g)$Betweenness,
  Closeness = V(g)$Closeness,
  PageRank = V(g)$PageRank,
  Eigenvector = V(g)$Eigenvector,
  Module = as.integer(V(g)$Module),
  HubScore = V(g)$HubScore,
  stringsAsFactors = FALSE
) %>% arrange(desc(HubScore))
write.csv(metrics, file.path(ppi_dir, "PPI_Hub_Ranking.csv"), row.names = FALSE)

module_summary <- metrics %>%
  group_by(Module) %>%
  summarise(
    Nodes = n(),
    MeanLogFC = mean(logFC, na.rm = TRUE),
    UpGenes = sum(logFC > 0, na.rm = TRUE),
    DownGenes = sum(logFC < 0, na.rm = TRUE),
    TopHub = Symbol[which.max(HubScore)],
    .groups = "drop"
  ) %>%
  arrange(desc(Nodes))
module_summary$TopFunction <- NA_character_
module_summary$TopFunctionFDR <- NA_real_

module_enrichment <- list()
if (requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  for (module_id in module_summary$Module[module_summary$Nodes >= 3]) {
    module_genes <- metrics$Symbol[metrics$Module == module_id]
    ego <- tryCatch(
      clusterProfiler::enrichGO(
        gene = module_genes,
        universe = unique(metrics$Symbol),
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL",
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 1,
        qvalueCutoff = 1,
        readable = FALSE
      ),
      error = function(e) NULL
    )
    if (!is.null(ego) && nrow(ego@result)) {
      x <- ego@result
      x$Module <- module_id
      module_enrichment[[as.character(module_id)]] <- x
      best <- x[order(x$p.adjust), , drop = FALSE][1, ]
      if (is.finite(best$p.adjust) && best$p.adjust < 0.05) {
        module_summary$TopFunction[module_summary$Module == module_id] <- best$Description
        module_summary$TopFunctionFDR[module_summary$Module == module_id] <- best$p.adjust
      }
    }
  }
}
module_summary$TopFunction[is.na(module_summary$TopFunction)] <- "No significant annotation"
write.csv(module_summary, file.path(ppi_dir, "PPI_Module_Summary.csv"), row.names = FALSE)
if (length(module_enrichment)) {
  write.csv(
    bind_rows(module_enrichment),
    file.path(ppi_dir, "PPI_Module_Enrichment_GO_BP.csv"),
    row.names = FALSE
  )
}

colors <- project_colors()
if (plot_enabled("ppi.module_overview")) {
  p <- ggplot(module_summary, aes(
    reorder(paste0("Module ", Module), Nodes),
    Nodes,
    fill = MeanLogFC
  )) +
    geom_col(width = 0.72) +
    coord_flip() +
    scale_fill_gradient2(low = colors[["down"]], mid = "white", high = colors[["up"]], midpoint = 0) +
    labs(
      title = "PPI functional modules",
      subtitle = "Louvain modules; color shows mean log2FC",
      x = NULL, y = "Genes", fill = "Mean log2FC"
    ) +
    project_theme()
  out <- file.path(ppi_dir, "PPI_Module_Overview.png")
  save_project_plot(p, out, width = 8, height = max(4.5, nrow(module_summary) * 0.45 + 1.5))
  write_figure_note(
    out,
    "PPI functional modules",
    "STRING Cytoscape edge and node tables",
    "Louvain community detection on confidence-weighted STRING interactions",
    "Modules with at least three genes are annotated by GO BP",
    "STRING + Gene Ontology"
  )
}

if (plot_enabled("ppi.module_networks")) {
  module_dir <- file.path(ppi_dir, "module_networks")
  dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
  for (module_id in module_summary$Module[module_summary$Nodes >= 5]) {
    selected <- metrics %>%
      filter(Module == module_id) %>%
      arrange(desc(HubScore)) %>%
      slice_head(n = 30)
    mg <- induced_subgraph(g, selected$Symbol)
    V(mg)$HubScore <- metrics$HubScore[match(V(mg)$name, metrics$Symbol)]
    V(mg)$logFC <- metrics$logFC[match(V(mg)$name, metrics$Symbol)]
    V(mg)$label <- ifelse(rank(-V(mg)$HubScore, ties.method = "first") <= 5, V(mg)$name, "")
    p <- ggraph(mg, layout = "fr") +
      geom_edge_link(aes(alpha = Score), color = "grey65", show.legend = FALSE) +
      geom_node_point(aes(size = HubScore, color = logFC)) +
      scale_color_gradient2(low = colors[["down"]], mid = "white", high = colors[["up"]], midpoint = 0) +
      scale_size_continuous(range = c(3, 9)) +
      geom_node_text(aes(label = label), repel = TRUE, size = 3.2, fontface = "bold") +
      theme_void() +
      labs(
        title = paste0("PPI Module ", module_id),
        subtitle = module_summary$TopFunction[module_summary$Module == module_id],
        color = "log2FC", size = "Hub score"
      )
    out <- file.path(module_dir, paste0("PPI_Module_", module_id, ".png"))
    save_project_plot(p, out, width = 7, height = 6)
    write_figure_note(
      out,
      paste0("PPI Module ", module_id),
      "STRING Cytoscape edge and node tables",
      "Louvain module; composite hub score uses degree, betweenness, PageRank, eigenvector centrality and |log2FC|",
      "Top five composite hubs are labelled",
      "STRING"
    )
  }
}

