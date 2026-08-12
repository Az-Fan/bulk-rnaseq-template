# Module 05 — frozen association-network views

run_network <- function(context) {
  decision <- module_decision(context$config, "network")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "05_Network", decision$status, decision$reason %||% ""))
  }
  resource_id <- context$config$resources$network %||% NULL
  if (is.null(resource_id)) stop("An enabled network module requires resources.network")
  path <- resource_path(context, resource_id)
  edges <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else read_project_table(path)
  required <- c("from_symbol", "to_symbol")
  if (length(setdiff(required, names(edges)))) stop("Network resource lacks from_symbol/to_symbol")
  profiles <- threshold_profiles(context$config)
  results <- list()
  for (contrast_id in names(context$state$differential)) {
    de <- context$state$differential[[contrast_id]]$result
    primary <- names(profiles)[[1L]]
    genes <- unique(de$gene_symbol[de[[paste0("significant_", primary)]]])
    genes <- genes[!is.na(genes) & nzchar(genes)]
    sub_edges <- edges[edges$from_symbol %in% genes & edges$to_symbol %in% genes, , drop = FALSE]
    root <- comparison_module_dir(context, contrast_id, "05_Network")
    tables <- file.path(root, "Tables"); figures <- file.path(root, "Figures")
    dir.create(tables, recursive = TRUE, showWarnings = FALSE); dir.create(figures, recursive = TRUE, showWarnings = FALSE)
    write_result_table(sub_edges, file.path(tables, "Network_Edges.tsv.gz"))
    if (!nrow(sub_edges)) {
      write_result_table(data.frame(status = "complete_no_edges", detail = "No frozen resource edges joined primary DE genes"),
                         file.path(tables, "Network_Status.tsv"))
      results[[contrast_id]] <- list(edges = sub_edges, nodes = data.frame())
      next
    }
    graph <- igraph::graph_from_data_frame(sub_edges, directed = FALSE)
    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
    communities <- igraph::cluster_louvain(graph, weights = if ("combined_score" %in% names(sub_edges)) igraph::E(graph)$combined_score else NULL)
    nodes <- data.frame(gene_symbol = igraph::V(graph)$name, degree = igraph::degree(graph),
                        betweenness = igraph::betweenness(graph, normalized = TRUE),
                        module = as.integer(igraph::membership(communities)), stringsAsFactors = FALSE)
    nodes <- merge(nodes, de[, c("gene_symbol", "log2FoldChange_ashr", "padj")], by = "gene_symbol", all.x = TRUE)
    nodes <- nodes[order(-nodes$degree, -nodes$betweenness), ]
    write_result_table(nodes, file.path(tables, "Network_Nodes_and_Hubs.tsv"))
    set.seed(as.integer(context$config$analysis$random_seed %||% 104729L))
    layout <- igraph::layout_with_fr(graph, niter = 1000)
    layout_table <- data.frame(gene_symbol = igraph::V(graph)$name, x = layout[, 1], y = layout[, 2])
    write_result_table(layout_table, file.path(tables, "Network_Layout.tsv"))
    edge_data <- igraph::as_data_frame(graph, what = "edges")
    edge_data$x <- layout_table$x[match(edge_data$from, layout_table$gene_symbol)]
    edge_data$y <- layout_table$y[match(edge_data$from, layout_table$gene_symbol)]
    edge_data$xend <- layout_table$x[match(edge_data$to, layout_table$gene_symbol)]
    edge_data$yend <- layout_table$y[match(edge_data$to, layout_table$gene_symbol)]
    node_data <- merge(layout_table, nodes, by = "gene_symbol", all.x = TRUE)
    labels <- head(nodes$gene_symbol, 12L)
    plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = edge_data, ggplot2::aes(x, y, xend = xend, yend = yend),
                            linewidth = 0.25, alpha = 0.22, colour = "#70777E") +
      ggplot2::geom_point(data = node_data, ggplot2::aes(x, y, size = degree, fill = log2FoldChange_ashr),
                          shape = 21, colour = "white", stroke = 0.35, alpha = 0.92) +
      ggrepel::geom_text_repel(data = node_data[node_data$gene_symbol %in% labels, ],
                               ggplot2::aes(x, y, label = gene_symbol), size = 2.4,
                               max.overlaps = Inf, seed = 104729, min.segment.length = 0) +
      ggplot2::scale_fill_gradient2(low = direction_palette[["Down"]], mid = "white", high = direction_palette[["Up"]], midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(2, 8)) + ggplot2::coord_equal() +
      ggplot2::labs(title = paste0(contrast_id, " — frozen association network"),
                    subtitle = "Node size: degree; colour: shrunken LFC; labels: highest-degree nodes",
                    fill = "ashr LFC", size = "Degree",
                    caption = "Layout uses a fixed seed; nodes and edges are unchanged. Network centrality does not establish causality.") +
      theme_publication() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                                          axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    save_publication_figure(plot, file.path(figures, "PPI_Module_Network"), 178, 150)
    results[[contrast_id]] <- list(edges = sub_edges, nodes = nodes, layout = layout_table)
  }
  context$state$network <- results
  record_module_status(context, "05_Network", "complete",
                       sprintf("Frozen network %s evaluated for %d contrasts", resource_id, length(results)))
}

