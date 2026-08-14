# Module 05 — legacy-compatible frozen STRING association network

run_network <- function(context) {
  decision <- module_decision(context$config, "network")
  if (!identical(decision$status, "enabled")) {
    return(record_module_status(context, "05_Network", decision$status, decision$reason %||% ""))
  }
  resource_id <- context$config$resources$network %||% NULL
  if (is.null(resource_id)) stop("An enabled network module requires resources.network")
  path <- resource_path(context, resource_id)
  edges <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else read_project_table(path)
  required <- c("from_symbol", "to_symbol", "combined_score")
  if (length(setdiff(required, names(edges)))) stop("STRING resource lacks from_symbol/to_symbol/combined_score")
  profiles <- threshold_profiles(context$config)
  network_config <- context$config$network %||% list()
  network_profile <- as.character(network_config$deg_profile %||% names(profiles)[[1L]])
  if (!network_profile %in% names(profiles)) stop("Unknown network.deg_profile: ", network_profile)
  max_input_genes <- as.integer(network_config$max_input_genes %||% 200L)
  minimum_string_score <- as.integer(network_config$minimum_string_score %||% 400L)
  module_fdr <- as.numeric(network_config$module_enrichment_fdr %||% 0.05)
  network_label_count <- as.integer(network_config$network_label_count %||% 30L)
  results <- list()
  for (contrast_id in names(context$state$differential)) {
    de <- context$state$differential[[contrast_id]]$result
    primary <- network_profile
    sig <- de[de[[paste0("significant_", primary)]] & !is.na(de$gene_symbol) & nzchar(de$gene_symbol), , drop = FALSE]
    sig <- sig[order(-abs(sig$log2FoldChange_ashr)), , drop = FALSE]
    sig <- sig[!duplicated(sig$gene_symbol), , drop = FALSE]
    target <- head(sig, min(max_input_genes, nrow(sig)))
    genes <- target$gene_symbol
    sub_edges <- edges[edges$from_symbol %in% genes & edges$to_symbol %in% genes &
                         edges$combined_score >= minimum_string_score, , drop = FALSE]
    root <- comparison_module_dir(context, contrast_id, "05_Network")
    tables <- file.path(root, "PPI", "Tables"); figures <- file.path(root, "PPI", "Figures")
    dir.create(file.path(tables, "cytoscape_input"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(figures, "module_networks"), recursive = TRUE, showWarnings = FALSE)
    write_result_table(target, file.path(tables, "PPI_Input_Selected_DEGs.tsv"))
    if (identical(max_input_genes, 200L)) {
      write_result_table(target, file.path(tables, "PPI_Input_Top200_DEGs.tsv"))
    }
    if (!nrow(sub_edges)) {
      write_result_table(data.frame(status = "complete_no_edges",
                                    detail = sprintf("No frozen STRING edges joined the selected top-%d DE genes at score >= %d",
                                                     max_input_genes, minimum_string_score)),
                         file.path(tables, "Network_Status.tsv"))
      results[[contrast_id]] <- list(edges = sub_edges, nodes = data.frame())
      next
    }
    graph <- igraph::graph_from_data_frame(
      sub_edges[, c("from_symbol", "to_symbol", "combined_score")], directed = FALSE
    )
    graph <- igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE,
                              edge.attr.comb = list(combined_score = "max"))
    weights <- igraph::E(graph)$combined_score
    weights[!is.finite(weights) | weights <= 0] <- 1
    names_graph <- igraph::V(graph)$name
    lfc <- target$log2FoldChange_ashr[match(names_graph, target$gene_symbol)]
    degree <- igraph::degree(graph)
    betweenness <- igraph::betweenness(graph, directed = FALSE, weights = 1 / weights, normalized = TRUE)
    closeness <- igraph::closeness(graph, weights = 1 / weights, normalized = TRUE)
    pagerank <- igraph::page_rank(graph, directed = FALSE, weights = weights)$vector
    eigenvector <- igraph::eigen_centrality(graph, directed = FALSE, weights = weights)$vector
    modules <- if (igraph::vcount(graph) >= 3L && igraph::ecount(graph) >= 2L) {
      igraph::membership(igraph::cluster_louvain(graph, weights = weights))
    } else rep(1L, igraph::vcount(graph))
    zscore <- function(x) {
      if (length(unique(x[is.finite(x)])) <= 1L) return(rep(0, length(x)))
      as.numeric(scale(x))
    }
    hub_score <- rowMeans(cbind(zscore(log1p(degree)), zscore(log1p(betweenness)),
                                zscore(pagerank), zscore(eigenvector), zscore(abs(lfc))), na.rm = TRUE)
    nodes <- data.frame(Symbol = names_graph, logFC = lfc, Degree = degree,
                        Betweenness = betweenness, Closeness = closeness,
                        PageRank = pagerank, Eigenvector = eigenvector,
                        Module = as.integer(modules), HubScore = hub_score,
                        stringsAsFactors = FALSE)
    nodes <- nodes[order(-nodes$HubScore), , drop = FALSE]
    edge_data <- igraph::as_data_frame(graph, what = "edges")
    cytoscape_edges <- data.frame(Source = edge_data$from, Target = edge_data$to,
                                  Score = edge_data$combined_score)
    write_result_table(nodes, file.path(tables, "PPI_Hub_Ranking.tsv"))
    write_result_table(nodes[, c("Symbol", "logFC", "Degree")], file.path(tables, "cytoscape_input", "Node_Attributes.tsv"))
    write_result_table(cytoscape_edges, file.path(tables, "cytoscape_input", "Network_Edges.tsv"))

    module_summary <- do.call(rbind, lapply(sort(unique(nodes$Module)), function(module_id) {
      x <- nodes[nodes$Module == module_id, , drop = FALSE]
      data.frame(Module = module_id, Nodes = nrow(x), MeanLogFC = mean(x$logFC, na.rm = TRUE),
                 UpGenes = sum(x$logFC > 0, na.rm = TRUE), DownGenes = sum(x$logFC < 0, na.rm = TRUE),
                 TopHub = x$Symbol[which.max(x$HubScore)], stringsAsFactors = FALSE)
    }))
    module_summary <- module_summary[order(-module_summary$Nodes), , drop = FALSE]
    module_summary$TopFunction <- "No significant annotation"
    module_summary$TopFunctionFDR <- NA_real_
    enrichment_list <- list()
    universe <- unique(nodes$Symbol)
    for (module_id in module_summary$Module[module_summary$Nodes >= 3L]) {
      module_genes <- nodes$Symbol[nodes$Module == module_id]
      ego <- tryCatch(clusterProfiler::enrichGO(
        gene = module_genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        keyType = "SYMBOL", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 1,
        qvalueCutoff = 1, readable = FALSE
      ), error = function(e) NULL)
      if (!is.null(ego) && nrow(ego@result)) {
        tab <- as.data.frame(ego); tab$Module <- module_id
        enrichment_list[[as.character(module_id)]] <- tab
        best <- tab[order(tab$p.adjust), , drop = FALSE][1, ]
        if (is.finite(best$p.adjust) && best$p.adjust < module_fdr) {
          module_summary$TopFunction[module_summary$Module == module_id] <- best$Description
          module_summary$TopFunctionFDR[module_summary$Module == module_id] <- best$p.adjust
        }
      }
    }
    write_result_table(module_summary, file.path(tables, "PPI_Module_Summary.tsv"))
    if (length(enrichment_list)) {
      write_result_table(do.call(rbind, enrichment_list), file.path(tables, "PPI_Module_Enrichment_GO_BP.tsv"))
    } else {
      write_result_table(data.frame(), file.path(tables, "PPI_Module_Enrichment_GO_BP.tsv"))
    }

    # Fixed layout is an explicitly controlled visual operation; graph membership
    # and all centrality statistics above do not depend on this layout.
    set.seed(as.integer(context$config$analysis$random_seed %||% 104729L))
    layout <- igraph::layout_with_fr(graph, niter = 1000)
    layout_table <- data.frame(Symbol = names_graph, x = layout[, 1], y = layout[, 2])
    write_result_table(layout_table, file.path(tables, "Network_Layout.tsv"))
    edge_plot <- cytoscape_edges
    edge_plot$x <- layout_table$x[match(edge_plot$Source, layout_table$Symbol)]
    edge_plot$y <- layout_table$y[match(edge_plot$Source, layout_table$Symbol)]
    edge_plot$xend <- layout_table$x[match(edge_plot$Target, layout_table$Symbol)]
    edge_plot$yend <- layout_table$y[match(edge_plot$Target, layout_table$Symbol)]
    node_plot <- merge(layout_table, nodes, by = "Symbol", all.x = TRUE)
    labels <- head(nodes$Symbol, min(network_label_count, nrow(nodes)))
    network_plot <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = edge_plot, ggplot2::aes(x, y, xend = xend, yend = yend),
                            linewidth = 0.3, alpha = 0.35, colour = "#70777E") +
      ggplot2::geom_point(data = node_plot, ggplot2::aes(x, y, size = Degree, fill = logFC),
                          shape = 21, colour = "white", stroke = 0.35, alpha = 0.95) +
      ggrepel::geom_text_repel(data = node_plot[node_plot$Symbol %in% labels, ],
                               ggplot2::aes(x, y, label = Symbol), size = 2.5,
                               max.overlaps = Inf, seed = 104729, min.segment.length = 0) +
      ggplot2::scale_fill_gradient2(low = direction_palette[["Down"]], mid = "white",
                                    high = direction_palette[["Up"]], midpoint = 0) +
      ggplot2::scale_size_continuous(range = c(2.5, 8)) + ggplot2::coord_equal() +
      ggplot2::labs(title = paste0(contrast_id, " — STRING PPI hub network"),
                    subtitle = sprintf("Top-%d %s DE input; STRING score >= %d; node size = degree; colour = ashr LFC",
                                       max_input_genes, network_profile, minimum_string_score),
                    fill = "ashr LFC", size = "Degree",
                    caption = "STRING edges are associations; hub status is not evidence of causality.") +
      theme_publication() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                                          axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
    save_publication_figure(network_plot, file.path(figures, "PPI_Hub_Network"), 178, 150)

    module_plot <- ggplot2::ggplot(module_summary,
      ggplot2::aes(stats::reorder(paste0("Module ", Module), Nodes), Nodes, fill = MeanLogFC)) +
      ggplot2::geom_col(width = 0.72) + ggplot2::coord_flip() +
      ggplot2::scale_fill_gradient2(low = direction_palette[["Down"]], mid = "white",
                                    high = direction_palette[["Up"]], midpoint = 0) +
      ggplot2::labs(title = "PPI functional modules", subtitle = "Louvain modules on frozen confidence-weighted STRING interactions",
                    x = NULL, y = "Genes", fill = "Mean ashr LFC") + theme_publication()
    save_publication_figure(module_plot, file.path(figures, "PPI_Module_Overview"), 150,
                            max(100, nrow(module_summary) * 11 + 40))

    for (module_id in module_summary$Module[module_summary$Nodes >= 5L]) {
      selected <- head(nodes[nodes$Module == module_id, , drop = FALSE], 30L)
      module_nodes <- node_plot[node_plot$Symbol %in% selected$Symbol, , drop = FALSE]
      module_edges <- edge_plot[edge_plot$Source %in% selected$Symbol & edge_plot$Target %in% selected$Symbol, , drop = FALSE]
      module_nodes$label <- ifelse(module_nodes$Symbol %in% head(selected$Symbol, 5L), module_nodes$Symbol, "")
      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = module_edges, ggplot2::aes(x, y, xend = xend, yend = yend),
                              colour = "#777D84", alpha = 0.4, linewidth = 0.3) +
        ggplot2::geom_point(data = module_nodes, ggplot2::aes(x, y, size = HubScore, fill = logFC),
                            shape = 21, colour = "white") +
        ggrepel::geom_text_repel(data = module_nodes[module_nodes$label != "", ],
                                 ggplot2::aes(x, y, label = label), seed = 104729, size = 2.7) +
        ggplot2::scale_fill_gradient2(low = direction_palette[["Down"]], mid = "white",
                                      high = direction_palette[["Up"]], midpoint = 0) +
        ggplot2::coord_equal() + ggplot2::labs(title = paste0("PPI Module ", module_id),
          subtitle = module_summary$TopFunction[module_summary$Module == module_id], fill = "ashr LFC", size = "Hub score") +
        theme_publication() + ggplot2::theme(axis.line = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
                                            axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank())
      save_publication_figure(p, file.path(figures, "module_networks", paste0("PPI_Module_", module_id)), 155, 135)
    }
    results[[contrast_id]] <- list(edges = cytoscape_edges, nodes = nodes, modules = module_summary, layout = layout_table)
  }
  context$state$network <- results
  record_module_status(context, "05_Network", "complete",
                       sprintf("STRING network %s evaluated for %d contrasts (%s; top %d; score >= %d)",
                               resource_id, length(results), network_profile, max_input_genes, minimum_string_score))
}
